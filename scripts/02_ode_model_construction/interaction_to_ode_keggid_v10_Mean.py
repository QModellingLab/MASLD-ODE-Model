#!/usr/bin/env python3
"""
Interaction to ODE Model Generator v9 + SBML L2V4 Output
=========================================================
New in v9 (upgraded from v8):

[Inhibition implementation fix (2026-03-13)]
  v8 issue: gene_inhibitors / pathway_inhibitors were parsed correctly, but
           generate_ode_system() never read them, so all inhibitory relationships were dropped.
  v9 fix:
    1. Added an inhibitors parameter to _act_eq()
       Inhibition factor: each inhibitor is multiplied by ksp^n / (ksp^n + inh_active^n)
       -> factor = 1 when inhibitor = 0 (no effect); factor = 0.5 when inhibitor = ksp
       -> multiple inhibitors are multiplied together (independent competitive inhibition)
    2. generate_ode_system() now passes gene_inhibitors / pathway_inhibitors through as well
    3. SBML generation now also includes the inhibitor modifier
  3 affected edges (PPrel inhibition explicitly annotated in hsa04932):
    AKT3  -> GSK3A  (inhibition)  -> Akt inhibits GSK_3 activation
    MAPK8 -> IRS1   (inhibition)  -> JNK1_2 inhibits IRS_1_2 activation (insulin resistance)
    SOCS3 -> IRS1   (inhibition)  -> SOCS3 inhibits IRS_1_2 activation (inflammation -> insulin resistance)
  Output filename changed from v8 -> v9

New in v8 (upgraded from v7):

[KEGG display name integration (2026-03-11)]
  v7 issue: ODE variable names used HGNC gene symbols (TNFRSF1A, NFKB1, AKT3),
           inconsistent with the KEGG hsa04932 diagram labels (TNFR1, NF-κB, Akt)
  v8 fix:
    1. Read the KEGG_display_name column from the v3 nodes file
    2. Built an HGNC -> sanitized KEGG variable-name mapping table
    3. ODE variables now use the KEGG display name (special characters converted to valid Python identifiers)
       α->a, β->b, γ->g, κ->kB, -->_, />_, space->removed
       e.g., NF-κB -> NF_kB, JNK1/2 -> JNK1_2, TGF-β1 -> TGF_b1
    4. Added a kegg_display_name column to gene_metadata
    5. Falls back automatically to the HGNC symbol when no v3 file is present (backward compatible)

New in v7 (upgraded from v6):

[Core change: initial-value design modeled on the previous two papers]
  v6 issue: gene initial values used absolute CPM -> condition differences <10% -> vanish after Hill-function saturation
  v7 fix:
    (1) Gene initial value = condition_expr / Normal_expr (ratio, Normal=1.0)
        -> amplifies condition differences so condition-specific dynamics remain visible in the ODE trajectory
    (2) Pathway output P_* initial inactive value fixed at 100
        -> fully consistent with the BIBE 2022 and Nutrition & Metabolism 2024 papers
        -> allows T50 to extend across tens of hours, revealing early/mid/late temporal differences

[Output model change: 3 ratio conditions x 5 IV types = 15 models]
  v6: Normal / Obese / NAFL / NASH (4 conditions, absolute values)
  v7: Obese_vs_Normal / NAFL_vs_Normal / NASH_vs_Normal (3 ratios)
      -> the Normal condition in every model is implicitly the 1.0 baseline

[Source of the ratio calculation]
  gene_iv_ratio = get_condition_iv(iv_type, condition) / get_condition_iv(iv_type, 'Normal')
  - Uses the same IV-type data sources as v6 (raw/normalized/log2_norm/log2_cpm/batch_corr)
  - Ratio values are NOT clipped (uniform no-clip policy; Zhu et al. 2019 argues against hard |log2FC| thresholds)

[Naming convention]
  ode_model_{iv_type}_{ratio_condition}_v7.py
  ratio_condition: Obese_vs_Normal | NAFL_vs_Normal | NASH_vs_Normal

Preserved from v6:
  - v3/v4 format auto-detection
  - SBML L2V4 generation
  - All node naming strategies
  - Backward-compatible imports
"""

import pandas as pd
import numpy as np
from collections import defaultdict
import re
import os

# =============================================================================
# Configuration
# =============================================================================

IV_TYPE_CONFIG = {
    'raw': {
        'label':    'Raw Counts',
        'file_key': 'v3',            # v3 unified file
        'sheet':    'Raw_Counts',
        'cols':     {'Normal': 'Normal_raw', 'Obese': 'Obese_raw',
                     'NAFL':   'NAFL_raw',   'NASH':  'NASH_raw'},
        'transform': None,
        'clip_min':  1e-6,
    },
    'normalized': {
        'label':    'Size-Factor Normalized',
        'file_key': 'v3',            # v3 unified file
        'sheet':    'Normalized_Counts',
        'cols':     {'Normal': 'Normal_norm', 'Obese': 'Obese_norm',
                     'NAFL':   'NAFL_norm',   'NASH':  'NASH_norm'},
        'transform': None,
        'clip_min':  1e-6,
    },
    'log2_norm': {
        'label':    'Log2-Norm -> Linear (2^x - 1)',
        'file_key': 'v3',            # v3 unified file
        'sheet':    'Log2_Norm',     # was Log2_Normalized in v1
        'cols':     {'Normal': 'Normal_log2norm', 'Obese': 'Obese_log2norm',
                     'NAFL':   'NAFL_log2norm',   'NASH':  'NASH_log2norm'},
        'transform': 'antilog2_minus1',
        'clip_min':  1e-6,
    },
    'log2_cpm': {
        'label':    'Log2-CPM -> Linear (2^x - 0.5), no BC',
        'file_key': 'v3',            # v3 unified file
        'sheet':    'Log2_CPM',      # was Log2_CPM_before_BC in v2
        'cols':     {'Normal': 'Normal_log2cpm', 'Obese': 'Obese_log2cpm',
                     'NAFL':   'NAFL_log2cpm',   'NASH':  'NASH_log2cpm'},
        'transform': 'antilog2_minus05',
        'clip_min':  1e-6,
    },
    'batch_corr': {
        'label':    'Batch-Corrected Linear [RECOMMENDED]',
        'file_key': 'v3',            # v3 unified file
        'sheet':    'ODE_Initial_Values',
        'cols':     {'Normal': 'Normal', 'Obese': 'Obese',
                     'NAFL':   'NAFL',   'NASH':  'NASH'},
        'transform': None,
        'clip_min':  1e-6,
    },
    'pydeseq2': {
        'label':    'pyDESeq2 ratio (2^log2FC) [UNIFIED]',
        'file_key': 'deg',           # reads from DEG Excel directly
        'sheet':    'Merged_3grp',
        'cols':     {'Gene': 'Gene',
                     'Obese': 'ratio_Obese',
                     'NAFL':  'ratio_NAFL',
                     'NASH':  'ratio_NASH'},
        'transform': None,           # ratios are already in linear scale
        'clip_min':  1e-6,
    },
}

# v7 conditions: ratio-based (condition / Normal)
ALL_IV_TYPES       = ['raw', 'normalized', 'log2_norm', 'log2_cpm', 'batch_corr',
                      'pydeseq2']
ALL_RATIO_CONDITIONS = ['Obese_vs_Normal', 'NAFL_vs_Normal', 'NASH_vs_Normal']
RATIO_TO_DENOM     = 'Normal'   # denominator is always Normal
RATIO_NUMERATOR    = {          # map ratio condition → numerator raw condition
    'Obese_vs_Normal': 'Obese',
    'NAFL_vs_Normal':  'NAFL',
    'NASH_vs_Normal':  'NASH',
}

# Ratio clip: DISABLED (uniform no-clip policy). Zhu et al. 2019 argues against hard |log2FC| thresholds.
# Bounds set effectively unbounded so true 2^log2FC values are preserved for all nodes and all cohorts.
RATIO_CLIP_MIN = 1e-6     # effectively no lower clip (uniform no-clip policy)
RATIO_CLIP_MAX = 1e6      # effectively no upper clip (uniform no-clip policy)

# Pathway output initial inactive (fixed = 100, per prior papers)
PATHWAY_INITIAL_INACTIVE = 100.0


def apply_iv_transform(value, transform, clip_min=1e-6):
    """Apply back-transform and floor-clip to a raw IV value."""
    try:
        val = float(value)
    except (TypeError, ValueError):
        val = 1.0
    if transform == 'antilog2_minus1':
        val = 2.0 ** val - 1.0
    elif transform == 'antilog2_minus05':
        val = 2.0 ** val - 0.5
    return max(val, clip_min)


# =============================================================================
# Main Converter Class
# =============================================================================

class InteractionToODEConverter_V7:
    """
    Converts interaction network to v8 ODE models.

    Design:
      - 5 IV types × 3 ratio conditions = 15 model files
      - Gene initial value = condition_expr / Normal_expr (ratio, Normal=1.0)
      - Pathway P_* initial inactive = 100 (fixed)

    Usage:
      conv = InteractionToODEConverter_V7(
          interactions_file='result_v4_symbol_interactions.xlsx',
          nodes_file='Node_Initial_Levels_4grp_v2_batchcorr.xlsx',
          nodes_file_v1='Node_Initial_Levels_4grp.xlsx',
          ratio_condition='NASH_vs_Normal',
          iv_type='batch_corr',
      )
      conv.load_interactions().load_nodes().parse_interactions()
      conv.generate_python_model('ode_model_batch_corr_NASH_vs_Normal_v8.py')
    """

    def __init__(self, interactions_file, nodes_file=None, nodes_file_v1=None,
                 nodes_file_v3=None, deg_file=None, node_name_table=None,
                 default_params=None, node_naming='symbol',
                 ratio_condition='NASH_vs_Normal', iv_type='batch_corr',
                 averaging_method='mean'):
        if iv_type not in IV_TYPE_CONFIG:
            raise ValueError(f"iv_type must be one of {list(IV_TYPE_CONFIG.keys())}, got '{iv_type}'")
        if ratio_condition not in ALL_RATIO_CONDITIONS:
            raise ValueError(
                f"ratio_condition must be one of {ALL_RATIO_CONDITIONS}, "
                f"got '{ratio_condition}'")

        self.interactions_file = interactions_file
        self.nodes_file_v2     = nodes_file
        self.nodes_file_v1     = nodes_file_v1
        self.nodes_file_v3     = nodes_file_v3
        self.deg_file          = deg_file
        self.node_name_table   = node_name_table
        self.node_naming       = node_naming
        self.ratio_condition   = ratio_condition
        self.iv_type           = iv_type
        if averaging_method not in ('mean', 'median'):
            raise ValueError(f"averaging_method must be 'mean' or 'median', got '{averaging_method}'")
        self.averaging_method  = averaging_method
        self.format_version    = None
        self.input_mode        = None

        # Resolved numerator condition (e.g. 'NASH')
        self.num_condition = RATIO_NUMERATOR[ratio_condition]

        self.interactions        = None
        self.nodes_data          = None
        self.genes               = set()
        self.pathways            = set()
        self.gene_metadata       = {}
        self.pathway_metadata    = {}
        self.gene_activators     = defaultdict(list)
        self.gene_inhibitors     = defaultdict(list)
        self.pathway_activators  = defaultdict(list)
        self.pathway_inhibitors  = defaultdict(list)
        self.initial_ratios      = {}   # gene_node -> ratio value
        self.kegg_display_map    = {}   # hgnc_symbol -> sanitized KEGG var name (v8)

        self.default_params = default_params or {
            'Vmax': 2.0, 'ksp': 2.0, 'n': 2.0, 'kcat': 2.0
        }

    # ------------------------------------------------------------------
    # Static helpers (identical to v6)
    # ------------------------------------------------------------------
    @staticmethod
    def sanitize_symbol(symbol):
        name = re.sub(r'[^a-zA-Z0-9_]', '_', str(symbol))
        if name and name[0].isdigit():
            name = 'gene_' + name
        return name

    @staticmethod
    def sanitize_kegg_display(kegg_name):
        """Convert KEGG display name → valid Python identifier.

        Rules (applied in order):
          1. Greek letters: α→a, β→b, γ→g, κ→kB
          2. Remove spaces
          3. Non-alphanumeric (except _) → _
          4. Collapse repeated _ and strip leading/trailing _
          5. Digit start → prefix 'gene_'

        Examples:
          NF-κB    → NF_kB     JNK1/2  → JNK1_2
          TGF-β1   → TGF_b1    IRS-1/2 → IRS_1_2
          PPAR-α   → PPAR_a    C/EBPα  → C_EBPa
          SREBP-1c → SREBP_1c  Cx I    → CxI
          eIF2α    → eIF2a     IKKβ    → IKKb
        """
        name = str(kegg_name)
        greek = [('κB','kB'), ('κ','k'), ('α','a'), ('β','b'), ('γ','g'),
                 ('δ','d'), ('ε','e'), ('ζ','z')]
        for gk, lat in greek:
            name = name.replace(gk, lat)
        name = name.replace(' ', '')                       # remove spaces
        name = re.sub(r'[^a-zA-Z0-9_]', '_', name)       # special → _
        name = re.sub(r'_+', '_', name).strip('_')        # collapse _ and strip
        if name and name[0].isdigit():
            name = 'gene_' + name
        return name or 'unknown'

    @staticmethod
    def sanitize_pathway_name(pathway_name):
        name = str(pathway_name).replace('TITLE:', '').strip()
        name = re.sub(r'[^a-zA-Z0-9_]', '_', name)
        name = re.sub(r'_+', '_', name).strip('_')
        if name and name[0].isdigit():
            name = 'pathway_' + name
        return name

    @staticmethod
    def sanitize_kegg_id(kegg_id):
        name = str(kegg_id).replace(':', '_')
        name = re.sub(r'[^a-zA-Z0-9_]', '_', name)
        if name and name[0].isdigit():
            name = 'gene_' + name
        return name

    # ------------------------------------------------------------------
    # Format detection (identical to v6)
    # ------------------------------------------------------------------
    def detect_format_version(self, df):
        cols = df.columns.tolist()
        if 'source_node' in cols and 'target_node' in cols:
            if 'gene1_symbol' in cols or 'species1_gene1' in cols:
                return 'v3'
            return 'v4'
        if 'gene1_symbol' in cols or 'species1_gene1' in cols:
            return 'v3'
        raise ValueError(f"Cannot determine format. Columns: {cols}")

    def detect_input_mode_v3(self, df):
        cols = df.columns.tolist()
        if 'gene1_symbol' in cols and 'gene2_symbol' in cols:
            return 'symbol'
        if 'species1_gene1' in cols and 'species1_gene2' in cols:
            return 'individual'
        raise ValueError(f"Unknown v3 format columns: {cols}")

    # ------------------------------------------------------------------
    # Load interactions
    # ------------------------------------------------------------------
    def load_interactions(self):
        print("Loading interaction data...")
        if not os.path.exists(self.interactions_file):
            raise FileNotFoundError(f"Not found: {self.interactions_file}")

        ext = os.path.splitext(self.interactions_file)[1].lower()
        if ext == '.xlsx':
            df = pd.read_excel(self.interactions_file, engine='openpyxl')
        elif ext == '.csv':
            df = pd.read_csv(self.interactions_file, encoding='utf-8-sig')
        else:
            df = pd.read_csv(self.interactions_file, encoding='utf-8-sig', sep='\t')

        print(f"  Loaded {len(df)} rows | Columns: {', '.join(df.columns.tolist())}")
        self.format_version = self.detect_format_version(df)
        print(f"  Format version: {self.format_version}")

        if self.format_version == 'v3':
            self.input_mode = self.detect_input_mode_v3(df)
            print(f"  v3 input mode: {self.input_mode}")

        self.interactions = df
        print(f"  Valid interactions: {len(self.interactions)}")
        return self

    # ------------------------------------------------------------------
    # Load nodes — v7 computes ratio = condition / Normal
    # ------------------------------------------------------------------
    def load_nodes(self):
        """
        Load gene expression values and compute ratio = num_condition / Normal.
        Both Normal and num_condition columns are read from the same sheet.
        Ratio is clipped to [RATIO_CLIP_MIN, RATIO_CLIP_MAX].
        Normal's ratio is implicitly 1.0.
        """
        if self.format_version != 'v4':
            print("  Skipping nodes file (not v4 format)")
            return self

        cfg       = IV_TYPE_CONFIG[self.iv_type]
        file_key  = cfg['file_key']

        # ==============================================================
        # pyDESeq2 mode: read ratios directly from DEG Excel
        # ==============================================================
        if file_key == 'deg':
            return self._load_nodes_pydeseq2(cfg)

        # ==============================================================
        # Original mode: read from Node_Initial_Levels_4grp Excel
        # ==============================================================
        sheet     = cfg['sheet']
        col_map   = cfg['cols']
        transform = cfg['transform']
        clip_min  = cfg['clip_min']

        # v3 is the unified file (all 5 IV types)
        if file_key == 'v3':
            if self.nodes_file_v3 and os.path.exists(self.nodes_file_v3):
                nodes_path = self.nodes_file_v3
            elif self.nodes_file_v2 and os.path.exists(self.nodes_file_v2):
                nodes_path = self.nodes_file_v2
                print(f"  [WARN] nodes_file_v3 not found, falling back to v2")
            else:
                nodes_path = self.nodes_file_v1
                print(f"  [WARN] nodes_file_v3/v2 not found, falling back to v1")
        elif file_key == 'v2':
            nodes_path = self.nodes_file_v2
        else:
            nodes_path = self.nodes_file_v1

        print(f"\nLoading nodes for ratio | iv_type={self.iv_type} ({cfg['label']})")
        print(f"  File            : {nodes_path}")
        print(f"  Sheet           : {sheet}")
        print(f"  Ratio condition : {self.ratio_condition}  "
              f"({self.num_condition} / Normal)")

        if not os.path.exists(nodes_path):
            raise FileNotFoundError(f"Nodes file not found: {nodes_path}")

        df = pd.read_excel(nodes_path, sheet_name=sheet, engine='openpyxl')
        print(f"  Rows  : {len(df)} | Columns: {', '.join(df.columns.tolist())}")

        if 'symbol' not in df.columns and 'node_name' in df.columns:
            df['symbol'] = df['node_name']
            print("  Using 'node_name' as 'symbol' (fallback)")

        # Build KEGG display name map (v8): hgnc_symbol → sanitized Python var
        if 'KEGG_display_name' in df.columns:
            n_mapped = 0
            for _, row in df.iterrows():
                hgnc = str(row.get('node_name', row.get('symbol', '')))
                kegg = str(row.get('KEGG_display_name', ''))
                if hgnc and kegg and kegg != 'nan' and kegg != hgnc:
                    self.kegg_display_map[hgnc] = self.sanitize_kegg_display(kegg)
                    n_mapped += 1
            print(f"  KEGG display map: {n_mapped} nodes with display name")
        else:
            print("  KEGG_display_name column not found — using HGNC symbols as-is")

        # Resolve Normal column
        col_normal = col_map.get('Normal')
        if col_normal not in df.columns:
            for fb in ['Normal', 'Normal_bc', 'Normal_raw', 'Normal_norm',
                       'Normal_log2', 'Normal_log2cpm']:
                if fb in df.columns:
                    col_normal = fb; break

        # Resolve numerator column
        col_num = col_map.get(self.num_condition)
        if col_num not in df.columns:
            for fb in [self.num_condition,
                       f'{self.num_condition}_bc',
                       f'{self.num_condition}_raw',
                       f'{self.num_condition}_norm',
                       f'{self.num_condition}_log2',
                       f'{self.num_condition}_log2cpm']:
                if fb in df.columns:
                    col_num = fb; break

        if col_normal is None or col_normal not in df.columns:
            raise ValueError(f"Cannot find Normal column in sheet='{sheet}'")
        if col_num is None or col_num not in df.columns:
            raise ValueError(
                f"Cannot find column for '{self.num_condition}' in sheet='{sheet}'")

        print(f"  Normal col : '{col_normal}' | Numerator col : '{col_num}'")
        self.nodes_data = df

        loaded = ratio_min = ratio_max = 0
        ratio_vals = []
        for _, row in df.iterrows():
            symbol = row.get('symbol', row.get('node_name', ''))
            if not symbol or pd.isna(symbol):
                continue
            hgnc_raw  = str(row.get('node_name', symbol))
            # v8: prefer KEGG display name as key; fallback to HGNC sanitized
            node_name = self.kegg_display_map.get(
                hgnc_raw, self.sanitize_symbol(str(symbol)))

            raw_normal = row.get(col_normal, 1.0)
            raw_num    = row.get(col_num,    1.0)
            if pd.isna(raw_normal): raw_normal = 1.0
            if pd.isna(raw_num):    raw_num    = 1.0

            val_normal = apply_iv_transform(raw_normal, transform, clip_min)
            val_num    = apply_iv_transform(raw_num,    transform, clip_min)

            # Ratio = condition / Normal, clipped
            ratio = np.clip(val_num / val_normal, RATIO_CLIP_MIN, RATIO_CLIP_MAX)
            self.initial_ratios[node_name] = float(ratio)
            ratio_vals.append(float(ratio))
            loaded += 1

        if ratio_vals:
            print(f"  Loaded {loaded} ratio values | "
                  f"min={min(ratio_vals):.3f} median={np.median(ratio_vals):.3f} "
                  f"max={max(ratio_vals):.3f}")
            n_up   = sum(1 for v in ratio_vals if v > 1.05)
            n_down = sum(1 for v in ratio_vals if v < 0.95)
            n_nc   = loaded - n_up - n_down
            print(f"  Up (>1.05): {n_up} | Down (<0.95): {n_down} | ~NC: {n_nc}")
        return self

    # ------------------------------------------------------------------
    # pyDESeq2 mode: read ratios from DEG Excel + node mapping from v4
    # ------------------------------------------------------------------

    # KEGG display names for hsa04932 nodes (HGNC → ODE variable name)
    # DEPRECATED in v10: use node_name_table Excel instead.
    # Kept as fallback when node_name_table is not provided.
    _KEGG_DISPLAY = {
        'ADIPOQ':'ACDC','PRKAG2':'AMPK','FOS':'AP_1','MAP3K5':'ASK1',
        'AKT3':'Akt','BAX':'Bax','BID':'Bid','BCL2L11':'Bim',
        'DDIT3':'CHOP','CEBPA':'C_EBPa','CDC42':'Cdc42','MLXIP':'ChREBP',
        'NDUFC2-KCTD14':'CxI','SDHA':'CxII','UQCR11':'CxIII','COX6B2':'CxIV',
        'CYCS':'Cytc','FAS':'Fas','FASLG':'FasL','GSK3A':'GSK_3',
        'IKBKB':'IKKb','IL1A':'IL_1','IL6':'IL_6','IL6R':'IL_6R',
        'CXCL8':'IL_8','ERN1':'IRE1a','IRS1':'IRS_1_2','MAPK8':'JNK1_2',
        'NR1H3':'LXR_a','PKLR':'L_PK','MAP3K11':'MLK3','NFKB1':'NF_kB',
        'LEPR':'ObR','EIF2AK3':'PERK','P3R3URF-PIK3R3':'PI3K',
        'PPARA':'PPAR_a','PPARG':'PPAR_g','RXRA':'RXR','RAC1':'Rac1',
        'SREBF1':'SREBP_1c','TGFB1':'TGF_b1','TNFRSF1A':'TNFR1',
        'TNF':'TNFa','ADIPOR1':'adipoR','EIF2S1':'eIF2a','MAPK14':'p38',
    }

    def _load_nodes_pydeseq2(self, cfg):
        """
        Load ODE initial ratios from pyDESeq2 results (DEG Excel).

        v10: If node_name_table is provided, use multi-gene averaging
             (consistent with N&M 2024 methodology).
             If not provided, fall back to v9 single-gene behavior.

        node_name_table Excel columns:
          source_node, kegg_display_name, ode_variable_name,
          n_genes, kegg_gene_ids, gene_symbols, entrez_ids, ensembl_ids
        """
        # --- 1. Resolve DEG file ---
        deg_path = self.deg_file
        if not deg_path or not os.path.exists(deg_path):
            raise FileNotFoundError(
                f"DEG file not found: {deg_path}\n"
                f"  For iv_type='pydeseq2', pass deg_file= to constructor.")

        ratio_col = cfg['cols'].get(self.num_condition)
        if not ratio_col:
            raise ValueError(
                f"No ratio column for '{self.num_condition}' in pydeseq2 config")

        use_nnt = (self.node_name_table is not None
                   and os.path.exists(self.node_name_table))
        avg_fn = np.median if self.averaging_method == 'median' else np.mean
        mode_label = (f'multi-gene {self.averaging_method} (v10)'
                      if use_nnt else 'single-gene (v9 fallback)')

        print(f"\nLoading nodes for ratio | iv_type=pydeseq2 ({cfg['label']})")
        print(f"  DEG file  : {deg_path}")
        print(f"  Ratio col : {ratio_col}  (condition: {self.num_condition})")
        print(f"  Mode      : {mode_label}")
        print(f"  Averaging : {self.averaging_method}  (NaN → 1.0 imputation)")
        if use_nnt:
            print(f"  NNT file  : {self.node_name_table}")

        # --- 2. Read DEG Excel → gene→ratio dict (all 3 conditions) ---
        deg_df = pd.read_excel(deg_path, sheet_name=cfg['sheet'], engine='openpyxl')
        gene_col = cfg['cols']['Gene']

        # Build per-condition ratio lookups
        ratio_lookup = {}
        for cond, col_name in cfg['cols'].items():
            if cond == 'Gene':
                continue
            ratio_lookup[cond] = dict(zip(
                deg_df[gene_col].astype(str),
                deg_df[col_name]))
        print(f"  DEG genes : {len(deg_df)} rows loaded")

        # Primary lookup for current condition
        gene_to_ratio = ratio_lookup.get(self.num_condition, {})

        # --- 3. Build kegg_display_map and gene lists ---
        if use_nnt:
            # ===== v10: Read from node_name_table =====
            nnt_df = pd.read_excel(self.node_name_table, engine='openpyxl')
            required_cols = ['source_node', 'ode_variable_name', 'gene_symbols']
            for c in required_cols:
                if c not in nnt_df.columns:
                    raise ValueError(
                        f"node_name_table missing column: '{c}'\n"
                        f"  Required: {required_cols}")

            # Build kegg_display_map from ode_variable_name column
            self.kegg_display_map = {}
            node_gene_lists = {}  # source_node → [gene_symbol, ...]
            for _, nrow in nnt_df.iterrows():
                src = str(nrow['source_node'])
                ode_var = str(nrow['ode_variable_name'])
                gene_syms_str = str(nrow.get('gene_symbols', ''))

                self.kegg_display_map[src] = ode_var

                # Parse gene list, skip MT-* and ?(xxx) entries
                genes = []
                for g in gene_syms_str.split(','):
                    g = g.strip()
                    if not g or g == 'nan':
                        continue
                    if g.startswith('MT-') or g.startswith('?') or g.endswith('*'):
                        continue
                    genes.append(g)
                node_gene_lists[src] = genes

            print(f"  NNT nodes : {len(nnt_df)} entries loaded")
            print(f"  Multi-gene: {sum(1 for v in node_gene_lists.values() if len(v) > 1)}"
                  f" nodes with >1 gene")
        else:
            # ===== v9 fallback: use _KEGG_DISPLAY =====
            self.kegg_display_map = {
                k: self.sanitize_kegg_display(v)
                for k, v in self._KEGG_DISPLAY.items()
            }
            node_gene_lists = {}
            print(f"  KEGG map  : {len(self._KEGG_DISPLAY)} built-in entries (v9 fallback)")

        # --- 4. Collect gene nodes from interaction file ---
        if self.interactions is None:
            raise RuntimeError("load_interactions() must be called before load_nodes()")

        all_gene_symbols = set()
        for _, row in self.interactions.iterrows():
            for col in ['source_node', 'target_node']:
                node = str(row.get(col, ''))
                if node and node != 'nan' and not node.startswith('P_'):
                    all_gene_symbols.add(node)
        print(f"  Gene nodes from interactions: {len(all_gene_symbols)}")

        # --- 5. Look up ratios (v10: multi-gene average) ---
        loaded = 0
        ratio_vals = []
        v4_records = []

        for symbol in sorted(all_gene_symbols):
            node_name = self.kegg_display_map.get(
                symbol, self.sanitize_symbol(symbol))

            # Get gene list for this node
            gene_list = node_gene_lists.get(symbol, [symbol])

            # Compute ratio for current condition (NaN → 1.0 imputation)
            cond_ratios = []
            for g in gene_list:
                r = gene_to_ratio.get(g)
                if r is None or pd.isna(r):
                    r = 1.0  # NaN imputation: assume no change
                cond_ratios.append(float(np.clip(r, RATIO_CLIP_MIN, RATIO_CLIP_MAX)))
            ratio = float(avg_fn(cond_ratios)) if cond_ratios else 1.0
            ratio = float(np.clip(ratio, RATIO_CLIP_MIN, RATIO_CLIP_MAX))

            self.initial_ratios[node_name] = ratio
            ratio_vals.append(ratio)
            loaded += 1

            # Collect all 3 conditions for output Excel
            rec = {
                'node_name': node_name,
                'symbol': symbol,
                'n_genes_in_node': len(gene_list),
                'gene_list': ','.join(gene_list),
            }
            for cond in ['Obese', 'NAFL', 'NASH']:
                cond_lookup = ratio_lookup.get(cond, {})
                cond_vals = []
                for g in gene_list:
                    r = cond_lookup.get(g)
                    if r is None or pd.isna(r):
                        r = 1.0  # NaN imputation
                    cond_vals.append(float(np.clip(r, RATIO_CLIP_MIN, RATIO_CLIP_MAX)))
                avg_r = float(avg_fn(cond_vals)) if cond_vals else 1.0
                rec[f'ratio_{cond}'] = float(np.clip(avg_r, RATIO_CLIP_MIN, RATIO_CLIP_MAX))
            v4_records.append(rec)

        if ratio_vals:
            print(f"  Loaded {loaded} ratio values | "
                  f"min={min(ratio_vals):.3f} median={np.median(ratio_vals):.3f} "
                  f"max={max(ratio_vals):.3f}")
            n_up   = sum(1 for v in ratio_vals if v > 1.05)
            n_down = sum(1 for v in ratio_vals if v < 0.95)
            n_nc   = loaded - n_up - n_down
            print(f"  Up (>1.05): {n_up} | Down (<0.95): {n_down} | ~NC: {n_nc}")

        # --- 6. Generate Node_Initial_Levels_pydeseq2_mean.xlsx ---
        try:
            v4_df = pd.DataFrame(v4_records)
            v4_out = os.path.join(os.path.dirname(deg_path),
                                  'Node_Initial_Levels_pydeseq2_mean.xlsx')
            with pd.ExcelWriter(v4_out, engine='openpyxl') as writer:
                v4_df.to_excel(writer, sheet_name='ODE_Ratios', index=False)
                meta = pd.DataFrame([
                    ['Source', os.path.basename(deg_path)],
                    ['Method', 'pyDESeq2 ratio = 2^(log2FoldChange)'],
                    ['Averaging', f'{self.averaging_method} (NaN→1.0 imputation)'],
                    ['Clip_range', f'[{RATIO_CLIP_MIN}, {RATIO_CLIP_MAX}]'],
                    ['Normal_ratio', '1.0 (implicit)'],
                    ['P_*_inactive', str(PATHWAY_INITIAL_INACTIVE)],
                    ['NNT_file', str(self.node_name_table) if use_nnt else 'N/A'],
                    ['Script', 'interaction_to_ode_keggid_v10_Mean.py'],
                    ['Date', pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')],
                ], columns=['Parameter', 'Value'])
                meta.to_excel(writer, sheet_name='Metadata', index=False)
            print(f"  Output    : {v4_out}")
        except Exception as e:
            print(f"  [WARN] Output Excel failed: {e}")

        return self

    # ------------------------------------------------------------------
    # Parse interactions (identical to v6)
    # ------------------------------------------------------------------
    def parse_interactions(self):
        print("\nParsing interactions...")
        if self.format_version == 'v3':
            if self.input_mode == 'symbol':
                self._parse_v3_symbol()
            else:
                self._parse_v3_individual()
        elif self.format_version == 'v4':
            self._parse_v4()
        return self

    def get_kegg_var(self, hgnc_symbol):
        """Return sanitized KEGG display name for ODE variable, or fallback to HGNC."""
        return self.kegg_display_map.get(str(hgnc_symbol),
                                         self.sanitize_symbol(str(hgnc_symbol)))

    def _parse_v4(self):
        print("  v4 format...")
        n_gene = n_path = 0
        for _, row in self.interactions.iterrows():
            src   = row.get('source_node', '')
            tgt   = row.get('target_node', '')
            itype = str(row.get('interaction_type', '')).lower()
            subs  = str(row.get('subtypes', '')).lower()

            if not src or not tgt or src == 'nan' or tgt == 'nan' \
                    or pd.isna(src) or pd.isna(tgt):
                continue

            # Pathway output: target starts with P_ (interaction_type may be GARel or None)
            if str(tgt).startswith('P_') or itype == 'garel':
                pnode = self.sanitize_pathway_name(tgt)
                snode = self.get_kegg_var(src)   # v8: KEGG display name
                self.pathways.add(pnode)
                self.genes.add(snode)
                self.pathway_metadata[pnode] = {
                    'pathway_name':  tgt,
                    'type':          'pathway_output',
                    'initial_level': PATHWAY_INITIAL_INACTIVE,
                }
                self.gene_metadata.setdefault(snode, {
                    'symbol':            src,
                    'kegg_display_name': self.kegg_display_map.get(src, src),
                    'initial_ratio':     self.initial_ratios.get(snode, 1.0),
                })
                if 'inhibition' in subs:
                    self.pathway_inhibitors[pnode].append(snode)
                else:
                    self.pathway_activators[pnode].append(snode)
                n_path += 1
            else:
                snode = self.get_kegg_var(src)   # v8: KEGG display name
                tnode = self.get_kegg_var(tgt)   # v8: KEGG display name
                self.genes.update([snode, tnode])
                self.gene_metadata.setdefault(snode, {
                    'symbol':            src,
                    'kegg_display_name': self.kegg_display_map.get(src, src),
                    'initial_ratio':     self.initial_ratios.get(snode, 1.0),
                })
                self.gene_metadata.setdefault(tnode, {
                    'symbol':            tgt,
                    'kegg_display_name': self.kegg_display_map.get(tgt, tgt),
                    'initial_ratio':     self.initial_ratios.get(tnode, 1.0),
                })
                if 'inhibition' in subs or 'repression' in subs:
                    self.gene_inhibitors[tnode].append(snode)
                else:
                    self.gene_activators[tnode].append(snode)
                n_gene += 1

        for d in [self.gene_activators, self.gene_inhibitors,
                  self.pathway_activators, self.pathway_inhibitors]:
            for k in d:
                d[k] = list(set(d[k]))

        print(f"  Genes: {len(self.genes)} | Pathways: {len(self.pathways)}")
        print(f"  Gene interactions: {n_gene} | Pathway interactions: {n_path}")

    def _parse_v3_symbol(self):
        print("  v3 symbol mode...")
        n_gene = n_path = 0
        for _, row in self.interactions.iterrows():
            src   = row.get('gene1_symbol', '')
            tgt   = row.get('gene2_symbol', '')
            itype = str(row.get('interaction_type', '')).lower()
            subs  = str(row.get('subtypes', '')).lower()
            init  = row.get('initial_level', None)
            if pd.isna(init): init = None

            if not src or not tgt or src == 'nan' or tgt == 'nan' \
                    or pd.isna(src) or pd.isna(tgt):
                continue

            if itype == 'garel':
                pnode = self.sanitize_pathway_name(str(tgt).strip())
                snode = self.get_kegg_var(src)   # v8: KEGG display name
                self.pathways.add(pnode)
                self.genes.add(snode)
                self.pathway_metadata[pnode] = {
                    'pathway_name':  str(tgt).strip(),
                    'type':          'pathway_output',
                    'initial_level': PATHWAY_INITIAL_INACTIVE,
                }
                self.gene_metadata.setdefault(snode, {
                    'symbol':            src,
                    'kegg_display_name': self.kegg_display_map.get(src, src),
                    'initial_ratio':     self.initial_ratios.get(snode, 1.0),
                })
                if 'inhibition' in subs:
                    self.pathway_inhibitors[pnode].append(snode)
                else:
                    self.pathway_activators[pnode].append(snode)
                n_path += 1
            else:
                snode = self.get_kegg_var(src)   # v8: KEGG display name
                tnode = self.get_kegg_var(tgt)   # v8: KEGG display name
                self.genes.update([snode, tnode])
                self.gene_metadata.setdefault(snode, {
                    'symbol':            src,
                    'kegg_display_name': self.kegg_display_map.get(src, src),
                    'initial_ratio':     self.initial_ratios.get(snode, 1.0),
                })
                self.gene_metadata.setdefault(tnode, {
                    'symbol':            tgt,
                    'kegg_display_name': self.kegg_display_map.get(tgt, tgt),
                    'initial_ratio':     self.initial_ratios.get(tnode, 1.0),
                })
                if 'inhibition' in subs or 'repression' in subs:
                    self.gene_inhibitors[tnode].append(snode)
                else:
                    self.gene_activators[tnode].append(snode)
                n_gene += 1

        for d in [self.gene_activators, self.gene_inhibitors,
                  self.pathway_activators, self.pathway_inhibitors]:
            for k in d:
                d[k] = list(set(d[k]))

        print(f"  Genes: {len(self.genes)} | Pathways: {len(self.pathways)}")
        print(f"  Gene interactions: {n_gene} | Pathway interactions: {n_path}")

    def _parse_v3_individual(self):
        print("  [WARN] v3 individual mode not implemented.")

    # ------------------------------------------------------------------
    # ODE equation generation (v9: inhibition implemented)
    # ------------------------------------------------------------------
    def _act_eq(self, node, regulators, inhibitors=None):
        ksp  = self.default_params['ksp']
        n    = self.default_params['n']
        kcat = self.default_params['kcat']
        Vmax = self.default_params['Vmax']
        sub  = f"{node}_inactive**{n}"
        den  = f"({ksp}**{n} + {node}_inactive**{n})"

        # Structural inhibition edges: Hill-function-based inhibition factor
        # Per Tseng (2024) N&M 21:65, Eq.(1)(2):
        #   enzyme-to-enzyme interactions (activation or inhibition) use Hill kinetics
        #   inh_factor = ksp^n / (ksp^n + inh_active^n)
        #   inh=0 → factor=1 (no inhibition); inh=ksp → factor=0.5 (half-inhibition)
        # NOTE: drug interventions use ki multiplier (Eq.4/5), handled separately
        #       by masld_generate_intervention_models_v4.py
        inh_factor = ""
        if inhibitors:
            inh_terms = [f"({ksp}**{n} / ({ksp}**{n} + {inh}_active**{n}))"
                         for inh in sorted(inhibitors)]
            inh_factor = " * " + " * ".join(inh_terms)

        if not regulators:
            return f"(({Vmax} * {sub}) / {den}){inh_factor}"
        terms = [f"({kcat} * {r}_active * {sub}) / {den}" for r in regulators]
        base  = " + ".join(f"({t})" for t in terms)
        if inh_factor:
            return f"({base}){inh_factor}"
        return base

    def _inact_eq(self, node):
        ksp  = self.default_params['ksp']
        n    = self.default_params['n']
        Vmax = self.default_params['Vmax']
        return (f"({Vmax} * {node}_active**{n}) / "
                f"({ksp}**{n} + {node}_active**{n})")

    def generate_ode_system(self):
        print("\nGenerating ODE system...")
        eqs = {}
        for node in sorted(self.genes):
            acts = self.gene_activators.get(node, [])
            inhs = self.gene_inhibitors.get(node, [])
            ae   = self._act_eq(node, acts, inhs)
            ie   = self._inact_eq(node)
            eqs[f"{node}_inactive"] = f"-({ae}) + ({ie})"
            eqs[f"{node}_active"]   = f"({ae}) - ({ie})"
        for node in sorted(self.pathways):
            acts = self.pathway_activators.get(node, [])
            inhs = self.pathway_inhibitors.get(node, [])
            ae   = self._act_eq(node, acts, inhs)
            ie   = self._inact_eq(node)
            eqs[f"{node}_inactive"] = f"-({ae}) + ({ie})"
            eqs[f"{node}_active"]   = f"({ae}) - ({ie})"
        return eqs

    # ------------------------------------------------------------------
    # Generate Python ODE model file
    # ------------------------------------------------------------------
    def generate_python_model(self, output_file=None):
        if output_file is None:
            output_file = (f'ode_model_{self.iv_type}_'
                           f'{self.ratio_condition}_v9.py')
        print(f"\nGenerating Python ODE model -> {output_file}")

        eqs       = self.generate_ode_system()
        cfg_label = IV_TYPE_CONFIG[self.iv_type]['label']

        # --- STATE_VARS ---
        sv_lines = []
        for gene in sorted(self.genes):
            meta  = self.gene_metadata.get(gene, {})
            sym   = meta.get('symbol', gene)
            kdisp = meta.get('kegg_display_name', gene)
            label = f"{kdisp} ({sym})" if kdisp != gene and kdisp != sym else sym
            sv_lines.append(f"    '{gene}_inactive',  # {label}")
            sv_lines.append(f"    '{gene}_active',    # {label}")
        for pw in sorted(self.pathways):
            pname = self.pathway_metadata.get(pw, {}).get('pathway_name', pw)
            sv_lines.append(f"    '{pw}_inactive',  # Pathway: {pname}")
            sv_lines.append(f"    '{pw}_active',    # Pathway: {pname}")

        # --- ODE body ---
        unpack   = []
        eq_lines = ["    dydt = np.zeros_like(y)", ""]
        idx = 0
        for gene in sorted(self.genes):
            meta  = self.gene_metadata.get(gene, {})
            sym   = meta.get('symbol', gene)
            kdisp = meta.get('kegg_display_name', gene)
            label = f"{kdisp} (HGNC: {sym})" if kdisp != gene and kdisp != sym else sym
            unpack   += [f"    # {label}",
                         f"    {gene}_inactive = y[{idx}]",
                         f"    {gene}_active   = y[{idx+1}]", ""]
            eq_lines += [f"    # {label}",
                         f"    dydt[{idx}]   = {eqs[gene+'_inactive']}",
                         f"    dydt[{idx+1}] = {eqs[gene+'_active']}", ""]
            idx += 2
        for pw in sorted(self.pathways):
            pname = self.pathway_metadata.get(pw, {}).get('pathway_name', pw)
            unpack   += [f"    # Pathway: {pname}",
                         f"    {pw}_inactive = y[{idx}]",
                         f"    {pw}_active   = y[{idx+1}]", ""]
            eq_lines += [f"    # Pathway: {pname}",
                         f"    dydt[{idx}]   = {eqs[pw+'_inactive']}",
                         f"    dydt[{idx+1}] = {eqs[pw+'_active']}", ""]
            idx += 2
        eq_lines.append("    return dydt")

        # --- Metadata strings ---
        gene_meta_str    = "\n".join(
            f"    '{g}': {self.gene_metadata.get(g, {})},  "
            f"# KEGG: {self.gene_metadata.get(g,{}).get('kegg_display_name',g)}"
            f" | ratio={self.gene_metadata.get(g,{}).get('initial_ratio',1.0):.4f}"
            for g in sorted(self.genes))
        pathway_meta_str = "\n".join(
            f"    '{p}': {self.pathway_metadata.get(p, {})},  "
            f"# initial_inactive={PATHWAY_INITIAL_INACTIVE}"
            for p in sorted(self.pathways))

        _intf = self.interactions_file.replace('\\', '/')
        _nv1  = (self.nodes_file_v1 or 'N/A').replace('\\', '/')
        _nv3  = (getattr(self, 'nodes_file_v3', None) or 'N/A').replace('\\', '/')

        code = f'''#!/usr/bin/env python3
"""
ODE Model v8 | iv_type={self.iv_type} | ratio_condition={self.ratio_condition}
{'='*72}
Label           : {cfg_label}
Ratio condition : {self.ratio_condition}
  Numerator     : {self.num_condition}
  Denominator   : {RATIO_TO_DENOM} (Normal, implicit = 1.0)
  Gene IV       : {self.num_condition}_expr / Normal_expr  (ratio-based)
  Pathway IV    : {PATHWAY_INITIAL_INACTIVE} (fixed, per prior papers)
Node naming     : KEGG hsa04932 display names (v8)
Format          : {self.format_version}
Generated       : {pd.Timestamp.now()}
Interactions    : {_intf}
Nodes v3        : {_nv3}
Genes           : {len(self.genes)}
Pathways        : {len(self.pathways)}
StateVars       : {(len(self.genes) + len(self.pathways)) * 2}
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd

MODEL_INFO = {{
    'iv_type':         '{self.iv_type}',
    'iv_label':        '{cfg_label}',
    'ratio_condition': '{self.ratio_condition}',
    'num_condition':   '{self.num_condition}',
    'denom_condition': '{RATIO_TO_DENOM}',
    'version':         'v8',
    'n_genes':         {len(self.genes)},
    'n_pathways':      {len(self.pathways)},
    'pathway_initial_inactive': {PATHWAY_INITIAL_INACTIVE},
    'gene_iv_method':  'ratio = condition_expr / Normal_expr',
}}

PARAMS = {{
    'Vmax': {self.default_params['Vmax']},
    'ksp' : {self.default_params['ksp']},
    'n'   : {self.default_params['n']},
    'kcat': {self.default_params['kcat']},
}}

STATE_VARS = [
{chr(10).join(sv_lines)}
]

GENE_METADATA = {{
{gene_meta_str}
}}

PATHWAY_METADATA = {{
{pathway_meta_str}
}}

# =============================================================================
# Initial conditions (v7 design):
#   Gene nodes   : inactive = ratio (condition/Normal), active = 0
#   Pathway nodes: inactive = {PATHWAY_INITIAL_INACTIVE} (fixed), active = 0
# =============================================================================
Y0 = np.zeros(len(STATE_VARS))
for _i, _var in enumerate(STATE_VARS):
    if '_inactive' in _var:
        _node = _var.replace('_inactive', '')
        if _node in PATHWAY_METADATA:
            Y0[_i] = PATHWAY_METADATA[_node].get('initial_level',
                                                  {PATHWAY_INITIAL_INACTIVE})
        elif _node in GENE_METADATA:
            Y0[_i] = GENE_METADATA[_node].get('initial_ratio', 1.0)
        else:
            Y0[_i] = 1.0


def ode_system(y, t, params=None):
    """ODE system equations."""
    if params is None:
        params = PARAMS
    Vmax = params['Vmax']
    ksp  = params['ksp']
    n    = params['n']
    kcat = params['kcat']

{chr(10).join(unpack)}
{chr(10).join(eq_lines)}


def simulate(t_max=200, n_points=2001, y0=None, params=None):
    """Run ODE simulation. Returns (t, y)."""
    if y0 is None:     y0     = Y0
    if params is None: params = PARAMS
    t = np.linspace(0, t_max, n_points)
    y = odeint(ode_system, y0, t, args=(params,), mxstep=10000)
    return t, y


def get_pathway_steady_states(y, last_frac=0.1):
    """Return {{var_name: steady_state_value}} for all active pathway variables."""
    n = max(1, int(len(y) * last_frac))
    return {{STATE_VARS[i]: float(np.mean(y[-n:, i]))
             for i in range(len(STATE_VARS))
             if '_active' in STATE_VARS[i] and
                STATE_VARS[i].replace('_active', '') in PATHWAY_METADATA}}


def get_gene_steady_states(y, last_frac=0.1):
    """Return {{gene_symbol: (inactive_ss, active_ss)}}."""
    n = max(1, int(len(y) * last_frac))
    result = {{}}
    for g, meta in GENE_METADATA.items():
        sym = meta.get('symbol', g)
        try:
            i_in  = STATE_VARS.index(g + '_inactive')
            i_act = STATE_VARS.index(g + '_active')
            result[sym] = (float(np.mean(y[-n:, i_in])),
                           float(np.mean(y[-n:, i_act])))
        except ValueError:
            pass
    return result


def get_pathway_timeseries(y):
    """Return {{pathway_name: active_curve}} for all pathway output nodes."""
    result = {{}}
    for pw, meta in PATHWAY_METADATA.items():
        pname = meta.get('pathway_name', pw)
        try:
            i = STATE_VARS.index(pw + '_active')
            result[pname] = y[:, i]
        except ValueError:
            pass
    return result


if __name__ == "__main__":
    print("Model Info:", MODEL_INFO)
    print(f"Genes: {{len(GENE_METADATA)}} | Pathways: {{len(PATHWAY_METADATA)}}")
    t, y = simulate()
    ss = get_pathway_steady_states(y)
    print("Pathway steady states (active):")
    for k, v in sorted(ss.items()):
        print(f"  {{k}}: {{v:.4f}}")
'''
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(code)
        print(f"  Saved: {output_file}")
        return output_file

    # ------------------------------------------------------------------
    # Generate SBML (preserved from v6, with ratio initial values)
    # ------------------------------------------------------------------
    def generate_sbml_model(self, output_file=None):
        if output_file is None:
            output_file = (f'ode_model_{self.iv_type}_'
                           f'{self.ratio_condition}_v8.xml')
        print(f"\nGenerating SBML L2V4 -> {output_file}")
        try:
            import libsbml
        except ImportError:
            print("  [ERROR] libsbml not installed.")
            return None

        doc   = libsbml.SBMLDocument(2, 4)
        model = doc.createModel()
        model.setId(f"ODE_V9_{self.iv_type}_{self.ratio_condition}")
        model.setName(f"ODE v9 {self.iv_type} {self.ratio_condition}")

        comp = model.createCompartment()
        comp.setId("cytosol"); comp.setSize(1.0); comp.setConstant(True)

        for pname, pval in self.default_params.items():
            p = model.createParameter()
            p.setId(pname); p.setValue(float(pval)); p.setConstant(True)

        for gene in sorted(self.genes):
            meta  = self.gene_metadata.get(gene, {})
            ratio = float(meta.get('initial_ratio', 1.0))
            sym   = meta.get('symbol', gene)
            for suffix, conc in [('_inactive', ratio), ('_active', 0.0)]:
                sp = model.createSpecies()
                sp.setId(f"{gene}{suffix}"); sp.setName(f"{sym}{suffix}")
                sp.setCompartment("cytosol"); sp.setInitialConcentration(conc)
                sp.setConstant(False); sp.setBoundaryCondition(False)

        for pw in sorted(self.pathways):
            pname = self.pathway_metadata.get(pw, {}).get('pathway_name', pw)
            for suffix, conc in [('_inactive', PATHWAY_INITIAL_INACTIVE),
                                  ('_active', 0.0)]:
                sp = model.createSpecies()
                sp.setId(f"{pw}{suffix}"); sp.setName(f"{pname[:30]}{suffix}")
                sp.setCompartment("cytosol"); sp.setInitialConcentration(conc)
                sp.setConstant(False); sp.setBoundaryCondition(False)

        rid = 0
        for node in sorted(self.genes) + sorted(self.pathways):
            is_gene = node in self.genes
            acts = (self.gene_activators if is_gene
                    else self.pathway_activators).get(node, [])
            rid += 1
            r = model.createReaction()
            r.setId(f"R{rid}_act_{node}"); r.setReversible(False)
            r.createReactant().setSpecies(f"{node}_inactive")
            r.createProduct().setSpecies(f"{node}_active")
            for a in acts:
                r.createModifier().setSpecies(f"{a}_active")
            inhs_sbml = (self.gene_inhibitors if is_gene
                         else self.pathway_inhibitors).get(node, [])
            for inh in inhs_sbml:
                r.createModifier().setSpecies(f"{inh}_active")
            kl = r.createKineticLaw()
            inh_sbml = ""
            if inhs_sbml:
                inh_parts = [f"(ksp^n / (ksp^n + {inh}_active^n))"
                             for inh in sorted(inhs_sbml)]
                inh_sbml = " * " + " * ".join(inh_parts)
            if not acts:
                kl.setFormula(
                    f"((Vmax * {node}_inactive^n) / (ksp^n + {node}_inactive^n)){inh_sbml}")
            else:
                terms = [f"(kcat * {a}_active * {node}_inactive^n) / "
                         f"(ksp^n + {node}_inactive^n)" for a in acts]
                formula = " + ".join(f"({t})" for t in terms)
                if inh_sbml:
                    formula = f"({formula}){inh_sbml}"
                kl.setFormula(formula)

            rid += 1
            r2 = model.createReaction()
            r2.setId(f"R{rid}_inact_{node}"); r2.setReversible(False)
            r2.createReactant().setSpecies(f"{node}_active")
            r2.createProduct().setSpecies(f"{node}_inactive")
            r2.createKineticLaw().setFormula(
                f"(Vmax * {node}_active^n) / (ksp^n + {node}_active^n)")

        libsbml.writeSBMLToFile(doc, output_file)
        print(f"  Saved: {output_file} | Species: {model.getNumSpecies()} "
              f"| Reactions: {model.getNumReactions()}")
        return output_file


# ---------------------------------------------------------------------------
# Backward-compatibility aliases
# ---------------------------------------------------------------------------
InteractionToODEConverter           = InteractionToODEConverter_V7
InteractionToODEConverter_V6        = InteractionToODEConverter_V7
InteractionToODEConverter_KEGGID_V5 = InteractionToODEConverter_V7


if __name__ == "__main__":
    print("Library module. Use run_complete_workflow_v7.py (or v8) to generate models.")
    print(f"\nSupported IV types       : {ALL_IV_TYPES}")
    print(f"Supported ratio conditions: {ALL_RATIO_CONDITIONS}")
    print(f"Total model files        : "
          f"{len(ALL_IV_TYPES)} x {len(ALL_RATIO_CONDITIONS)} = "
          f"{len(ALL_IV_TYPES)*len(ALL_RATIO_CONDITIONS)}")
