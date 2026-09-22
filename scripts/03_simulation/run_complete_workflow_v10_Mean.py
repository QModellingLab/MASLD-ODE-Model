#!/usr/bin/env python3
"""
Complete Workflow - Version 10 (pyDESeq2 + multi-gene averaging)
================================================================
Default: pyDESeq2 ratio × 3 conditions = 3 ODE model files.
v10: Uses node_name_table for multi-gene averaged initial values
     (consistent with N&M 2024 methodology).
Legacy IV types (raw/normalized/log2_cpm/batch_corr) still available via --iv-types.

Input files (pydeseq2 mode):
  - result_v4/result_v4_symbol_interactions.xlsx  (network topology)
  - DEG_pydeseq2_3comparisons.xlsx                (pyDESeq2 ratios)
  - node_name_table_hsa04932.xlsx                  (node→gene mapping, v10)

Output:
  - models_pydeseq2_mean/ode_model_pydeseq2_{condition}_v10_mean.py  (3 files)
  - Node_Initial_Levels_pydeseq2_mean.xlsx                (summary record)

Usage:
  python run_complete_workflow_v10.py                         # interactive menu
  python run_complete_workflow_v10.py --step 3               # generate 3 models
  python run_complete_workflow_v10.py --step 3 --ratio-conds NASH_vs_Normal
"""

import os
import sys

if sys.platform == 'win32':
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    if hasattr(sys.stderr, 'reconfigure'):
        sys.stderr.reconfigure(encoding='utf-8')
    try:
        import ctypes
        ctypes.windll.kernel32.SetConsoleOutputCP(65001)
    except Exception:
        pass

# =============================================================================
# Configuration Parameters  <- EDIT THESE
# =============================================================================

KGML1          = 'hsa04932.xml'
KGML2          = 'panu04932.xml'
OUTPUT_PREFIX  = 'result_v4'
WORK_DIR       = '.'
NODES_FILE_V1  = 'Node_Initial_Levels_4grp.xlsx'          # legacy v1 (fallback)
NODES_FILE_V3  = 'Node_Initial_Levels_4grp_v3.xlsx'          # v3 unified (all 5 IV types)
DEG_FILE       = 'DEG_pydeseq2_3comparisons.xlsx'            # pyDESeq2 ratios
NNT_FILE       = 'node_name_table_hsa04932.xlsx'             # v10: node→gene mapping
AVERAGING_METHOD = 'mean'  # 'mean' or 'median' for multi-gene nodes

IV_TYPES       = ['pydeseq2']
RATIO_CONDITIONS = ['Obese_vs_Normal', 'NAFL_vs_Normal', 'NASH_vs_Normal']
MODELS_SUBDIR  = 'models_pydeseq2_mean'
GENERATE_SBML  = False

# =============================================================================
# Labels
# =============================================================================

IV_LABELS = {
    'pydeseq2':   'pyDESeq2 ratio (2^log2FC) [UNIFIED]',
    # Legacy (still usable via --iv-types flag):
    'raw':        'Raw Counts',
    'normalized': 'Size-Factor Normalized',
    'log2_norm':  'Log2-Norm -> Linear (2^x-1)',
    'log2_cpm':   'Log2-CPM -> Linear (2^x-0.5), no BC',
    'batch_corr': 'Batch-Corrected Linear',
}

RATIO_LABELS = {
    'Obese_vs_Normal': 'Obese / Normal (early obesity)',
    'NAFL_vs_Normal':  'NAFL  / Normal (steatosis)',
    'NASH_vs_Normal':  'NASH  / Normal (steatohepatitis)',
}

IV_NEEDS_V3 = {'raw', 'normalized', 'log2_norm', 'log2_cpm', 'batch_corr'}  # need v3 for node expression
IV_NEEDS_V1 = {'raw', 'normalized', 'log2_norm'}   # fallback if v3 not present
# pydeseq2 does NOT need v3/v1 — uses built-in KEGG mapping + DEG Excel

# =============================================================================
# Step functions
# =============================================================================

def step1_run_mapper():
    print("=" * 80)
    print("Step 1: Run KEGG Mapper v4")
    print("=" * 80)
    try:
        from kegg_mapper_union_intersection_v4 import run_workflow
        result = run_workflow(
            species1_xml=KGML1, species2_xml=KGML2,
            output_prefix=OUTPUT_PREFIX, node_mode='symbol', show_summary=True)
        if result == 0:
            print(f"\n[OK] Generated: {OUTPUT_PREFIX}_symbol_interactions.xlsx")
            return (f"{OUTPUT_PREFIX}_symbol_interactions.xlsx",
                    f"{OUTPUT_PREFIX}_symbol_nodes.xlsx")
        else:
            print("[ERROR] Mapper failed"); return None, None
    except ImportError as e:
        print(f"[ERROR] kegg_mapper module not found: {e}"); return None, None


def step2_wait_for_editing():
    print("=" * 80)
    print("Step 2: Edit Node Initial Levels (Optional)")
    print("=" * 80)
    print(f"\n  Open {NODES_FILE_V3} in Excel if needed.")
    input("\nPress Enter to continue...")


def step3_generate_models(iv_types=None, ratio_conditions=None,
                           work_dir=None, models_subdir=None,
                           generate_sbml=None):
    """
    Generate ODE model files v10 (multi-gene averaged IVs).

    Output: {models_subdir}/ode_model_{iv_type}_{ratio_condition}_v10.py
    """
    if iv_types        is None: iv_types        = IV_TYPES
    if ratio_conditions is None: ratio_conditions = RATIO_CONDITIONS
    if work_dir        is None: work_dir        = WORK_DIR
    if models_subdir   is None: models_subdir   = MODELS_SUBDIR
    if generate_sbml   is None: generate_sbml   = GENERATE_SBML

    print("=" * 80)
    print("Step 3: Generate ODE Models v10 (multi-gene averaged IVs + node_name_table)")
    print(f"  IV types         : {iv_types}")
    print(f"  Ratio conditions : {ratio_conditions}")
    print(f"  SBML             : {'Yes' if generate_sbml else 'No'}")
    print("=" * 80)

    interactions_file = os.path.join(work_dir,
                                      f'{OUTPUT_PREFIX}_symbol_interactions.xlsx')
    nodes_v1 = os.path.join(work_dir, NODES_FILE_V1)
    nodes_v3 = os.path.join(work_dir, NODES_FILE_V3)

    if not os.path.exists(interactions_file):
        print(f"[ERROR] Not found: {interactions_file}")
        return []

    need_v3 = any(iv in IV_NEEDS_V3 for iv in iv_types)
    if need_v3 and not os.path.exists(nodes_v3):
        if not os.path.exists(nodes_v1):
            print(f"[ERROR] Neither v3 nor v1 nodes file found: {nodes_v3}")
            return []
        else:
            print(f"[WARN] {NODES_FILE_V3} not found, will use v1 fallback for applicable IV types")

    # DEG file for pydeseq2 IV type
    deg_path = os.path.join(work_dir, DEG_FILE)
    if 'pydeseq2' in iv_types and not os.path.exists(deg_path):
        print(f"[ERROR] DEG file not found for pydeseq2 IV type: {deg_path}")
        return []

    # v10: node_name_table for multi-gene averaging
    nnt_path = os.path.join(work_dir, NNT_FILE)
    if 'pydeseq2' in iv_types and not os.path.exists(nnt_path):
        print(f"[WARN] node_name_table not found: {nnt_path}")
        print(f"       Will use v9 single-gene fallback.")
        nnt_path = None

    models_dir = os.path.join(work_dir, models_subdir)
    os.makedirs(models_dir, exist_ok=True)
    print(f"\nOutput directory: {models_dir}\n")

    try:
        from interaction_to_ode_keggid_v10_Mean import InteractionToODEConverter_V7
    except ImportError:
        print("[ERROR] interaction_to_ode_keggid_v10_Mean.py not found")
        return []

    generated = []
    failed    = []
    total     = len(iv_types) * len(ratio_conditions)
    count     = 0

    for iv_type in iv_types:
        for rc in ratio_conditions:
            count += 1
            print(f"\n[{count}/{total}] {iv_type} × {rc}")
            print(f"  IV    : {IV_LABELS.get(iv_type, iv_type)}")
            print(f"  Ratio : {RATIO_LABELS.get(rc, rc)}")
            print("-" * 60)

            try:
                conv = InteractionToODEConverter_V7(
                    interactions_file=interactions_file,
                    nodes_file_v3=nodes_v3 if os.path.exists(nodes_v3) else None,
                    nodes_file_v1=nodes_v1 if os.path.exists(nodes_v1) else None,
                    deg_file=deg_path if os.path.exists(deg_path) else None,
                    node_name_table=nnt_path,
                    ratio_condition=rc,
                    iv_type=iv_type,
                    averaging_method=AVERAGING_METHOD,
                )
                conv.load_interactions()
                conv.load_nodes()
                conv.parse_interactions()

                py_out = os.path.join(models_dir,
                                       f'ode_model_{iv_type}_{rc}_v10_mean.py')
                conv.generate_python_model(py_out)
                generated.append(('py', iv_type, rc, py_out))

                if generate_sbml:
                    xml_out = os.path.join(models_dir,
                                            f'ode_model_{iv_type}_{rc}_v10_mean.xml')
                    result = conv.generate_sbml_model(xml_out)
                    if result:
                        generated.append(('xml', iv_type, rc, xml_out))

                print(f"  [OK] {iv_type} × {rc}")

            except Exception as e:
                import traceback
                print(f"  [ERROR] {iv_type} × {rc}: {e}")
                traceback.print_exc()
                failed.append((iv_type, rc, str(e)))

    # Summary
    print("\n" + "=" * 80)
    print("Step 3 Completed  (v10 — multi-gene averaged IVs + node_name_table)")
    print("=" * 80)

    py_files  = [r for r in generated if r[0] == 'py']
    xml_files = [r for r in generated if r[0] == 'xml']
    print(f"\nGenerated: {len(py_files)} Python models"
          + (f", {len(xml_files)} SBML models" if xml_files else ""))

    if failed:
        print(f"\n[WARN] Failed: {len(failed)}")
        for iv, rc, err in failed:
            print(f"  {iv} × {rc}: {err[:80]}")

    # Model matrix
    all_rcs  = ['Obese_vs_Normal', 'NAFL_vs_Normal', 'NASH_vs_Normal']
    header   = f"{'IV Type':>12}  " + "  ".join(f"{c:17}" for c in all_rcs)
    print(f"\nModel matrix (OK/ERR/-):\n{header}")
    for iv in iv_types:
        row = f"{iv:>12}  "
        for rc in all_rcs:
            if rc not in ratio_conditions:
                row += f"{'     -     ':19}"
            elif any(f[1] == iv and f[2] == rc for f in failed):
                row += f"{'   ERROR   ':19}"
            elif any(r[1] == iv and r[2] == rc for r in py_files):
                row += f"{'    OK     ':19}"
            else:
                row += f"{'     -     ':19}"
        print(row)

    if py_files:
        print(f"\nNext step: evaluate model biological plausibility:")
        print(f"  python evaluate_bio_plausibility_v10.py --model-dir {models_dir}")

    return generated


# =============================================================================
# Settings / Interactive menu
# =============================================================================

def show_settings():
    print("\n[Current Settings]")
    print(f"  Work dir         : {WORK_DIR}")
    print(f"  DEG file         : {DEG_FILE}")
    print(f"  Node name table  : {NNT_FILE}")
    print(f"  Averaging method : {AVERAGING_METHOD}")
    print(f"  IV types         : {IV_TYPES}")
    print(f"  Ratio conditions : {RATIO_CONDITIONS}")
    print(f"  Models subdir    : {MODELS_SUBDIR}")
    print(f"  KGML1 (Step 1)   : {KGML1}")
    print(f"  KGML2 (Step 1)   : {KGML2}")
    print(f"  Output prefix    : {OUTPUT_PREFIX}")


def interactive_menu():
    global GENERATE_SBML, IV_TYPES, RATIO_CONDITIONS, OUTPUT_PREFIX
    global KGML1, KGML2, MODELS_SUBDIR

    print("\n" + "=" * 80)
    print("Complete Workflow v10: pyDESeq2 ratio × 3 conditions ODE Models")
    print("  Gene IV  = pyDESeq2 2^(log2FC) ratio (Normal = 1.0)")
    print("  Path IV  = 100 (fixed, per prior BIBE 2022 / N&M 2024 papers)")
    print("=" * 80)

    while True:
        n_models = len(IV_TYPES) * len(RATIO_CONDITIONS)
        print("\n" + "-" * 50)
        print("  1. Run Step 1 (KEGG Mapper)")
        print(f"  2. Generate {n_models} models (pydeseq2 × 3 conditions)")
        print("  3. Show current settings")
        print("  4. Exit")
        print(f"\n  Current: IVs={IV_TYPES} | Conditions={RATIO_CONDITIONS}")

        choice = input("\nEnter option (1-4): ").strip()

        if choice == '1':
            step1_run_mapper()
        elif choice == '2':
            step3_generate_models()
        elif choice == '3':
            show_settings()
        elif choice == '4':
            print("Goodbye!"); break
        else:
            print("[ERROR] Invalid option")


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(
        description='Complete Workflow v10: pyDESeq2 ratio × 3 conditions ODE models',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive menu
  python run_complete_workflow_v10.py

  # Generate 3 pydeseq2 models (default)
  python run_complete_workflow_v10.py --step 3

  # Generate only NASH_vs_Normal model
  python run_complete_workflow_v10.py --step 3 --ratio-conds NASH_vs_Normal

  # Legacy: use other IV types (requires Node_Initial_Levels_4grp_v3.xlsx)
  python run_complete_workflow_v10.py --step 3 --iv-types batch_corr log2_cpm
        """)
    parser.add_argument('--step', choices=['1', '3', 'all', 'interactive'],
                        default='interactive')
    parser.add_argument('--kgml1',         default=None)
    parser.add_argument('--kgml2',         default=None)
    parser.add_argument('--output-prefix', default=None)
    parser.add_argument('--work-dir',      default=None)
    parser.add_argument('--nodes-v1',      default=None)
    parser.add_argument('--nodes-v3',      default=None,
                        help='Path to unified v3 nodes file (legacy IV types only)')
    parser.add_argument('--iv-types', nargs='+',
                        choices=['raw','normalized','log2_norm',
                                 'log2_cpm','batch_corr','pydeseq2'],
                        default=None)
    parser.add_argument('--ratio-conds', nargs='+',
                        choices=['Obese_vs_Normal','NAFL_vs_Normal','NASH_vs_Normal'],
                        default=None,
                        dest='ratio_conds')
    parser.add_argument('--models-subdir', default=None)
    parser.add_argument('--no-sbml', action='store_true')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if args.kgml1:          KGML1          = args.kgml1
    if args.kgml2:          KGML2          = args.kgml2
    if args.output_prefix:  OUTPUT_PREFIX  = args.output_prefix
    if args.work_dir:       WORK_DIR       = args.work_dir
    if args.nodes_v1:       NODES_FILE_V1  = args.nodes_v1
    if args.nodes_v3:       NODES_FILE_V3  = args.nodes_v3
    if args.models_subdir:  MODELS_SUBDIR  = args.models_subdir
    if args.no_sbml:        GENERATE_SBML  = False

    iv_types    = args.iv_types   if args.iv_types   else IV_TYPES
    ratio_conds = args.ratio_conds if args.ratio_conds else RATIO_CONDITIONS

    try:
        if args.step == 'interactive':
            interactive_menu()
        elif args.step == '1':
            step1_run_mapper()
        elif args.step == '3':
            step3_generate_models(iv_types=iv_types,
                                   ratio_conditions=ratio_conds,
                                   work_dir=WORK_DIR,
                                   models_subdir=MODELS_SUBDIR,
                                   generate_sbml=GENERATE_SBML)
        elif args.step == 'all':
            intf, nf = step1_run_mapper()
            if intf:
                step2_wait_for_editing()
                step3_generate_models(iv_types=iv_types,
                                       ratio_conditions=ratio_conds,
                                       work_dir=WORK_DIR,
                                       models_subdir=MODELS_SUBDIR,
                                       generate_sbml=GENERATE_SBML)
    except KeyboardInterrupt:
        print("\n\nCancelled by user.")
    except Exception as e:
        import traceback
        print(f"\n[ERROR] Fatal error: {e}")
        traceback.print_exc()
