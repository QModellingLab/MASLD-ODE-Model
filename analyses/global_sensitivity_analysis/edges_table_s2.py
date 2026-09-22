#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary Table S2 (all 106 regulatory interactions of the hsa04932 ODE
network), transcribed with Source/Target converted from gene-symbol naming
(as used in Table S2 itself) to the ODE node display names used in
GENE_METADATA / STATE_VARS (as used in Table S1 and the model source files).
P_* pathway-output targets are already ODE node names and are left as-is.

Used by step2_topology_comparison.py to build the directed interaction graph
independently of the Morris global-sensitivity results.
"""

# gene-symbol (as it appears in Table S2) -> ODE node display name (Table S1)
SYMBOL_TO_NODE = {
    'ADIPOQ': 'ACDC', 'ADIPOR1': 'adipoR', 'AKT3': 'Akt', 'ATF4': 'ATF4',
    'BAX': 'Bax', 'BCL2L11': 'Bim', 'BID': 'Bid', 'CASP3': 'CASP3',
    'CASP7': 'CASP7', 'CASP8': 'CASP8', 'CDC42': 'Cdc42', 'CEBPA': 'C_EBPa',
    'COX6B2': 'CxIV', 'CXCL8': 'IL_8', 'CYCS': 'Cytc', 'CYP2E1': 'CYP2E1',
    'DDIT3': 'CHOP', 'EIF2AK3': 'PERK', 'EIF2S1': 'eIF2a', 'ERN1': 'IRE1a',
    'FAS': 'Fas', 'FASLG': 'FasL', 'FOS': 'AP_1', 'GSK3A': 'GSK_3',
    'IKBKB': 'IKKb', 'IL1A': 'IL_1', 'IL6': 'IL_6', 'IL6R': 'IL_6R',
    'INS': 'INS', 'INSR': 'INSR', 'IRS1': 'IRS_1_2', 'ITCH': 'ITCH',
    'LEP': 'LEP', 'LEPR': 'ObR', 'MAP3K11': 'MLK3', 'MAP3K5': 'ASK1',
    'MAPK14': 'p38', 'MAPK8': 'JNK1_2', 'MLXIP': 'ChREBP',
    'NDUFC2-KCTD14': 'CxI', 'NFKB1': 'NF_kB', 'NR1H3': 'LXR_a',
    'P3R3URF-PIK3R3': 'PI3K', 'PKLR': 'L_PK', 'PPARA': 'PPAR_a',
    'PPARG': 'PPAR_g', 'PRKAG2': 'AMPK', 'RAC1': 'Rac1', 'RXRA': 'RXR',
    'SDHA': 'CxII', 'SOCS3': 'SOCS3', 'SREBF1': 'SREBP_1c',
    'TGFB1': 'TGF_b1', 'TNF': 'TNFa', 'TNFRSF1A': 'TNFR1', 'TRAF2': 'TRAF2',
    'UQCR11': 'CxIII', 'XBP1': 'XBP1',
}

# (source_symbol, target_symbol_or_P*, relation) -- relation kept for
# reference but not used by the topology comparison (which treats all
# edges as unweighted directed edges regardless of activation/inhibition).
_RAW_EDGES = [
    ('IL6', 'IL6R'), ('TNF', 'TNFRSF1A'), ('INS', 'INSR'), ('LEP', 'LEPR'),
    ('ADIPOQ', 'ADIPOR1'), ('TNF', 'MAPK8'), ('FASLG', 'FAS'),
    ('IL6R', 'SOCS3'), ('TNFRSF1A', 'NFKB1'), ('NFKB1', 'IL6'),
    ('INSR', 'IRS1'), ('IRS1', 'P3R3URF-PIK3R3'), ('P3R3URF-PIK3R3', 'AKT3'),
    ('AKT3', 'GSK3A'), ('GSK3A', 'INS'), ('INSR', 'NR1H3'),
    ('NR1H3', 'RXRA'), ('RXRA', 'SREBF1'), ('MLXIP', 'PKLR'),
    ('LEPR', 'PRKAG2'), ('ADIPOR1', 'PRKAG2'), ('PRKAG2', 'MAPK14'),
    ('MAPK14', 'PPARA'), ('CDC42', 'MAP3K11'), ('RAC1', 'MAP3K11'),
    ('MAP3K11', 'MAPK8'), ('MAPK8', 'ITCH'), ('ITCH', 'CASP8'),
    ('FAS', 'CASP8'), ('CASP8', 'BID'), ('BID', 'BAX'), ('CYCS', 'CASP3'),
    ('CYCS', 'CASP7'), ('ERN1', 'XBP1'), ('ERN1', 'TRAF2'),
    ('TRAF2', 'MAP3K5'), ('MAP3K5', 'MAPK8'), ('TRAF2', 'IKBKB'),
    ('MAPK8', 'IRS1'), ('MAPK8', 'FOS'), ('IKBKB', 'NFKB1'),
    ('FOS', 'IL1A'), ('FOS', 'IL6'), ('FOS', 'TNF'), ('NFKB1', 'IL1A'),
    ('NFKB1', 'TNF'), ('EIF2AK3', 'EIF2S1'), ('CYP2E1', 'MAPK8'),
    ('CYP2E1', 'IKBKB'), ('CYP2E1', 'FASLG'), ('CYP2E1', 'TNF'),
    ('CYP2E1', 'CXCL8'), ('CYP2E1', 'TGFB1'), ('EIF2S1', 'ATF4'),
    ('ATF4', 'DDIT3'), ('DDIT3', 'BCL2L11'), ('BCL2L11', 'BAX'),
    ('BAX', 'CYCS'), ('NDUFC2-KCTD14', 'MAPK8'), ('SDHA', 'MAPK8'),
    ('UQCR11', 'MAPK8'), ('COX6B2', 'MAPK8'), ('NDUFC2-KCTD14', 'IKBKB'),
    ('SDHA', 'IKBKB'), ('UQCR11', 'IKBKB'), ('COX6B2', 'IKBKB'),
    ('NDUFC2-KCTD14', 'FASLG'), ('SDHA', 'FASLG'), ('UQCR11', 'FASLG'),
    ('COX6B2', 'FASLG'), ('NDUFC2-KCTD14', 'TNF'), ('SDHA', 'TNF'),
    ('UQCR11', 'TNF'), ('COX6B2', 'TNF'), ('NDUFC2-KCTD14', 'CXCL8'),
    ('SDHA', 'CXCL8'), ('UQCR11', 'CXCL8'), ('COX6B2', 'CXCL8'),
    ('NDUFC2-KCTD14', 'TGFB1'), ('SDHA', 'TGFB1'), ('UQCR11', 'TGFB1'),
    ('COX6B2', 'TGFB1'), ('SOCS3', 'SREBF1'), ('XBP1', 'CEBPA'),
    ('RXRA', 'PPARG'), ('SOCS3', 'IRS1'),
    ('GSK3A', 'P_Hyperinsulinemia'),
    ('SREBF1', 'P_De_novo_fatty_acid_synthesis'),
    ('MLXIP', 'P_De_novo_fatty_acid_synthesis'),
    ('XBP1', 'P_De_novo_fatty_acid_synthesis'),
    ('PRKAG2', 'P_Improvement_of_NAFLD'), ('PPARA', 'P_Improvement_of_NAFLD'),
    ('PPARG', 'P_Development_of_NAFLD'), ('CEBPA', 'P_Adipogenesis'),
    ('MAPK8', 'P_HCC_proliferation'), ('MAPK8', 'P_Apoptosis'),
    ('IL1A', 'P_Development_of_steatohepatitis'),
    ('IL6', 'P_Development_of_steatohepatitis'),
    ('TNF', 'P_Development_of_steatohepatitis'),
    ('FASLG', 'P_Cell_death'), ('TNF', 'P_Cell_death'),
    ('CXCL8', 'P_Inflammation'), ('TGFB1', 'P_Inflammation'),
    ('TGFB1', 'P_Fibrosis'), ('CASP3', 'P_Hepatocyte_injury'),
    ('CASP7', 'P_Hepatocyte_injury'),
]

assert len(_RAW_EDGES) == 106, f'expected 106 edges, got {len(_RAW_EDGES)}'


def _convert(name):
    return SYMBOL_TO_NODE.get(name, name)  # P_* names pass through unchanged


EDGES = [(_convert(s), _convert(t)) for s, t in _RAW_EDGES]
