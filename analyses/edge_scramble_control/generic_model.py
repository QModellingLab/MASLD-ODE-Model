"""
generic_model.py
=================
Topology-driven ODE model for the MASLD hsa04932 network, built to exactly
reproduce the hardcoded ode_model_pydeseq2_{NASH,NAFL}_vs_Normal_v10_mean.py
model files -- so that we can then perform degree-preserving edge-scrambling
(Maslov & Sneppen, 2002, Science) on the 103 activating edges while holding
node initialization and kinetics fixed (response to Reviewer 3, Major
Comment 1).

This module is pure topology + kinetics logic. It does not read or modify
your hardcoded model files; it only needs their STATE_VARS, Y0, and PARAMS
at runtime (imported by the run scripts).
"""
import re
import random
import numpy as np

# ---------------------------------------------------------------------
# Gene symbol -> ODE node name mapping.
# Extracted programmatically (once) from the STATE_VARS comments in the
# hardcoded model files, e.g. "'ADIPOQ_inactive',  # ADIPOQ" -> ACDC.
# Kept here as a static table so this module has no dependency on the
# hardcoded files at import time.
# ---------------------------------------------------------------------
GENE_TO_NODE = {
    "ADIPOQ": "ACDC", "PRKAG2": "AMPK", "FOS": "AP_1", "MAP3K5": "ASK1", "ATF4": "ATF4",
    "AKT3": "Akt", "BAX": "Bax", "BID": "Bid", "BCL2L11": "Bim", "CASP3": "CASP3",
    "CASP7": "CASP7", "CASP8": "CASP8", "DDIT3": "CHOP", "CYP2E1": "CYP2E1",
    "CEBPA": "C_EBPa", "CDC42": "Cdc42", "MLXIP": "ChREBP", "NDUFC2-KCTD14": "CxI",
    "SDHA": "CxII", "UQCR11": "CxIII", "COX6B2": "CxIV", "CYCS": "Cytc", "FAS": "Fas",
    "FASLG": "FasL", "GSK3A": "GSK_3", "IKBKB": "IKKb", "IL1A": "IL_1", "IL6": "IL_6",
    "IL6R": "IL_6R", "CXCL8": "IL_8", "INS": "INS", "INSR": "INSR", "ERN1": "IRE1a",
    "IRS1": "IRS_1_2", "ITCH": "ITCH", "MAPK8": "JNK1_2", "LEP": "LEP", "NR1H3": "LXR_a",
    "PKLR": "L_PK", "MAP3K11": "MLK3", "NFKB1": "NF_kB", "LEPR": "ObR",
    "EIF2AK3": "PERK", "P3R3URF-PIK3R3": "PI3K", "PPARA": "PPAR_a", "PPARG": "PPAR_g",
    "RXRA": "RXR", "RAC1": "Rac1", "SOCS3": "SOCS3", "SREBF1": "SREBP_1c",
    "TGFB1": "TGF_b1", "TNFRSF1A": "TNFR1", "TNF": "TNFa", "TRAF2": "TRAF2",
    "XBP1": "XBP1", "ADIPOR1": "adipoR", "EIF2S1": "eIF2a", "MAPK14": "p38",
}

# ---------------------------------------------------------------------
# Supplementary Table S2: all 106 regulatory interactions (source, target,
# effect). "act" = activation, "inh" = inhibition. Gene symbols here match
# Table S2 / KEGG hsa04932; converted to ODE node names via GENE_TO_NODE
# (P_* names are already ODE node names and pass through unchanged).
# ---------------------------------------------------------------------
RAW_EDGES = [
    ("IL6","IL6R","act"), ("TNF","TNFRSF1A","act"), ("INS","INSR","act"),
    ("LEP","LEPR","act"), ("ADIPOQ","ADIPOR1","act"), ("TNF","MAPK8","act"),
    ("FASLG","FAS","act"), ("IL6R","SOCS3","act"), ("TNFRSF1A","NFKB1","act"),
    ("NFKB1","IL6","act"), ("INSR","IRS1","act"), ("IRS1","P3R3URF-PIK3R3","act"),
    ("P3R3URF-PIK3R3","AKT3","act"), ("AKT3","GSK3A","inh"), ("GSK3A","INS","act"),
    ("INSR","NR1H3","act"), ("NR1H3","RXRA","act"), ("RXRA","SREBF1","act"),
    ("MLXIP","PKLR","act"), ("LEPR","PRKAG2","act"), ("ADIPOR1","PRKAG2","act"),
    ("PRKAG2","MAPK14","act"), ("MAPK14","PPARA","act"), ("CDC42","MAP3K11","act"),
    ("RAC1","MAP3K11","act"), ("MAP3K11","MAPK8","act"), ("MAPK8","ITCH","act"),
    ("ITCH","CASP8","act"), ("FAS","CASP8","act"), ("CASP8","BID","act"),
    ("BID","BAX","act"), ("CYCS","CASP3","act"), ("CYCS","CASP7","act"),
    ("ERN1","XBP1","act"), ("ERN1","TRAF2","act"), ("TRAF2","MAP3K5","act"),
    ("MAP3K5","MAPK8","act"), ("TRAF2","IKBKB","act"), ("MAPK8","IRS1","inh"),
    ("MAPK8","FOS","act"), ("IKBKB","NFKB1","act"), ("FOS","IL1A","act"),
    ("FOS","IL6","act"), ("FOS","TNF","act"), ("NFKB1","IL1A","act"),
    ("NFKB1","TNF","act"), ("EIF2AK3","EIF2S1","act"), ("CYP2E1","MAPK8","act"),
    ("CYP2E1","IKBKB","act"), ("CYP2E1","FASLG","act"), ("CYP2E1","TNF","act"),
    ("CYP2E1","CXCL8","act"), ("CYP2E1","TGFB1","act"), ("EIF2S1","ATF4","act"),
    ("ATF4","DDIT3","act"), ("DDIT3","BCL2L11","act"), ("BCL2L11","BAX","act"),
    ("BAX","CYCS","act"), ("NDUFC2-KCTD14","MAPK8","act"), ("SDHA","MAPK8","act"),
    ("UQCR11","MAPK8","act"), ("COX6B2","MAPK8","act"), ("NDUFC2-KCTD14","IKBKB","act"),
    ("SDHA","IKBKB","act"), ("UQCR11","IKBKB","act"), ("COX6B2","IKBKB","act"),
    ("NDUFC2-KCTD14","FASLG","act"), ("SDHA","FASLG","act"), ("UQCR11","FASLG","act"),
    ("COX6B2","FASLG","act"), ("NDUFC2-KCTD14","TNF","act"), ("SDHA","TNF","act"),
    ("UQCR11","TNF","act"), ("COX6B2","TNF","act"), ("NDUFC2-KCTD14","CXCL8","act"),
    ("SDHA","CXCL8","act"), ("UQCR11","CXCL8","act"), ("COX6B2","CXCL8","act"),
    ("NDUFC2-KCTD14","TGFB1","act"), ("SDHA","TGFB1","act"), ("UQCR11","TGFB1","act"),
    ("COX6B2","TGFB1","act"), ("SOCS3","SREBF1","act"), ("XBP1","CEBPA","act"),
    ("RXRA","PPARG","act"), ("SOCS3","IRS1","inh"),
    ("GSK3A","P_Hyperinsulinemia","act"), ("SREBF1","P_De_novo_fatty_acid_synthesis","act"),
    ("MLXIP","P_De_novo_fatty_acid_synthesis","act"), ("XBP1","P_De_novo_fatty_acid_synthesis","act"),
    ("PRKAG2","P_Improvement_of_NAFLD","act"), ("PPARA","P_Improvement_of_NAFLD","act"),
    ("PPARG","P_Development_of_NAFLD","act"), ("CEBPA","P_Adipogenesis","act"),
    ("MAPK8","P_HCC_proliferation","act"), ("MAPK8","P_Apoptosis","act"),
    ("IL1A","P_Development_of_steatohepatitis","act"), ("IL6","P_Development_of_steatohepatitis","act"),
    ("TNF","P_Development_of_steatohepatitis","act"), ("FASLG","P_Cell_death","act"),
    ("TNF","P_Cell_death","act"), ("CXCL8","P_Inflammation","act"),
    ("TGFB1","P_Inflammation","act"), ("TGFB1","P_Fibrosis","act"),
    ("CASP3","P_Hepatocyte_injury","act"), ("CASP7","P_Hepatocyte_injury","act"),
]
assert len(RAW_EDGES) == 106, f"expected 106 edges, got {len(RAW_EDGES)}"


def _to_node(x):
    return x if x.startswith('P_') else GENE_TO_NODE[x]


EDGES = [(_to_node(s), _to_node(t), e) for s, t, e in RAW_EDGES]
ACT_EDGES = [(s, t) for s, t, e in EDGES if e == 'act']
INH_EDGES = [(s, t) for s, t, e in EDGES if e == 'inh']
assert len(ACT_EDGES) == 103 and len(INH_EDGES) == 3

# All 70 ODE node names (58 molecular + 12 pathway-output), in the same
# order used by the hardcoded STATE_VARS list. Provided here as a static
# fallback; run scripts normally derive node order from the hardcoded
# model's own STATE_VARS instead of this list, to guarantee consistency.
ALL_MOLECULAR_NODES = sorted(set(GENE_TO_NODE.values()))
ALL_P_NODES = sorted(set(t for _, t, _ in EDGES if t.startswith('P_')))

SILYMARIN_TARGETS = ['CASP3', 'CASP7', 'CASP8', 'CYP2E1', 'IL_8', 'TNFa', 'NF_kB', 'TGF_b1']
CORE_OUTPUTS = ['P_Hepatocyte_injury', 'P_Cell_death', 'P_Inflammation']


def build_generic_ode(state_vars, edge_list_act, edge_list_inh, ki_multipliers=None):
    """Return an ode_system(y, t, params) function built purely from the
    given edge lists (source, target) -- topology-driven, not hardcoded.

    state_vars: the STATE_VARS list from the hardcoded model (defines node
        order / indexing so results line up 1:1 with the original files).
    edge_list_act / edge_list_inh: (source_node, target_node) tuples using
        ODE node names (as in ACT_EDGES / INH_EDGES above).
    ki_multipliers: optional {node_name: ki} dict; multiplies that node's
        forward (activation) flux by ki, reproducing the silymarin
        intervention (Eq. 4) without any source-code regex hacking.
    """
    idx = {name: i for i, name in enumerate(state_vars)}
    n_states = len(state_vars)
    node_names = [s[:-9] for s in state_vars if s.endswith('_inactive')]
    ki_multipliers = ki_multipliers or {}

    act_reg = {node: [] for node in node_names}
    for s, t in edge_list_act:
        act_reg[t].append(s)
    inh_reg = {node: [] for node in node_names}
    for s, t in edge_list_inh:
        inh_reg[t].append(s)

    def ode_system(y, t, params=None):
        if params is None:
            params = {'Vmax': 2.0, 'ksp': 2.0, 'n': 2.0, 'kcat': 2.0}
        Vmax, ksp, n, kcat = params['Vmax'], params['ksp'], params['n'], params['kcat']
        dydt = np.zeros(n_states)
        for node in node_names:
            i_idx = idx[node + '_inactive']
            a_idx = idx[node + '_active']
            inactive = y[i_idx]
            active = y[a_idx]
            regs = act_reg[node]
            if regs:
                fwd = sum(kcat * y[idx[r + '_active']] * inactive**n / (ksp**n + inactive**n)
                          for r in regs)
            else:
                fwd = Vmax * inactive**n / (ksp**n + inactive**n)
            for r in inh_reg[node]:
                fwd = fwd * (ksp**n / (ksp**n + y[idx[r + '_active']]**n))
            if node in ki_multipliers:
                fwd = fwd * ki_multipliers[node]
            bwd = Vmax * active**n / (ksp**n + active**n)
            dydt[i_idx] = -fwd + bwd
            dydt[a_idx] = fwd - bwd
        return dydt

    return ode_system


def degree_preserving_scramble(edge_list, n_swaps, rng):
    """Maslov-Sneppen double edge-swap on a directed edge list of (source,
    target) tuples. Preserves in-degree and out-degree of every node
    exactly (Maslov & Sneppen, 2002, Science). P_* nodes never appear as
    sources in the input edge list, so a double edge-swap -- which only
    ever recombines existing sources with existing targets from OTHER
    edges -- can never make a P_* node a source; "P_* sink nodes remain
    sinks" is therefore preserved automatically, with no extra logic
    needed.

    n_swaps: number of successful swaps to perform. Milo et al. (2002)
        recommend on the order of 10-100x the edge count for thorough
        randomization; the run scripts default to 20x.
    rng: a random.Random instance (pass one with a fixed seed for
        reproducibility).
    """
    edges = list(edge_list)
    edge_set = set(edges)
    n = len(edges)
    swapped = 0
    attempts = 0
    max_attempts = n_swaps * 50
    while swapped < n_swaps and attempts < max_attempts:
        attempts += 1
        i, j = rng.sample(range(n), 2)
        a, b = edges[i]
        c, d = edges[j]
        if a == c or b == d or a == d or c == b:
            continue
        new1, new2 = (a, d), (c, b)
        if new1 in edge_set or new2 in edge_set:
            continue
        edge_set.discard((a, b)); edge_set.discard((c, d))
        edge_set.add(new1); edge_set.add(new2)
        edges[i], edges[j] = new1, new2
        swapped += 1
    return edges, swapped


def erdos_renyi_edges(n_edges, source_pool, target_pool, rng):
    """Uniform-random directed edges (source_pool -> target_pool), no
    self-loops, no duplicate edges. Unlike degree_preserving_scramble,
    this does NOT constrain the degree sequence to match hsa04932 -- it
    is a genuinely different null-model class, used as a robustness check
    on whether findings depend on matching hsa04932's specific degree
    distribution (see README, "Two null models" section)."""
    edges = set()
    while len(edges) < n_edges:
        s = rng.choice(source_pool)
        t = rng.choice(target_pool)
        if s == t:
            continue
        edges.add((s, t))
    return list(edges)
