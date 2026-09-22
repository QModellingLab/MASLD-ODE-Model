"""
model_loader.py
================
Locates and loads the hardcoded NASH/NAFL model files, using the same
search-path convention as equilibrate_then_inhibit.py, so this analysis
folder can sit anywhere under the repo (e.g.
analyses/edge_scramble_control/) and still find the models automatically.
"""
import os
import importlib.util

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))  # .../analyses/edge_scramble_control -> repo root

MODEL_FILES = {
    'NASH': 'ode_model_pydeseq2_NASH_vs_Normal_v10_mean.py',
    'NAFL': 'ode_model_pydeseq2_NAFL_vs_Normal_v10_mean.py',
}

_SEARCH_DIRS = [
    os.path.join(REPO_ROOT, 'data', 'models'),
    os.path.join(REPO_ROOT, 'data', 'models_pydeseq2_mean'),
    os.path.join(SCRIPT_DIR, 'models_pydeseq2_mean'),
    SCRIPT_DIR,
]


def _resolve(fname):
    for d in _SEARCH_DIRS:
        p = os.path.join(d, fname)
        if os.path.exists(p):
            return p
    return None


def load_models():
    """Returns (nash_module, nafl_module). Raises with a clear message
    (listing every path searched) if either file cannot be found."""
    mods = {}
    for cond, fname in MODEL_FILES.items():
        path = _resolve(fname)
        if path is None:
            searched = '\n    '.join(_SEARCH_DIRS)
            raise FileNotFoundError(
                f"Could not find '{fname}'. Searched:\n    {searched}\n"
                f"Fix: copy the two ode_model_pydeseq2_*.py files into "
                f"this analysis folder, or into data/models/ under the "
                f"repo root."
            )
        spec = importlib.util.spec_from_file_location(f'{cond}_model', path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        mods[cond] = mod
        print(f"Loaded {cond}: {path}")
    return mods['NASH'], mods['NAFL']
