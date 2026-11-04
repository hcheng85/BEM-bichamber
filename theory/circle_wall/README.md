# Circular forward-ray wall theory validation

This directory validates `ANALYTICAL_NOTE.md` against the source-exact N=1,
wall-only specialization of `BEM_circular_parallel_instance.py`.  It introduces
no fitted parameters and does not claim to model collective switching.

Run from the repository root:

```bash
uv run --with numpy --with scipy --with matplotlib --with numba \
  python theory/circle_wall/run_validation.py
```

Use `--quick` only for a smoke test.  The committed tables and figures use the
full default independent-trajectory budgets. Exact settings and deterministic
seed formulas are also recorded in `analysis_config.json`.
