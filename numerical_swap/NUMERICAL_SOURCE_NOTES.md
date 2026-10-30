# Numerical source provenance and packaging-only changes

Input: `w_swap_robustness(1).py` and `plot_w_swap_robustness(1).py`, uploaded by the user.

The SWAP model, SDP equations, numerical dual-residual audit, solver settings, and plotted values have **not** been deliberately modified.

Packaging-only edits to `src/w_swap_robustness.py`:
1. Removed a hash-only dependency on a historical `reference_2026/numerical_swap/w3_marginal_selftest_fig2.py` file, which was not supplied. The original implementation did **not** import this missing file for the numerical computations; it read it solely to compute an extra SHA256 environment entry. The package continues hashing its own solver script.
2. Replaced the potentially unavailable `ncpol2sdpa.__version__` with package metadata lookup `importlib.metadata.version("ncpol2sdpa")` in provenance logging.
3. Updated only the top-level comment to accurately describe the self-contained calculation source.

`src/plot_w_swap_robustness.py` and `results/all_curves.csv` are exact byte-for-byte copies of the supplied input files.

The exact old/revised certificate verifiers were copied without editing.
