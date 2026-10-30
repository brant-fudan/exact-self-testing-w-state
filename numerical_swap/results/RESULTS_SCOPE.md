# Numerical result files

* `all_curves.csv` is the unmodified supplied 48-point summary.
* `{marginal,joint}_{box,white}/curve.csv` re-groups those original rows into four 12-point inputs for the supplied plotting script.
* `numerical_lower_bound` is a **floating-point, residual-adjusted** SWAP/NPA estimate, not an exact rigorous SDP certificate.
* Original per-point SDP logs, dual matrices, primal matrices and coefficient maps were not provided and are not part of this package. To re-create them, use the numerical solver and a compatible MOSEK environment.
* `diagnostics_pass=True` reflects the supplied numerical post-solve checks; it is not proof that all floating-point rounding effects have been bounded rigorously.
