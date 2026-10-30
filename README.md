# Exact self-testing of the three-qubit W state

This repository accompanies *Exact Self-Testing of the Three-Qubit W State from Marginal Correlations* (Physics Letters A manuscript PHYLA-D-26-02336, revision R1). The `main` branch contains the exact ideal-point certificate, the revised analytic robustness certificate, and finite-noise SWAP/NPA computation code and data.

The self-testing result concerns the complete ideal one- and two-body marginal behavior of three binary settings `Z`, `X`, and `D` per party. No three-body correlator is required by the marginal self-test. A separate two-setting scenario with a three-body correlator is included only as a numerical comparison.

## Contents

| Location | Purpose | Evidence |
|---|---|---|
| `certificate/` and `verifier/` | Original ideal-point level-2 certificate and standalone verifier | Exact algebra over `Q(sqrt(2))` |
| `robustness/` | Revised Gram certificate, branch-alignment domination and robustness constants | Exact algebra and rational positivity checks |
| `numerical_swap/src/` | Finite-level SWAP/NPA solver, plotting and archived-data checks | Numerical optimization |
| `numerical_swap/results/` | 48 archived fidelity estimates, grouped data and figures | Floating-point numerical estimates |
| `verification_runs/` | Logs from checks performed for the supplementary package | Recorded checks |
| `SOURCE_PROVENANCE.json` and `MANIFEST_SHA256.txt` | Source provenance and content-integrity manifest | File provenance |

The original `certificate/`, `verifier/`, `CITATION.cff`, and `SHA256SUMS.txt` remain available at their established paths. The revised certificate in `robustness/` is distinct from the original ideal-point certificate.

## Verify the exact certificates

Python 3.11 and SymPy 1.13.1 are the pinned dependencies for the exact checkers.

```bash
python -m pip install -r requirements.txt
python verifier/verify_certificate.py
python robustness/verification/verify_revision.py
```

The original checker verifies the 55-word basis, reduced affine constraints, the scalar bound `u_bar >= sqrt(2)`, Gram positivity and the ideal W nullspace. The revised checker also verifies the exact branch-alignment Gram domination and the constants appearing in the analytic robustness bound.

Run the combined checks from the repository root:

```bash
python verify_all.py
```

This command checks the two exact certificates and performs a consistency check of the archived numerical CSV. It does not re-solve the numerical SDPs.

## Numerical SWAP/NPA results

The marginal protocol uses all 36 one- and two-body moments of `Z`, `X`, and `D`. The two numerical noise models are (i) componentwise absolute-error boxes and (ii) a specified white-noise path. The comparison `joint` protocol uses a three-body correlator, so it makes different observational assumptions.

The 48 archived points are in `numerical_swap/results/all_curves.csv`, and the plotted figures are under `numerical_swap/results/figures/`. To check the archived table and reproduce the plots:

```bash
python numerical_swap/src/verify_archived_curves.py
python numerical_swap/src/plot_w_swap_robustness.py --results numerical_swap/results
```

For new numerical optimizations, install compatible versions of the packages in `numerical_swap/requirements_numerical.txt` and configure a MOSEK license. For example:

```bash
python numerical_swap/src/w_swap_robustness.py --protocol marginal --model box --output-dir numerical_swap/recomputed/marginal_box
python numerical_swap/src/w_swap_robustness.py --protocol marginal --model white --output-dir numerical_swap/recomputed/marginal_white
```

The solver writes per-point matrices, logs, and diagnostic files to the selected new output directory. The archived CSV consists of floating-point, dual-residual-adjusted estimates: it is not an exact or interval-certified SDP proof. Original point-by-point MOSEK matrices and logs were not included in the supplied archive; fresh runs can generate those files. See `numerical_swap/results/RESULTS_SCOPE.md` and `numerical_swap/NUMERICAL_SOURCE_NOTES.md` for details.

## Integrity and credit

To check the manifest from the repository root:

```bash
sha256sum -c MANIFEST_SHA256.txt
```

The original ideal-point software authorship and license information is preserved in `CITATION.cff` and `LICENSE`. The associated paper has its own author list. The manuscript's full mathematical proof, assumptions and exact formulas are in the article; the repository provides the machine-readable certificates and research code supporting those claims.
