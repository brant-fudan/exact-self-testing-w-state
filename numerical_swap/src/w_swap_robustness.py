"""Audited numerical SWAP bounds for the marginal and joint-data W protocols.

Finite-level numerical SWAP/NPA model and floating-point dual residual audit.
No external historical source files are needed for a fresh numerical solve.
Only reflection squares and different-party commutation are algebraic rules.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import importlib.metadata
import itertools
import json
import platform
import shlex
import sys
import time
from pathlib import Path

import ncpol2sdpa
import numpy as np
import sympy as sp
from ncpol2sdpa import SdpRelaxation, generate_operators
from scipy import sparse

PARTIES = ("A", "B", "C")
W_SUPPORT = ((0, 0, 1), (0, 1, 0), (1, 0, 0))
DEFAULT_EPSILONS = ("0.00001", "0.00005", "0.0001", "0.0002", "0.0005",
                    "0.001", "0.002", "0.003", "0.004", "0.006", "0.008", "0.01")


def operator_map(protocol):
    """Generate independent Hermitian reflections; D has no algebraic link to X/Z."""
    if protocol not in ("marginal", "joint"):
        raise ValueError("protocol must be marginal or joint")
    names = ("Z", "X", "D") if protocol == "marginal" else ("Z", "X")
    ops = {p: dict(zip(names, generate_operators(p, len(names), hermitian=True))) for p in PARTIES}
    return [o for p in PARTIES for o in ops[p].values()], ops


def reflection_substitutions(ops):
    r"""Impose O^2=I and [O_p,O_q]=0 only for p != q."""
    rules = {o**2: sp.S.One for p in PARTIES for o in ops[p].values()}
    for p, q in itertools.combinations(PARTIES, 2):
        for a in ops[p].values():
            for b in ops[q].values():
                rules[b * a] = a * b
    return rules


def observations(ops, protocol):
    r"""Return (name, polynomial, exact ideal mean) for the declared protocol.

    Marginal: 9 one-body plus 27 two-body moments of Z,X,D. Joint: three
    Z_i, three Z_i Z_j, three X_i X_j, and ONLY the normalized orbit
    (ZXX+XZX+XXZ)/3 = 2/3. Its box tolerance epsilon means a tolerance
    3 epsilon on the unnormalized sum. No other joint moments are fixed.
    """
    one = {"Z": sp.Rational(1, 3), "X": sp.S.Zero, "D": sp.sqrt(2) / 6}
    pair = {("Z", "Z"): -sp.Rational(1, 3), ("Z", "X"): sp.S.Zero,
            ("X", "Z"): sp.S.Zero, ("Z", "D"): -sp.sqrt(2) / 6,
            ("D", "Z"): -sp.sqrt(2) / 6, ("X", "X"): sp.Rational(2, 3),
            ("X", "D"): sp.sqrt(2) / 3, ("D", "X"): sp.sqrt(2) / 3,
            ("D", "D"): sp.Rational(1, 6)}
    result = []
    for p in PARTIES:
        for name in (("Z", "X", "D") if protocol == "marginal" else ("Z",)):
            result.append((f"{p}_{name}", ops[p][name], one[name]))
    for p, q in itertools.combinations(PARTIES, 2):
        for (a, b), value in pair.items():
            if protocol == "joint" and (a, b) not in (("Z", "Z"), ("X", "X")):
                continue
            result.append((f"{p}_{a}_{q}_{b}", ops[p][a] * ops[q][b], value))
    if protocol == "joint":
        orbit = sum(sp.prod(ops[p]["Z" if p == excited else "X"] for p in PARTIES)
                    for excited in PARTIES) / 3
        result.append(("C_W/3", sp.expand(orbit), sp.Rational(2, 3)))
    return result


def noise_constraints(ops, protocol, epsilon, model):
    r"""Box: |<O>-c|<=epsilon; white: <O>=(1-epsilon)c, as scalar moments.

    Evaluate irrational constants at the numerical SDP boundary. In the
    installed ncpol2sdpa parser, an unevaluated sqrt(2) can be classified as
    a scalar monomial and lose its value. Keep exact ideals in observations().
    """
    eps = sp.Rational(str(epsilon))
    if not 0 <= eps <= 1:
        raise ValueError("epsilon must be in [0,1]")
    eq, ineq = [], []
    for _, word, value in observations(ops, protocol):
        if not value.is_Rational:
            value = value.evalf(17)
        if model == "box":
            ineq.extend((eps + word - value, eps - word + value))
        elif model == "white":
            eq.append(word - (1 - eps) * value)
        else:
            raise ValueError("model must be box or white")
    return eq, ineq


def fidelity_objective(ops):
    r"""Return F_SWAP = (1/3) sum_{u,v in W support} K_v^dagger K_u.

    K_0=(I+Z)/2, K_1=X(I-Z)/2; K_u is the ordered tensor product over
    parties A,B,C. K_0^dagger K_0+K_1^dagger K_1=I for any reflections.
    The expectation is <W|rho_out|W>, NOT its square root.
    Same-party multiplication order is retained before and after expansion.
    """
    local = {p: [(1 + ops[p]["Z"]) / 2,
                 ops[p]["X"] * (1 - ops[p]["Z"]) / 2] for p in PARTIES}
    k = {bits: sp.expand(sp.prod(local[p][b] for p, b in zip(PARTIES, bits))) for bits in W_SUPPORT}
    return sp.expand(sum(sp.adjoint(k[v]) * k[u] for u in W_SUPPORT for v in W_SUPPORT) / 3)


def extra_monomials(ops, basis="base"):
    r"""Adjoint-closed local tensor basis {I,Z,X,XZ,ZX,(D)}^3.

    The extended basis adds XZX and ZXZ to each local set. Cartesian products
    contain all SWAP branch words and hence every objective matrix entry.
    ncpol2sdpa level 1 already supplies identity and single measurements.
    """
    automatic = {sp.S.One} | {o for p in PARTIES for o in ops[p].values()}
    choices = []
    for p in PARTIES:
        z, x = ops[p]["Z"], ops[p]["X"]
        local = [sp.S.One, z, x, x * z, z * x]
        if "D" in ops[p]:
            local.append(ops[p]["D"])
        if basis == "extended":
            local.extend((x * z * x, z * x * z))
        elif basis != "base":
            raise ValueError("basis must be base or extended")
        choices.append(local)
    return [sp.prod(words) for words in itertools.product(*choices) if sp.prod(words) not in automatic]


def build_relaxation(protocol, epsilon, model, basis="base"):
    """Build the outer moment relaxation without any ideal-only consequences."""
    variables, ops = operator_map(protocol)
    eq, ineq = noise_constraints(ops, protocol, epsilon, model)
    relaxation = SdpRelaxation(variables, verbose=0)
    relaxation.get_relaxation(1, objective=fidelity_objective(ops),
                             substitutions=reflection_substitutions(ops),
                             momentequalities=eq, momentinequalities=ineq,
                             extramonomials=extra_monomials(ops, basis))
    relaxation._w_ops = ops
    return relaxation


def _dual_vector(blocks):
    """Contract symmetric matrices with an upper-triangle-only affine SDP map."""
    return np.concatenate([(q * (2 * np.triu(np.ones_like(q), 1) + np.eye(len(q)))).ravel()
                           for q in blocks])


def audit_dual(coefficient_map, objective, constant, blocks, guard=1e-12):
    r"""Compute a floating-point safeguarded dual bound, NOT an exact certificate.

    For G(y)=F0+sum Fi yi >=0 and f=c0+c.y, any Q>=0 gives
    f >= c0-<F0,Q>-sum_i |c_i-<Fi,Q>|, since |yi|<=1 for unitary words.
    Q is shifted by max(0,-lambda_min(Q))+guard. No moment substitutions or
    equality elimination are used, so each variable is a pure word moment.
    Eigenvalues and residuals remain floating point, requiring NUMERICAL labels.
    """
    raw, corrected, eigs, shifts = [], [], [], []
    for block in blocks:
        q = np.asarray(block, dtype=np.float64)
        q = (q + q.T) / 2
        eig = float(np.linalg.eigvalsh(q).min())
        shift = max(0., -eig) + guard
        raw.append(q)
        corrected.append(q + shift * np.eye(len(q)))
        eigs.append(eig)
        shifts.append(shift)
    raw_coeff = np.asarray(coefficient_map.T @ _dual_vector(raw)).ravel()
    coeff = np.asarray(coefficient_map.T @ _dual_vector(corrected)).ravel()
    residual = objective - coeff[1:]
    bound = float(constant - coeff[0] - np.abs(residual).sum())
    return {"dual_witness_objective": float(constant - raw_coeff[0]),
            "dual_coefficient_residual_max": float(np.max(np.abs(objective - raw_coeff[1:]))),
            "corrected_dual_residual_l1": float(np.abs(residual).sum()),
            "min_dual_eigenvalue": min(eigs), "dual_diagonal_shift_max": max(shifts),
            "uncapped_corrected_bound": bound, "numerical_lower_bound": max(0., bound)}


def physical_witness_fidelity(epsilon, model):
    r"""Feasible upper sanity bound from rho_eta=(1-eta)|W><W|+eta I/8.

    All protocol means have |c_j|<=2/3, so eta=min(3 epsilon/2,1) fits
    the box. With ideal Pauli X/Z, F_SWAP=1-7 eta/8. This is a feasible
    realization, NOT the DI lower bound that the SDP computes.
    """
    eta = epsilon if model == "white" else min(1.5 * epsilon, 1.)
    return 1 - 7 * eta / 8


def solve_point(relaxation, epsilon, protocol, model, basis, directory, tolerance, threads):
    """Solve and archive primal moments, dual matrices, affine map, and diagnostics."""
    eps = float(epsilon)
    point = directory / f"eps_{str(epsilon).replace('.', 'p')}"
    point.mkdir(exist_ok=False)
    eq, ineq = noise_constraints(relaxation._w_ops, protocol, epsilon, model)
    relaxation.verbose = 0
    relaxation.process_constraints(momentequalities=eq, momentinequalities=ineq)
    parameters = {"iparam.num_threads": threads,
                  "dparam.intpnt_co_tol_pfeas": tolerance,
                  "dparam.intpnt_co_tol_dfeas": tolerance,
                  "dparam.intpnt_co_tol_rel_gap": tolerance}
    relaxation.verbose = 1
    started = time.perf_counter()
    with (point / "solver.log").open("w") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        relaxation.solve("mosek", solverparameters=parameters)
    elapsed = time.perf_counter() - started
    if str(relaxation.status) != "optimal":
        raise RuntimeError(f"{point}: non-optimal status {relaxation.status}")
    measured = observations(relaxation._w_ops, protocol)
    observed = np.array([float(relaxation[word]) for _, word, _ in measured])
    ideals = np.array([float(value) for _, _, value in measured])
    violation = (np.maximum(np.abs(observed - ideals) - eps, 0).max() if model == "box"
                 else np.abs(observed - (1 - eps) * ideals).max())
    moment_fidelity = float(relaxation[fidelity_objective(relaxation._w_ops)])
    objective = np.asarray(relaxation.obj_facvar, dtype=np.float64)
    coefficient_map = relaxation.F[:, :len(objective) + 1].tocsr()
    audit = audit_dual(coefficient_map, objective, relaxation.constant_term, relaxation.y_mat)
    row = {"evidence_class": "NUMERICAL", "protocol": protocol, "noise_model": model,
           "epsilon": eps, "basis": basis, "monomials": len(relaxation.monomial_sets[0]),
           "sdp_variables": relaxation.n_vars, "observations": len(measured),
           "status": str(relaxation.status), "solve_seconds": elapsed,
           "ncpol_primal": float(relaxation.primal), "ncpol_dual": float(relaxation.dual),
           "gap_abs": abs(float(relaxation.primal - relaxation.dual)),
           "constraint_violation_max": float(violation), "moment_fidelity": moment_fidelity,
           "min_moment_eigenvalue": min(float(np.linalg.eigvalsh((np.asarray(q) + np.asarray(q).T) / 2).min())
                                        for q in relaxation.x_mat),
           "physical_witness_fidelity": physical_witness_fidelity(eps, model), **audit}
    # The MOSEK adapter reverses SDP primal/dual roles; verify via actual matrices.
    row["dual_objective_match_residual"] = min(abs(audit["dual_witness_objective"] - float(v))
                                               for v in (relaxation.primal, relaxation.dual))
    row["moment_objective_match_residual"] = min(abs(moment_fidelity - float(v))
                                                 for v in (relaxation.primal, relaxation.dual))
    row["diagnostics_pass"] = bool(
        row["gap_abs"] <= 2e-5 and row["constraint_violation_max"] <= 2e-6
        and row["min_moment_eigenvalue"] >= -2e-6
        and row["dual_objective_match_residual"] <= 2e-6
        and row["moment_objective_match_residual"] <= 2e-6
        and row["dual_coefficient_residual_max"] <= 2e-6
        and row["numerical_lower_bound"] <= row["physical_witness_fidelity"] + 2e-5)
    arrays = {"objective": objective, "constant_term": np.asarray(relaxation.constant_term),
              "observed": observed, "ideal": ideals}
    for prefix, blocks in (("moment", relaxation.x_mat), ("dual", relaxation.y_mat)):
        arrays.update({f"{prefix}_{i}": np.asarray(q) for i, q in enumerate(blocks)})
    np.savez_compressed(point / "solution.npz", **arrays)
    sparse.save_npz(point / "coefficient_map.npz", coefficient_map)
    (point / "result.json").write_text(json.dumps(row, indent=2) + "\n")
    return row


def write_table(rows, path):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", choices=("marginal", "joint"), required=True)
    parser.add_argument("--model", choices=("box", "white"), required=True)
    parser.add_argument("--basis", choices=("base", "extended"), default="base")
    parser.add_argument("--epsilons", nargs="+", default=DEFAULT_EPSILONS)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    source_paths = [Path(__file__)]
    environment = {"python": sys.version, "executable": sys.executable,
                   "platform": platform.platform(), "hostname": platform.node(),
                   "device": "CPU", "dtype": "float64; complex128 matrix tests", "seed": 20261016,
                   "random_sampling_in_sdp": False, "ncpol2sdpa_module_version": importlib.metadata.version("ncpol2sdpa"),
                   "package_versions": {name: importlib.metadata.version(name) for name in
                                        ("numpy", "scipy", "sympy", "ncpol2sdpa", "Mosek", "matplotlib")},
                   "source_sha256": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                                     for p in source_paths}}
    _, ops = operator_map(args.protocol)
    parameters = {**vars(args), "output_dir": str(args.output_dir),
                  "observations": [{"name": n, "polynomial": str(p), "ideal": str(v)}
                                   for n, p, v in observations(ops, args.protocol)],
                  "propagate_ideal_relations": False, "objective": "<W|rho_out|W>",
                  "evidence_class": "NUMERICAL"}
    (args.output_dir / "environment.json").write_text(json.dumps(environment, indent=2) + "\n")
    (args.output_dir / "parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
    (args.output_dir / "command.txt").write_text(f"cwd={Path.cwd()}\n" +
        shlex.join([sys.executable, "src/w_swap_robustness.py", *sys.argv[1:]]) + "\n")
    print(f"Building {args.protocol} / {args.model} / {args.basis}", flush=True)
    start = time.perf_counter()
    relaxation = build_relaxation(args.protocol, args.epsilons[0], args.model, args.basis)
    basis_text = "\n".join(map(str, relaxation.monomial_sets[0])) + "\n"
    (args.output_dir / "basis.txt").write_text(basis_text)
    (args.output_dir / "basis.sha256").write_text(hashlib.sha256(basis_text.encode()).hexdigest() + "\n")
    print(f"basis={len(relaxation.monomial_sets[0])} variables={relaxation.n_vars} build_seconds={time.perf_counter()-start:.2f}", flush=True)
    rows = []
    for epsilon in args.epsilons:
        print(f"Solving epsilon={epsilon}", flush=True)
        row = solve_point(relaxation, epsilon, args.protocol, args.model, args.basis,
                          args.output_dir, args.tolerance, args.threads)
        rows.append(row)
        write_table(rows, args.output_dir / "curve.csv")
        print(json.dumps(row), flush=True)
        if not row["diagnostics_pass"]:
            raise RuntimeError("Point failed diagnostics; saved for inspection, not for plotting")
    if not all(row["diagnostics_pass"] for row in rows):
        raise RuntimeError("At least one point failed diagnostics; inspect raw outputs")
    print("SWAP_ROBUSTNESS_SWEEP_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
