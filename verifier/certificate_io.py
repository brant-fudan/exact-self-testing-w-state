r"""Decode the exact marginal-W certificate over :math:`\mathbb{Q}(\sqrt{2})`.

Each JSON scalar is stored as ``[a, b]`` and represents

.. math:: a + b\sqrt{2}, \qquad a,b\in\mathbb{Q}.

Only the language-neutral JSON format used by the published verifier is
supported. Certificate-generation and legacy conversion utilities are kept
out of this release.
"""
from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import sympy as sp

SQRT2 = sp.sqrt(2)
FORMAT_NAME = "marginal_w_exact_certificate_v1"


def pair_to_expr(pair: Iterable[str]) -> sp.Expr:
    r"""Decode ``[a, b]`` as the exact scalar :math:`a+b\sqrt{2}`."""
    a, b = pair
    return sp.Rational(a) + sp.Rational(b) * SQRT2


def matrix_from_json(obj: dict[str, Any]) -> sp.MutableDenseMatrix:
    r"""Decode a matrix whose entries lie in :math:`\mathbb{Q}(\sqrt{2})`."""
    rows = int(obj["rows"])
    cols = int(obj["cols"])
    data = obj["data"]
    if len(data) != rows:
        raise ValueError(
            f"matrix row count mismatch: expected {rows}, got {len(data)}"
        )

    entries: list[sp.Expr] = []
    for row in data:
        if len(row) != cols:
            raise ValueError(
                f"matrix column count mismatch: expected {cols}, got {len(row)}"
            )
        entries.extend(pair_to_expr(value) for value in row)
    return sp.Matrix(rows, cols, entries)


def word_from_json(word: list[str]) -> tuple[str, str, str]:
    """Decode a three-party operator word."""
    if not (
        isinstance(word, list)
        and len(word) == 3
        and all(isinstance(local_word, str) for local_word in word)
    ):
        raise ValueError(f"invalid JSON word: {word!r}")
    return tuple(word)  # type: ignore[return-value]


def load_certificate(path: str | Path) -> dict[str, Any]:
    """Load the language-neutral exact certificate from ``path``."""
    certificate_path = Path(path)
    if certificate_path.suffix.lower() != ".json":
        raise ValueError("the published certificate must use the JSON format")

    with certificate_path.open("r", encoding="utf-8") as stream:
        raw = json.load(stream)
    if raw.get("format") != FORMAT_NAME:
        raise ValueError("unknown certificate JSON format")

    return {
        "DEN": int(raw["DEN"]),
        "P_den": int(raw["P_den"]),
        "P_int": [[int(value) for value in row] for row in raw["P_int"]],
        "B": [word_from_json(word) for word in raw["B"]],
        "reps": [word_from_json(word) for word in raw["reps"]],
        "A": matrix_from_json(raw["A"]),
        "b": matrix_from_json(raw["b"]),
        "c": matrix_from_json(raw["c"]),
        "Bmat": matrix_from_json(raw["Bmat"]),
        "lam": [pair_to_expr(value) for value in raw["lam"]],
        "S_upper": [pair_to_expr(value) for value in raw["S_upper"]],
        "pairs": [tuple(int(value) for value in pair) for pair in raw["pairs"]],
        "R": matrix_from_json(raw["R"]),
        "target_witnesses": {
            name: [pair_to_expr(value) for value in witness]
            for name, witness in raw["target_witnesses"].items()
        },
    }


def contains_float_atoms(obj: Any) -> bool:
    """Return whether a nested exact object contains a SymPy ``Float``."""
    if isinstance(obj, sp.Basic):
        return bool(obj.atoms(sp.Float))
    if isinstance(obj, sp.MatrixBase):
        return any(contains_float_atoms(value) for value in obj)
    if isinstance(obj, dict):
        return any(contains_float_atoms(value) for value in obj.values())
    if isinstance(obj, (list, tuple)):
        return any(contains_float_atoms(value) for value in obj)
    return False
