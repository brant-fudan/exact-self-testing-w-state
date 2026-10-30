#!/usr/bin/env python3
"""Check exact ideal/revised certificates and the archived numerical data.

This command does not re-solve the numerical SWAP/NPA optimizations.
"""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
STEPS = [
    ("original ideal-point certificate", ROOT / "verifier/verify_certificate.py"),
    ("revised exact certificate and robustness bounds", ROOT / "robustness/verification/verify_revision.py"),
    ("archived numerical curve consistency", ROOT / "numerical_swap/src/verify_archived_curves.py"),
]


def main():
    for description, script in STEPS:
        print(f"\n=== {description} ===", flush=True)
        subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)
    print("\nALL PACKAGE CHECKS PASSED (numerical SDPs not re-solved).")


if __name__ == "__main__":
    main()
