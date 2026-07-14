# Exact self-testing of the three-qubit W state

This repository contains the exact certificate and standalone verifier for the
computer-assisted step in the strict marginal self-test of the three-qubit
W state. The associated theorem uses only the complete ideal one- and two-body
marginal behavior of the local observables `Z`, `X`, and
`D = (Z + X) / sqrt(2)`.

The verifier establishes the exact level-2 moment-matrix bound

```text
u_bar >= sqrt(2)
```

over `Q(sqrt(2))`. It does not contain certificate-discovery code, numerical
SDP searches, manuscript files, the incomplete certificate-free route, or the
separate joint-data `[2,2,2]` theorem.

## Contents

- `certificate/w_exact_certificate.json`: exact certificate data over
  `Q(sqrt(2))`;
- `verifier/verify_certificate.py`: standalone full verifier;
- `verifier/certificate_io.py`: decoder for the language-neutral JSON format;
- `SHA256SUMS.txt`: integrity manifest for the release files.

## Requirements

- Python 3.11;
- SymPy 1.13.1.

Install the pinned dependency with

```bash
python -m pip install -r requirements.txt
```

MOSEK and other SDP solvers are not required.

## Verification

From the repository root, run

```bash
python verifier/verify_certificate.py
```

The verifier independently reconstructs and checks:

1. the absence of floating-point atoms;
2. the 55-word level-2 basis and 110 real adjoint-orbit coordinates;
3. equivalence of the stored and generated affine constraint spaces;
4. the Gram-adjoint map;
5. the exact dual identity and objective value over `Q(sqrt(2))`;
6. positive definiteness of the Gram matrix by rational interval arithmetic
   and Gershgorin's theorem;
7. the ideal W-state nullspace;
8. containment of the target relation vectors in the certificate column
   space.

A successful run ends with

```text
FULL CERTIFICATE VERIFIED.
```

## Scope

The certificate supplies one scalar lower bound in the proof. The remaining
state-dependent rigidity argument and local-isometry extraction are analytic.
The repository does not claim a finite-noise robustness theorem or a global
self-test of the measurement operators.

## License

The verification code is released under the MIT License.
