#!/usr/bin/env python3
"""Check the completeness/declared diagnostics of archived NUMERICAL SWAP data.

This does not re-solve the SDP, recover missing raw witnesses or provide
rigorous interval certification of floating-point lower bounds.
"""
import argparse
import csv
from collections import Counter
from pathlib import Path

MODELS = ('box', 'white')
PROTOCOLS = ('marginal', 'joint')
FIELDS = {'evidence_class', 'protocol', 'noise_model', 'epsilon', 'monomials', 'observations',
          'status', 'diagnostics_pass', 'numerical_lower_bound', 'gap_abs',
          'constraint_violation_max', 'min_moment_eigenvalue',
          'physical_witness_fidelity'}

def verify(path):
    with path.open(newline='') as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 48 or not rows or not FIELDS.issubset(rows[0]):
        raise ValueError('Expected 48 rows with complete audit fields')
    counts = Counter((r['protocol'], r['noise_model']) for r in rows)
    if any(counts[(p,m)] != 12 for p in PROTOCOLS for m in MODELS):
        raise ValueError(f'Unexpected protocol/noise-model groups: {counts}')
    eps_sets = []
    for protocol in PROTOCOLS:
        for model in MODELS:
            group=[r for r in rows if (r['protocol'],r['noise_model'])==(protocol,model)]
            eps=[float(r['epsilon']) for r in group]
            if len(set(eps)) != 12:
                raise ValueError(f'Duplicate epsilon point: {(protocol,model)}')
            eps_sets.append(set(eps))
            for r in group:
                if r['status'] != 'optimal' or r['diagnostics_pass'] != 'True':
                    raise ValueError(f'Unaccepted solver diagnostic: {r}')
                if r['evidence_class'] != 'NUMERICAL':
                    raise ValueError('Expected NUMERICAL evidence class')
                if int(r['monomials']) != (216 if protocol=='marginal' else 125):
                    raise ValueError('Unexpected basis size')
                if int(r['observations']) != (36 if protocol=='marginal' else 10):
                    raise ValueError('Unexpected declared observed moment count')
                b=float(r['numerical_lower_bound'])
                if not 0 <= b <= float(r['physical_witness_fidelity']) + 2e-5:
                    raise ValueError('Lower bound inconsistent with upper physical witness')
                if float(r['gap_abs']) > 2e-5 or float(r['constraint_violation_max']) > 2e-6:
                    raise ValueError('Gap or constraint diagnostics failed')
                if float(r['min_moment_eigenvalue']) < -2e-6:
                    raise ValueError('PSD diagnostic failed')
    if any(e!=eps_sets[0] for e in eps_sets[1:]):
        raise ValueError('Epsilon grid differs between groups')
    print(f'ARCHIVED NUMERICAL RESULTS CHECKED: {len(rows)} points, four groups; '
          'all declared diagnostics pass.')
    print('NOTE: these are numerical checks, not an exact numerical SDP certificate.')

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,default=Path(__file__).resolve().parents[1]/'results'/'all_curves.csv')
    args=parser.parse_args()
    verify(args.results)
