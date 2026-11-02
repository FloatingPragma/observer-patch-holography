"""Reproducible arithmetic and source custody for the analytic bound."""
import hashlib
import json
from pathlib import Path
from .bounds import constants, example, power_exponents, F

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENTS = (
    'model/__init__.py',
    'model/centered_source_ground/README.md',
    'model/continuum_free_limit/README.md',
    'model/strong_coupling_refinement/README.md',
    'model/continuum_composite_channel/README.md',
)


def bindings():
    paths = list(PARENTS) + [str(p.relative_to(ROOT)) for p in sorted(HERE.glob('*.py'))]
    paths.append(str((HERE/'README.md').relative_to(ROOT)))
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def build():
    before = bindings()
    rows = []
    for alpha in (F(0), F(1, 5), F(1, 4), F(2, 7), F(1, 3)):
        powers, passed = power_exponents(alpha)
        rows.append({'coupling_power': str(alpha),
                     'powers': {k: str(v) for k, v in powers.items()},
                     'sufficient_free_condition': passed})
    result = {
        'schema': 'oph.ballistic_free_boundary.v1',
        'scope': 'Analytic sufficient free-continuum class for the centered source; exact arithmetic checks are not a numerical proof of a continuum limit.',
        'source_bindings': before,
        'constants': constants(), 'power_controls': rows,
        'finite_bounds': [example(b) for b in (2, 4, 8, 16, 64)],
        'theorem': {
            'clock': 'same centered K=H-(U/2)N',
            'spacing': 'c=v/a with fixed v>0',
            'sufficient_coupling_class': '0<=U(c)=o(c^(2/7))',
            'full_field_limit': 'free eight-valley Weyl vacuum hierarchy',
            'excitation_energy': 'v*|p|',
            'rest_masses_in_this_free_class': 'zero',
            'sharp_boundary_at_power_2_over_7_claimed': False,
            'interactions_above_boundary_proved': False,
            'physical_species_identified': False,
            'protected_spin_limit_excluded': False,
        },
    }
    if bindings() != before:
        raise RuntimeError('source changed during arithmetic replay')
    return result


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def validate(value):
    if canonical(value) != canonical(build()):
        raise ValueError('arithmetic, source custody or theorem scope changed')
