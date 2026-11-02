"""Narrow source custody and exact arithmetic for the written theorem."""
import hashlib
import json
from pathlib import Path
from . import algebra, bounds

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENTS = (
    'model/__init__.py',
    'model/protected_spin_limit/README.md',
    'model/spin_goldstone_sector/README.md',
    'model/spin_particle_scaling/README.md',
)
FRAGMENT = 'paper/tex_fragments/SPIN_CANONICAL_FIELD.tex'


def bindings():
    paths = list(PARENTS)+[str(p.relative_to(ROOT)) for p in sorted(HERE.glob('*.py'))]
    paths += [str((HERE/'README.md').relative_to(ROOT)), FRAGMENT]
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def build():
    before = bindings()
    result = {
        'schema': 'oph.spin_canonical_field.v1',
        'source_bindings': before,
        'direct_exact_checks': algebra.exact_checks(),
        'constants': bounds.constants(),
        'transport_examples': [bounds.scaling_example(n) for n in (2, 4, 8, 16, 64)],
        'theorem': {
            'source': 'same cubic antiferromagnetic spin-1/2 H_J',
            'reference': 'same maximal staggered-field-selected zero-field ground',
            'scaled_coordinates': 'X=p*x, T=J*p*t, p->0',
            'canonical_pair': 'Q=p^2 sum eps_x F(px) Sx^3; P=-p sum G(px) Sx^2',
            'symplectic_form': 'm_s integral(F G_prime-G F_prime)',
            'regular_local_positive_energy_Weyl_net': True,
            'state_is_assumed_ergodic_or_Gaussian': False,
            'two_independent_canonical_pairs_proved': False,
            'common_Wightman_field_domain_proved': False,
            'surviving_connected_interaction_proved': False,
            'sharp_mass_shell_proved': False,
            'Lorentz_covariance_proved': False,
            'physical_species_identified': False,
            'effective_ordered_state_preparation_rate_proved': False,
        },
        'inherited_custody': {
            'kind': 'Explicit mathematical inputs with current source-byte bindings',
            'prior_full_producers_rerun_here': False,
            'finite_arithmetic_proves_infinite_limit_by_itself': False,
        },
        'primary_import': {
            'url': 'https://arxiv.org/html/1712.09018v4',
            'statement': 'Theorem 2.1: field-selected m_s is the supremum over all infinite-volume zero-field ground states',
        },
    }
    if bindings() != before:
        raise RuntimeError('source changed during replay')
    return result


def validate(value):
    if canonical(value) != canonical(build()):
        raise ValueError('source, arithmetic or scientific scope changed')
