"""Narrow, current source bindings for the written continuum-energy proof."""
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from . import algebra, bounds

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENTS = ('model/__init__.py', 'model/spin_canonical_field/README.md',
           'model/spin_goldstone_sector/README.md',
           'model/spin_particle_scaling/README.md',
           'model/spin_shell_diagnostic/README.md')
FRAGMENT = 'paper/tex_fragments/SPIN_CANONICAL_ENERGY.tex'


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def bindings():
    paths = list(PARENTS) + [str(p.relative_to(ROOT)) for p in sorted(HERE.glob('*.py'))]
    paths += [str((HERE/'README.md').relative_to(ROOT)), FRAGMENT]
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in paths}


def build():
    before = bindings()
    m = Fraction(1, 1000)
    result = {
        'schema': 'oph.spin_canonical_energy.v1',
        'source_bindings': before,
        'direct_exact_checks': algebra.exact_checks(),
        'stiffness_enclosure': bounds.serialize(bounds.stiffness_interval(m)),
        'spectral_windows': bounds.serialize([
            bounds.window(channel, m, radius)
            for radius in (1, Fraction(1, 2), Fraction(1, 4))
            for channel in ('P', 'Q')]),
        'rotation_remainder_coefficient': str(bounds.rotation_remainder_coefficient()),
        'energy_escape_controls': bounds.serialize([bounds.escaping_energy(n) for n in (2, 16, 1024)]),
        'continuous_spectral_control': bounds.serialize(bounds.continuous_countercontrol(m, 1)),
        'theorem': {
            'source': 'same cubic spin-1/2 antiferromagnet H_J with each bond counted once',
            'reference': 'same infinitesimal-axis-1-staggered-field-selected zero-field ground',
            'scaled_clock': 'L_p=K/(Jp); X=px and T=Jpt',
            'actual_Stone_field_vectors_identified': True,
            'canonical_fields_have_no_zero_energy_component': True,
            'nonzero_positive_band_at_arbitrarily_small_momentum': True,
            'actual_continuum_dynamics_nonstatic_and_gapless': True,
            'kappa': '-omega(S0^1 Se^1+S0^3 Se^3), m^2/2<=kappa<=1/2',
            'exact_microscopic_twist_energy_limit': 'kappa/2 integral|grad G|^2+6 kappa integral F^2',
            'continuum_twist_energy_equality_proved': False,
            'strong_unbounded_field_convergence_assumed': False,
            'state_assumed_Gaussian': False,
            'sharp_mass_shell_proved': False,
            'Lorentz_covariance_proved': False,
            'surviving_connected_interaction_proved': False,
            'physical_species_or_particle_masses_identified': False,
        },
        'inherited_custody': {
            'kind': 'Current narrow mathematical source bindings; analytic proof in README',
            'prior_full_producers_rerun_here': False,
            'finite_checks_alone_prove_infinite_volume_result': False,
        },
        'primary_import': {
            'url': 'https://arxiv.org/html/2105.04970v2',
            'statement': 'Eq. (3.4): transverse inverse-energy bound 1/(2J e(k+Q)); finite positive-field real ground',
        },
    }
    if bindings() != before:
        raise RuntimeError('source changed during replay')
    return result


def validate(value):
    if canonical(value) != canonical(build()):
        raise ValueError('source, arithmetic or scientific scope changed')
