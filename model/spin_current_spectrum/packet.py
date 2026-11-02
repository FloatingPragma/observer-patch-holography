"""Narrow reproducible proof and exact-arithmetic packet."""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
from . import algebra, bounds

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
PARENTS=('tools/strict_json.py','model/__init__.py','model/spin_canonical_energy/README.md',
         'model/spin_canonical_field/README.md','model/spin_goldstone_sector/README.md')
FRAGMENT='paper/tex_fragments/SPIN_CURRENT_SPECTRUM.tex'
OUTPUT=HERE/'results.json'


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)


def bindings():
    names=list(PARENTS)+[str(p.relative_to(ROOT)) for p in sorted(HERE.glob('*.py'))]
    names += [str((HERE/'README.md').relative_to(ROOT)),FRAGMENT]
    return {n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}


def build():
    before=bindings()
    result={'schema':'oph.spin_current_spectrum.v1','source_bindings':before,
            'exact_finite_checks':algebra.exact_checks(),
            'soft_current_fractions':bounds.serialize([bounds.soft_fraction(F(1,1000),F(1,2),r) for r in (1,F(1,2))]),
            'fractional_tail_controls':bounds.serialize([bounds.fractional_tail(F(1,4),d,s) for d in (2,3,4) for s in (2,16,256)]),
            'energy_escape_controls':bounds.serialize([bounds.escaping_current(n) for n in (3,16,1024)]),
            'theorem':{'source':'unchanged cubic spin-1/2 antiferromagnet; each bond counted once',
                       'reference':'same infinitesimal-axis-1-staggered-field-selected zero-field ground',
                       'clock':'L_p=K/(Jp)',
                       'exact_current_continuity_and_field_torque':True,
                       'exact_current_inverse_moment':'J*kappa/2',
                       'nonzero_soft_current_fraction':True,
                       'fractional_energy_uniform_integrability_for_0_s_1':True,
                       'first_energy_uniform_integrability_proved':False,
                       'actual_ground_hard_current_tail_vanishes_proved':False,
                       'bare_twist_is_relaxed_stiffness':False,
                       'sharp_mass_shell_proved':False,
                       'entire_response_concentration_necessary_for_any_sharp_subbranch':False,
                       'physical_species_or_particle_masses_identified':False},
            'inherited_custody':{'kind':'current narrow mathematical source bindings; complete derivation in README',
                                 'all_prior_producers_rerun_here':False,
                                 'finite_controls_prove_infinite_ordered_state':False},
            'primary_sources':{'current_sum_rule':'https://arxiv.org/pdf/cond-mat/9311020',
                               'inherited_susceptibility':'https://arxiv.org/html/2105.04970v2'}}
    if bindings()!=before:raise RuntimeError('source changed during replay')
    return result


def validate(value):
    if canonical(value)!=canonical(build()):raise ValueError('arithmetic, source or scientific scope changed')
