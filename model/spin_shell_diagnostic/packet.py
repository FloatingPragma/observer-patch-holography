"""A new diagnostic record; it does not promote the missing particle premise."""
import argparse
import hashlib
import json
from pathlib import Path
from .algebra import local_certificate,cubic_trace_certificate
from .criteria import ratio_moments,uniform_continuum_control
from .experiment import finite_control
from fractions import Fraction as F

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'
INPUTS=(
    'model/__init__.py',
    'model/spin_goldstone_sector/__init__.py',
    'model/protected_spin_limit/README.md',
    'model/protected_spin_limit/algebra.py',
    'model/spin_goldstone_sector/README.md',
    'model/spin_goldstone_sector/algebra.py',
    'model/spin_goldstone_sector/experiment.py',
    'model/spin_particle_scaling/README.md',
    'model/spin_particle_scaling/bounds.py',
    'model/requirements.txt',
    'tools/strict_json.py',
)


def bindings():
    names=set(INPUTS)
    names.update(p.relative_to(ROOT).as_posix() for p in HERE.iterdir()
                 if p.suffix in ('.py','.md'))
    return {n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in sorted(names)}


def canonical(value):return json.dumps(value,sort_keys=True,allow_nan=False,separators=(',',':'))


def report():
    before=bindings()
    two=ratio_moments([F(1),F(2)],[F(1,2),F(1,2)])
    value={'schema':'oph.spin_shell_diagnostic.v1',
            'scope':'EXACT_SOURCE_IDENTITIES_AND_OPEN_MASS_SHELL_DIAGNOSTIC',
            'mass_shell_established':False,'sharp_lattice_pole_established':False,
            'whole_response_concentration_necessary_for_any_sharp_branch':False,
            'physical_species_identified':False,'new_finite_attachment_claimed':False,
            'new_checks_pass':True,'source_bindings':before,
            'inherited_scope':{'old_full_receipt_closures_replayed':False,
                               'direct_mathematical_inputs_pinned':True,
                               'source_or_repair_law_changed':False},
            'exact_local_algebra':local_certificate(),
            'exact_tracial_force':cubic_trace_certificate(),
            'broad_band_criterion':{'energy_lower_per_Jp':'1/16000',
                                    'energy_upper_per_Jp':'64000000',
                                    'ratio_upper_rational':'1408000000/7',
                                    'lattice_error_coefficient':'c/12',
                                    'test_scope':'entire selected response, or an independently justified component with proved nonvanishing weight',
                                    'required_common_ratio_mean':'one unknown c>0',
                                    'required_ratio_defect_limit':'0',
                                    'speed_identified':False},
            'two_ray_ratio_control':{k:str(v) for k,v in two.items()},
            'continuous_ward_compatible_control':uniform_continuum_control(),
            'finite_source_control':finite_control(),
            'remaining_estimate':'A sufficient route is common-speed concentration of the entire broad-window response, or of an independently justified spectral component with proved nonvanishing weight. A sharp component may coexist with a continuum; whole-response concentration is not necessary for its existence.',
            'primary_source':'https://arxiv.org/html/2105.04970v2',
            'primary_theorem_scope':'Theorem2.1 bounds the filtered mean energy between positive multiples of momentum; it does not establish a vanishing spectral width.'}
    if bindings()!=before:
        raise ArithmeticError('source changed during spin-shell computation')
    return value


def validate(value,expected=None):
    expected=report() if expected is None else expected
    if canonical(value)!=canonical(expected):raise ValueError('spin shell diagnostic does not replay')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    args=parser.parse_args();value=report()
    if args.write:args.write.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')
    elif args.check:
        from tools import strict_json
        validate(strict_json.load(args.check),value)
    else:print(json.dumps(value,sort_keys=True,indent=2,allow_nan=False))

if __name__=='__main__':main()
