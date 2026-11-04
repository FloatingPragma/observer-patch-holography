"""Source-bound rational proof packet for uniform moving two-body scattering."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import bounds, vacuum

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)


def serialize(value):
    if type(value) is F:return str(value)
    if type(value) is dict:return {k:serialize(v) for k,v in value.items()}
    if type(value) in (list,tuple):return [serialize(v) for v in value]
    return value


def bindings():
    names=integration_source_paths(ROOT,{
        'domain_and_proof':'model/moving_contact_scattering/PROOF.md',
        'source_receipts':['model/band_edge_contact/results.json','model/moving_bound_pair/results.json'],
        'independent_checks':['model/moving_contact_scattering/test_scattering.py',
                              'model/moving_contact_scattering/test_vacuum.py']})
    names.add('model/moving_contact_scattering/README.md')
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


_IMPORTED_BINDINGS=bindings()


def build():
    before=bindings()
    if before!=_IMPORTED_BINDINGS:raise ValueError('source changed after import; use a fresh process')
    value=serialize({'schema':'oph.moving_contact_scattering.v1','source_bindings':before,
        'coarea':bounds.coarea_certificate(),
        'additive_Fock_extension':vacuum.certificate(),
        'joined_finite_witness':bounds.joined_certificate(F(1,16384),1,1,1,1),
        'tolerance_schedules':[bounds.spacing_certificate(1,1,1,eps) for eps in (F(1,4),F(1,100),F(1,10000))],
        'scope':{'particle_number_sector':2,'spatial_dimensions':3,
                 'additive_generator_particle_number_sectors':'all nonnegative integers',
                 'fixed_N_self_adjointness_and_finite_spacing_stability':True,
                 'empty_Fock_ground_vacuum_excluded_for_additive_reflected_generator':True,
                 'pair_branch_exact_affine_Lorentz_shell_excluded':True,
                 'uniform_nonzero_total_momentum_on_shell_limit':True,
                 'uniform_threshold_endpoint_included':True,
                 'same_source_coupling_and_global_clock_for_scattering_and_binding':True,
                 'energy_reflection_and_subtraction_declared':True,
                 'all_eight_edge_valleys_retained':True,
                 'physical_species_identified':0,
                 'relativistic_many_body_QFT_constructed':False,
                 'positive_energy_many_body_vacuum_constructed':False,
                 'full_momentum_space_wave_operator_embedding_constructed':False,
                 'native_scattering_detector_executed':False,
                 'gamma_derived_from_source_without_family_choice':False,
                 'entire_completion_gate_accepted':False,
                 'historical_receipt_producers_all_replayed':False}})
    if before!=bindings():raise ValueError('source changed during proof replay')
    return value


def validate(value,replay=None):
    actual=build() if replay is None else replay
    if canonical(value)!=canonical(actual):raise ValueError('scattering proof, source or scope differs from replay')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    modes.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    args=parser.parse_args()
    result=build()
    if args.check:validate(strict_json.load(args.check),result)
    else:args.write.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print('PASS uniform moving two-body scattering and bound-band common-clock join')
