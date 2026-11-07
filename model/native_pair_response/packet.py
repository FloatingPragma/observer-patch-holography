"""Replay same-occupation native spatial responses and quadratic fields."""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from math import comb
from pathlib import Path

from tools import strict_json
from model.native_correlation_sectors import packet as correlation_packet
from model.native_sector_readout import packet as tree_packet
from model.observer_time_fields import continuum_dynamics
from . import response as r,records,smearing as w

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCE_FILES=tuple(sorted(set((*continuum_dynamics.SOURCE_PATHS,
    'model/observer_time_fields/REFINEMENT.md','model/native_local_response/PROOF.md',
    'model/native_local_response/bounds.py','model/native_spatial_scaling/PROOF.md',
    'tools/strict_json.py','model/requirements.txt',
    *(f'model/native_pair_response/{name}' for name in (
        '__init__.py','response.py','records.py','smearing.py','packet.py',
        'test_response.py','test_smearing.py','test_records_packet.py','PROOF.md','README.md'))))))


def source_bindings():
    paths=set(SOURCE_FILES)|set(correlation_packet.source_bindings())|set(tree_packet.source_bindings())
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths) if not p.startswith('docs/')}


def rational_tree(value):
    if type(value) is F:return str(value)
    if type(value) in (tuple,list):return [rational_tree(x) for x in value]
    if type(value) is dict:return {k:rational_tree(v) for k,v in value.items()}
    return value


def finite_controls():
    family=r.s.source.Family((1,)*6,((0,1),(1,2),(2,3),(3,4),(4,5)),(1,)*5)
    probes=(((0,1),(2,3)),((1,2),(4,5)))
    matrix=r.response_matrix(family,probes,3);rows=[]
    for q in (2,3,4):
        variance=4*r.s.norm_factor(6,2,q)
        for T in (0,1,2):
            for i in range(2):
                for j in range(2):
                    actual=records.direct_cross(family,q,T,probes[j],probes[i])
                    expected=variance*matrix[T][i][j]
                    if actual!=expected:raise ArithmeticError('complete native cross records differ from response matrix')
                    rows.append({'q':q,'attempts':T,'future':i,'past':j,'complete_score':str(actual)})
    path=r.s.source.Family((1,)*4,((0,1),(1,2),(2,3)),(1,)*3);edgepair=((0,1),(2,3))
    native=r.response_matrix(path,(edgepair,),1)[1][0][0]
    killed=r.response_matrix(path,(edgepair,),1,killed=True)[1][0][0]
    if native-killed!=F(1,12):raise ArithmeticError('finite exclusion interaction control changed')
    return {'all_occupation_cross_controls':rows,
            'path4_one_attempt':{'native':str(native),'killed_independent':str(killed),
                                'difference':str(native-killed),'collision_retained':True}}


def smearings():
    f=(F(2),F(1),F(),F(-1),F(-1),F(-1));g=(F(1),F(),F(-2),F(1),F(),F())
    h=w.projection(f,g);rows=[]
    for q in (2,3,4):
        chi=F(q*(6-q),6*5);states=tuple(combinations(range(6),q));values=[]
        for x in states:
            field=w.field(f,g,q,x);included=sum(v for pair,v in h.items() if set(pair)<=set(x))/(6*chi)
            if field!=included:raise ArithmeticError('exact primitive quadratic subtraction failed')
            values.append(field)
        if sum(values):raise ArithmeticError('primitive field has nonzero canonical mean')
        for i in range(6):
            if sum(v for x,v in zip(states,values,strict=True) if i in x):
                raise ArithmeticError('primitive field retained a one-body component')
        fourth=sum(sum(f[i] for i in x)**4 for x in states)/(len(states)*(6*chi)**2)
        variance=sum(x*x for x in f)/6
        rows.append({'q':q,'configurations':len(states),**w.correction_bound(f,g,q),
                     'finite_density_fourth_moment':fourth,'Gaussian_fourth_moment':3*variance**2})
    return rational_tree({'six_port_projection_rows':rows,
        'uniform_eighth_moment_at_mean_one_third':2*3**4*24,
        'bounded_profile_subGaussian_variance_upper':'3*C^2/2',
        'moment_convergence_uses_uniform_integrability':True,
        'continuum_claim':'Wick quadratic of the inherited stationary diffusive Gaussian field',
        'microscopic_seam_probe_given_same_continuum_limit':False,
        'unbounded_L2_or_pointwise_product_claimed':False})


def operational_access():
    family,probes=r.federation_menu(0);rows=[]
    for q in (2,80):
        value=records.snapshot_program(0,q,2)
        # The decoder receives only records; these are independent known-input controls.
        initial=sum(1<<(i*(240//q)) for i in range(q))
        coins=(0,1)+(0,)*239
        result=records.trace_snapshots(value,initial,((0,1),(0,1)),coins)
        if result['decoded']['body_masks'][-1]!=initial:
            raise ArithmeticError('control body should move one unit away and return it')
        rows.append({'q':q,'program_digest':result['program_digest'],'body_attempts':2,
                     'detector_attempts':239,'complete_record_count':len(result['records']),
                     'exported_records':len(result['archive']),'FIFO_slots':len(result['memory']),
                     'word_probability':str(result['weight']),
                     'original_kernel_and_record_only_snapshots_agree':True,
                     'query_status':result['status']})
    cross=records.program(0,80,2,probes[0],probes[2])
    stopped=records.source.execute(cross,{(sum(1<<i for i in range(80)),0,sum(1<<i for i in range(80)),0):F(1)},max_branches=1)
    if stopped.status!='RESOURCE_STOP' or stopped.attempted_steps!=0:
        raise ArithmeticError('expanded host-stop control changed')
    return {'snapshot_controls':rows,'cross_program_digest':records.source.digest(cross),
            'cross_attempts':cross['declared_resources']['attempts_upper'],
            'expanded_control_status':stopped.status,'expanded_control_attempts':stopped.attempted_steps,
            'complete_pair_slots_reused':625,'initial_census_supplied_to_decoder':False,
            'global_snapshot_tape_given_local_cone_error':False}


def produce():
    before=source_bindings();levels=[]
    for level in (0,1,2):
        family,probes=r.federation_menu(level);native=r.response_matrix(family,probes,4)
        killed=r.response_matrix(family,probes,4,killed=True);errors=[]
        for j,past in enumerate(probes):
            error=r.comparison_error(family,past,probes[0],4)
            if abs(native[4][0][j]-killed[4][0][j])>error['normalized_cross_error_upper']:
                raise ArithmeticError('signed cross-response exceeded complete collision bound')
            errors.append(error)
        q=len(family.weights)//3
        levels.append(rational_tree({'level':level,'ports':len(family.weights),'q':q,
            'seams':len(family.seams),'settings':probes,'variance':4*r.s.norm_factor(len(family.weights),2,q),
            'native_matrix_by_attempt':native,'killed_matrix_by_attempt':killed,
            'collision_error_at_four_attempts':errors,
            'local_boundary_at_four_attempts':r.local_error(family,q,probes[0],probes[2],4,1)}))
    result={'schema':'oph.native_pair_response.v1','source_bindings':before,
            'spatial_matrices':levels,'independent_controls':finite_controls(),
            'smooth_primitive_fields':smearings(),'native_record_access':operational_access(),
            'scope':{'same_working_q_reference_and_native_attempt_clock':True,
                     'auxiliary_killed_comparison_is_native_repair_law':False,
                     'collision_survivors_postselected':False,'source_caps_changed':False,
                     'physical_clock_Lorentzian_QFT_or_particle_masses_derived':False,
                     'Gaussian_collective_limit_is_interacting_continuum':False}}
    if source_bindings()!=before:raise ArithmeticError('sources changed during pair response replay')
    return result


def canonical(value):return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path)))!=canonical(produce()):
        raise ValueError('native pair-response receipt differs from complete source replay')
    return True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:OUTPUT.write_text(canonical(produce()));print('Native pair-response receipt written.')
    else:verify();print('Native pair-response receipt verified.')
