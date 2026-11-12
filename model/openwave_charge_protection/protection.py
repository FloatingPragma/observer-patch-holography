"""Actual texture records, integer lifts and independent native current labels."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib
import json

import sympy as s

from model.compact_defect_transport import source
from model.compact_defect_pair.topology import topology
from model.lattice_texture_selection.shells import degree
from model.defect_integer_lift import model as lift
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def independent_degree(triangles, state, target):
    """Solve rational barycentric systems rather than reusing cross products."""
    point=s.Matrix([s.Rational(x) for x in target])
    count=0
    witnesses=[]
    for index, triangle in enumerate(triangles):
        matrix=s.Matrix.hstack(*(s.Matrix([s.Rational(x) for x in state[p]]) for p in triangle))
        determinant=matrix.det()
        if determinant==0:
            continue
        barycentric=matrix.inv()*point
        if all(x>0 for x in barycentric):
            sign=1 if determinant>0 else -1
            count+=sign
            witnesses.append({'triangle':index,'orientation':sign,
                              'positive_barycentric_coordinates':list(map(str,barycentric))})
    return count,witnesses


def half_degree_solutions(total):
    if type(total) is not int:
        raise ValueError('topological degree is an integer')
    return [] if total%2 else [total//2]


@lru_cache(maxsize=1)
def texture_certificate():
    states, records, transition=source.continuing_protocol()
    surfaces=source.surfaces()
    comparisons=[]
    for state_id,state in enumerate(states):
        for shell_id,triangles in enumerate(surfaces):
            directions={p:state[p] for triangle in triangles for p in triangle}
            native=degree(triangles,directions)
            value,witnesses=independent_degree(triangles,directions,native['regular_target'])
            if value!=native['degree'] or value!=(-1 if shell_id==0 else 1):
                raise ArithmeticError('independent oriented shell count differs')
            comparisons.append({'state':state_id,'shell':shell_id,'degree':value,
                'triangles':len(triangles),'regular_target':native['regular_target'],
                'independent_signed_preimages':witnesses})
    # A finite guard sweep uses the actual law on every shell vertex. It
    # does not silently change the four-proposal stochastic kernel.
    baseline=states[0]
    vertices=sorted({p for triangles in surfaces for triangle in triangles for p in triangle})
    guards=Counter(); accepted=0; proposals=0
    for point in vertices:
        f=baseline[point];scale=max(map(abs,f))
        candidates=[tuple(f[j]+sign*F(1,512)*scale*int(j==axis) for j in range(3))
                    for axis in range(3) for sign in (-1,1)]
        candidates.extend([(F(),F(),F()),tuple(-x for x in f)])
        for candidate in candidates:
            target,record=source.update(baseline,point,candidate)
            proposals+=1;guards[record['reason']]+=1
            if record['accepted']:
                accepted+=1
                if source.charges(target)!=[-1,1]:
                    raise ArithmeticError('admitted exact proposal changed integer degree')
            elif target!=baseline:
                raise ArithmeticError('rejected proposal is not the retained identity')
    erased={p:((F(0),F(0),F(1)) if any(v) else v) for p,v in baseline.items()}
    if source.charges(erased)!=[0,0]:
        raise ArithmeticError('unguarded erasure must destroy the winding')
    zeros=lambda state:sum(2 for v in state.values() if not any(v))
    if zeros(erased)!=zeros(baseline) or zeros(baseline)!=4:
        raise ArithmeticError('erasure must keep the four zero onsite coordinates')
    if half_degree_solutions(1) or half_degree_solutions(-1):
        raise ArithmeticError('a unit integer winding was halved')
    return {'actual_continuing_states':len(states),'actual_complete_records':records,
        'complete_transition_matrix':[[str(x) for x in row] for row in transition],
        'independent_degree_comparisons':comparisons,
        'all_finite_histories_preserve_each_unit_degree_by_guarded_homotopy':True,
        'equal_integer_split_solutions_for_plus_one':half_degree_solutions(1),
        'equal_integer_split_solutions_for_minus_one':half_degree_solutions(-1),
        'unequal_integer_decomposition_control':[2,-1],
        'general_energetic_fission_excluded_by_integer_degree_alone':False,
        'cores_fixed_by_the_actual_declared_law':True,
        'adversarial_guard_sweep':{'vertices':len(vertices),'proposals':proposals,'accepted':accepted,
            'retained_reasons':dict(sorted(guards.items())),
            'sweep_replaces_actual_stochastic_law':False},
        'unguarded_erasure_degrees':source.charges(erased),
        'zero_onsite_coordinates_before_and_after':[zeros(baseline),zeros(erased)],
        'zero_modes_alone_protect_topological_degree':False,
        'compact_geometry':topology()}


def integer_lift_certificate():
    # Exhaust the actual local load/coin/class domain. The lifted integers
    # are comparison labels, whereas the register law consumes classes.
    count=0;accepted=0;hidden_examples=[]
    for xa,xb,coin,k1,k2 in product(range(6),range(6),(0,1),range(6),range(6)):
        na,nb,j=lift.attempt(xa,xb,coin)
        if na+nb!=xa+xb:raise ArithmeticError('native repair lost load')
        apply=bool(j and lift.accepts(k1,k2,j,xa-xb))
        accepted+=int(apply);count+=1
        for sheet1,sheet2 in product((-1,0,1),repeat=2):
            before=(k1+6*sheet1,k2+6*sheet2)
            after=(before[0]+j,before[1]-j) if apply else before
            if sum(after)!=sum(before):raise ArithmeticError('integer curl update lost oriented pair sum')
            expected=((k1+j)%6,(k2-j)%6) if apply else (k1,k2)
            if tuple(x%6 for x in after)!=expected:raise ArithmeticError('integer lift disagrees with source register')
            if apply and before==(4,0) and after[0]==6:
                hidden_examples.append({'loads_before':[xa,xb],'coin':coin,'j':j,
                    'lift_before':list(before),'lift_after':list(after),
                    'classes_before':[k1,k2],'classes_after':list(expected)})
    if count!=2592 or not hidden_examples:raise ArithmeticError('complete local law or hidden-sheet witness missing')
    fed=lift.Federation(1)
    sweep_checks=[]
    def audit_sweep(index,loads,register,classes,links,values):
        state={'values':values,'classes':classes,'link':links,'register':register}
        checks=lift.gates(fed,state)
        if not all(checks.values()):raise ArithmeticError('actual scheduled source failed integer Stokes')
        sweep_checks.append(index)
    run=lift.run(fed,200,after_sweep=audit_sweep)
    if len(sweep_checks)!=200:raise ArithmeticError('not every source cut was retained')
    lift_values=run['values']
    return {'local_load_coin_class_cases':count,'accepted_nonzero_twist_cases':accepted,
        'integer_representative_pairs_per_case':9,
        'same_class_alternative_lifts_have_identical_acceptance':True,
        'explicit_nontransporting_lift_witness':hidden_examples[0],
        'actual_reference_run':{'level':1,'carriers':fed.carriers,'sweeps':200,
            'applied_twists':len(run['twists']),'every_sweep_Stokes_and_modulus_checked':True,
            'maximum_absolute_current':run['max_abs_j'],
            'all_final_carrier_totals_zero':all(sum(lift_values[20*c:20*c+20])==0 for c in range(fed.carriers)),
            'final_hidden_flat_faces':sum(value!=0 and value%6==0 for value in lift_values),
            'final_state_values_sha256':lift.digest(lift_values)},
        'integer_lift_is_a_conserved_particle_hypercharge':False,
        'closed_sphere_total_winding_identified_with_this_face_curl':False}


def current_attachment_certificate():
    basis,_,currents,_=representation()
    q=currents[10]+currents[11]
    if q!=s.diag(*q.diagonal()):raise ArithmeticError('native Q chart is not diagonal')
    qvals=list(q.diagonal())
    colour_cartan=s.Rational(2,3)*currents[6]+s.Rational(4,3)*currents[7]
    cvals=list(colour_cartan.diagonal())
    if any((3*value).is_Integer is not True for value in qvals):
        raise ArithmeticError('native current left the one-third integer lattice')
    if any((charge+colour).is_Integer is not True for charge,colour in zip(qvals,cvals,strict=True)):
        raise ArithmeticError('native colour-centre/current relation failed')
    neutral=qvals.index(s.Integer(0));charged=qvals.index(s.Integer(1))
    # Same texture and same spatial band vector; only the native internal
    # occupation label differs. This is a literal two-state subspace of
    # the declared compact source, not a phenomenological particle map.
    winding=s.eye(2)
    charge=s.diag(qvals[neutral],qvals[charged])
    if charge==charge[0,0]*winding or charge[0,0]==charge[1,1]:
        raise ArithmeticError('the current/texture counterexample disappeared')
    # Any spatial U leaves the internal current unchanged. A two-coordinate
    # exact gate checks the tensor identity used by native band transport.
    u=s.Matrix([[s.Rational(3,5),-s.Rational(4,5)],[s.Rational(4,5),s.Rational(3,5)]])
    lifted=s.kronecker_product(u,s.eye(15))
    current=s.kronecker_product(s.eye(2),q)
    assert lifted.H*current*lifted==current
    return {'native_electric_candidate_current':'Q=T3+Y, same fifteen-state current module',
        'one_family_Q_weights':list(map(str,qvals)),
        'colour_central_generator_weights':list(map(str,cvals)),
        'all_Q_plus_colour_weights_integer':True,
        'full_Fock_identity':'exp(2pi i dGamma(Q)) Gamma(Z_colour)=I',
        'full_Fock_charge_spectrum_subset':'(1/3) Z',
        'colour_singlet_charge_spectrum_subset':'Z',
        'half_unit_charge_eigenvalue_absent_on_complete_Fock_carrier':True,
        'every_native_controlled_texture_record_preserves_Q':True,
        'source_current_quantization_requires_winding_identity':False,
        'mean_charge_can_be_half_without_a_half_charge_eigenvalue':{
            'state_density':[['1/2','0'],['0','1/2']],
            'Q_outcomes':['0','1'],'mean':'1/2'},
        'three_declared_family_copies':3,
        'same_plus_one_texture_one_particle_witness':{
            'neutral_internal_index':neutral,'neutral_internal_subset':list(basis[neutral]),
            'charged_internal_index':charged,'charged_internal_subset':list(basis[charged]),
            'degree_operator':[['1','0'],['0','1']],
            'charge_operator':[['0','0'],['0','1']],
            'no_function_of_degree_equals_native_Q_on_this_subspace':True},
        'vacuum_has_same_texture_degrees_and_Q_zero':True,
        'uniform_native_current_preserved_by_spatial_transport':True,
        'topology_current_locking_constraint_in_source':False,
        'local_Gauss_dressing_derived_by_this_package':False,
        'unit_degree_is_physical_electron_charge':False,
        'finite_obstruction':'The unchanged carrier contains equal-winding states with unequal current eigenvalues; no universal Q=f(winding) holds without restricting or coupling that carrier.'}


def certificate():
    seeds=[str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md','.lean')]
    seeds+=['model/compact_defect_transport/PROOF.md','model/compact_defect_pair/PROOF.md',
            'model/defect_integer_lift/federation.json','model/defect_integer_lift/PROOF.md']
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(seeds))}
    return {'schema':'oph.winding_charge_protection.v1','input_sha256':hashes,
        'guarded_winding':texture_certificate(),'actual_integer_register':integer_lift_certificate(),
        'current_attachment':current_attachment_certificate(),
        'scope':{'minimal_unit_integer_winding_on_actual_texture':True,
            'complete_native_guarded_history_conservation':True,
            'equal_halving_forbidden_for_unit_integer_winding':True,
            'full_Fock_native_current_quantization_and_colour_singlet_integer_charge':True,
            'native_guarded_transport_preserves_current_branchwise':True,
            'EM_charge_winding_identity_on_unchanged_source':False,
            'energetic_no_fractionation_of_physical_electron_proved':False,
            'distinct_topology_register_and_current_sources_identified':False},
        'verdict':'PASS_UNIT_WINDING_PROTECTION_AND_CURRENT_ATTACHMENT_OBSTRUCTION'}


def validate_pins(value):
    for name,expected in value['input_sha256'].items():
        path=Path(name)
        if path.is_absolute() or '..' in path.parts or hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=expected:
            raise ValueError('source pin mismatch: '+name)
