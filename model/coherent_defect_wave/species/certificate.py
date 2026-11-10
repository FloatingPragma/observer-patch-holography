"""Exact occupation extension and finite species-selection boundaries."""
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import combinations, product
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from ..model import SparseState
from .fock import FockObserver, FockState, coordinates, mode
from . import spatial_fock
from . import joint_readout

ROOT = Path(__file__).resolve().parents[3]


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def lift_one(state):
    return FockState({(g, (mode(p,j,r,a),)): value
                      for (g,p,j,r,a),value in state.amplitudes}, state.norm_denominator)


def replace_internal(state, target):
    return SparseState({(g,p,target,r,a): value
                        for (g,p,j,r,a),value in state.amplitudes},state.norm_denominator)


def create(state, index):
    result = {}
    for occupied,value in state.items():
        if index not in occupied:
            output = tuple(sorted((*occupied,index)))
            sign = (-1)**sum(i < index for i in occupied)
            result[output] = result.get(output,0)+sign*value
    return {k:v for k,v in result.items() if v}


def annihilate(state, index):
    result = {}
    for occupied,value in state.items():
        if index in occupied:
            output = tuple(i for i in occupied if i != index)
            sign = (-1)**sum(i < index for i in occupied)
            result[output] = result.get(output,0)+sign*value
    return {k:v for k,v in result.items() if v}


def add(left,right):
    result = dict(left)
    for key,value in right.items():
        result[key] = result.get(key,0)+value
    return {k:v for k,v in result.items() if v}


def algebra_certificate():
    # All occupations of four chosen modes, with all 16 creation/annihilation pairs.
    # The proof supplies the same identity with an arbitrary spectator set.
    states = [tuple(i for i in range(4) if mask >> i & 1) for mask in range(16)]
    for occupied,i,j in product(states,range(4),range(4)):
        state = {occupied:1}
        require(add(annihilate(create(state,j),i),create(annihilate(state,i),j))
                == (state if i == j else {}), 'mixed CAR identity')
        require(add(create(create(state,j),i),create(create(state,i),j)) == {},
                'creation CAR identity')
    def pair(state,i,j):
        return create(create(state,j),i)
    require(pair(pair({():1},2,3),0,1) == pair(pair({():1},0,1),2,3),
            'disjoint even auxiliary composites commute')
    return {'tested_fock_modes':4,'tested_occupation_basis_states':16,
            'mixed_CAR_checks':256,'creation_CAR_checks':256,
            'full_7200_mode_CAR':'exterior-algebra proof; executable checks are finite witnesses',
            'even_auxiliary_composites_exchange_sign':1,
            'physical_fermion_statistics':'declared exterior wave architecture',
            'spin_statistics_theorem_instantiated':False}


def representation_certificate():
    rows=[]
    for degree in (2,4):
        for occupied in combinations(range(5),degree):
            n3=sum(i < 3 for i in occupied); n2=degree-n3
            color = tuple(sum(((1,0),(-1,1),(0,-1))[i][a]
                              for i in occupied if i < 3) for a in range(2))
            twice_t3 = int(3 in occupied)-int(4 in occupied)
            q = -2*n3+3*n2
            rows.append({'id':'aux-even-'+''.join(map(str,occupied)),
                         'auxiliary_occupation':list(occupied),'n3':n3,'n2':n2,
                         'color_cartan':list(color),'twice_T3':twice_t3,
                         'six_Y':q,'electric_Q':str(F(twice_t3,2)+F(q,6))})
    grouped=Counter((row['n3'],row['n2'],row['six_Y']) for row in rows)
    expected={(2,0,-4):3,(1,1,1):6,(0,2,6):1,(3,1,-3):2,(2,2,2):3}
    require(grouped == expected,'exact 15-state auxiliary exterior weights')
    require(sum(row['six_Y'] for row in rows)==0,'hypercharge trace')
    require(sum(row['six_Y']**3 for row in rows)==0,'cubic hypercharge trace')
    return {'declared_domain':'Lambda^2 C^5 plus Lambda^4 C^5, auxiliary internal register',
            'architecture_relation':'compatible representation control; no source-to-species identification',
            'rows':rows,'dimension':15,'trace_six_Y':0,'trace_cube_six_Y':0,
            'multiplets':[
                {'SU3':'3bar','SU2':'1','Y':'-2/3','dimension':3,'comparison':'up singlet conjugate'},
                {'SU3':'3','SU2':'2','Y':'1/6','dimension':6,'comparison':'quark doublet'},
                {'SU3':'1','SU2':'1','Y':'1','dimension':1,'comparison':'charged-lepton singlet conjugate'},
                {'SU3':'1','SU2':'2','Y':'-1/2','dimension':2,'comparison':'lepton doublet'},
                {'SU3':'3bar','SU2':'1','Y':'1/3','dimension':3,'comparison':'down singlet conjugate'}],
            'physical_spin_supplied':False,'three_generations_supplied':False,
            'field_dynamics_supplied':False}


def same_source_certificate(observer):
    checks=0
    for g,p,j,r,a in product((1,3), (0,5,6),range(5),(0,1),(0,1,59)):
        state=observer.base.basis_state(g,p,j,r,a)
        require(observer.step(lift_one(state))==lift_one(observer.base.step(state)),
                'one-occupation restriction equals the literal joint clock')
        checks+=1
    # Matrix units on the internal triplet and doublet intertwine the full clock.
    copies=0
    for g,p,r,a in product((1,3),(0,5,6),(0,1),(0,59)):
        for block in ((0,1,2),(3,4)):
            state=observer.base.basis_state(g,p,block[0],r,a)
            out=observer.base.step(state)
            for j in block:
                require(observer.base.step(replace_internal(state,j))==replace_internal(out,j),
                        'exact internal spectator intertwiner')
                copies+=1
    triplet=observer.base.basis_state(1,0,0,0,0)
    left=replace_internal(observer.base.step(triplet),3).as_dict()
    right=observer.base.step(replace_internal(triplet,3)).as_dict()
    cross_norm=sum((left.get(key,0)-right.get(key,0))**2 for key in left.keys()|right.keys())
    require(cross_norm==F(421,225),'cross-block matrix unit fails to commute')
    return {'one_occupation_executed_basis_witnesses':checks,
            'same_source_identity':'U_F restricted to one wave equals U by exterior degree-one identity',
            'internal_matrix_unit_witnesses':copies,
            'exact_decomposition':'I3 tensor U_triplet direct_sum I2 tensor U_doublet',
            'multiplicity_dimensions':[3,2], 'reduced_joint_dimension':14400,
            'reduced_positive_Gauss_dimension':7200,
            'pure_internal_unitary_commutant':'U(3) x U(2)',
            'pure_internal_tracefree_antihermitian_dimension':12,
            'cross_block_commutator_witness_norm_squared':str(cross_norm),
            'same_operator_internal_generation_mass_splitting_possible':False,
            'reason':'every spectral projector is identical on the spectator copies',
            'arbitrary_detector_can_break_spectator_symmetry':True}


def occupation_certificate(observer):
    indices=(mode(5,3,0,0),mode(5,4,0,0))
    raw=observer.basis_state(1,indices)
    paired=observer.gamma(raw)
    state=FockState({**raw.as_dict(),**paired.as_dict()},2)
    require(observer.gamma(state)==state,'normalized positive-Gauss preparation')
    output=observer.step(state)
    require(state.norm_squared()==output.norm_squared()==1,'exact two-wave normalization')
    require(observer.gamma(output)==output,'executed two-wave Gauss preservation')
    require(observer.controls(observer.controls(output),inverse=True)==output,
            'ordered control inverse')
    numbers={len(occupied) for (g,occupied),_ in output.amplitudes}
    charges={sum((2 if coordinates(i)[1]<3 else -3) for i in occupied)
             for (g,occupied),_ in output.amplitudes}
    require(numbers=={2} and charges=={-6},'fixed occupation and charge')
    require(all(observer.base.law.region_charge(g)==1 for (g,_),_ in output.amplitudes),
            'odd regional boundary flux')
    branches=observer.face_measurement(output)
    require(sum(branch.norm_squared() for _,branch in branches)==1,'complete five-face readout')
    require(all(observer.gamma(branch)==branch for _,branch in branches),'each face branch constrained')
    probabilities={str(face):str(branch.norm_squared()) for face,branch in branches}
    expected={'0':'7417/7500','2':'449/202500','5':'1/202500',
              '16':'112/50625','17':'1343/202500'}
    require(probabilities==expected,'fixed two-wave exact detector distribution')
    alternative=FockObserver('presence')
    other=alternative.step(state)
    other_prob={str(face):str(branch.norm_squared()) for face,branch in alternative.face_measurement(other)}
    require(other.norm_squared()==1 and alternative.gamma(other)==other,'alternative constrained law')
    require(other_prob!=probabilities,'one-wave agreement does not select the pair law')
    # Direct two-wave controlled move separates the laws without a large inversion.
    require(observer.controls(raw)==raw,'even endpoint count leaves parity control idle')
    require(alternative.controls(raw)!=raw,'presence control moves the same two-wave input')
    vacua=[]
    for g in observer.base.link_words:
        vacuum=observer.basis_state(g)
        require(observer.step(vacuum)==vacuum,'empty-wave invariant vector')
        vacua.append(g)
    pairs={tuple(sorted((g,g^observer.base.law.vertex_stars[5]))) for g in vacua}
    require(len(pairs)==5,'five positive Gauss matter-vacuum vectors')
    return {'full_Hilbert_dimension':'10 * 2^7200','wave_modes':7200,
            'shared_link_registers':1,'architecture_control':'odd endpoint occupation',
            'prepared_modes':[list(coordinates(i)) for i in indices],
            'prepared_link_words':[1,1^observer.base.law.vertex_stars[5]],
            'preparation':'normalized positive-Gauss projection of two-wave basis input',
            'output_support':len(output.amplitudes),'output_norm_squared':'1',
            'complete_face_probabilities':probabilities,
            'alternative_presence_complete_face_probabilities':other_prob,
            'Gauss_preserved':True,'complete_face_branches_Gauss_preserved':True,
            'wave_number':2,'total_Q5':-6,'regional_boundary_charge':1,
            'fixed_empty_wave_basis_vectors':10,'fixed_positive_Gauss_empty_wave_dimension':5,
            'empty_wave_is_zero_flux_physical_vacuum':False,
            'number_changing_or_pair_creation_terms':False,
            'relativistic_scattering_identification':False}


def input_hashes():
    own={str(p.relative_to(ROOT)) for p in Path(__file__).parent.glob('*.py')
         if not p.name.startswith('test_')}
    sources=own|{
        'model/requirements.txt','model/particle_common_contract/provenance.py',
        'evidence/native_geometric_source_20260925/spec.json',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
        'Lean/Screen/A2HolonomyBridge.lean','Lean/Screen/FermionSectorAssembly.lean',
        'Lean/Geometry/CanonicalLorentzModule.lean','Lean/Geometry/MassShellKinematics.lean',
        'Lean/QFT/InheritanceMatrix.lean'}
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(import_closure(sources))}


def certificate():
    observer=FockObserver()
    return {'program':'coherent-defect-wave-shared-link-CAR-v1',
            'architecture_inputs':['quantum Hilbert semantics','exterior wave Fock space',
                                   'odd endpoint occupation transport control'],
            'empirical_particle_inputs':[], 'CAR':algebra_certificate(),
            'auxiliary_weights':representation_certificate(),
            'same_source_and_internal_reduction':same_source_certificate(observer),
            'two_occupation_execution':occupation_certificate(observer),
            'spatial_Fock_join':spatial_fock.certificate(),
            'joint_translation_phase_readout':joint_readout.certificate(),
            'physical_scope':{'SM_species_identified':False,'physical_spin_derived':False,
                              'relativistic_dispersion_derived':False,'rest_masses_derived':False,
                              'dynamical_nonabelian_gauge_fields_derived':False},
            'input_sha256':input_hashes()}
