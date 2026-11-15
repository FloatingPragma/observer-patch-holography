"""Energy and complete spectral records of the one frozen local walk.

No Hamiltonian is supplied to source execution. A logarithm branch is an
explicit interpretation of the unit-circle records.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
WALK = ROOT / 'openwave-tuning/auditor/native_walk_bridge'
APPARATUS = ROOT / 'openwave-tuning/auditor/common_walk_particle'
FIELD = ROOT / 'openwave-tuning/auditor/native_current_walk'
sys.path[:0] = [str(WALK), str(APPARATUS), str(FIELD), str(HERE.parent / 'clock'),
               str(ROOT / 'reverse-engineering-reality')]
import common_walk as common
import compile_walk as compiled
import current_walk as field
import instrument_contract as apparatus
from native_phase_probe import actual_native_pointer
from clock_contract import effective_time
from model.observer.source_scalar import compiler
from model.certified_reference import certificate as cert

BASE_ID = common.MODEL_ID


def admitted_state(state):
    value = np.asarray(state, complex)
    if value.shape != (128, 16) or not np.isfinite(value).all():
        raise ValueError('finite full128-by16 state required')
    if abs(np.vdot(value, value).real - 1) > 1e-10:
        raise ValueError('normalized source input required; no normalization is applied')
    return value


def inverse_walk(value):
    built = compiled.controls()
    swapped = np.asarray(value).reshape(-1)[built['swap_permutation']].reshape(128, 16)
    return swapped @ built['reflection'].conj()


def positive_witness(value):
    """D(W)=I-(W+W*)/2, a conserved observable, not a time generator."""
    value = np.asarray(value)
    return value - (common.apply_walk(value) + inverse_walk(value)) / 2


def quasienergy_labels(bits, h, tau0=F(1)):
    if type(bits) is not int or not 1 <= bits <= 5:
        raise ValueError('one through five probe bits required')
    dt = effective_time(1, h, tau0)
    N = 1 << bits
    return [{'bin': k, 'angle_pi_multiple': str(F(2*k, N) if 2*k < N else F(2*k, N)-2),
             'energy_pi_multiple_per_clock_unit': str(-(F(2*k, N) if 2*k < N else F(2*k, N)-2)/dt)}
            for k in range(N)]


def qft_word(bits):
    """Inverse Fourier transform on little-endian phase bits, full phase."""
    if type(bits) is not int or not 1 <= bits <= 5:
        raise ValueError('one through five phase bits required')
    word, template = [], []
    for target in reversed(range(bits)):
        word += compiler.remap(compiled.controls()['primitive_words']['hadamard'], [target])
        template += [('I', (target,), cert._pi(-1, 2)), ('H', (target,), cert._pi(1, 2))]
        for control in range(target):
            multiple = F(1, 2**(target-control))
            local = compiler.compile_diagonal([0., 0., 0., math.pi*float(multiple)])
            word += compiler.remap(local, [control, target])
            for mask in range(4):
                label = ''.join('Z' if mask >> bit & 1 else 'I' for bit in range(2))
                local_template = cert.pauli(label, cert.PiMultiple(multiple*(-1)**mask.bit_count()/4))
                template += [(kind, tuple((control, target)[q] for q in targets), angle)
                             for kind, targets, angle in local_template]
    for a in range(bits//2):
        b = bits-1-a
        for c, t in ((a,b),(b,a),(a,b)):
            word += compiler.cnot(c,t)
            template += cert.cnot(c,t)
    return compiler.inverse(word), cert.inverse(template)


def phase_estimation(state, bits):
    """Complete Kraus outputs K_k=1/N sum_n exp(-2pi ikn/N) W^n.

    The returned vectors include every source/event/padding coordinate.
    Source pulse errors are bounded separately; no output is conditioned.
    """
    value = admitted_state(state)
    quasienergy_labels(bits, F(1))
    N = 1 << bits
    history = []
    for n in range(N):
        history.append(value.copy())
        if n+1 < N:
            value = common.apply_walk(value)
    output = np.fft.fft(np.stack(history), axis=0)/N
    probabilities = (abs(output)**2).sum(axis=(1,2))
    assert abs(probabilities.sum()-1) < 2e-10
    return output, probabilities


def qpe_contract(bits):
    N = 1 << bits
    word, template = qft_word(bits)
    end, fixed, variables = cert.consume(word, 0, template)
    assert end == len(word) and not variables
    qft_error = fixed + cert.controls_bound(word, {})
    source = json.loads((ROOT/'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json').read_text())
    probe = json.loads((APPARATUS/'INSTRUMENT_CONTRACT.json').read_text())
    source_error = F(source['whole_W_operator_error_exact'])
    prep_error = F(source['whole_source_preparation_error_exact'])
    wrapper = F(probe['controlled_W_wrapper_error_exact'])
    pointer_error = apparatus.pointer_bound()
    hword = compiled.controls()['primitive_words']['hadamard']
    _, hfixed, hv = cert.consume(hword, 0, [('I',(0,),cert._pi(-1,2)),('H',(0,),cert._pi(1,2))])
    assert not hv
    hadamard_error = hfixed+cert.controls_bound(hword,{})
    resets = apparatus.reset_contract(23+bits, bits)
    total = (prep_error+(N-1)*(source_error+wrapper)+qft_error
             +bits*(hadamard_error+pointer_error)+F(resets['timeout_upper_exact']))
    return {'bits':bits,'all_phase_bins':N,
            'controlled_W_executions':N-1,'native_pulses_excluding_source_preparation':
                (N-1)*probe['controlled_W_native_pulses']+len(word)+bits*(len(hword)+6),
            'inverse_QFT_word':word,'inverse_QFT_error_exact':str(qft_error),
            'total_variation_bound_exact':str(total),'total_variation_bound':math.nextafter(float(total),math.inf),
            'reset_contract':resets,
            'macro_program': {'phase_bits':list(range(23,23+bits)),
                'prepare':'native blank resets, common source preparation, native Hadamard on each phase bit',
                'controlled_powers':[{'control':23+j,'controlled_common_W_repetitions':1<<j} for j in range(bits)],
                'analysis':'remap the emitted inverse_QFT_word onto phase bits; read each with actual native C/E pointer',
                'fresh_pointer_auxiliaries_required':True},
            'scope':'Frozen base source only, with all native reset/export failures. Full field or Fock extension requires control of its entire augmented word and fresh disjoint wires.'}


@lru_cache(None)
def physical_spectral_data():
    n = 96
    A = np.kron(np.eye(n),np.ones((12,1))/np.sqrt(12))
    neigh,_ = common.geometry_tables()
    swap = np.arange(12*n)
    for v in range(n):
        for d in range(6):
            swap[12*v+2*d+1]=12*neigh[v,d]+2*d+1
    SA=A[swap]
    T=A.T@SA
    values,vectors=np.linalg.eigh(T)
    moving=values<1-1e-12
    theta=np.arccos(values[moving])
    Af=A@vectors[:,moving]
    SAf=SA@vectors[:,moving]
    plus=(Af-SAf*np.exp(1j*theta)[None,:])/np.sqrt(2*(1-values[moving]**2))[None,:]
    minus=plus.conj()
    symmetric=(A+SA)@vectors/np.sqrt(2*(1+values))[None,:]
    Pminus=(np.eye(12*n)+np.eye(12*n)[swap])/2-symmetric@symmetric.T
    H=(minus*theta[None,:])@minus.conj().T-(plus*theta[None,:])@plus.conj().T+math.pi*Pminus
    def act(x):return (2*A@(A.T@x)-x)[swap]
    def inverse(x):
        x=x[swap]
        return 2*A@(A.T@x)-x
    assert np.max(abs(H-H.conj().T))<2e-13
    assert np.max(abs(act(H)-inverse(H.conj().T).conj().T))<2e-12
    assert np.max(abs(H@plus+plus*theta[None,:]))<2e-12
    return A,swap,values,vectors,theta,plus,minus,Pminus,H


def continuity(state):
    state=admitted_state(state)
    Dstate=positive_witness(state)
    Wstate=common.apply_walk(state)
    WDstate=common.apply_walk(Dstate)
    before=(abs(state)**2).sum(axis=1)
    after=(abs(Wstate)**2).sum(axis=1)
    ebefore=np.real((state.conj()*Dstate).sum(axis=1))
    eafter=np.real((Wstate.conj()*WDstate).sum(axis=1))
    reflection=compiled.controls()['reflection']
    reflected=state@reflection.T
    refD=Dstate@reflection.T
    neighbours,edge_ids=common.geometry_tables()
    flux=np.zeros(288);eflux=np.zeros(288)
    B=np.zeros((128,288))
    for a in range(96):
        for d,b in enumerate(neighbours[a]):
            if a>=b:continue
            e=int(edge_ids[a,d]);pa,pb=common.encode_position(a),common.encode_position(int(b));slot=2*d+1
            B[pa,e]=1;B[pb,e]=-1
            flux[e]=abs(reflected[pa,slot])**2-abs(reflected[pb,slot])**2
            eflux[e]=np.real(reflected[pa,slot].conjugate()*refD[pa,slot]-reflected[pb,slot].conjugate()*refD[pb,slot])
    p_error=float(np.max(abs(after-before+B@flux)))
    e_error=float(np.max(abs(eafter-ebefore+B@eflux)))
    assert max(p_error,e_error)<2e-12
    assert abs(ebefore.sum()-eafter.sum())<2e-12
    return {'position_continuity_error':p_error,'D_density_continuity_error':e_error,
            'total_D_before':float(ebefore.sum()),'total_D_after':float(eafter.sum()),
            'D_density_minimum':float(ebefore.min()),'all128_position_outputs':after.tolist(),
            'all288_native_mean_current_records':flux.tolist(),
            'all288_D_flux_bilinears':eflux.tolist(),
            'scope':'D density and its bilinear flux are a conserved diagnostic. They are not an independently identified physical stress-energy tensor.'}


def run():
    files=[Path(__file__),WALK/'common_walk.py',WALK/'compile_walk.py',
           WALK/'COMMON_WALK_RECEIPT.json',APPARATUS/'MANIFEST.json',APPARATUS/'INSTRUMENT_CONTRACT.json',
           FIELD/'current_walk.py',FIELD/'MANIFEST.json',HERE.parent/'clock/clock_contract.py',
           HERE.parent/'clock/MANIFEST.json',ROOT/'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json',
           Path(cert.__file__),Path(compiler.__file__),APPARATUS/'native_phase_probe.py',
           HERE.parent/'manybody/MANIFEST.json',HERE.parent/'manybody/manybody.py',
           HERE.parent/'manybody/mode_field.py',HERE.parent/'manybody/AMENDMENT.json']
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    state=common.initial_cell_state()
    pointer_laws,_=actual_native_pointer()
    fine_pointer={}
    for law,branches in pointer_laws.items():
        fine_pointer[law]={}
        for bit in (0,1):
            rows=[]
            for metadata,kraus in branches:
                probability=sum(float(np.vdot(K[:,bit],K[:,bit]).real) for K in kraus)
                rows.append({'native_fine_record':metadata,'probability_given_bit':probability})
            assert abs(sum(row['probability_given_bit'] for row in rows)-1)<1e-12
            fine_pointer[law][str(bit)]=rows
    qpe=[]
    for bits in (2,3,4):
        N=1<<bits
        word,_=qft_word(bits)
        actual=compiler.execute(word,bits)
        target=np.exp(-2j*np.pi*np.outer(np.arange(N),np.arange(N))/N)/np.sqrt(N)
        error=float(np.linalg.norm(actual-target,2))
        assert error<1e-11
        out,probs=phase_estimation(state,bits)
        contract=qpe_contract(bits)
        contract['native_QFT_all_columns_operator_error']=error
        contract['all_phase_probabilities']=probs.tolist()
        contract['all_output_norm_squared']=float(probs.sum())
        contract['quasienergy_labels']=quasienergy_labels(bits,F(1,2))
        contract['Kraus_formula']='K_k=(1/N)sum_n exp(-2pi*i*k*n/N) W^n; outputs retain the full padded source state and any unchanged reference'
        contract['circular_phase_tail_bound']='Pr(distance(theta,2pi*k/N)>=delta)<=min(1,1/(N*sin(delta/2)^2)) for0<delta<=pi; no branch-unwrapped bound across pi is inferred'
        qpe.append(contract)
    A,swap,values,vectors,theta,plus,minus,Pminus,H=physical_spectral_data()
    seed=np.zeros(96);seed[:12]=1/np.sqrt(12)
    weights=abs(vectors.T@seed)**2
    positive_generator_return=np.sum(weights*np.exp(-1j*np.arccos(np.clip(values,-1,1))))
    actual_return=np.vdot(state,common.apply_walk(state))
    assert abs(positive_generator_return.imag)>.05 and abs(actual_return.imag)<1e-12
    # Local holding-slot difference is an exact -1 mode. Two such fermionic
    # modes have additive energy2pi but total Floquet phase1.
    flat=np.zeros((128,16),complex);flat[0,0]=1/np.sqrt(2);flat[0,2]=-1/np.sqrt(2)
    assert np.max(abs(common.apply_walk(flat)+flat))<1e-12
    second_flat=np.zeros((128,16),complex);second_flat[1,0]=1/np.sqrt(2);second_flat[1,2]=-1/np.sqrt(2)
    assert np.max(abs(common.apply_walk(second_flat)+second_flat))<1e-12
    orbital=np.column_stack((flat.reshape(-1),second_flat.reshape(-1)))
    evolved=np.column_stack((common.apply_walk(flat).reshape(-1),common.apply_walk(second_flat).reshape(-1)))
    fock_overlap=np.linalg.det(orbital.conj().T@evolved)
    assert abs(fock_overlap-1)<1e-12
    locality=[]
    g=common.graph()
    from collections import deque
    for origin in (0,1,17):
        distance={origin:0};todo=deque([origin])
        while todo:
            x=todo.popleft()
            for y in g.neighbours[x]:
                if y not in distance:distance[y]=distance[x]+1;todo.append(y)
        for target,d in distance.items():
            if d>2:
                block=H[12*origin:12*origin+12,12*target:12*target+12]
                locality.append((float(np.max(abs(block))),origin,target,d))
    far=max(locality)
    assert far[0]>1e-5
    histories=[]
    for _ in range(5):
        histories.append(continuity(state));state=common.apply_walk(state)
    rng=np.random.default_rng(1729)
    arbitrary=rng.normal(size=(128,16))+1j*rng.normal(size=(128,16));arbitrary/=np.linalg.norm(arbitrary)
    histories.append(continuity(arbitrary))
    # Same augmented field evolution; only the retained character is changed.
    initial=common.initial_cell_state()
    field_rows=[]
    for character in (0,1,2,3):
        joint=initial[:,:,None]*field.ideal_character(character)
        moment=np.vdot(joint,field.apply_current_walk(joint))
        observed=1-moment.real
        exact=F(1,12)+F(1,72)*(1-(1 if character==0 else F(1,2) if character==1 else F(-1,2) if character==2 else -1))
        assert abs(observed-float(exact))<1e-12
        field_rows.append({'character':character,'D_full_current_walk':float(observed),'exact':str(exact),
                           'one_step_X_probe_probabilities':[(1+moment.real)/2,(1-moment.real)/2]})
    # Finite-box gaps from existing development graphs, with no reserved run.
    gaps=[]
    from model.native_spatial_refinement import source as spatial
    for size in (2,3,4):
        eigen=[]
        for x in range(size):
            for y in range(size):
                for z in range(size):
                    eigen.extend(np.linalg.eigvalsh(spatial.symbol(2*np.pi*np.array([x,y,z])/size)))
        nonzero=np.array(eigen)[np.array(eigen)>1e-10]
        gap=float(np.arccos(1-nonzero.min()/12))
        assert gap<=math.pi**2/(6*size)+1e-12
        gaps.append({'size':size,'least_nonzero_phase':gap,'gap_is_not_rest_mass':True,
                     'analytic_infinite_volume_upper':math.pi**2/(6*size)})
    rejected=0
    for call in (lambda:phase_estimation(np.zeros((128,16)),2),
                 lambda:quasienergy_labels(3,F(0)),
                 lambda:quasienergy_labels(3,F(1),F(-1)),
                 lambda:quasienergy_labels(True,F(1)),
                 lambda:phase_estimation(np.full((128,16),np.nan),2)):
        try:call()
        except ValueError:rejected+=1
        else:raise AssertionError('malformed source or clock admitted')
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    result={'schema':'oph.common.local12.energy.v1','candidate_id':BASE_ID,
            'checks':{'same_W_and_clock_only':True,'complete_native_QPE_laws':True,
                      'all_QFT_columns_and_scalar_phases':True,'full_source_error_and_reset_composition':True,
                      'signed_log_conserved_and_branch_explicit':True,'positive_generator_substitution_rejected':True,
                      'local_native_probability_and_D_continuity':True,'local_log_assumption_rejected':True,
                      'field_intervention_uses_same_augmented_W':True,'finite_box_gap_not_mass':True,
                      'Fock_additive_principal_log_alias_control':True,'invalid_source_and_clock_inputs_rejected':rejected==5,
                      'source_pins_unchanged':True},
            'clock':'Delta=tau0*h, exactly the frozen completed-word clock; no physical seconds or hbar calibration',
            'phase_estimation':qpe,'continuity_records':histories,
            'native_fine_pointer_laws':fine_pointer,
            'QPE_full_fine_law':'For bin k with bits b_j and fixed law C or E, multiply its unnormalized source Kraus output by the tensor product of the displayed native pointer branches for b_j. Every zero record and first reset/export failure remains in the complete instrument.',
            'signed_log':{'convention':'W=exp(-i Delta H); Arg in[-pi,pi), so -1 assigned+pi/Delta',
                          'physical_negative_modes':95,'physical_zero_modes':194,
                          'physical_positive_modes':95+768,
                          'arbitrary_reference_ground_shift_is_not_source_selected':True,
                          'largest_log_entry_beyond_two_graph_edges':{'magnitude':far[0],'source':far[1],'target':far[2],'distance':far[3]}},
            'positive_witness':'D=I-(W+W*)/2>=0, commuting withW; its density/flux and native moments are derived, but D is not the time generator',
            'positive_generator_control':{'actual_W_return':[float(actual_return.real),float(actual_return.imag)],
                'exp_minus_i_abs_log_return':[float(positive_generator_return.real),float(positive_generator_return.imag)],
                'native_Y_probe_detects_difference':True},
            'field_intervention':{'rows':field_rows,'character1_minus0_D_exact':'1/144',
                'interpretation':'A declared field-character preparation intervention changes the conserved total-step witness. This is not an autonomous energy transfer or a physical force law.',
                'apparatus_requirement':'Control the full current-augmented W, including field shifts; bare controlledW and its error budget cannot be substituted. Full affine preparation retains its explicit integration obligation.'},
            'finite_box_gaps':gaps,
            'rest_gap_statement':'At k=0 the acoustic phase is exactly0. Fixed optical phases acos(2/3),acos(1/2),acos(1/3),pi scale as1/h underDelta=tau0*h; they do not produce a finite positive continuum rest mass without a declared parameter/refinement change.',
            'manybody_scope':{'fermionic_amendment':'native-local12-fermion-modes-v1-proposed',
                'vacuum':'Gamma(W)|empty>=|empty>; empty state is not the ground of signed dGamma(H) when negative modes are admitted',
                'fixed_charge_reference':'Closed torus requires sumG=Q. A filled branch changesQ/background; energy and Gauss references must be changed together.',
                'alias_control':'Two occupied distinct localized -1 modes: dGamma(H) energy2pi/Delta but Gamma(W) phase1 and principal many-body quasienergy0.',
                'actual_two_orbital_determinant_overlap':[float(fock_overlap.real),float(fock_overlap.imag)],
                'supplied_statistics':True,'no_private_H_in_source_execution':True},
            'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'ENERGY_CONTRACT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'QPE_bounds':[r['total_variation_bound'] for r in qpe],
                      'finite_box_gaps':gaps,'field_D_shift':'1/144','long_range_log':result['signed_log']['largest_log_entry_beyond_two_graph_edges']}))


if __name__=='__main__':run()
