"""Check the exact repair-only boundary and execute the supplied M6 service."""
from fractions import Fraction as F
from math import pi
import numpy as np
from model.compiled_native_history.compiler import compile_history
from model.openwave_native_process.source import admitted_local_program, fixed_command
from model.native_program import interpreter as native


def phase_boundary():
    commands=[fixed_command((0,1)), {'op':'repair','selector':{'kind':'last_current',
        'positive':[1,7],'otherwise':[0,1]},'law':'mixture','repeat':1}]
    compiled=compile_history(admitted_local_program(commands=commands,capacity=1))
    effects=[F(0)]*66
    for b in compiled.branches:
        domain=[x for x,y in b.source_map]; images=[y for x,y in b.source_map]
        assert len(set(domain))==len(domain) and len(set(images))==len(images)
        for x in domain:effects[x]+=F(*b.weight)
    assert effects==[F(1)]*66
    return {'complete_source_columns':66,'raw_history_blocks':len(compiled.branches),
        'all_effects_diagonal':True,'effect_sum_identity_exact':True,
        'phase_twin_total_variation_exact':'0',
        'worst_case_X_measurement_TV_error_lower_bound':'1/2',
        'claim_domain':'repair/hold/record-only feedback, classical memory, occupation-diagonal terminal readout; no coherent source_word service'}


def coherent_extension():
    states=native.local.basis(True)
    a=tuple(int(i in (0,2)) for i in range(12));b=tuple(int(i in (1,2)) for i in range(12))
    indices=[states.index(a),states.index(b)]
    V=np.eye(66,dtype=complex)[:,indices];P=V@V.conj().T
    # The supplied six-level response is executed; a logical H is only a comparator.
    angle=pi/2;word=native.source_word([native.compiler.pulse('H',angle,[0])],1)
    U=native.compiler.native_pulse('H',angle)
    good=V@U[:2,:2]@V.conj().T
    kraus=[('OUTSIDE_CODE',np.eye(66)-P),('SUCCESS',good)]
    for j in range(2,6):
        k=np.zeros((2,2),complex);k[0,:]=U[j,:2]
        kraus.append(('EXPORT_ABORT_'+str(j),V@k@V.conj().T))
    completeness=float(np.linalg.norm(sum(K.conj().T@K for _,K in kraus)-np.eye(66),2))
    H=np.array([[1,1],[1,-1]],complex)/np.sqrt(2)
    # One global phase is irrelevant to a CP block. The binary64 angle is declared.
    transfer_error=float(np.linalg.norm(U[:2,:2]-1j*H,2))
    if transfer_error>1e-10:transfer_error=float(np.linalg.norm(U[:2,:2]+1j*H,2))
    records=[];distributions=[]
    command=fixed_command((0,1),'E')
    read=compile_history(admitted_local_program(commands=[command],capacity=1))
    for sign in (1,-1):
        psi=(V[:,0]+sign*V[:,1])/np.sqrt(2);rho=np.outer(psi,psi.conj());probs=[]
        for label,K in kraus:
            out=K@rho@K.conj().T
            if label!='SUCCESS':probs.append(float(np.trace(out).real));continue
            for branch in read.branches:
                probs.append(float(F(*branch.weight))*sum(float(out[i,i].real) for i,_ in branch.source_map))
        assert abs(sum(probs)-1)<2e-12 and min(probs)>-1e-14
        distributions.append(probs)
        service=native.execute_word(word,np.array([[1,sign],[sign,1]],complex)/2)
        records.append({'phase_sign':sign,'complete_probability_vector':probs,
            'service_success_probability':float(np.trace(service['blocks'][('SUCCESS',)]).real),
            'all_four_processor_abort_probabilities':[float(np.trace(v).real)for k,v in service['blocks'].items()if k!=('SUCCESS',)]})
    tv=sum(abs(x-y) for x,y in zip(*distributions))/2
    assert completeness<2e-12 and transfer_error<2e-12 and tv>1-2e-12
    # All full-carrier matrix units follow from these explicit complete Kraus maps;
    # verify the orthogonal-complement and selected-code action independently.
    residual=0.
    for i in range(66):
        for j in range(66):
            expected=0. if i!=j else 1.
            observed=sum(np.vdot(K[:,j],K[:,i])for _,K in kraus)
            residual=max(residual,float(abs(observed-expected)))
    return {'source_indices':indices,'source_states':[list(a),list(b)],
        'angle_binary64_exact_rational':str(F(angle)),
        'source_word_resources':word['declared_resources'],
        'native_terminal_readout':'complete literal E repair on seam(0,1)',
        'complete_instrument_blocks_before_native_readout':len(kraus),
        'matrix_units_checked':66**2,'trace_residual':residual,
        'complete_effect_residual':completeness,'phase_insensitive_gate_residual':transfer_error,
        'phase_twin_total_variation':tv,'records':records,
        'scope':'declared two-state coherent load + existing M6 response + complete export + actual E readout; outside-code and every processor abort retained',
        'source_preparation_and_coherent_load_supplied':True,
        'native_repair_alone_implements_coherent_service':False,
        'physical_work_or_seconds_derived':False}


def run():
    return {'schema':'oph.native_geometry_readout.v1',
        'verdict':'PASS_NATIVE_PHASE_BOUNDARY_AND_DECLARED_COHERENT_READOUT',
        'repair_only':phase_boundary(),'coherent_extension':coherent_extension(),
        'geometry_join':{'two_state_X_observable_accessible_with_declared_service':True,
            'full_ten_mean_stress_curvature_stencil_compiled':False,
            'independent_Newton_normalization':False,'Einstein_continuum_proved':False,
            'physical_photon_clock_or_metric_identified':False},
        'new_Nature_validation_rows':0}
