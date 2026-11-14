"""Actual regional electric-flux instruments on the g=1 common source."""
from bridge import *
from itertools import product
from functools import lru_cache
import math

def phase_word(values,probe):
    """Controlled exp(i pi Phi_R), using exact source Gauss eigenvalues."""
    for mask in range(8):
        coefficient=sum((F(v)*(-1)**((mask&j).bit_count()) for j,v in enumerate(values)),F())/8
        for control,c in (('I',F(1,2)),('Z',F(-1,2))):
            label=''.join('Z' if mask>>j&1 else 'I' for j in range(3))+control
            yield from words.remap(native.pauli(label,native.angle(coefficient*c)),(0,1,2,probe))

def apparatus(regions,clock_setting):
    count=len(regions);clock=3+count;word=[]
    for j,region in enumerate(regions):
        probe=3+j
        word+=list(words.remap(native.primitive_templates()['hadamard'],(probe,)))
        word+=list(phase_word(region['flux_eigenvalues_by_cell'],probe))
        word+=list(words.remap(native.primitive_templates()['hadamard'],(probe,)))
    if clock_setting=='Y':word.append(('Z',(clock,),native.angle(F(1,4))))
    word+=list(words.remap(native.primitive_templates()['hadamard'],(clock,)))
    return word,clock

def source_diagnostics():
    # Pinned g=1, actual mass, field reference and all source branches.
    state=source.seed(plane=True);reports=[]
    for macros in (0,1,2):
        records=source.records(state,electric_moments=False)
        probs=np.zeros(8)
        for row in records['positions']:probs[row['occupied_root_cells'][0]]=row['probability']
        assert np.max(abs(probs-np.ones(8)/8))<1e-10
        reports.append({'macros':macros,'retained_branches':len(state.amps),'all8_cell_probabilities':probs.tolist(),
            'native_source_mass':'1','spacing':'1/8','coupling':True,'electric_charge':True,'free_field':True,
            'maximum_transverse_momentum_norm':float(max(np.linalg.norm(state.p,axis=1))),
            'physical_source_time':str(F(macros,8)),'field_lattice_time':state.time,
            'all_spin_probabilities':records['all_spin_basis_probabilities'],'normalization':records['norm_squared']})
        if macros<2:state=source.step(state,h=.125,mass=1.,coupling=True,free=True,electric_charge=True)
    assert reports[1]['maximum_transverse_momentum_norm']>1
    return reports

def first_failures(tape,bits):
    """Disjoint first-export effects, preserving pulse and source-column IDs."""
    dim=1<<bits;M=np.eye(dim,dtype=complex);total=np.zeros((dim,dim),complex);rows=[]
    for number,event in enumerate(tape):
        native.verify_pulse(event)
        response=common.response(tuple(event['controls_exact']));d=1<<len(event['targets'])
        ix=common.indices(bits,tuple(event['targets']));old=M[ix]
        failed=np.einsum('ab,bcr->acr',response[d:,:d],old).reshape(-1,dim)
        effect=failed.conj().T@failed;total+=effect
        rows.append({'first_export_pulse':number,'targets':event['targets'],
                     'probability_on_every_input_basis_column':np.diag(effect).real.tolist(),
                     'effect_trace':float(np.trace(effect).real)})
        M[ix]=np.einsum('ab,bcr->acr',response[:d,:d],old)
    residual=float(np.linalg.norm(M.conj().T@M+total-np.eye(dim),2));assert residual<3e-10
    return {'first_export_records':rows,'failure_effect_real':total.real.tolist(),
            'failure_effect_imag':total.imag.tolist(),'complete_Gram_residual':residual},total

def run():
    pins=verify();topology=read(HERE/'TOPOLOGY.json');regions=topology['regions']
    # Reuse an independently executed native pointer implementation, whose
    # controls are the same frozen field25 rational primitive source.
    import native_tomography
    branches,pointer_tape,pgram,pbound=native_tomography.pointer()
    pointer_failures,pfailed=first_failures(pointer_tape,3)
    pointer_failure=pfailed[np.ix_([0,4],[0,4])]
    fine_effects={law:[{'record':meta,'effect_real':E.real.tolist(),'effect_imag':E.imag.tolist()} for meta,E in rows] for law,rows in branches.items()}
    common_result=read(common.FIELD/'RESULT.json');prefix={row['macros']:row for row in common_result['source_contracts']}
    diagnostics=source_diagnostics();programs={'native_READY_pointer':pointer_tape};failure_records={'native_READY_pointer':pointer_failures};rows=[];maxgram=0.;maxlaw=0.;fine_total=0
    plus_tape=events(words.remap(native.primitive_templates()['hadamard'],(0,)))
    H,_,_,hbound=execute(plus_tape,1,np.eye(2));plus=H[:,0]
    experiments=[([name], [region]) for name,region in regions.items()]
    experiments.append((['daughter_A','daughter_B'],[regions['daughter_A'],regions['daughter_B']]))
    for names,selected in experiments:
        count=len(selected);bits=4+count;clock=3+count;dim=1<<bits
        for setting in ('X','Y'):
            word,actualclock=apparatus(selected,setting);assert actualclock==clock
            tape=events(word);U,failure,gram,wordbound=execute(tape,bits,np.eye(dim));maxgram=max(maxgram,gram)
            label='+'.join(names)+'_'+setting;programs[label]=tape
            failure_records[label],repeated_failure=first_failures(tape,bits)
            assert np.linalg.norm(repeated_failure-failure,2)<1e-14
            # The ideal apparatus is a complete affine integer flux readout
            # for every source cell and every supplied clock input.
            for cell in range(8):
                expected_bits=[abs(r['flux_eigenvalues_by_cell'][cell]) for r in selected]
                for clock_input in (0,1):
                    input_code=cell+(clock_input<<clock)
                    out=U[:,input_code]
                    for output,amplitude in enumerate(out):
                        if any(((output>>(3+j))&1)!=b for j,b in enumerate(expected_bits)):
                            assert abs(amplitude)<2e-11
            extra_registers=3*count+2 # each charge probe+two work, plus clock work
            resets=native.reset_repetition_plan(extra_registers,F(1,10**9))
            meterbound=wordbound+(count+1)*pbound+F(resets['complete_timeout_union_bound_exact'])
            for macros in (1,2):
                # Supply the actual ideal source cell marginal and its h/8
                # clock; the frozen all-source bound prices the coupled prefix.
                tick=events([('Z',(0,),native.angle(0,0,F(-1,128)))]*macros)
                C,Cfail,cgram,cbound=execute(tick,1,np.eye(2));cv=C@plus
                # All eight source cells are retained as separate columns.
                initial=np.zeros((dim,8),complex)
                for cell in range(8):
                    initial[cell,cell]=cv[0]/math.sqrt(8);initial[cell+(1<<clock),cell]=cv[1]/math.sqrt(8)
                out=U@initial
                # Partial trace over actual cell labels, keeping the charge
                # probes and clock. This includes every probe/padding output.
                measured=1<<(count+1);rho=np.zeros((measured,measured),complex)
                for cell in range(8):
                    block=out[cell::8,:];rho+=block@block.conj().T
                records={};forbidden=0.
                for laws in product(('C','E'),repeat=count+1):
                    choices=[branches[law] for law in laws];fine=[];coarse=np.zeros(measured)
                    for entries in product(*choices):
                        effect=np.ones((1,1),complex);code=0;metadata=[]
                        for j,(meta,E) in enumerate(entries):
                            effect=np.kron(E,effect);code|=meta['bit']<<j;metadata.append(meta)
                        probability=float(np.trace(effect@rho).real);assert probability>-1e-12
                        fine.append({'records':metadata,'probability':probability});coarse[code]+=probability;fine_total+=1
                    assert abs(sum(coarse)-np.trace(rho).real)<3e-11
                    expected=np.zeros(measured);clockmean=(math.cos if setting=='X'else math.sin)(macros/64)
                    for cell in range(8):
                        chargecode=sum(abs(r['flux_eigenvalues_by_cell'][cell])<<j for j,r in enumerate(selected))
                        expected[chargecode]+=F(1,8)*(1+clockmean)/2
                        expected[chargecode+(1<<count)]+=F(1,8)*(1-clockmean)/2
                    err=float(np.max(abs(coarse-expected)));assert err<3e-10;maxlaw=max(maxlaw,err)
                    if count==2:
                        forbidden=max(forbidden,sum(coarse[x] for x in range(measured) if x&3==3));assert forbidden<2e-11
                    first_pointer_failures=[]
                    for failed_pointer in range(count+1):
                        effect=np.ones((1,1),complex)
                        for j in range(count+1):
                            local=(sum((E for _,E in choices[j]),np.zeros((2,2),complex)) if j<failed_pointer
                                   else pointer_failure if j==failed_pointer else np.eye(2))
                            effect=np.kron(local,effect)
                        first_pointer_failures.append(float(np.trace(effect@rho).real))
                    assert abs(sum(coarse)+sum(first_pointer_failures)-np.trace(rho).real)<3e-11
                    records[''.join(laws)]={'all_joint_binary_probabilities':coarse.tolist(),'fine_records':fine,
                        'first_pointer_export_probability_by_pointer':first_pointer_failures,
                        'source_charge_decode':['Q_R=-(2 G_R-1)*probe_bit' for r in selected]}
                complete=min(F(1),F(prefix[macros]['complete_source_record_TV_bound_exact'])+meterbound)
                rows.append({'regions':names,'clock_setting':setting,'macros':macros,'native_meter_pulses':len(tape),
                    'laws':records,'forbidden_both_daughters_occupied_probability':forbidden,
                    'all_source_prefix_bound_exact':prefix[macros]['complete_source_record_TV_bound_exact'],
                    'additional_meter_bound_exact':str(meterbound),'complete_common_source_record_TV_bound_exact':str(complete),
                    'complete_common_source_record_TV_bound':float(complete),'extra_blank_reset_contract':resets,
                    'component_law':'Before actual reset-success weights. Compose the complete frozen coupled source and rank-one native resets, retaining first timeouts and source/pointer export failures.'})
    # Omitting the declared reference flips its physical flux bit.
    ref=regions['reference_singleton'];wrong=dict(ref)
    wrong['flux_eigenvalues_by_cell']=[v-1 for v in ref['flux_eigenvalues_by_cell']]
    wrongword,_=apparatus([wrong],'X');wrongU,_,_,_=execute(events(wrongword),5,np.eye(32))
    correctword,_=apparatus([ref],'X');correctU,_,_,_=execute(events(correctword),5,np.eye(32))
    wrong_reference=float(np.linalg.norm(wrongU-correctU,2));assert wrong_reference>1.9
    # A 2pi field-flux Fourier loop is identity on every admitted integer
    # sector. The half-charge mutation has a -1 loop on the occupied sector.
    particle=regions['particle_singleton'];integer_loop=dict(particle)
    integer_loop['flux_eigenvalues_by_cell']=[2*v for v in particle['flux_eigenvalues_by_cell']]
    loop_word,_=apparatus([integer_loop],'X');half_word,_=apparatus([particle],'X')
    loop_tape=events(loop_word);half_tape=events(half_word)
    loopU,_,loopgram,loopbound=execute(loop_tape,5,np.eye(32))
    halfU,_,halfgram,halfbound=execute(half_tape,5,np.eye(32))
    # Actual source basis cell1 has physical charge1 inside this sphere.
    source_column=1
    loop_defect=float(sum(abs(loopU[x,source_column])**2 for x in range(32) if x&8))
    half_signal=float(sum(abs(halfU[x,source_column])**2 for x in range(32) if x&8))
    assert loop_defect<2e-12 and half_signal>1-2e-12
    programs['integer_2pi_loop']=loop_tape;programs['half_charge_2pi_mutant']=half_tape
    failure_records['integer_2pi_loop'],_=first_failures(loop_tape,5)
    failure_records['half_charge_2pi_mutant'],_=first_failures(half_tape,5)
    assert verify()==pins
    save('NATIVE_PROGRAMS.json',programs)
    save('NATIVE_FIRST_FAILURES.json',failure_records)
    save('NATIVE_RECORDS.json',{'source_diagnostics':diagnostics,'experiments':rows,'fine_pointer_effects':fine_effects})
    save('NATIVE_RESULT.json',{'schema':'oph.common_source.native_integer_flux.v1','all_checks_pass':True,
        'g':1,'source_preparation':'Existing uniform8 root-cell plane, spin00, all190 coexact vacuum and3 normalized harmonic Gaussians.',
        'source_clock':'Existing h=1/8 and h/8 phase tick after completed g=1 macro.',
        'source_prefix_native_whole_iterator_executed':False,'native_appended_meter_all_columns_executed':True,
        'complete_meter_experiments':len(rows),'fine_probabilities_evaluated':fine_total,'maximum_complete_Gram_residual':maxgram,
        'maximum_law_diagnostic_error':maxlaw,'pointer_Gram_residual':pgram,'wrong_reference_operator_separation':wrong_reference,
        'maximum_complete_common_source_TV_bound':max(r['complete_common_source_record_TV_bound'] for r in rows),
        'actual_joint_daughter_readout':True,'extra_history_register':False,'source_success_normalized':False,
        'native_2pi_integer_vs_half_control':{'integer_loop_wrong_bit_probability':loop_defect,
            'half_mutant_wrong_bit_probability':half_signal,'integer_loop_native_error_exact':str(loopbound),
            'half_mutant_native_error_exact':str(halfbound),'complete_Gram_residuals':[loopgram,halfgram],
            'scope':'Half mutation is a rejected source/curvature control, not an admitted preparation.'},
        'field_Fourier_equivalence':'exp(i t Phi_R)=exp(i t (G_R-n_R))=product_e exp(i t (B^T1_R)_e E_e) on the exact U Gauss fiber. This compiles the exact cancelled observable, not the rounded-Q edge implementation.',
        'new_state_measurement':'Only named terminal regional fluxes and source clock; no reconstruction of unmeasured intermediate currents.'})
if __name__=='__main__':run()
