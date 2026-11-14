"""Regenerate finite native controls, exact bounds, and scoped development replay."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import islice
from pathlib import Path
import hashlib
import json
import math
import sys
import numpy as np
from scipy.linalg import expm

from providers import ROOT, HERE, PINS, native, verify
import geometry as g
import native_oscillators as words
import transfer


def write(name,value):
    (HERE/name).write_text(json.dumps(value,indent=2,default=str)+'\n')


@lru_cache(None)
def source_generators():
    response=json.loads((ROOT/'reverse-engineering-reality/code/source_selection_model/response.json').read_text())
    def scalar(a,b):return float(F(a))+math.sqrt(5)*float(F(b))
    result=[]
    for blocks in response['generators']:
        matrix=np.zeros((6,6),complex)
        for k,rows in enumerate(blocks):
            for i,row in enumerate(rows):
                for j,cell in enumerate(row):
                    matrix[3*k+i,3*k+j]=scalar(*cell[:2])+1j*scalar(*cell[2:])
        result.append(matrix)
    assert len(result)==12
    return result


def replay(events,bits):
    """All code columns from independent frozen full six-level source matrices."""
    dim=1<<bits; state=np.eye(dim,dtype=complex); failure=np.zeros((dim,dim),complex)
    cache={}
    for event in events:
        native.verify_pulse(event)
        key=tuple(event['controls_exact'])
        if key not in cache:
            generator=sum(float(F(c))*G for c,G in zip(key,source_generators()))
            cache[key]=expm(generator)
        full=cache[key]; targets=event['targets']; local=1<<len(targets)
        mask=sum(1<<q for q in targets)
        bases=[x for x in range(dim) if not x&mask]
        indices=np.array([[base+sum(((value>>(len(targets)-1-j))&1)<<q
                                   for j,q in enumerate(targets)) for base in bases]
                          for value in range(local)])
        gathered=state[indices]
        rejected=np.einsum('ab,bcr->acr',full[local:,:local],gathered).reshape(-1,dim)
        failure+=rejected.conj().T@rejected
        state[indices]=np.einsum('ab,bcr->acr',full[:local,:local],gathered)
    gram=float(np.linalg.norm(state.conj().T@state+failure-np.eye(dim),2))
    assert gram<2e-10
    return state,failure,gram


def small_controls():
    bits=80; cap=2; wires=(0,1); control=((2,1),)
    ann=np.zeros((4,4),complex)
    for n in range(1,cap+1):ann[n-1,n]=math.sqrt(n)
    X=ann+ann.conj().T; Y=1j*(ann.conj().T-ann)
    rows=[]; tapes={}
    for quadrature,H in (('X',X),('Y',Y)):
        for predicate in (0,1):
            template=list(words.translation_word(cap,F(2,7),4,bits,wires,((2,predicate),),quadrature))
            events=list(words.pulses(template,bits)); actual,failure,gram=replay(events,3)
            target=np.eye(8,dtype=complex)
            target[4*predicate:4*predicate+4,4*predicate:4*predicate+4]=expm(2j*H/7)
            discrepancy=float(np.linalg.norm(actual-target,2))
            native_bound=sum(F(e['full_response_operator_error_exact']) for e in events)
            bound=words.translation_error(cap,F(2,7),4,bits)+native_bound
            assert discrepancy<float(bound)+1e-10
            # Both cap padding columns and every false predicate are defined.
            assert np.linalg.norm(actual[:,3]-np.eye(8)[:,3])<1e-10
            assert np.linalg.norm(actual[:,7]-np.eye(8)[:,7])<1e-10
            name=f'{quadrature}_control{predicate}'
            tapes[name]=events
            rows.append({'name':name,'pulses':len(events),'all_columns':8,
                         'operator_error':discrepancy,'bound_exact':str(bound),
                         'native_rounding_bound_exact':str(native_bound),
                         'complete_Gram_residual':gram,'all_export_failure_Gram_trace':float(np.trace(failure).real)})
    # The same test includes a vacuum scalar, local padding, a noncommuting
    # displacement and a literal compressed harmonic p^2 free interval.
    template=list(words.translation_word(cap,F(1,5),8,bits,wires))
    template+=list(words.diagonal_linear_word(cap,F(-11,30),F(-11,60),wires))
    template+=list(words.harmonic_free_word(cap,F(1,3),16,bits,wires))
    events=list(words.pulses(template,bits)); actual,failure,gram=replay(events,2)
    bigann=np.zeros((6,6),complex)
    for n in range(1,6):bigann[n-1,n]=math.sqrt(n)
    momentum=1j*(bigann.conj().T-bigann)/math.sqrt(2)
    free=np.zeros((4,4),complex);free[:3,:3]=(momentum@momentum)[:3,:3]/4
    diagonal=np.diag([-11*(n+.5)/30 if n<=cap else 0 for n in range(4)])
    target=expm(-1j*free/3)@expm(1j*diagonal)@expm(1j*X/5)
    error=float(np.linalg.norm(actual-target,2))
    native_bound=sum(F(e['full_response_operator_error_exact']) for e in events)
    bound=words.translation_error(cap,F(1,5),8,bits)+words.harmonic_free_error(cap,F(1,3),16,bits)+native_bound
    assert error<float(bound)+1e-10
    assert np.linalg.norm(actual[:,3]-np.eye(4)[:,3])<1e-10
    tapes['combined_oscillator']=events
    rows.append({'name':'combined_oscillator','pulses':len(events),'all_columns':4,
                 'operator_error':error,'bound_exact':str(bound),'complete_Gram_residual':gram,
                 'all_export_failure_Gram_trace':float(np.trace(failure).real),
                 'vacuum_phase_retained':True,'compressed_full_p_squared':True})
    # Exact prefix Walsh sums must include every physical and padding column.
    for cap in range(1,17):
        for mask in range(1<<cap.bit_length()):
            assert words.signed_prefix_sum(cap,mask)==tuple(sum((n**power)*(-1)**((n&mask).bit_count()) for n in range(cap+1)) for power in (0,1))
    mutated=dict(tapes['combined_oscillator'][0]); mutated['controls_exact']=mutated['controls_exact'].copy()
    mutated['controls_exact'][0]=str(F(mutated['controls_exact'][0])+F(1,1000))
    try:native.verify_pulse(mutated)
    except ArithmeticError:pass
    else:raise AssertionError('changed control certificate admitted')
    write('SMALL_NATIVE_PROGRAMS.json',{'schema':'oph.native_field.local_words.v1','words':tapes,
          'full_first_failure_rule':'Every event has six-level native exp(sum twelve rational controls*source generator); retain each first complement row and all later-stop labels. No conditional renormalization.',
          'reference':'all-column operator replay implies the same bound after arbitrary reference tensoring'})
    return rows


def program_contract(macros,co,harm,tolerance=F(1,10**6)):
    p,L,S=transfer.parameters(macros)
    # Sum all 16 branch predicates on each axis and all retained modes.
    # This deliberately uses triangle sums, including zero coefficients.
    gamma=max(L,S,F(1))+1
    translation_factor=macros*48*(190*gamma*gamma*co*(co+1)+3*gamma*gamma*harm*(harm+1))/4
    kinetic_norm=g.A*((2*harm+1)+harm*(harm+1))/4
    splitting_numerator=translation_factor+macros*3*kinetic_norm**2/2
    slices=math.ceil(splitting_numerator/tolerance)+1
    resources=words.resource_bound(co,harm,slices,macros)
    pulses=resources['native_pulses_upper']
    maxcap=max(co,harm)
    omega=transfer.up(g.A*g.B*F(g.data()['spectral']['maximum_eigenvalue']))
    angle_bound=max(F(2),gamma*(maxcap+1),omega*(maxcap+1),F(maxcap+1)**2)
    precision=max(80,(pulses*tolerance.denominator//tolerance.numerator+1).bit_length()+100)
    while True:
        compiler=pulses*words.universal_pulse_error(precision,angle_bound)
        rounding=F()
        for axis,cell,sign in product(range(3),range(8),(-1,1)):
            for mode,cap in enumerate([co]*190+[harm]*3):
                lo,hi=g.coefficient(cell,axis,sign,mode,precision)
                assert abs((lo+hi)/2)<=gamma
                rounding+=(hi-lo)*transfer.up(cap)+gamma*cap/F(1<<precision)
        rounding*=macros
        rounding+=macros*(190*(co+1)+3*g.A*harm/4)/F(1<<precision)
        if compiler+rounding<tolerance:break
        precision*=2
    return {'macros':macros,'coexact_cap':co,'harmonic_cap':harm,
            'subdivisions':slices,'precision_bits':precision,
            'whole_program_Trotter_upper_exact':str(splitting_numerator/slices),
            'coefficient_rounding_upper_exact':str(rounding),
            'native_pulse_export_operator_upper_exact':str(compiler),
            'resources':resources,'full_program_enumerated_or_executed':False,
            'finite_iterator':'native_oscillators.pulses(macro_word(Rc,Rh,slices,bits),bits), repeated macros times',
            'native_scope':'one charged root-cell spinor, all8 initial root cells and arbitrary spectator; no two-fermion route word claimed'}


def edge_program_contract(edge, time, co, harm, tolerance=F(1,10**6)):
    caps=[co]*190+[harm]*3
    omega=transfer.up(g.A*g.B*F(g.data()['spectral']['maximum_eigenvalue']))
    gamma=abs(F(time))*transfer.up((1+F(g.data()['spectral']['frame_Gram_error_upper']))*omega/(2*g.A))+1
    splitting=sum(gamma*gamma*cap*(cap+1)/4 for cap in caps)
    slices=max(1,math.ceil(splitting/tolerance))
    terms=2*slices*sum(cap*(1<<words.register_bits(cap))//2 for cap in caps)+16
    pulses=(8*(max(map(words.register_bits,caps))+1)+1)*terms+20
    precision=max(80,(pulses*tolerance.denominator//tolerance.numerator+1).bit_length()+100)
    angle=gamma*(max(caps)+1)+1
    while pulses*words.universal_pulse_error(precision,angle)>tolerance:
        precision*=2
    rounding=F()
    for mode,cap in enumerate(caps):
        lo,hi=g.edge_coefficient(edge,mode,precision)
        lo,hi=sorted((F(time)*lo,F(time)*hi))
        assert abs((lo+hi)/2)<=gamma
        rounding+=(hi-lo)*transfer.up(cap)+gamma*cap/F(1<<precision)
    # The pointer is two native CNOTs, with one logical flip preparing code2.
    pointer_template=list(words.pointer_word(2,0,1))
    pointer_bound=sum(F(native.certified_pulse(k,t,a,precision)['full_response_operator_error_exact']) for k,t,a in pointer_template)
    return {'edge':edge,'Fourier_parameter_exact':str(time),'subdivisions':slices,'precision_bits':precision,
            'native_pulses_upper':pulses,'Trotter_upper_exact':str(splitting/slices),
            'coefficient_rounding_upper_exact':str(rounding),
            'native_export_operator_upper_exact':str(pulses*words.universal_pulse_error(precision,angle)),
            'pointer_native_pulses':8,'pointer_native_operator_error_exact':str(pointer_bound),
            'probe_and_pointer_blank_reset':native.reset_repetition_plan(3,tolerance),
            'actual_record_rule':'After pointer CNOT(2,0), CNOT(0,1), use pinned native coin/current C or E grammar. Keep every fine label, first export complement and reset timeout.',
            'complete_iterator':'pulses(edge_ramsey_word(edge,time,layout(Rc,Rh),subdivisions,precision_bits,setting),precision_bits)',
            'field_source_error':'add edge_transfer source-state/measurement norm bound; all bounded terminal TV laws contract it',
            'sharp_bin_joint_measurement_claimed':False}


def pointer_control():
    sys.path.insert(0,str(ROOT/'openwave-tuning/auditor/common_walk_particle'))
    import native_phase_probe
    instruments,run=native_phase_probe.actual_native_pointer()
    tape=list(words.pulses(words.pointer_word(2,0,1),80))
    replayed,failure,gram=replay(tape,3)
    target=np.zeros((8,8),complex)
    for x in range(8):
        y=x^2
        y^=((y>>2)&1)
        y^=(y&1)<<1
        target[y,x]=1
    assert np.linalg.norm(replayed-target,2)<1e-10
    result={}
    for law,branches in instruments.items():
        effects=[sum((K.conj().T@K for K in matrices),np.zeros((2,2),complex)) for meta,matrices in branches]
        assert np.linalg.norm(sum(effects)-np.eye(2),2)<1e-11
        result[law]=[{'record':meta,'effect_real':effect.real.tolist(),'effect_imag':effect.imag.tolist()} for (meta,matrices),effect in zip(branches,effects)]
    return {'native_fine_effects':result,'full_source_pointer_export_failure_trace':float(np.trace(run.failure_gram).real),
            'rational_ready_pointer_word':tape,'rational_all8column_operator_residual':float(np.linalg.norm(replayed-target,2)),
            'rational_complete_Gram_residual':gram,
            'pointer_word':'CNOT(system,work0); CNOT(work0,work1), native READY code2',
            'outcome_formula':'probability = trace[(I_field tensor effect_fine)rho_after_Ramsey]; retain native source failure outcomes separately',
            'binary_Fourier_law':'P(bit0)=(1+Re chi(t))/2 for X; (1+Im chi(t))/2 for Y'}


def complete_edge_bound(row):
    field_comparison=sum(F(row[name]['state_error_upper_exact']) for name in ('coexact_box','harmonic_box','spectral_transfer'))
    native_and_reset=F(row['complete_source_record_TV_bound_exact'])-field_comparison
    probe=row['edge_Fourier_native_program']
    total=native_and_reset+F(row['edge_Fourier_source_transfer']['source_state_and_measurement_state_error_upper_exact'])
    total+=sum(F(probe[key]) for key in ('Trotter_upper_exact','coefficient_rounding_upper_exact',
                'native_export_operator_upper_exact','pointer_native_operator_error_exact'))
    total+=F(probe['probe_and_pointer_blank_reset']['complete_timeout_union_bound_exact'])
    row['complete_edge_Fourier_record_TV_bound_exact']=str(min(F(1),total))
    row['complete_edge_Fourier_record_TV_bound_upper']=float(min(F(1),total))
    return total


def freeze_outputs():
    verify()
    local=[p for p in HERE.iterdir() if p.is_file() and p.name not in ('MANIFEST.json',) and p.suffix not in ('.log',)]
    write('MANIFEST.json',{'schema':'oph.native_field.freeze.v1','files_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(local)},
           'producer_files':PINS,'official_rows_changed':[]})


def main():
    verify(); print('frozen producer pins verified',flush=True)
    controls=small_controls(); pointer=pointer_control(); print('complete small M6 words and C/E pointer checked',flush=True)
    geometry=g.data(); print('exact actual geometry and harmonic basis checked',flush=True)
    rows=[]
    for macros in (1,2):
        co,harm=transfer.choose_caps(macros)
        co_bound=transfer.coexact_box(macros,co); h_bound=transfer.harmonic_box(macros,harm)
        spectral=transfer.geometry_transfer(macros)
        native_contract=program_contract(macros,co,harm)
        prep=words.preparation_contract(co,harm)
        prep['matter_and_clock_preparation']={'matter':'spin00 and uniform all8 root cells',
            'native_pulses':8,'clock':'plus state; h=1/8, one tick -h Z/16 after every completed macro',
            'additional_blank_resets':6,
            'additional_reset_contract':native.reset_repetition_plan(6,F(1,10**8)),
            'native_word':[native.certified_pulse(k,t,a,native_contract['precision_bits'])
                for k,t,a in words.plane_preparation_word(words.layout(co,harm))]}
        preparation_error=sum(F(e['full_response_operator_error_exact']) for e in prep['matter_and_clock_preparation']['native_word'])
        total=sum(F(x) for x in (co_bound['state_error_upper_exact'],h_bound['state_error_upper_exact'],
                  spectral['state_error_upper_exact'],native_contract['whole_program_Trotter_upper_exact'],
                  native_contract['coefficient_rounding_upper_exact'],native_contract['native_pulse_export_operator_upper_exact'],
                  prep['native_reset']['complete_timeout_union_bound_exact']))
        total+=preparation_error+F(prep['matter_and_clock_preparation']['additional_reset_contract']['complete_timeout_union_bound_exact'])
        edge_contract=edge_program_contract(0,F(1,8),co,harm)
        rows.append({'macros':macros,'coexact_box':co_bound,'harmonic_box':h_bound,
                     'spectral_transfer':spectral,'native_program':native_contract,'preparation':prep,
                     'complete_source_record_TV_bound_exact':str(total),
                     'complete_source_record_TV_bound_upper':float(total),
                     'edge_Fourier_source_transfer':transfer.edge_transfer(macros,co,harm,F(1,8)),
                     'edge_Fourier_native_program':edge_contract})
        complete_edge_bound(rows[-1])
        print(json.dumps({'macros':macros,'caps':[co,harm],'bound':float(total),'bits':native_contract['precision_bits']}),flush=True)
    currents=[]
    for axis,cell,sign in product(range(3),range(8),(-1,1)):
        row=g.current(cell,axis,sign)
        assert len(row['coexact_rational_numerators'])==190 and len(row['harmonic_rational_numerators'])==3
        currents.append({'cell':cell,'axis':axis,'sign':sign,**row,
                         'coefficients_at80bits':[[str(x) for x in g.coefficient(cell,axis,sign,m,80)] for m in range(193)]})
    write('COORDINATES.json',{'currents':currents,'harmonic_orthogonal_rational_vectors':geometry['harmonics'],
          'harmonic_squared_norms':geometry['harmonic_norms'],'charge_sectors':geometry['sectors'],
          'coexact_frame':'frozen Q; exact physical frame U=polar(KQ Lambda^-1), with priced transfer',
          'edge_observable':'P_q(edge)+U(edge,:)p_co+H(edge,:)p_h; implemented characteristic uses Q with certified frame error',
          'all288_edges_all190plus3_modes':True})
    # Prefix witness is an actual emitted family member, not an ideal matrix.
    info=rows[0]['native_program']; word=words.macro_word(info['coexact_cap'],info['harmonic_cap'],info['subdivisions'],info['precision_bits'])
    prefix=list(islice(words.pulses(word,info['precision_bits']),32))
    assert all(native.verify_pulse(e) for e in prefix)
    write('FULL_PROGRAM_PREFIX.json',{'emitted_pulses':prefix,'prefix_only':True,'complete_iterator_contract':info})
    verify()
    static=words.static_holding_contract(F(8))
    write('STATIC_HOLDING.json',static)
    result={'schema':'oph.native_interacting_field.connection.v1',
            'status':'PASS_SCOPED_CONSTRUCTIVE_TENSOR_FIELD_INTERFACE',
            'small_complete_native_replays':controls,'source_contracts':rows,'complete_native_pointer':pointer,
            'actual_current_count':len(currents),'retained_modes':193,
            'microscopic_coefficients':{'a':'1/2','b':'3','empirical_coupling_installed':False},
            'complete_record_rule':'retain every reset timeout, all M6 first export complements, every padding result, and all native C/E pointer fine records; no selected-success normalization',
            'physical_time':'completed-word clock with phase n*h/8, h=1/8 for routed examples; cell field evolves time1 per macro; no seconds or elapsed pulse-duration theorem',
            'boundaries':['one-particle size2 routed source only; no spatial interacting refinement',
                          'edge instrument supplies a bounded Fourier characteristic, not complete joint sharp-edge bins',
                          'finite words and analytic whole-program bounds; no full193-mode hardware or dense simulation',
                          'native coherent load/export and reset services are retained declared architecture services',
                          'no stationary interacting species or observed coupling claim'],
            'static_holding_contract':'STATIC_HOLDING.json',
            'official_rows_changed':[],'physical_rows_validated':[], 'producer_input_files':len(PINS)}
    write('RESULT.json',result)
    freeze_outputs()
    print(json.dumps({'status':result['status'],'small_word_count':len(controls),'retained_modes':193}),flush=True)


if __name__=='__main__':
    from itertools import product
    main()
