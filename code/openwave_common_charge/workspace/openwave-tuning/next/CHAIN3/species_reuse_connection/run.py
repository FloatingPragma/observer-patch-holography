"""Execute the register-to-mass-slot connection and retain complete records."""
from fractions import Fraction as F
from itertools import product, islice
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import sys
import numpy as np
from scipy.linalg import expm
import source as s

HERE = Path(__file__).resolve().parent


def save(name, value):
    (HERE/name).write_text(json.dumps(value, indent=2, default=str)+'\n')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix(value):
    return {'real':value.real.tolist(), 'imag':value.imag.tolist()}


def providers():
    expected = '719036fe46769e68202731f765045dc3d75ebc8fe6e901e85319ce435776779f'
    assert sha(s.FIELD/'MANIFEST.json') == expected
    manifest = json.loads((s.FIELD/'MANIFEST.json').read_text())
    for relative, value in manifest['files_sha256'].items():
        assert sha(s.FIELD/relative) == value
    for relative, value in manifest['producer_files'].items():
        assert sha(s.FIELD/'inputs'/relative) == value
    return manifest


def register_certificate(data, boxes):
    weak = []
    for N in (1,2,8,32,s.ROUNDS):
        row=s.complete_register_instrument(data,N)
        for key in ('success_effect','analyzer_complement_effect','all_first_failure_effect'):
            row[key]=matrix(row[key])
        row['bounds']=s.weak_limit_bounds(data['P'],N)
        weak.append(row)
    actual=s.complete_register_instrument(data,s.ROUNDS)['success_effect']
    omega=np.exp(2j*np.pi/3)
    Fourier=np.array([[omega**(-j*k)/math.sqrt(3) for k in range(3)]for j in range(3)])
    target=Fourier@np.diag(list(map(float,boxes['mass_midpoints'])))@Fourier.conj().T
    err=float(np.linalg.norm(actual-target,2));assert err<2e-11
    z=(1-float(data['P'])+float(data['P'])*np.exp(1j/s.ROUNDS))**s.ROUNDS
    normalized=z/abs(z)
    wrong=data['matrices']['0']+normalized*data['matrices']['+']+normalized.conjugate()*data['matrices']['-']
    wrong_effect=wrong.conj().T@wrong
    delta=float(np.linalg.norm(actual-wrong_effect,2));assert delta>1e-9
    # Changing the source record changes the effect. This prevents a supplied
    # fixed list of ratios from passing while the executed source is ignored.
    plain=s.register.compression(False)
    plain={key:(A+math.sqrt(2)*B)/(96*math.sqrt(3))for key,(A,B)in plain.items()}
    plainM=plain['0']+z*plain['+']+z.conjugate()*plain['-']
    plain_delta=float(np.linalg.norm(actual-plainM.conj().T@plainM,2));assert plain_delta>.01
    return {'executed_tick_dimension':432,'integer_tick_Gram':'16 I',
            'integer_compression_parts':{k:{'A':A.tolist(),'B':B.tolist()}for k,(A,B)in data['parts'].items()},
            'exclusive_count':str(data['P']),'joint_weights':{k:str(v)for k,v in data['weights'].items()},
            'weak_rounds':weak,'selected_rounds':s.ROUNDS,
            'round_selection':'declared dyadic accuracy 2^-20 for unit-kick subdivision; no mass target consulted',
            'finite_weak_amplitude_rectangle':[v.strings() for v in boxes['weak_amplitude']],
            'lambda_intervals':[v.strings()for v in boxes['lambda']],
            'effect_intervals':[v.strings()for v in boxes['mass']],
            'mass_midpoints_exact':list(map(str,boxes['mass_midpoints'])),
            'actual_instrument_effect_vs_interval_midpoint_operator_error':err,
            'wrong_unit_modulus_return_effect_difference':delta,
            'plain_record_effect_difference':plain_delta,
            'weak_all_return_not_normalized':True}


def native_mass_certificate(masses):
    spec=importlib.util.spec_from_file_location('frozen_field_run',s.FIELD/'run.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    template=list(s.flavor_mass_word(masses))
    events=list(s.field_words.pulses(template,s.PRECISION))
    actual,failure,gram=module.replay(events,3)
    values=list(map(float,masses))+[0.]
    ideal=np.diag([np.exp(-1j/8*(1 if not (x&1)else -1)*values[x>>1])for x in range(8)])
    residual=float(np.linalg.norm(actual-ideal,2));assert residual<1e-11
    allowance=sum(F(e['full_response_operator_error_exact'])for e in events)
    assert allowance<F(1,10**20)
    # The flavor padding label is an identity mass block. It is retained as
    # compiler padding, never identified as a fourth physical family.
    assert np.linalg.norm(actual[:,6:8]-np.eye(8)[:,6:8])<1e-11
    firstfail=float(np.trace(failure).real)
    programs={'mass_slot':events};records=[]
    for setting in ('X','Y'):
        for clock_setting in ('X','Y'):
            gates=[]
            for wire in (0,3):gates+=list(s.field_words.remap(s.native.primitive_templates()['hadamard'],(wire,)))
            gates+=template
            gates.append(('Z',(3,),s.native.angle(0,0,F(-1,128))))
            for wire,basis in ((0,setting),(3,clock_setting)):
                if basis=='Y':gates.append(('Z',(wire,),s.native.angle(F(1,4))))
                gates+=list(s.field_words.remap(s.native.primitive_templates()['hadamard'],(wire,)))
            pulses=list(s.field_words.pulses(gates,s.PRECISION))
            U,fail,gr=module.replay(pulses,4)
            programs[setting+clock_setting]=pulses
            pointer_template=list(s.field_words.pointer_word(0,4,5))+list(s.field_words.pointer_word(3,6,7))
            pointer_pulses=list(s.field_words.pulses(pointer_template,s.PRECISION))
            programs[setting+clock_setting+'_both_native_pointers']=pulses+pointer_pulses
            for k in range(3):
                initial=2*k
                state=U[:,initial]
                q=np.zeros((2,2),float)
                for x,amp in enumerate(state):
                    assert (x>>1)&3==k or abs(amp)<1e-11
                    q[x&1,(x>>3)&1]+=abs(amp)**2
                # Same finite-clock source law. Mass uses phase2hm, clock
                # uses phaseh/8. Both are measured with terminal apparatus.
                h=1/8
                mu=(math.cos if setting=='X'else math.sin)(2*h*values[k])
                ck=(math.cos if clock_setting=='X'else math.sin)(h/8)
                expected=np.outer([(1+mu)/2,(1-mu)/2],[(1+ck)/2,(1-ck)/2])
                err=float(np.max(abs(q-expected)));assert err<2e-11,(setting,clock_setting,k,q,expected)
                failprob=float(fail[initial,initial].real)
                assert abs(q.sum()+failprob-1)<2e-11
                records.append({'family_character':k,'mass_basis':setting,'clock_basis':clock_setting,
                    'all_four_joint_binary_probabilities':q.tolist(),'all_native_first_export_failures':failprob,
                    'total':float(q.sum()+failprob),'ideal_probability_max_error':err,
                    'native_pulses':len(pulses),'source_time_exact':'1/8',
                    'both_native_pointer_pulses':len(pointer_pulses),
                    'both_native_pointer_operator_error_upper_exact':str(sum(F(e['full_response_operator_error_exact'])for e in pointer_pulses)),
                    'energy_from_phase':'atan2(<Y_mass>,<X_mass>)/(2h)',
                    'clock_from_phase':'8 atan2(<Y_clock>,<X_clock>) = h'})
    # The frozen native terminal C/E instrument includes every fine label.
    receipt=json.loads((s.FIELD/'RESULT.json').read_text())
    pointer=receipt['pointer'] if 'pointer'in receipt else None
    if pointer is None:
        pointer=next(v for v in receipt.values()if isinstance(v,dict)and'native_fine_effects'in v)
    fine=[]
    for row in records:
        probabilities=np.array(row['all_four_joint_binary_probabilities'])
        for law, branches in pointer['native_fine_effects'].items():
            for clock_law, clock_branches in pointer['native_fine_effects'].items():
                labels=[]
                for branch in branches:
                    effect=np.array(branch['effect_real'])+1j*np.array(branch['effect_imag'])
                    assert np.linalg.norm(effect-np.diag(np.diag(effect)),2)<2e-11
                    for clock_branch in clock_branches:
                        clock_effect=np.array(clock_branch['effect_real'])+1j*np.array(clock_branch['effect_imag'])
                        assert np.linalg.norm(clock_effect-np.diag(np.diag(clock_effect)),2)<2e-11
                        value=float(np.sum(probabilities*np.outer(np.diag(effect).real,np.diag(clock_effect).real)))
                        labels.append({'mass_record':branch['record'],'clock_record':clock_branch['record'],'probability':value})
                assert abs(sum(x['probability']for x in labels)-probabilities.sum())<2e-11
                fine.append({'family_character':row['family_character'],'mass_basis':row['mass_basis'],
                     'clock_basis':row['clock_basis'],'mass_law':law,'clock_law':clock_law,'fine_labels':labels,
                     'pointer_application':'two disjoint copies of the frozen all-column native pointer, with explicit16pulse joint word',
                     'unconditional_source_export_failure':row['all_native_first_export_failures']})
    reset=s.native.reset_repetition_plan(8,F(1,10**8))
    # Compare the actual generator suffix without enumerating the huge
    # provider word. The function delegates that suffix unchanged.
    flavored=s.flavored_macro(1,2,1,s.PRECISION,masses)
    info=s.field_words.layout(1,2)
    prefix=list(s.flavor_mass_word(masses,flavor=(info['memory_bits'],info['memory_bits']+1),beta=info['spin'][0]))
    assert list(islice(flavored,len(prefix)))==prefix
    old=iter(s.field_words.macro_word(1,2,1,s.PRECISION));next(old)
    assert list(islice(flavored,48))==list(islice(old,48))
    save('NATIVE_PROGRAMS.json',programs)
    save('NATIVE_RECORDS.json',{'joint_source_clock':records,'native_pointer_fine_records':fine,
           'pointer_contract':pointer,'blank_reset':reset,
           'complete_bound_rule':'Add declared reset timeout and all first native export failures, retaining their labels. No success conditioning.'})
    return {'mass_slot_pulses':len(events),'all_eight_columns_checked':True,
            'mass_operator_residual':residual,'native_operator_error_upper_exact':str(allowance),
            'complete_Gram_residual':gram,'first_export_failure_Gram_trace':firstfail,
            'joint_source_clock_instruments':len(records),'fine_pointer_instruments':len(fine),
            'blank_reset':reset,'shared_clock_ticks_after_source_word':True,
            'mass_and_clock_both_have_native_C_E_fine_readout':True,
            'arbitrary_flavor_quantum_reference_preserved':True,
            'actual_provider_suffix_first48operations_equal':True}


def full_native_spectrum(masses):
    from model.common_source.chiral_candidate import kinetic
    upstream=json.loads((s.RER/'model/openwave_full_higgs_floor/results.json').read_text())
    for relative,expected in upstream['source_sha256'].items():
        assert sha(s.RER/relative)==expected,relative
    theorem=upstream['exact_native_symbol']
    D=np.zeros((15,15),int);charge3=np.zeros(15,int)
    for pair in theorem['all_native_pairs']:
        D[pair['left'],pair['right']]=int(pair['coefficient'])
        charge3[pair['left']]=int(3*F(pair['left_charge']))
        charge3[pair['right']]=int(3*F(pair['right_charge']))
    C=D+D.T
    assert np.array_equal(D@D,np.zeros((15,15),int))
    assert np.array_equal(C@C,np.diag((charge3!=0).astype(int)))
    assert np.array_equal(np.diag(charge3)@C+C@np.diag(charge3),np.zeros((15,15),int))
    lepton=next(p for p in theorem['all_native_pairs']if p['channel']=='charged_lepton')
    basis=np.zeros((15,2),complex);basis[lepton['left'],0]=int(lepton['coefficient']);basis[lepton['right'],1]=1
    small=np.zeros((8,4),complex);small[0,0]=small[2,1]=small[7,2]=1;small[5,3]=-1
    frames=[]
    for k in range(3):
        u=np.kron(basis,np.eye(3)[:,k:k+1]);spin=np.kron(np.eye(2),u)
        full=np.block([[spin,np.zeros_like(spin)],[np.zeros_like(spin),spin.conj()]])
        frames.append(full@small)
    def dirac_hamiltonian(p,m):
        return -np.kron(np.diag([1,-1]),sum((x*a for x,a in zip(p,kinetic.pauli_matrices())),np.zeros((2,2),complex)))+m*np.kron(np.array([[0,1],[1,0]]),np.eye(2))
    values=np.array(list(map(float,masses)));rows=[];worst=intertwine=0.
    Q=np.kron(np.eye(2),np.kron(np.diag(charge3/3),np.eye(3)))
    charge=np.block([[Q,np.zeros_like(Q)],[np.zeros_like(Q),-Q.T]])
    minima=[0,0,0]
    for address in product(range(4),repeat=3):
        p=np.array([(0,1,0,-1)[i]for i in address]);w=sum(1-(1,0,-1,0)[i]for i in address)
        H=kinetic.bdg_symbol(p,np.kron(C,np.diag(values+w)),speed=1.)
        eig=np.linalg.eigvalsh(H)
        positive=np.r_[np.repeat(np.sqrt(p@p+(values+w)**2),28),np.repeat(np.sqrt(p@p),6)]
        err=float(np.max(abs(eig-np.sort(np.r_[-positive,positive]))));assert err<3e-12
        assert np.linalg.norm(H@charge-charge@H)<2e-12
        for k,U in enumerate(frames):
            delta=float(np.linalg.norm(H@U-U@dirac_hamiltonian(p,values[k]+w)))
            assert delta<2e-12;intertwine=max(intertwine,delta)
            if not np.any(p)and w==0:minima[k]+=1
        rows.append({'momentum_address':address,'all180_Nambu_energies':eig.tolist(),'Wilson_shift':w,'spectral_error':err})
        worst=max(worst,err)
    assert minima==[1,1,1]
    save('FULL90_SPECTRUM.json',{'exact_arbitrary_mass_theorem':theorem,'all64momenta':rows,
        'physical_CAR_modes':5760,'all_Nambu_rows':11520,'charged_minima_by_family':minima,
        'neutral_zero_CAR_modes_retained':48,'neutral_reference_degeneracy':'2^48',
        'charged_channels_total':21,'source_amendment':'same native Higgs-Wilson action with executed finite register coefficients substituted; this is not the routed source',
        'upstream237sourcepins_verified':len(upstream['source_sha256']),
        'not_claimed':['native field25 acts on all90coordinates','old routed QW equals this Hamiltonian','interacting charged mass poles','target-free observed spectrum']})
    return {'full_Nambu_rows':11520,'full_native_eigenvalue_max_error':worst,
            'native_frame_intertwiner_max_error':intertwine,'charged_minima_by_family':minima,
            'all_other_states_retained':True,'common_routed_source_transfer_claimed':False}


def provenance():
    paths=set()
    for module in list(sys.modules.values()):
        filename=getattr(module,'__file__',None)
        if not filename:continue
        p=Path(filename).resolve()
        if p.suffix=='.pyc':
            try:p=Path(importlib.util.source_from_cache(str(p)))
            except ValueError:pass
        if p.is_file()and p.is_relative_to(s.ROOT)and (s.ROOT!=s.WORKSPACE or not p.is_relative_to(HERE)):paths.add(p)
    manifest=json.loads((s.FIELD/'MANIFEST.json').read_text())
    paths.add((s.FIELD/'MANIFEST.json').resolve())
    paths.update((s.FIELD/p).resolve()for p in manifest['files_sha256'])
    paths.update((s.FIELD/'inputs'/p).resolve()for p in manifest['producer_files'])
    for rel in ('model/lepton_register_round/PROOF.md','model/lepton_register_round/certificate.py',
                'model/openwave_single_valley_leptons/PROOF.md','model/openwave_single_valley_leptons/spectrum.py',
                'model/lepton_mass_compiler/PROOF.md','model/lepton_cycle/experimental_inputs.json',
                'model/openwave_full_higgs_floor/results.json'):
        paths.add((s.RER/rel).resolve())
    upstream=json.loads((s.RER/'model/openwave_full_higgs_floor/results.json').read_text())
    paths.update((s.RER/p).resolve()for p in upstream['source_sha256'])
    return {str(p.relative_to(s.ROOT)):sha(p)for p in sorted(paths)}


def field_join_budgets(boxes):
    provider=json.loads((s.FIELD/'RESULT.json').read_text())
    extra_reset=s.native.reset_repetition_plan(2,F(1,10**8))
    eta=max((b.hi-b.lo)/2 for b in boxes['mass'])
    rows=[]
    for contract in provider['source_contracts']:
        count=contract['macros'];bits=contract['native_program']['precision_bits']
        events=list(s.field_words.pulses(s.flavor_mass_word(boxes['mass_midpoints']),bits))
        mass_error=count*sum(F(e['full_response_operator_error_exact'])for e in events)
        coefficient=count*F(1,8)*eta
        reset=F(extra_reset['complete_timeout_union_bound_exact'])
        bound=F(contract['complete_source_record_TV_bound_exact'])+mass_error+coefficient+reset
        rows.append({'macros':count,'provider_record_TV_bound_exact':contract['complete_source_record_TV_bound_exact'],
                     'extra_flavor_blank_reset':extra_reset,'native_mass_word_error_exact':str(mass_error),
                     'finite_effect_coefficient_error_exact':str(coefficient),
                     'complete_flavored_record_TV_bound_exact':str(min(F(1),bound)),
                     'complete_flavored_record_TV_bound_upper':float(bound),
                     'source_precision_bits':bits,'mass_pulses_per_macro':len(events),
                     'provider_mass_word_error_not_subtracted':True,
                     'complete_full_source_executed':False})
    return rows


def main():
    providers();print('frozen field provider verified',flush=True)
    data=s.executed_register();boxes=s.enclosures(data)
    register=register_certificate(data,boxes);save('REGISTER.json',register)
    print('actual432register and complete weak instruments passed',flush=True)
    native=native_mass_certificate(boxes['mass_midpoints']);print('native flavor mass and common clock passed',flush=True)
    census=full_native_spectrum(boxes['mass_midpoints']);print('full90native coefficient substitution passed',flush=True)
    join=field_join_budgets(boxes)
    m=boxes['mass_midpoints'];anchor=F('0.00051099895069');unit=anchor/m[1]
    energies=[{'character':k,'mass_source_exact':str(m[k]),'GeV_exact':str(m[k]*unit),'GeV':float(m[k]*unit),
            'role':'calibration input'if k==1 else'retrospective conditional comparison'}for k in (1,2,0)]
    pins=provenance();save('RUNTIME_INPUTS.json',{'files_sha256':pins})
    result={'schema':'oph.species.common_mass_connection.v1','status':'PASS_SCOPED_FINITE_REGISTER_AND_NATIVE_MASS_SLOT',
        'register':register,'native_mass_and_clock':native,'full90_substitution':census,
        'field_provider_manifest_sha256':sha(s.FIELD/'MANIFEST.json'),
        'flavored_field_complete_bounds':join,
        'native_field_join':{'source':'task25 declared one-charge193mode spectral-code source',
             'flavor_qubits_added':2,'flavor_physical_labels':[0,1,2],'padding_label':3,
             'field_current_and_clock_are_identity_on_flavor':True,
             'word':'source.flavored_macro replaces exactly first scalar mass word; every following source operation is unchanged',
             'Gauss':'mass is diagonal on flavor and spin, hence commutes with charge/cell/field Gauss constraint',
             'retained_bounds':'task25 finite-macro branch/tail counts are unchanged for arbitrary controlled mass phases; add finite-register coefficient and native mass-word errors',
             'full_field_word_execution':False,'full_source_spatial_limit':False},
        'charged_energies_GeV':energies,'source_energy_unit_GeV_exact':str(unit),
        'calibration_ancestry':['electron energy anchor','retrospective analyzer balance and oriented read','exclusive count interpreted as unit-kick phase','squared analyzer effect interpreted as source mass'],
        'physical_full_rows_validated':[],'empirical_targets_used_to_select_rounds':False,
        'same_construction_all_species_established':False,
        'runtime_self_estimation_of_mass_from_sampled_counts':False,
        'missing_connections':['native full90 field/current compiler','interacting physical lepton identification','unfitted spectrum criterion','neutrino/color/decay measurement sector'],
        'official_rows_changed':[],'source_sha256':pins}
    save('RESULT.json',result)
    print(json.dumps({'status':result['status'],'GeV':energies,'provider_files':len(pins)},default=str),flush=True)


if __name__=='__main__':main()
