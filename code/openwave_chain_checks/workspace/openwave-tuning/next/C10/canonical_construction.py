"""Reproduce one selected C10 construction/protocol and its source contract."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
from scipy.linalg import expm
from joint_source import JointNativeSource,BASE
from spectral_readout import execute as execute_spectral_readout
HERE=Path(__file__).resolve().parent


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def execute():
    joint=HERE/'JOINT_SOURCE.json';record=json.loads(joint.read_text())
    for relative,pin in record['source_sha256'].items():
        if digest(BASE/relative)!=pin:raise ValueError('joint source custody drift: '+relative)
    archive=HERE/record['state_archive']['filename']
    if digest(archive)!=record['state_archive']['sha256']:raise ValueError('joint archive custody drift')
    model=JointNativeSource();data=np.load(archive);h=model.hamiltonian;rho=data['global_allN_reference']
    if not np.array_equal(h,data['hamiltonian']):raise ValueError('canonical model differs from archived joint H')
    # The one selected protocol: H-native time .37, then probability .29 of MIX.
    time=.37;probability=.29;u=expm(-1j*time*h);coherent=u@rho@u.conj().T
    committed,ledger=model.field.commit(coherent,'MIX');post=(1-probability)*coherent+probability*committed
    outcomes=[{**r,'probability':probability*r['probability'],'unnormalized_total_energy':probability*r['unnormalized_total_energy'],'unnormalized_field_energy':probability*r['unnormalized_field_energy']} for r in ledger]
    outcomes.append({'status':'no_commit','probability':1-probability,'unnormalized_total_energy':float((1-probability)*np.trace(h@coherent).real),'unnormalized_field_energy':float((1-probability)*np.trace(model.field.field@coherent).real)})
    occupation=np.array([[(mask>>z)&1 for z in range(6)] for mask,_ in model.states]);position=np.diag(post).real@occupation
    charge_records=[];region={0,3}
    q=np.array([sum((mask>>z)&1 for z in region) for mask,_ in model.states])
    for value in range(7):charge_records.append({'charge':value,'probability':float(np.diag(post).real[q==value].sum())})
    p=np.diag(occupation[:,0]);current=1j*(h@p-p@h)
    zero=model.field.__class__(model.source,kappa=0,face_coefficient=0)
    gauss=all(model.source.gauss(s)==(1,1,1,0,0,0) for s in model.states)
    reference=record['references']['global_allN'];before=float(np.trace(h@rho).real);after=float(np.trace(h@post).real)
    spectral_inputs={key:data[key] for key in data.files if key!='hamiltonian' and data[key].shape==(model.n,model.n)}
    spectral_inputs['selected_protocol_poststate']=post
    spectral_result,spectral_arrays=execute_spectral_readout(model,spectral_inputs)
    sources=[Path(__file__),HERE/'spectral_readout.py',HERE/'joint_source.py',HERE/'charged_reservoir.py',HERE/'modular_field_adapter.py',joint,archive]
    result={'schema':'oph.c10.canonical_construction.v1','construction_id':'triangle-CAR-allN-Z6-periodic-half-record-native-MIX-v1','checks':{'same_joint_H_exact':True,'canonical_global_reference_used':bool(abs(before-reference['reference_energy'])<1e-12),'reference_coherent_stationarity':bool(np.max(abs(coherent-rho))<1e-12),'zero_field_is_exact_native_baseline':bool(np.array_equal(zero.hamiltonian,model.source.hamiltonian)),'fixed_background_Gauss':bool(gauss),'selected_protocol_trace':bool(abs(np.trace(post)-1)<1e-12),'selected_native_complete':bool(abs(sum(x['probability'] for x in outcomes)-1)<1e-12),'whole_region_charge_PVM_complete':bool(abs(sum(x['probability'] for x in charge_records)-1)<1e-12),'selected_instrument_energy_work':bool(abs(sum(x['unnormalized_total_energy'] for x in outcomes)-after)<1e-12)},'selected':{'H':'JOINT_SOURCE.json parameters.source_H; joint_source_states.npz hamiltonian','reference':'JOINT_SOURCE.json references.global_allN; joint_source_states.npz global_allN_reference','clock':'dimensionless source hopparameter; no physicalclock conversion','time':time,'event_probability':probability,'event_rate_interpretation':'one supplied Bernoulli split-step; no derivedPoissonrate or uniquecontinuumgenerator','instrument':'MIX: uniform seam, fair C/E and faircoin; no individualFcollapse','boundary_energy':0,'reservoir_energy':.23},'readouts':{'complete_selected_native_outcomes':outcomes,'all_six_position_occupations':position.tolist(),'whole_region':[0,3],'complete_whole_region_charge_PVM':charge_records,'same_H_port0_current_expectation_before':float(np.trace(current@rho).real),'same_H_port0_current_expectation_after':float(np.trace(current@post).real),'source_energy_before':before,'source_energy_after':after,'mean_selected_native_work':after-before,'field_energy_before':float(np.trace(model.field.field@rho).real),'field_energy_after':float(np.trace(model.field.field@post).real),'complete_final_number_energy_PVM':model.energy_number_readout(post)},'controls':{'fixed_N3_reference_equals_global_here':bool(np.max(abs(data['fixed_N3_reference']-rho))<1e-12),'zero_field_relation':'exact all384state H equality after kappa=g=0, not transfer of interacting spectral gaps','pure_C_and_E':'retained in JOINT_SOURCE.json references.global_allN.native_instruments as controls','C9':'bare/character results are separate controls; no pure C9 gap imported'},'scope':'One finite Hamiltonian/reference/hopclock/native instrument plus explicit boundarydressed reservoir apparatus. Occupation/record/current/energy readouts are declared complete mathematical apparatus. No physicalspecies, spin, vacuum, GeV mass, Maxwell radiation, Lorentz or gravity identification.','full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):digest(p) for p in sources}}
    result['spectral_readout_instrument']=spectral_result
    result['checks']['grouped_Luders_instrument_and_continuation']=all(spectral_result['checks'].values())
    return result,{'selected_protocol_poststate':post,**spectral_arrays}


def contract():
    targets={'hamiltonian_and_all_number_spectra':('JOINT_SOURCE.json','actual_all_number_spectra'),'reference':('JOINT_SOURCE.json','references.global_allN'),'boundary_dressed_addition_removal':('JOINT_SOURCE.json','references.global_allN.complete_reservoir_response'),'same_reference_C_E_MIX_work':('JOINT_SOURCE.json','references.global_allN.native_instruments'),'charged_pulse_field_energy':('JOINT_SOURCE.json','references.global_allN.field_energy_before_reservoir / field_energy_after_reservoir'),'selected_protocol_and_readouts':('CANONICAL_CONSTRUCTION.json','readouts'),'baseline_controls':('CANONICAL_CONSTRUCTION.json','controls'),'complete_grouped_Luders_states_and_continuation':('CANONICAL_CONSTRUCTION.json','spectral_readout_instrument')}
    files=[HERE/n for n in ['canonical_construction.py','spectral_readout.py','joint_source.py','charged_reservoir.py','modular_field_adapter.py','JOINT_SOURCE.json','joint_source_states.npz','CANONICAL_CONSTRUCTION.json','canonical_states.npz']]
    return {'schema':'oph.c10.canonical_source_contract.v1','construction_id':'triangle-CAR-allN-Z6-periodic-half-record-native-MIX-v1','one_source':{'states':'all six-orbital CAR occupations with six shared Z6half-records;384states','H':'same actual CARHm+periodic-half-record potential+cyclic face term','reference':'global interacting all-number ground projector; unique N3 in this case','native_clock':'dimensionless half-hopparameter','native_instrument':'selectedMIX on actual fullsource; C/E standalonecontrols','reservoir':'same fullsource, complete7modeCARswap, integerB=3−N; zero boundaryextraenergy, .23nr'},'input_ancestry':{'CAR_statistics_and_phases':'C9configuration-level CAR extension with explicit transport phases; not derivedphysicalspin','Gauss_preparation':'existing affine modularfiber/tree arithmetic; sourcebackground111000 fixed','field':'L05 half-record quadratic normalform adapted to explicitly periodicrepresentatives; faceW samecycle; no hardwallidentitytransfer','reservoir':'existing mass_instrument CARannihilator/swapmethod with explicitboundarydressing','ideal_preparations':'global spectralreference and fair reservoir occupation supplied','numerical_execution':'finite float64 dense/source-number block algebra; no interval-certified trajectories'},'completed_observables':{name:{'file':str((HERE/file).relative_to(BASE)),'json_path':key} for name,(file,key) in targets.items()},'boundaries':['fixedGaussbackground not automatically chargeconjugationinvariant','zero boundaryenergy is supplied ledgerchoice, not completephysicalpreparationcost','full384source retained; no uniform/character compression after periodicHE','allnumberPVM uses actualoccupation labels; no unjustified genericmodularintegerdecoder','reference/clock/readouts lack physicalspecies/field/spacetimeidentification','separate component experiments are controls, not measurements transferred into this interacting construction'],'full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):digest(p) for p in files}}


def compare(a,b,path=()):
    if isinstance(a,dict):return isinstance(b,dict) and a.keys()==b.keys() and all(compare(a[k],b[k],path+(k,)) for k in a)
    if isinstance(a,list):return isinstance(b,list) and len(a)==len(b) and all(compare(x,y,path+(str(i),)) for i,(x,y) in enumerate(zip(a,b)))
    if type(a) is bool:return type(b) is bool and a==b
    if type(a) is int:return type(b) is int and a==b
    if type(a) is float:
        if type(b) is not float:return False
        measured=bool(path and path[0] in ('readouts','spectral_readout_instrument'))
        declared=bool(path and path[-1] in ('native_time','native_probability','energy_group_tolerance'))
        return bool(np.isclose(a,b,atol=1e-10,rtol=1e-10)) if measured and not declared else a==b
    return type(a) is type(b) and a==b

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);mode=parser.add_mutually_exclusive_group(required=True);mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    result,arrays=execute()
    if not all(result['checks'].values()):raise ValueError('canonical construction failed')
    path=HERE/'CANONICAL_CONSTRUCTION.json';archive=HERE/'canonical_states.npz'
    if args.write:
        np.savez_compressed(archive,**arrays);path.write_text(json.dumps(result,indent=2)+'\n');(HERE/'CANONICAL_SOURCE_CONTRACT.json').write_text(json.dumps(contract(),indent=2)+'\n')
    else:
        if not compare(result,json.loads(path.read_text())):raise ValueError('canonical scientific replay differs')
        saved=np.load(archive)
        if any(not np.allclose(value,saved[key],atol=1e-10,rtol=1e-10) for key,value in arrays.items()):raise ValueError('canonical archived poststate differs')
        if not compare(contract(),json.loads((HERE/'CANONICAL_SOURCE_CONTRACT.json').read_text())):raise ValueError('canonical source contract drift')
    print('PASS: one canonical interacting source, selected native protocol, complete readouts and source contract')
