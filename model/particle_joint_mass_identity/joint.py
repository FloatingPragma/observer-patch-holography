"""Commuting source observables, complete probes and finite resolution."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json
import math

import numpy as np
from scipy import linalg as la
import sympy as s

from model.common_source import calibrated_matter as matter
from model.common_source.chiral_candidate import action
from model.conditional_particle_spectrum import source,readout
from model.coherent_defect_wave.mass_instrument import model as apparatus
from model.lepton_mass_compiler import analyzer,compiler
from model.native_flavor_orbit import intervals as iv
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    paths.update(('model/conditional_particle_spectrum/results.json','model/conditional_particle_spectrum/PROOF.md',
                  'model/lepton_mass_compiler/results.json','model/lepton_mass_compiler/PROOF.md',
                  'model/lepton_cycle/experimental_inputs.json','requirements.txt','model/requirements.txt'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


def authenticate():
    for name in ('conditional_particle_spectrum','lepton_mass_compiler'):
        value=json.loads((ROOT/'model'/name/'results.json').read_text())
        for path,digest in value.get('source_bindings',value.get('input_sha256',{})).items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('changed source: '+path)


def rational(value):
    value=s.simplify(value)
    if value.is_Rational is not True:raise ValueError('rational square required')
    return F(int(value.p),int(value.q))


def mass_box(mass):
    return iv.sqrt_rational(rational(mass*mass),digits=76)


def energy_box(mass):
    return iv.sqrt_rational(1+rational(mass*mass),digits=76)


def measured_groups(modes,values,masses):
    groups={}
    for index,(mode,observables) in enumerate(zip(modes,values,strict=True)):
        mass=masses[mode['species']];energy=s.sqrt(1+mass*mass)
        key=tuple(observables)+(str(energy),)
        if key not in groups:
            groups[key]={'charge':observables[0],'color_C2':observables[1],
                'helicity_twice':observables[2],'momentum_z_GeV':observables[3],
                'energy_GeV_exact':str(energy),'energy_GeV_interval':energy_box(mass).strings(),
                'rest_coefficient_GeV_interval':mass_box(mass).strings(),
                'mass_squared_GeV2_exact':str(s.simplify(mass*mass)),
                'mode_indices':[],'candidate_species':set(),'orientations':set()}
        row=groups[key];row['mode_indices'].append(index)
        row['candidate_species'].add(mode['species']);row['orientations'].add(mode['orientation'])
    result=[]
    for row in groups.values():
        row['candidate_species']=sorted(row['candidate_species']);row['orientations']=sorted(row['orientations'])
        row['projector_rank']=len(row['mode_indices']);result.append(row)
    return result


@lru_cache(maxsize=1)
def fermion_joint():
    f=source.fermions();frame=f['named_positive_frame'];h=f['native_pair_h'];d=f['native_pair_delta']
    hn=np.block([[h,d],[d.conj().T,-h.T]])
    charge=source.numeric(action.electric_charge());color=source.numeric(action.identification_operators()['color_C2'])
    qcar=np.kron(np.eye(4),charge);ccar=np.kron(np.eye(4),color)
    szcar=np.kron(np.eye(2),np.kron(np.diag([.5,-.5]),np.eye(45)))
    pcar=np.diag(np.r_[np.ones(90),-np.ones(90)])
    operators={'charge':la.block_diag(qcar,-qcar.T),'color_C2':la.block_diag(ccar,ccar.T),
               'spin_z':la.block_diag(szcar,-szcar.T),'momentum_z':la.block_diag(pcar,-pcar.T),
               'energy':hn}
    operators['helicity']=operators['spin_z']@operators['momentum_z']
    values=[]
    for mode in f['modes']:
        spin=mode.get('spin_twice')
        helicity=mode['helicity_twice'] if spin is None else spin*mode['momentum_GeV'][2]
        values.append((mode['electric_charge'],'4/3' if mode['species'] in ('up','charm','top','down','strange','bottom') else '0',
                       helicity,mode['momentum_GeV'][2]))
    errors={}
    for i,(name,operator) in enumerate(operators.items()):
        errors['Hermitian_'+name]=float(la.norm(operator-operator.conj().T))
        for other,second in list(operators.items())[:i]:
            errors['commutator_'+name+'_'+other]=float(la.norm(operator@second-second@operator))
    for name,column in (('charge',0),('color_C2',1),('momentum_z',3),('helicity',2)):
        diagonal=np.array([float(F(row[column]))/(2 if name=='helicity' else 1) for row in values])
        errors['actual_joint_eigenvector_'+name]=float(la.norm(operators[name]@frame-frame*diagonal))
    identity=frame@frame.conj().T
    tau=np.block([[np.zeros((180,180)),np.eye(180)],[np.eye(180),np.zeros((180,180))]])
    errors['complete_Nambu_positive_and_negative']=float(la.norm(identity+tau@identity.conj()@tau-np.eye(360)))
    errors['positive_projector_isometry']=float(la.norm(frame.conj().T@frame-np.eye(180)))
    if max(errors.values())>3e-10:raise ArithmeticError('actual commuting fermion observable failed')
    # Exact identities controlling the unresolved tiny-mass numerical block.
    mass, p=s.symbols('m p',real=True);block=s.Matrix([[p,mass],[mass,-p]])
    if block*block!=(p*p+mass*mass)*s.eye(2):raise ArithmeticError('exact source mass pair failed')
    color_exact=action.identification_operators()['color_C2'];q_exact=action.electric_charge()
    if not matter._zero(color_exact*f['K']-f['K']*color_exact) or not matter._zero(q_exact.T*f['K']+f['K']*q_exact):
        raise ArithmeticError('exact source charge or color compatibility failed')
    masses={row['species']:s.sympify(row['rest_coefficient_GeV']['exact']) for row in f['targets']}
    groups=measured_groups(f['modes'],values,masses)
    if len(groups)!=84 or sum(g['projector_rank'] for g in groups)!=180:
        raise ArithmeticError('complete fermion joint eigenspaces differ')
    return {'groups':groups,'targets':f['targets'],'mode_count':180,'Nambu_dimension':360,
            'errors':errors,'exact_mass_pair_identity':True,'exact_source_charge_and_color_identities':True,
            'ideal_energy_eigenspaces_from_exact_native_mass_frame':True,
            'binary64_neutrino_energy_splittings_resolved':False,
            'all_negative_Nambu_partners_retained_without_species_doubling':True}


@lru_cache(maxsize=1)
def adjoint_color_control():
    """Exact normalized SU3 adjoint Casimir on the full eight-color fiber."""
    z=s.S.Zero;i=s.I;r=s.sqrt(3)
    raw=[s.Matrix([[0,1,0],[1,0,0],[0,0,0]]),s.Matrix([[0,-i,0],[i,0,0],[0,0,0]]),
         s.diag(1,-1,0),s.Matrix([[0,0,1],[0,0,0],[1,0,0]]),s.Matrix([[0,0,-i],[0,0,0],[i,0,0]]),
         s.Matrix([[0,0,0],[0,0,1],[0,1,0]]),s.Matrix([[0,0,0],[0,0,-i],[0,i,0]]),s.diag(1,1,-2)/r]
    t=[x/2 for x in raw]
    adj=[s.Matrix(8,8,lambda c,b:s.simplify(2*s.trace(t[c]*(a*t[b]-t[b]*a)))) for a in t]
    casimir=sum((a*a for a in adj),s.zeros(8)).applyfunc(s.simplify)
    if casimir!=3*s.eye(8) or any(a.H!=a for a in adj):raise ArithmeticError('adjoint color Casimir failed')
    return {'fiber_dimension':8,'generator_trace_normalization':'Tr_fund(T_a T_b)=delta_ab/2',
            'quadratic_Casimir_exact':'3 I8','canonical_basis_role':'declared orthonormal adjoint coordinates of the existing eight equal-frequency color gauge modes',
            'independent_color_component_observable_measured':False}


@lru_cache(maxsize=1)
def boson_joint():
    b=source.bosons();adjoint_color_control();n=len(b['modes']);frame=np.zeros((n,n),complex)
    spin_z=np.zeros((n,n),complex);q=np.diag([float(F(m['electric_charge'])) for m in b['modes']])
    p=np.diag([m['momentum_GeV'][2] for m in b['modes']]);h=np.diag([m['energy_GeV']['numeric'] for m in b['modes']])
    c=np.diag([3 if m['species']=='gluon' else 0 for m in b['modes']]);values=[None]*n
    for group in b['groups']:
        bycolor=defaultdict(list)
        for j in group['modes']:bycolor[b['modes'][j]['color_component']].append(j)
        for indices in bycolor.values():
            size=len(indices)
            if size==1:
                rotation=np.eye(1);spin=np.zeros((1,1));eigen=[0]
            elif size==2:
                rotation=np.array([[1,1],[-1j,1j]])/np.sqrt(2);spin=np.array([[0,-1j],[1j,0]]);eigen=[-1,1]
            elif size==3:
                rotation=np.array([[1/np.sqrt(2),0,1/np.sqrt(2)],[-1j/np.sqrt(2),0,1j/np.sqrt(2)],[0,1,0]])
                spin=np.array([[0,-1j,0],[1j,0,0],[0,0,0]]);eigen=[-1,0,1]
            else:raise ArithmeticError('wrong polarization count')
            frame[np.ix_(indices,indices)]=rotation;spin_z[np.ix_(indices,indices)]=spin
            for index,sz in zip(indices,eigen,strict=True):
                mode=b['modes'][index]
                values[index]=(mode['electric_charge'],'3' if mode['species']=='gluon' else '0',2*sz*mode['momentum_GeV'][2],mode['momentum_GeV'][2])
    helicity=spin_z@p;errors={'polarization_frame_unitarity':float(la.norm(frame.conj().T@frame-np.eye(n)))}
    ops={'charge':q,'color_C2':c,'helicity':helicity,'momentum_z':p,'energy':h}
    for i,(name,operator) in enumerate(ops.items()):
        for other,second in list(ops.items())[:i]:errors['commutator_'+name+'_'+other]=float(la.norm(operator@second-second@operator))
    for name,column in (('charge',0),('color_C2',1),('momentum_z',3),('helicity',2)):
        diagonal=np.array([float(F(row[column]))/(2 if name=='helicity' else 1) for row in values])
        errors['joint_eigenvector_'+name]=float(la.norm(ops[name]@frame-frame*diagonal))
    if max(errors.values())>2e-12:raise ArithmeticError('boson joint polarizations failed')
    groups=measured_groups(b['modes'],values,b['masses'])
    if len(groups)!=28 or sum(g['projector_rank'] for g in groups)!=56:raise ArithmeticError('boson joint rank differs')
    return {'groups':groups,'targets':b['targets'],'mode_count':56,'errors':errors,
            'adjoint_color':adjoint_color_control(),
            'same_canonical_gauge_mass_and_scalar_Hessian':b['checks'],
            'massless_longitudinal_modes_added':False,'helicity_is_defined_at_both_nonzero_momenta':True}


def coarse_classes(groups,half_width=None):
    classes=defaultdict(list)
    for index,row in enumerate(groups):
        key=(row['charge'],row['color_C2'],row['helicity_twice'],row['momentum_z_GeV'])
        if half_width is not None:
            width=2*F(half_width);lo,hi=map(F,row['energy_GeV_interval'])
            a,b=lo//width,hi//width
            if a!=b:raise ArithmeticError('energy enclosure crosses declared bin boundary')
            key+=(int(a),)
        classes[key].append(index)
    return [{'quantum_numbers':list(key),'joint_rows':indices,
             'candidate_species':sorted({species for i in indices for species in groups[i]['candidate_species']}),
             'rank':sum(groups[i]['projector_rank'] for i in indices)} for key,indices in classes.items()]


@lru_cache(maxsize=1)
def complete_joint_records():
    f,b=fermion_joint(),boson_joint();groups=[]
    for family,data in (('fermion',f),('boson',b)):
        for row in data['groups']:groups.append(dict(row,source_sector=family,joint_index=len(groups)))
    sectors=[{'joint_index':-1,'kind':'vacuum','rank_exact':'1'}]
    sectors +=[{'joint_index':row['joint_index'],'kind':'one_excitation','rank_exact':str(row['projector_rank'])} for row in groups]
    dimension=(1<<180)*1653
    sectors.append({'joint_index':112,'kind':'complete_multi_excitation_complement','rank_exact':str(dimension-237)})
    pointers=('fermion_create','fermion_remove','boson_create','boson_remove','boson_no_operation')
    rows=[]
    for pointer in pointers:
        for sector in sectors:
            probability=F(0);index=sector['joint_index']
            if sector['kind']=='one_excitation':
                group=groups[index]
                if pointer==group['source_sector']+'_create':probability=F(group['projector_rank'],360 if group['source_sector']=='fermion' else 116)
            elif pointer=='boson_no_operation' and sector['kind']=='vacuum':probability=F(1,58)
            rows.append({'pointer':pointer,'joint_index':index,'projector_kind':sector['kind'],'probability_exact':str(probability)})
    if sum(F(r['probability_exact']) for r in rows)!=1 or len(rows)!=570:
        raise ArithmeticError('complete unconditioned joint law differs')
    states,creators,no_op,cap_error,closure_error=readout.boson_operators()
    # Direct source columns of every capped creator, keeping the quantum label.
    errors=[]
    for j,a in enumerate(creators):
        column=a[:,0].toarray().ravel()/np.sqrt(116)
        expected=np.zeros(len(states));expected[j+1]=1/np.sqrt(116)
        errors.append(float(la.norm(column-expected)))
    if max(errors)>1e-14 or not math.isclose(no_op[0]/2,1/58,abs_tol=1e-16):
        raise ArithmeticError('actual capped-vacuum probe differs')
    no_energy=coarse_classes(groups);resolution=coarse_classes(groups,F(1,10**6))
    fine=coarse_classes(groups,F(1,10**26))
    color_error=sum(2*F(row['projector_rank']-1,360 if row['source_sector']=='fermion' else 116)
                    for row in groups if row['projector_rank']>1)
    return {'source_Hilbert_dimension_exact':str(dimension),'source_preparation':'same native Bogoliubov vacuum tensor same capped oscillator vacuum',
            'source_reference_and_clock':'the exact parent quadratic source, quasiparticle vacuum, dual charge/color reference and declared GeV-conjugate evolution parameter',
            'joint_groups':groups,'complete_projectors':sectors,'complete_records':rows,
            'positive_records':sum(F(row['probability_exact'])>0 for row in rows),
            'zero_records':sum(F(row['probability_exact'])==0 for row in rows),
            'probability_sum':'1','all_input_instrument_identity':'one-half of the full CAR instrument plus one-half of the complete capped boson instrument, followed by a complete orthogonal projector family',
            'probe_pointers_do_not_contain_species_or_family_labels':True,
            'reference_dimension':474,'reference_includes_all_original_mode_coordinates':True,
            'joint_color_components_unresolved_and_coherent':True,
            'nonvacuum_many_body_inputs_retained_in_complete_complement':True,
            'all_fermions_filled_control':{'pointer_probabilities':['0','1/2','14/29','0','1/58'],'all_outputs_in_multi_excitation_complement':True},
            'source_execution_errors':{'capped_source_columns':max(errors),'capped_all1653_column_identity':closure_error,'capped_CAR_sum_identity':cap_error},
            'without_energy_readout':{'class_count':len(no_energy),'classes':no_energy},
            'finite_energy_binning':{'half_bin_width_GeV':'1/1000000','class_count':len(resolution),'classes':resolution,
                'neutrino_mass_species_resolved':False,'binning_is_declared_ideal_detector_service':True},
            'fine_energy_binning':{'half_bin_width_GeV':'1/100000000000000000000000000','class_count':len(fine),
                'ideal_exact_source_species_resolved':True,'hardware_resolution_certified':False},
            'negative_control_color_reference_dephasing':{'unhalved_complete_output_trace_distance_exact':str(color_error),
                'original_law_preserved':False}}


@lru_cache(maxsize=1)
def compiled_lepton_joint():
    a=[apparatus.annihilator(i,4) for i in range(4)]
    beta=s.kronecker_product(s.Matrix([[0,1],[1,0]]),s.eye(2))
    h=sum((beta[i,j]*a[i].H*a[j] for i in range(4) for j in range(4)),s.zeros(16))+2*s.eye(16)
    q=s.diag(*(2-word.bit_count() for word in range(16)))
    spin=s.diag(*(sum((1 if j%2==0 else -1)*s.Rational(1,2) for j in range(4) if word>>j&1) for word in range(16)))
    if h*q!=q*h or h*spin!=spin*h:raise ArithmeticError('compiled rest observables do not commute')
    negative=(s.eye(4)-beta)/2;rho=s.zeros(16)
    for left in combinations(range(4),2):
        for right in combinations(range(4),2):rho[sum(1<<j for j in left),sum(1<<j for j in right)]=negative.extract(left,right).det()
    effects=[]
    for energy,projector in enumerate(apparatus.energy_projectors(h,s.Integer(1))):
        for charge in range(-2,3):
            for twice_spin in range(-2,3):
                diagonal=s.diag(*(int(q[j,j]==charge and 2*spin[j,j]==twice_spin) for j in range(16)))
                effect=projector*diagonal
                if effect!=s.zeros(16):effects.append((energy,charge,twice_spin,effect))
    if sum((p for *_,p in effects),s.zeros(16))!=s.eye(16) or len(effects)!=15:
        raise ArithmeticError('complete spin-charge-energy rest effects differ')
    rows=[];coarse=defaultdict(lambda:s.S.Zero)
    for initial in (1,0):
        for flag,state in enumerate(apparatus.reservoir_branches(rho,initial)):
            for energy,charge,twice_spin,effect in effects:
                probability=s.trace(effect*state)/2
                if probability.is_nonnegative is not True:raise ArithmeticError('negative rest record')
                row={'setting':'add' if initial else 'remove','reservoir_out':flag,'native_electric_charge':charge,
                     'spin_z_twice':twice_spin,'source_energy_in_same_m':energy,'probability_exact':str(probability)}
                rows.append(row);coarse[(row['setting'],flag,charge,energy)]+=probability
    old=compiler.exact_rest_instrument()['complete_records']
    if any(coarse[(r['setting'],r['reservoir_out'],r['native_electric_charge'],r['energy_in_units_of_same_m'])]!=s.Rational(r['probability_exact']) for r in old):
        raise ArithmeticError('spin refinement changed an original mass record')
    precision=compiler.ramsey_readout();intervals=analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2']
    ratio_rows=precision['ratio_bounds'];anchor=precision['electron_anchor_GeV_exact']
    families=[]
    for name,k in (('electron',1),('muon',2),('tau',0)):
        box=[anchor,anchor] if name=='electron' else next(r['electron_anchored_GeV_interval'] for r in ratio_rows if r['label']==name)
        families.append({'species':name,'character':k,'source_rest_energy_interval':intervals[k],
                         'conditional_GeV_interval':box,'role':'empirical_anchor' if name=='electron' else 'compiled_free_source_ratio',
                         'spin_charge_energy_records_per_family':60})
    return {'source_Fock_dimension_per_family':16,'with_original_reservoir_dimension':32,
            'source_color_C2_exact':'0','source_momentum_exact':[0,0,0],
            'spin_axis_defined':True,'helicity_at_rest_claimed':False,
            'joint_projector_count':15,'joint_projector_ranks':[int(s.trace(p)) for *_,p in effects],
            'complete_records_per_family':rows,'complete_record_count_per_family':60,
            'original36_mass_records_recovered_by_forgetting_spin':True,
            'families':families,'empirical_anchor_GeV_exact':anchor,
            'finite_compiler_unhalved_law_error_exact':precision['compiler']['complete_unhalved_trace_error_exact'],
            'finite_compiler_rounds_executed':False,'empirical_anchor_uncertainty_included':False,
            'source_phase_and_family_dictionary_retain_parent_ancestry':True,
            'observed_interacting_lepton_identification_proved':False}


def calculate():
    pins=bindings();authenticate();f=fermion_joint();b=boson_joint();joint=complete_joint_records()
    targets=[]
    for row in f['targets']+b['targets']:
        mass=s.sympify(row['rest_coefficient_GeV']['exact'])
        targets.append({'species':row['species'],'conditional_rest_coefficient_GeV_interval':mass_box(mass).strings(),
            'mass_squared_GeV2_exact':str(s.simplify(mass*mass)),
            'input_convention':row['mass_convention'],'physical_matching':row['missing_physical_attachment'],
            'empirical_inputs':row.get('empirical_inputs',['matched gauge/scalar parameters and supplied running Higgs background']),
            'native_prediction':False})
    if len(targets)!=17:raise ArithmeticError('target inventory changed')
    result={'schema':'oph.particle_joint_mass_identity.v1','source_bindings':pins,'all_checks_pass':True,
            'calibrated_fermion_joint':f,'calibrated_boson_joint':b,'complete_calibrated_joint_instrument':joint,
            'conditional_GeV_targets':targets,'compiled_lepton_joint':compiled_lepton_joint(),
            'scope':{'actual_same_source_joint_observable_constructed':True,
                'all_preparation_probe_failure_and_complement_records_retained':True,
                'arbitrary_untouched_reference_instrument_identity':True,
                'physical_Standard_Model_particle_dictionary_derived':False,
                'calibrated_masses_are_independent_predictions':False,
                'full_many_body_color_and_energy_spectrum_measured_in_complement':False,
                'sterile_neutrino_or_extra_particle_species_inferred':False,
                'finite_bin_hardware_resolution_certified':False,
                'compiled_lepton_and_calibrated_action_mass_conventions_merged':False}}
    if bindings()!=pins:raise ArithmeticError('joint-readout source changed during execution')
    return result


def equal(a,b):
    if isinstance(b,dict):return isinstance(a,dict) and set(a)==set(b) and all(equal(a[k],v) for k,v in b.items())
    if isinstance(b,list):return isinstance(a,list) and len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b,strict=True))
    if isinstance(b,float):return type(a) in (int,float) and math.isfinite(a) and abs(a-b)<=2e-12+2e-9*abs(b)
    return type(a) is type(b) and a==b


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=p.parse_args()
    result=calculate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    elif not equal(json.loads(path.read_text()),result):raise ArithmeticError('joint particle receipt differs')
    print('Native joint quantum-number, mass, reference and complete-outcome records verified.')
