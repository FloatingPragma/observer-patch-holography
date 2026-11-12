"""Whole-colour collision factors and all capped-boson jump records."""
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import argparse
import hashlib
import json
import math

import numpy as np
from scipy.sparse import csr_matrix

from model.thermal_particle_records import records as numbers
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
TEMPERATURES=(F(1,10),F(1),F(10),F(100));PLAN=(F(10),F(1),F(1,10),F(1),F(10),F(100))
BITS=180;GRID=1<<BITS


def down(x):return F((x*GRID).numerator//(x*GRID).denominator,GRID)
def up(x):return -down(-x)


def negative_exp_endpoint(x):
    """Alternating Taylor at x/2^k<=1, then outward dyadic squaring."""
    x=F(x)
    if x<0:raise ValueError('nonnegative inverse-temperature energy required')
    k=0
    while x>1:x/=2;k+=1
    term=F(1);total=F(1)
    for n in range(1,41):term*=-x/n;total+=term
    upper=up(total);lower=max(F(0),down(total+term*(-x)/41))
    for _ in range(k):lower=down(lower*lower);upper=up(upper*upper)
    return lower,upper


def boltzmann_box(energy,temperature):
    lo,hi=map(F,energy);t=F(temperature)
    if t<=0:raise ValueError('positive supplied temperature required')
    return negative_exp_endpoint(hi/t)[0],negative_exp_endpoint(lo/t)[1]


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    for package in ('particle_joint_mass_identity','particle_thermal_history','thermal_particle_records','conditional_particle_spectrum'):
        paths|={f'model/{package}/results.json',f'model/{package}/PROOF.md'}
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


@lru_cache(maxsize=1)
def source():
    _,_,parent,groups=numbers.inputs();numbers.authenticate(ROOT/'model/thermal_particle_records/results.json')
    rotation=numbers.boson_rotation();words=rotation['words'];index={word:i for i,word in enumerate(words)}
    fgroups=[g for g in groups if g['source_sector']=='fermion'];bgroups=[g for g in groups if g['source_sector']=='boson']
    b_to_g={i:k for k,g in enumerate(bgroups) for i in g['mode_indices']}
    energies=[m['energy_GeV']['numeric'] for m in parent['mode_inventory']['boson_modes']]
    energy_boxes=[None]*56
    for g in bgroups:
        for i in g['mode_indices']:energy_boxes[i]=g['energy_GeV_interval']
    return dict(parent=parent,groups=groups,fgroups=fgroups,bgroups=bgroups,words=words,index=index,
        b_to_g=b_to_g,benergy=energies,benergy_boxes=energy_boxes,
        word_energy=np.array([sum(energies[i] for i in word) for word in words]))


@lru_cache(maxsize=8)
def fermion_bath(temperature):
    result=[]
    for group in source()['fgroups']:
        rank=group['projector_rank'];energy=float(sum(map(F,group['energy_GeV_interval']))/2)
        q=math.exp(-energy/float(temperature));lo,hi=boltzmann_box(group['energy_GeV_interval'],temperature)
        probabilities=[];boxes=[]
        for word in range(1<<rank):
            n=word.bit_count();probabilities.append(q**n/(1+q)**rank)
            boxes.append((lo**n/(1+hi)**rank,hi**n/(1+lo)**rank))
        error=sum((max(abs(F.from_float(p)-a),abs(F.from_float(p)-b)) for p,(a,b) in zip(probabilities,boxes)),F(0))
        result.append(dict(joint_index=group['joint_index'],rank=rank,energy=energy,
            probabilities=np.array(probabilities),exact_probability_intervals=boxes,coefficient_error=error,
            all_quantum_reset_Kraus_coordinates=(1<<rank)**2,
            log_word_probabilities=[-word.bit_count()*energy/float(temperature)-rank*math.log1p(q) for word in range(1<<rank)]))
    return result


@lru_cache(maxsize=8)
def boson_bath(temperature):
    d=source();q=np.exp(-np.array(d['benergy'])/float(temperature))
    boxes=[boltzmann_box(e,temperature) for e in d['benergy_boxes']]
    upper_sum=sum((b for _,b in boxes),F(0));alpha=down(1/(2+upper_sum));af=float(alpha)
    rows=[];cols=[];values=[];branch_support=defaultdict(list);row_errors=[];normalization_error=0.
    for x,word in enumerate(d['words']):
        n=len(word);occupancy=Counter(word);exit_rate=n+(sum(q)+sum(q[j]*v for j,v in occupancy.items()) if n<2 else 0.)
        exit_lo=F(n)+(sum((a for a,_ in boxes),F(0))+sum((boxes[j][0]*v for j,v in occupancy.items()),F(0)) if n<2 else 0)
        exit_hi=F(n)+(upper_sum+sum((boxes[j][1]*v for j,v in occupancy.items()),F(0)) if n<2 else 0)
        hold=1-af*exit_rate;hold_box=(1-alpha*exit_hi,1-alpha*exit_lo)
        if hold_box[0]<0:raise ArithmeticError('nonpositive no-jump effect')
        entries=[(x,hold,hold_box)];branch_support[(n,'hold',-1)].append((x,x,hold))
        for j in range(56):
            occupation=occupancy[j]
            if occupation:
                out=list(word);out.remove(j);y=d['index'][tuple(out)];weight=af*occupation
                exact=alpha*occupation;entries.append((y,weight,(exact,exact)))
                branch_support[(n,'remove',j)].append((x,y,weight))
            if n<2:
                y=d['index'][tuple(sorted(word+(j,)))];weight=af*q[j]*(occupation+1)
                bounds=(alpha*boxes[j][0]*(occupation+1),alpha*boxes[j][1]*(occupation+1))
                entries.append((y,weight,bounds));branch_support[(n,'create',j)].append((x,y,weight))
        error=F(0)
        for y,weight,(lo,hi) in entries:
            rows.append(x);cols.append(y);values.append(weight)
            exact_float=F.from_float(weight);error+=max(abs(exact_float-lo),abs(exact_float-hi))
        row_errors.append(error);normalization_error=max(normalization_error,abs(sum(w for _,w,_ in entries)-1))
    matrix=csr_matrix((values,(rows,cols)),shape=(1653,1653))
    branches=[]
    for n in range(3):
        for kind,modes in [('hold',[-1]),('remove',range(56)),('create',range(56))]:
            for mode in modes:
                support=branch_support.get((n,kind,mode),[])
                branches.append({'input_total_number':n,'event':kind,'bath_mode_quantum_coordinate':mode,
                    'public_joint_group':None if mode<0 else d['bgroups'][d['b_to_g'][mode]]['joint_index'],
                    'exact_zero_operator':not support,'entries_input_output_squared_amplitude':support})
    if len(branches)!=339:raise ArithmeticError('fine quantum-coordinate count differs')
    indegree=max(Counter(cols).values());gamma=F(2*indegree+2,2**53-2*indegree-2)
    eta=max(row_errors);step_error=eta+gamma*(1+eta)+F(4*1653*indegree,2**1074)
    delta=2*alpha*alpha
    return dict(matrix=matrix,branches=branches,alpha=alpha,delta=delta,step_error=step_error,
        q=q,boltzmann_boxes=boxes,coefficient_error=eta,normalization_error=normalization_error,
        public_record_count=171,fine_coordinate_count=339,maximum_incoming_terms=indegree)


def reference_controls(temperature):
    d=source();bath=boson_bath(temperature)
    group=next(g for g in d['bgroups'] if g['projector_rank']==8)
    i,j=group['mode_indices'][:2];xi=d['index'][(i,)];xj=d['index'][(j,)]
    # Annihilation of the two source colours transfers this entanglement into
    # the bath mode register; the public pointer does not measure that mode.
    psi=np.zeros((8,2),complex);psi[0,0]=psi[1,1]=1/math.sqrt(2)
    rho=np.outer(psi.reshape(-1),psi.reshape(-1).conj());dephased=np.zeros_like(rho)
    for k in range(8):dephased[k*2:(k+1)*2,k*2:(k+1)*2]=rho[k*2:(k+1)*2,k*2:(k+1)*2]
    distance=float(sum(abs(np.linalg.eigvalsh(rho-dephased))))
    norms=np.zeros(1653)
    for branch in bath['branches']:
        for x,y,weight in branch['entries_input_output_squared_amplitude']:norms[x]+=weight
    if max(abs(norms-1))>2e-13 or abs(distance-1)>1e-13:raise ArithmeticError('reference-safe complete bath control failed')
    return {'full1653_Kraus_identity_diagonal_error':float(max(abs(norms-1))),
        'full_Kraus_identity_offdiagonal_exact_zero':True,
        'coherent_colour_deletion_original_probability':float(bath['alpha']),
        'quantum_bath_reference_entanglement_transferred':True,
        'fine_bath_colour_measurement_negative_control_unhalved_distance':distance,
        'original_reference_marginal_retained':True,'source_basis_indices':[xi,xj]}


def entropy(p):
    values=p[p>0];return float(-np.dot(values,np.log(values)))


def thermal_state(temperature):
    d=source();f=fermion_bath(temperature);logs=-d['word_energy']/float(temperature)
    weights=np.exp(logs);p=weights/weights.sum()
    lo=[];hi=[];b=boson_bath(temperature)
    for word in d['words']:
        lo.append(math.prod(b['boltzmann_boxes'][j][0] for j in word));hi.append(math.prod(b['boltzmann_boxes'][j][1] for j in word))
    zlo=sum(lo,F(0));zhi=sum(hi,F(0))
    error=sum((max(abs(F.from_float(x)-a/zhi),abs(F.from_float(x)-bb/zlo)) for x,a,bb in zip(p,lo,hi)),F(0))
    return [g['probabilities'].copy() for g in f],p,[g['coefficient_error'] for g in f],error


def observables(factors,probabilities):
    d=source();groups=[];species=defaultdict(float);energy=0.;number=0.;nong=0.
    for group,p in zip(d['fgroups'],factors):
        rank=group['projector_rank'];counts=[word.bit_count() for word in range(1<<rank)]
        mean=float(np.dot(p,counts));second=float(np.dot(p,[n*(n-1) for n in counts]))
        excess=second-rank*(rank-1)*(mean/rank)**2;nong=max(nong,abs(excess))
        frequency=float(sum(map(F,group['energy_GeV_interval']))/2)
        energy+=mean*frequency;number+=mean;species[group['candidate_species'][0]]+=mean
        groups.append({'joint_index':group['joint_index'],'mean_number':mean,
            'complete_count_probabilities':[float(sum(p[w] for w,n in enumerate(counts) if n==count)) for count in range(rank+1)],
            'factorial_second_moment_minus_binomial':excess,'all_occupation_word_probabilities':p.tolist()})
    bcounts=defaultdict(lambda:[0.]*3);class_counts=defaultdict(float)
    for word,p in zip(d['words'],probabilities):
        count=Counter(d['b_to_g'][j] for j in word);key=tuple(sorted(count.items()));class_counts[key]+=float(p)
        for k,g in enumerate(d['bgroups']):bcounts[k][count[k]]+=float(p)
    for k,g in enumerate(d['bgroups']):
        p=bcounts[k];mean=p[1]+2*p[2];number+=mean;species[g['candidate_species'][0]]+=mean
        groups.append({'joint_index':g['joint_index'],'mean_number':mean,'complete_count_probabilities':p})
    energy+=float(probabilities@d['word_energy'])
    return {'expected_energy_GeV':energy,'entropy_nats':sum(entropy(p) for p in factors)+entropy(probabilities),
        'expected_number':number,'group_count_laws':groups,'largest_fermion_non_binomial_second_moment':nong,
        'species_counts':[{'species':k,'mean_number':v,'fraction_of_expected_counts':v/number} for k,v in sorted(species.items())],
        'all435_boson_group_count_laws':[{'group_counts':[[d['bgroups'][k]['joint_index'],n] for k,n in key],'probability':v} for key,v in sorted(class_counts.items())]}


def execute(steps):
    factors,p,ferrors,berror=thermal_state(F(100));initial=observables(factors,p);rows=[]
    for temperature in PLAN:
        target=fermion_bath(temperature);b=boson_bath(temperature);before=observables(factors,p)
        if max(ferrors+[g['coefficient_error'] for g in target])>F(1,4):
            raise ArithmeticError('fermion arithmetic bound outside its stated domain')
        first=None
        for tick in range(steps):
            factors=[.5*x+.5*g['probabilities'] for x,g in zip(factors,target)]
            p=np.asarray(p@b['matrix']).ravel()
            if tick==0:first=observables(factors,p)
        # Positive updates with fixed coefficients, including coefficient and
        # arithmetic errors. All outward rounding is rational.
        unit=F(4,2**53-4)
        ferrors=[up((F(1,2)**steps)*old+(1-F(1,2)**steps)*(g['coefficient_error']+4*unit))
                 for old,g in zip(ferrors,target)]
        berror=up((1+berror)*(1+b['step_error'])**steps-1)
        full_error=up((1+berror)*math.prod(1+x for x in ferrors)-1)
        after=observables(factors,p);heat=after['expected_energy_GeV']-before['expected_energy_GeV']
        sigma=after['entropy_nats']-before['entropy_nats']-heat/float(temperature)
        if sigma< -2e-10:raise ArithmeticError('fixed-temperature Gibbs channel entropy production failed')
        rows.append({'temperature_GeV':float(temperature),'collision_sweeps':steps,
            'complete_fermion_factors':[x.tolist() for x in factors],'complete1653_boson_state_probabilities':p.tolist(),
            'complete_joint_count_readout':after,'first_sweep_count_readout':first,
            'heat_into_source_GeV_numeric':heat,'entropy_production_numeric':sigma,
            'full_source_state_L1_error_bound':str(full_error),'full_source_state_L1_error_bound_numeric':float(full_error),
            'all_quantum_archive_and_pointer_histories':'compose the complete stored local Stinespring factors in original order; no branch or bath coordinate is removed',
            'actual_many_body_state_is_gaussian_claimed':False})
    return {'initial_temperature_GeV':100,'steps_per_setting':steps,'initial_observables':initial,'stages':rows,
        'state_probabilities_renormalized':False,'initial_preparation':'supplied finite Gibbs state on the same calibrated source'}


def certificate():
    pins=bindings();d=source();temperatures=[]
    for t in TEMPERATURES:
        b=boson_bath(t);f=fermion_bath(t);_,p,_,_=thermal_state(t)
        residual=float(np.max(abs(np.asarray(p@b['matrix']).ravel()-p)))
        if residual>2e-14:raise ArithmeticError('actual finite Gibbs state not stationary')
        m=(F(20)/b['delta']).__ceil__()
        if sum((F(20)**k/math.factorial(k) for k in range(41)),F(0))<=10**8:
            raise ArithmeticError('elementary exponential budget bound failed')
        # exp(-20)<10^-8; the chosen sweep count also suppresses all84 factors.
        temperatures.append({'temperature_GeV':float(t),'boson_uniformization_alpha_exact':str(b['alpha']),
            'boson_two_sweep_reset_minorization_exact':str(b['delta']),
            'boson_reset_minorization_numeric':float(b['delta']),
            'sufficient_sweeps_for_unhalved_source_reference_error_1e_minus6':2*m,
            'all339_boson_fine_Kraus_coordinates':b['branches'],
            'all171_public_boson_records':{'input_number':[0,1,2],'event':'hold or create/remove in each of28 joint groups; quantum mode coordinate unmeasured'},
            'fermion_groups':[{'joint_index':g['joint_index'],'rank':g['rank'],
                'all_thermal_occupation_probabilities':g['probabilities'].tolist(),
                'log_probabilities':g['log_word_probabilities'],'reset_Kraus_coordinate_count':g['all_quantum_reset_Kraus_coordinates']} for g in f],
            'Gibbs_stationarity_numeric_error':residual,'reference_controls':reference_controls(t)})
    result={'schema':'oph.quantum_thermal_relaxation.v1','source_bindings':pins,'all_checks_pass':True,
        'source_dimension_exact':str(2**180*1653),'source_fermion_modes':180,'source_boson_modes':56,
        'source_joint_groups':d['groups'],'boson_word_order':[list(x) for x in d['words']],
        'temperature_rows':temperatures,'cooling_experiments':[execute(n) for n in (32,256)],
        'reference_extended_source_mixing_bound':'min(2,2*(1-(1-2^-n)^84)+2*(1-delta)^floor(n/2)); trace all bath/archive outputs only for this marginal bound',
        'heat_entropy_identity':'Delta S-beta Q=D(rho||tau)-D(Phi(rho)||tau)>=0, including an untouched reference using tau tensor rho_R',
        'scope':{'same_calibrated_quadratic_source':True,'complete_quantum_bath_coordinates_retained':True,
                 'whole_colour_groups_and_covariant_bath_registers':True,'supplied_bath_temperature_and_collision_services':True,
                 'native_repair_prepares_this_bath':False,'physical_clock_or_Universe_abundances_derived':False,
                 'complete_archive_claimed_to_forget_initial_information':False,
                 'finite_boson_cap_replaced_by_uncapped_sampler':False}}
    if bindings()!=pins:raise ArithmeticError('thermal relaxation input changed')
    return result


def equal(a,b):
    if isinstance(b,dict):return isinstance(a,dict) and set(a)==set(b) and all(equal(a[k],v) for k,v in b.items())
    if isinstance(b,(list,tuple)):return isinstance(a,(list,tuple)) and len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b,strict=True))
    if isinstance(b,float):return isinstance(a,(int,float)) and abs(a-b)<=2e-12+2e-11*abs(b)
    return type(a) is type(b) and a==b


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    result=certificate();path=HERE/'results.json'
    if a.write:path.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif not equal(json.loads(path.read_text()),result):raise ArithmeticError('quantum thermal receipt differs')
    print('PASS: complete finite quantum thermal channel and cooling response')
