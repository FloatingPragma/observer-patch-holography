"""Same CAR mass kernel, exact pair factorization and interval pulse evolution."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import argparse
import hashlib
import json
import math

import numpy as np
import sympy as sp
from scipy.linalg import block_diag

from model.conditional_particle_spectrum import source as quadratic
from model.thermal_particle_records import records as joint
from model.native_effect_synthesis import intervals as iv
from model.quantum_thermal_relaxation.relaxation import negative_exp_endpoint
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
Q=iv.Q;ZERO=(0,0);ONE=(Q,Q)


def add(a,b):return iv.add(a,b)
def neg(a):return iv.neg(a)
def sub(a,b):return add(a,neg(b))
def mul(a,b):return iv.mul(a,b)
def exact(x):return iv.exact(F(x))
def divide(a,b):
    if b[0]<=0:raise ValueError('strictly positive divisor required')
    return mul(a,(Q*Q//b[1],iv.ceildiv(Q*Q,b[0])))
def square(a):
    lo=0 if a[0]<=0<=a[1] else min(a[0]**2,a[1]**2)
    return lo//Q,iv.ceildiv(max(a[0]**2,a[1]**2),Q)
def sqrt(a):
    if a[0]<0:raise ValueError('nonnegative root interval required')
    lo=math.isqrt(a[0]*Q);hi=math.isqrt(a[1]*Q)
    return lo,hi if hi*hi==a[1]*Q else hi+1
def mid(a):return float(F(a[0]+a[1],2*Q))
def encode(a):return [str(F(x,Q)) for x in a]
def cadd(a,b):return add(a[0],b[0]),add(a[1],b[1])
def cmul(a,b):return sub(mul(a[0],b[0]),mul(a[1],b[1])),add(mul(a[0],b[1]),mul(a[1],b[0]))
def conj(a):return a[0],neg(a[1])
def cnorm(a):return add(square(a[0]),square(a[1]))
def cmid(a):return complex(mid(a[0]),mid(a[1]))
def cencode(a):return {'real':encode(a[0]),'imaginary':encode(a[1])}
def matrix_product(a,b):
    return [[cadd(cmul(a[i][0],b[0][j]),cmul(a[i][1],b[1][j])) for j in range(2)] for i in range(2)]


def centered(matrix):
    result=[];max_error=0
    for row in matrix:
        out=[]
        for z in row:
            parts=[];error=0
            for lo,hi in z:
                c=(lo+hi)//2;parts.append((c,c));error+=max(c-lo,hi-c)
            max_error=max(max_error,error);out.append(tuple(parts))
        result.append(out)
    # Spectral norm <= 2*largest absolute entry, bounded by real+imag errors.
    return result,2*max_error


def stable_product(unitary_interval,centre,error):
    pulse,local=centered(unitary_interval)
    product,roundoff=centered(matrix_product(pulse,centre))
    # Exact pulse is unitary. Norm(centre)<=1+error/Q.
    return product,error+iv.ceildiv(local*(Q+error),Q)+roundoff


def sincos(angle):
    centre=F(angle[0]+angle[1],2*Q);radius=F(angle[1]-angle[0],2*Q);halves=0
    while abs(centre)>4:centre/=2;halves+=1
    sine,cosine=iv.sincos(centre)
    for _ in range(halves):
        sine,cosine=mul(exact(2),mul(sine,cosine)),sub(square(cosine),square(sine))
    allowance=iv.ceildiv(radius.numerator*Q,radius.denominator)
    return (sine[0]-allowance,sine[1]+allowance),(cosine[0]-allowance,cosine[1]+allowance)


def frame(mass,momentum):
    p=exact(momentum);energy=sqrt(add(square(p),square(mass)))
    cosine=sqrt(divide(add(energy,p),mul(exact(2),energy)))
    sine=divide(mass,sqrt(mul(exact(2),mul(energy,add(energy,p)))))
    return energy,cosine,sine


def ramp_scale(s):return 1+3*(3*s*s-2*s*s*s)


def pair_evolution(mass_box,duration,steps):
    mass=(exact(F(mass_box[0]))[0],exact(F(mass_box[1]))[1])
    e0,c0,s0=frame(mass,F(1));ef,cf,sf=frame(mass,F(1,4))
    cosine=add(mul(cf,c0),mul(sf,s0));sine=sub(mul(sf,c0),mul(cf,s0))
    change=[[(cosine,ZERO),(sine,ZERO)],[(neg(sine),ZERO),(cosine,ZERO)]]
    evolved=[[(ONE,ZERO),(ZERO,ZERO)],[(ZERO,ZERO),(ONE,ZERO)]];operator_error=0
    if steps:
        dt=F(duration,steps)
        for step in range(steps):
            p=1/ramp_scale(F(2*step+1,2*steps));ep=exact(p)
            energy=sqrt(add(square(mass),square(ep)));sin,cos=sincos(mul(energy,exact(dt)))
            scale=divide(sin,energy)
            hz=divide(add(square(mass),ep),e0);hx=divide(mul(mass,sub(ONE,ep)),e0)
            iz=neg(mul(scale,hz));ix=neg(mul(scale,hx))
            pulse=[[(cos,iz),(ZERO,ix)],[(ZERO,ix),(cos,neg(iz))]]
            evolved,operator_error=stable_product(pulse,evolved,operator_error)
    centre,operator_error=stable_product(change,evolved,operator_error)
    result=[[[ (part[0]-operator_error,part[1]+operator_error) for part in z] for z in row] for row in centre]
    if mass==ZERO:
        result[0][1]=(ZERO,ZERO);result[1][0]=(ZERO,ZERO)
    probability=cnorm(result[0][1])
    probability=(max(0,probability[0]),min(Q,probability[1]))
    if probability[0]>probability[1]:raise ArithmeticError('invalid production probability enclosure')
    numeric=np.array([[cmid(x) for x in row] for row in result])
    error=float(np.linalg.norm(numeric.conj().T@numeric-np.eye(2)))
    if error>1e-11:raise ArithmeticError('pair evolution midpoint lost unitarity')
    return {'matrix':result,'numeric':numeric,'probability':probability,'initial_energy':e0,'final_energy':ef,
            'unitarity_midpoint_error':error,'interval_operator_error':str(F(operator_error,Q))}


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    for package in ('conditional_particle_spectrum','particle_joint_mass_identity','thermal_particle_records'):
        paths|={f'model/{package}/results.json',f'model/{package}/PROOF.md'}
    paths|={'model/native_effect_synthesis/PROOF.md','model/quantum_thermal_relaxation/PROOF.md'}
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


@lru_cache(maxsize=1)
def source():
    _,_,parent,groups=joint.inputs();f=quadratic.fermions()
    h,delta=quadratic.kinetic.momentum_pair_source((0,0,1),kernel=f['mass_basis'],speed=1.)
    symbol=np.block([[h,delta],[delta.conj().T,-h.T]]);kinetic=np.diag(np.diag(symbol));off=symbol-kinetic
    change=block_diag(np.kron(np.eye(4),quadratic.numeric(f['W'])),np.kron(np.eye(4),quadratic.numeric(f['W'])).conj())
    mass_frame=change.conj().T@f['named_positive_frame'];positive=np.where(np.diag(kinetic).real>0)[0]
    location={k:int(positive[np.argmax(abs(mass_frame[positive,k]))]) for k in range(180)}
    inverse={v:k for k,v in location.items()}
    if len(inverse)!=180:raise ArithmeticError('named frame does not cover every positive mass-pair coordinate')
    # A unit structural coefficient chooses only the zero-mass mode pairing;
    # the actual mass and its entire evolution remain exactly zero.
    template=f['mass_basis'].copy().as_mutable()
    for i in range(45):
        if all(template[i,j]==0 for j in range(45)):template[i,i]=1
    ht,dt=quadratic.kinetic.momentum_pair_source((0,0,1),kernel=template,speed=1.)
    support=np.block([[ht,dt],[dt.conj().T,-ht.T]])-kinetic
    eps=sp.Matrix([[0,1],[-1,0]]);zero90=sp.zeros(90);zero180=sp.zeros(180)
    def exact_off(kernel):
        small=sp.kronecker_product(eps,kernel.conjugate())
        delta=zero90.row_join(small).col_join((-small.T).row_join(zero90))
        return zero180.row_join(delta).col_join(delta.H.row_join(zero180))
    literal=exact_off(f['mass_basis']);structural=exact_off(template)
    targets_by_species={r['species']:sp.sympify(r['rest_coefficient_GeV']['exact']) for r in f['targets']}
    pairs=[];used=set();phase_errors=[];mode_to_group={i:g['joint_index'] for g in groups if g['source_sector']=='fermion' for i in g['mode_indices']}
    for k,i in location.items():
        if k in used:continue
        targets=np.where(abs(support[:,i])>0)[0]
        if len(targets)!=1:raise ArithmeticError('mass-kernel component is not a literal two-coordinate block')
        j=int(targets[0]);partner=inverse[(j+180)%360]
        if partner==k or partner in used:raise ArithmeticError('physical pair partition fails')
        first=f['modes'][k];second=f['modes'][partner]
        if first['species']!=second['species']:raise ArithmeticError('mass kernel mixes paired species')
        mass=targets_by_species[first['species']]
        entry=literal[j,i];sign=int(sp.sign(structural[j,i]))
        if sign not in (-1,1) or sp.simplify(entry-sign*mass)!=0:
            raise ArithmeticError('literal mass entry has an unaccounted phase')
        ib=(j+180)%360;jb=(i+180)%360
        if sp.simplify(literal[jb,ib]+sign*mass)!=0:
            raise ArithmeticError('canonical partner does not have opposite pairing sign')
        omega=math.sqrt(1+float(mass)**2);c=math.sqrt((omega+1)/(2*omega));sn=float(mass)/math.sqrt(2*omega*(omega+1))
        expected=np.zeros(360,complex);expected[i]=c;expected[j]=sign*sn
        phase_errors.append(float(np.linalg.norm(mass_frame[:,k]-expected)))
        expected=np.zeros(360,complex);expected[ib]=c;expected[jb]=-sign*sn
        phase_errors.append(float(np.linalg.norm(mass_frame[:,partner]-expected)))
        used|={k,partner};pairs.append({'mode_indices':[k,partner],'joint_groups':[mode_to_group[k],mode_to_group[partner]],
            'species':first['species'],'mass_symbol_coordinates':[i,j],
            'canonical_pairing_sign':sign,'literal_mass_entry_exact':str(entry),
            'partner_mass_symbol_coordinates':[ib,jb],
            'exact_partner_entry_is_negative':True,
            'structural_massless_pairing_only':bool(abs(off[j,i])==0)})
    if len(pairs)!=90 or len(used)!=180:raise ArithmeticError('full CAR source not partitioned into90 pairs')
    if max(phase_errors)>2e-12:raise ArithmeticError('named source mode phase differs from the literal analytic block')
    controls=[{'canonical_named_frame_phase_error':max(phase_errors),'all90_pair_signs_exact':True}]
    for scale in (F(1),F(2),F(4)):
        p=F(1)/scale;h,d=quadratic.kinetic.momentum_pair_source((0,0,float(p)),kernel=f['K'],speed=1.)
        actual=np.block([[h,d],[d.conj().T,-h.T]])
        target=float(p)*kinetic+off;error=float(np.linalg.norm(change.conj().T@actual@change-target))
        if error>3e-10:raise ArithmeticError('actual mass kernel and expansion blocks differ')
        controls.append({'scale_factor':str(scale),'full360_Nambu_intertwiner_error':error})
    masses={row['species']:row['rest_coefficient_GeV']['exact'] for row in f['targets']}
    boxes={g['candidate_species'][0]:g['rest_coefficient_GeV_interval'] for g in groups if g['source_sector']=='fermion'}
    return dict(parent=parent,fermions=f,groups=[g for g in groups if g['source_sector']=='fermion'],pairs=pairs,masses=masses,boxes=boxes,controls=controls)


def thermal_occupation(energy,temperature):
    if temperature is None:return ZERO
    lo,hi=(F(x,Q) for x in energy)
    a=negative_exp_endpoint(hi/F(temperature))[0];b=negative_exp_endpoint(lo/F(temperature))[1]
    return exact(a/(1+b))[0],exact(b/(1+a))[1]


def pair_state(evolution,temperature,sign=1):
    if sign not in (-1,1):raise ValueError('canonical sign must be plus or minus one')
    f=thermal_occupation(evolution['initial_energy'],temperature);c=evolution['probability'];one_minus=sub(ONE,f)
    odd=mul(f,one_minus);paired=add(square(f),mul(sub(ONE,mul(exact(2),f)),c))
    empty=sub(ONE,add(mul(exact(2),odd),paired));weights=[empty,odd,odd,paired]
    fp=mid(f);cp=float(abs(evolution['numeric'][0,1])**2)
    values=[(1-fp)**2-(1-2*fp)*cp,fp*(1-fp),fp*(1-fp),fp**2+(1-2*fp)*cp]
    if any(not float(F(w[0],Q))-2e-15<=v<=float(F(w[1],Q))+2e-15 for w,v in zip(weights,values)):
        raise ArithmeticError('consistent pair diagnostic lies outside its enclosure')
    if min(values)<-1e-13 or abs(sum(values)-1)>1e-12:raise ArithmeticError('complete pair occupation law fails')
    s=evolution['matrix'];coherence=cmul(conj(s[0][0]),conj(s[0][1]));factor=mul(exact(sign),sub(ONE,mul(exact(2),f)));coherence=(mul(coherence[0],factor),mul(coherence[1],factor))
    return {'occupation00_10_01_11_probabilities':values,
        'occupation_probability_intervals':[encode(x) for x in weights],
        'even_vacuum_pair_coherence':cencode(coherence),'initial_occupation_interval':encode(f),
        'canonical_pairing_sign':sign,
        'pair_unitary_even_block':'D_s X S X D_s with D_s=diag(1,s), X=swap, and the stored positive-mass S; odd states have phase one',
        'exact_zero_odd_branches_retained':temperature is None}


def convolution(a,b):
    out=[ZERO for _ in range(len(a)+len(b)-1)]
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]=add(out[i+j],mul(x,y))
    return out


def number_law(pair_rows):
    total=[ONE];point_total=np.array([1.]);point_species=defaultdict(lambda:np.array([1.]));point_marginals={g['joint_index']:np.array([1.]) for g in source()['groups']};species=defaultdict(lambda:[ONE]);marginals={g['joint_index']:[ONE] for g in source()['groups']}
    factors=[]
    for pair in pair_rows:
        w=[(exact(F(x[0]))[0],exact(F(x[1]))[1]) for x in pair['state']['occupation_probability_intervals']]
        polynomial=[w[0],add(w[1],w[2]),w[3]]
        total=convolution(total,polynomial);species[pair['species']]=convolution(species[pair['species']],polynomial)
        a,b=pair['joint_groups'];v=pair['state']['occupation00_10_01_11_probabilities'];vp=[v[0],v[1]+v[2],v[3]]
        point_total=np.convolve(point_total,vp);point_species[pair['species']]=np.convolve(point_species[pair['species']],vp)
        if a==b:point_marginals[a]=np.convolve(point_marginals[a],vp)
        else:
            point_marginals[a]=np.convolve(point_marginals[a],[v[0]+v[2],v[1]+v[3]])
            point_marginals[b]=np.convolve(point_marginals[b],[v[0]+v[1],v[2]+v[3]])
        factors.append({'joint_groups':[a,b],'coefficient_intervals00_10_01_11':[encode(x) for x in w]})
        if a==b:marginals[a]=convolution(marginals[a],polynomial)
        else:
            marginals[a]=convolution(marginals[a],[add(w[0],w[2]),add(w[1],w[3])])
            marginals[b]=convolution(marginals[b],[add(w[0],w[1]),add(w[2],w[3])])
    values=point_total.tolist();mean=sum(i*x for i,x in enumerate(values))
    if min(values)<-1e-14 or abs(sum(values)-1)>2e-11 or not 0<=mean<=180:
        raise ArithmeticError('point count law is not a complete normalized probability law')
    widths=[float(F(hi-lo,Q)) for lo,hi in total]
    if max(widths)>1e-14:raise ArithmeticError('count probability enclosure is too wide for the declared diagnostic')
    if any(not float(F(lo,Q))-2e-14<=v<=float(F(hi,Q))+2e-14 for (lo,hi),v in zip(total,values)):
        raise ArithmeticError('point count law disagrees with the certified interval law')
    normalization=(sum(x[0] for x in total),sum(x[1] for x in total))
    if len(total)!=181 or not normalization[0]<=Q<=normalization[1]:raise ArithmeticError('complete181-outcome count law differs')
    return {'joint_count_generating_polynomial':'product over90 complete factors (p00+p10*z_g+p01*z_h+p11*z_g*z_h); includes every cross-group correlation',
        'all90_joint_count_factors':factors,'all181_total_number_probabilities':values,
        'all181_total_number_probability_intervals':[encode(x) for x in total],
        'maximum_total_number_probability_interval_width':max(widths),
        'point_law_normalization_error':abs(sum(values)-1),
        'point_law_is_binary64_diagnostic':True,
        'point_vs_interval_diagnostic_tolerance':2e-14,
        'diagnostic_tolerance_is_rigorous_point_rounding_bound':False,
        'normalization_interval':encode(normalization),'expected_total_number':mean,
        'expected_total_number_interval':encode((sum(i*x[0] for i,x in enumerate(total)),sum(i*x[1] for i,x in enumerate(total)))),
        'complete84_group_count_marginals':[{'joint_index':g,'probabilities':point_marginals[g].tolist(),'probability_intervals':[encode(x) for x in p]} for g,p in sorted(marginals.items())],
        'complete12_species_number_laws':[{'species':k,'probabilities':point_species[k].tolist(),
            'probability_intervals':[encode(x) for x in p],
            'mean_number_interval':encode((sum(i*x[0] for i,x in enumerate(p)),sum(i*x[1] for i,x in enumerate(p)))),
            'mean_number':float(point_species[k]@np.arange(len(point_species[k])))} for k,p in sorted(species.items())]}


def smooth_discretization_bound(mass_upper,duration,steps):
    if not steps:return None
    t=F(duration);h=t/steps;m=F(mass_upper);first=F(9,2)/t;second=F(117,2)/(t*t)
    integral=first*h*h/4
    # For the admitted resolutions integral<1, exp(integral)<3.
    if integral>=1:raise ValueError('midpoint bound requires the declared refined regime')
    delta=steps*(h**3*(second/24+2*m*first/12)+F(3,2)*integral**2)
    return str(min(F(2),delta))


def run(duration,steps):
    results={};species=[]
    for name,box in source()['boxes'].items():
        evolution=pair_evolution(box,duration,steps);results[name]=evolution
        species.append({'species':name,'mass_GeV_interval':box,
            'vacuum_pair_production_probability_interval':encode(evolution['probability']),
            'vacuum_pair_production_probability_midpoint':mid(evolution['probability']),
            'full_Bogoliubov_pair_matrix_intervals':[[cencode(x) for x in row] for row in evolution['matrix']],
            'matrix_midpoint_unitarity_error':evolution['unitarity_midpoint_error'],
            'certified_Bogoliubov_matrix_operator_error':evolution['interval_operator_error'],
            'single_pair_continuous_cubic_ramp_operator_error_bound':smooth_discretization_bound(box[1],duration,steps)})
    preparations=[]
    for temperature in (None,F(1)):
        pairs=[dict(pair,state=pair_state(results[pair['species']],temperature,pair['canonical_pairing_sign'])) for pair in source()['pairs']]
        preparations.append({'initial_temperature_GeV':None if temperature is None else float(temperature),
            'initial_state':'instantaneous vacuum' if temperature is None else 'supplied finite fermion Gibbs Gaussian',
            'all90_full_pair_states':pairs,'complete_number_readout':number_law(pairs),
            'complete_Fock_dimension_exact':str(2**180),'state_is_full_fermionic_Gaussian':True})
    if steps:
        errors=defaultdict(F);multiplicity=defaultdict(int)
        for pair in source()['pairs']:
            name=pair['species'];multiplicity[name]+=1
            errors[name]+=min(F(2),2*F(smooth_discretization_bound(source()['boxes'][name][1],duration,steps)))
        total_error=sum(errors.values(),F())
        for prep in preparations:
            law=prep['complete_number_readout'];lo,hi=map(F,law['expected_total_number_interval'])
            prep['continuous_ramp_mean_number_bounds']={
                'local_observable_error_rule':'sum over contributing pairs min(2,2*delta_pair), since0<=N_pair<=2; no full-source dimension multiplier',
                'total_absolute_error_bound':str(total_error),
                'total_exact_mean_interval':[str(max(F(0),lo-total_error)),str(min(F(180),hi+total_error))],
                'species':[{'species':row['species'],'absolute_error_bound':str(errors[row['species']]),
                    'exact_mean_interval':[str(max(F(0),F(row['mean_number_interval'][0])-errors[row['species']])),
                        str(min(F(2*multiplicity[row['species']]),F(row['mean_number_interval'][1])+errors[row['species']]))]}
                    for row in law['complete12_species_number_laws']]}
    massless=[x for x in species if F(x['mass_GeV_interval'][1])==0]
    if len(massless)!=1 or massless[0]['vacuum_pair_production_probability_interval']!=['0','0']:
        raise ArithmeticError('zero-mass production control is not exactly zero')
    return {'law':'instantaneous1-to4 quench' if not steps else 'midpoint piecewise constant approximation to a=1+3(3s^2-2s^3)',
        'duration_in_declared_GeV_inverse_parameter':str(duration),'segments':steps,'scale_initial':'1','scale_final':'4',
        'species_production':species,'preparations':preparations,
        'full_Fock_continuous_ramp_operator_error_bound':None if not steps else str(min(F(2),sum((F(smooth_discretization_bound(source()['boxes'][p['species']][1],duration,steps)) for p in source()['pairs']),F()))),
        'full_Fock_reference_extended_output_trace_error_bound':None if not steps else str(min(F(2),2*sum((F(smooth_discretization_bound(source()['boxes'][p['species']][1],duration,steps)) for p in source()['pairs']),F())))}


def certificate():
    pins=bindings();d=source();protocols=[run(F(0),0)]+[run(t,n) for t in (F(1),F(8)) for n in (64,256,1024)]
    fast=protocols[3]['preparations'][0]['continuous_ramp_mean_number_bounds']['total_exact_mean_interval']
    slow=protocols[6]['preparations'][0]['continuous_ramp_mean_number_bounds']['total_exact_mean_interval']
    quench=protocols[0]['preparations'][0]['complete_number_readout']['expected_total_number_interval']
    if not (F(slow[1])<F(fast[0]) and F(fast[1])<F(quench[0])):
        raise ArithmeticError('continuous vacuum production ordering is not certified')
    result={'schema':'oph.expanding_quadratic_particles.v1','source_bindings':pins,'all_checks_pass':True,
        'actual360_Nambu_source_controls':d['controls'],'all180_modes':[x for x in d['fermions']['modes']],
        'initial_joint_groups':d['groups'],
        'final_instantaneous_joint_groups':[dict(g,initial_comoving_momentum_z_GeV=g['momentum_z_GeV'],
            momentum_z_GeV=str(F(g['momentum_z_GeV'])/4),
            final_energy_GeV_interval=encode(frame((exact(F(d['boxes'][g['candidate_species'][0]][0]))[0],exact(F(d['boxes'][g['candidate_species'][0]][1]))[1]),F(1,4))[0]),
            projector_convention='same named field lineage in the final instantaneous positive-energy frame') for g in d['groups']],
        'exterior_reorder':{'pair_order_mode_indices':[i for pair in d['pairs'] for i in pair['mode_indices']],
            'occupation_word_phase':'(-1)^I, I=number of occupied inversions of this pair-order permutation relative to original mode order'},
        'protocols':protocols,
        'certified_continuous_vacuum_mean_ordering':{'duration8_interval':slow,'duration1_interval':fast,'quench_interval':quench,
            'strict_slow_less_fast_less_quench':True,'refinement_used_segments':1024},
        'scope':{'same_native_calibrated_mass_kernel_and_CAR':True,'all_twelve_fermion_species_and_conjugates':True,
                 'all_group_count_correlations_and_zero_branches_retained':True,'arbitrary_reference_unitary_lift':True,
                 'boson_production_calculated':False,'scale_factor_and_clock_are_supplied':True,
                 'Einstein_backreaction_or_physical_cosmic_time_derived':False,
                 'calibrated_mass_values_are_new_predictions':False,
                 'finite_interval_pulse_results_identified_with_exact_continuous_ramp':False}}
    if bindings()!=pins:raise ArithmeticError('expansion source changed during execution')
    return result


def main():
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=parser.parse_args()
    result=certificate();path=HERE/'results.json';text=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if args.write:path.write_text(text)
    elif path.read_text()!=text:raise ArithmeticError('complete expansion receipt differs')
    print('PASS: all twelve same-source fermion production and complete count laws')
