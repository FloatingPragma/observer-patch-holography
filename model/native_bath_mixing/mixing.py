"""Rational entry enclosures, eigenbasis error bounds and complete cooling lag."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math

import numpy as np
from scipy.linalg import eigh

from model.native_thermal_completion import bath
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
PARENT='model/native_thermal_completion/results.json'


def authenticate():
    value=json.loads((ROOT/PARENT).read_text())
    for path,digest in value['source_bindings'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:
            raise ValueError('changed bath source: '+path)
    return value


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    paths|={PARENT,'model/native_thermal_completion/PROOF.md','model/native_thermal_admission/PROOF.md',
            'Lean/Dynamics/StationaryLyapunov.lean','model/requirements.txt'}
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


def sqrt_box(value,bits=180):
    value=F(value)
    if value<0:raise ValueError('nonnegative square root required')
    denominator=1<<bits
    lower=math.isqrt(value.numerator*(denominator**2)//value.denominator)
    a=F(lower,denominator)
    return a,a if a*a==value else F(lower+1,denominator)


def gamma(operations):
    return F(operations,2**53-operations)


def float_fraction(x):return F.from_float(float(x))


def maximum_fraction(array):
    return max((abs(float_fraction(x)) for x in np.asarray(array).flat),default=F(0))


def eigensystem_enclosure(matrix,entry_error):
    """Certify the eigenvalues of the exact symmetric matrix behind matrix.

    Every numeric matrix product is bounded by the standard binary64 dot-product
    gamma bound; exact Fraction calculations enclose residual/Gram/polar errors.
    The returned interval is not merely an unvalidated eigensolver residual.
    """
    matrix=np.asarray(matrix)
    if matrix.ndim!=2 or matrix.shape[0]!=matrix.shape[1] or not np.isfinite(matrix).all():
        raise ValueError('finite square symmetric matrix required')
    if not np.array_equal(matrix,matrix.T) or F(entry_error)<0:
        raise ValueError('symmetric matrix and nonnegative entry enclosure required')
    n=len(matrix)
    if n==1:raise ValueError('nonconstant component required')
    values,vectors=eigh(matrix)
    if not np.isfinite(values).all() or not np.isfinite(vectors).all():
        raise ArithmeticError('nonfinite candidate eigensystem')
    positive=abs(vectors[vectors!=0])
    if len(positive) and float(min(positive))**2*max(float(min(abs(values[values!=0]))),1e-300)<np.finfo(float).tiny:
        raise ArithmeticError('eigenbasis products leave the normal binary64 domain')
    gram=vectors.T@vectors
    absolute_gram=abs(vectors).T@abs(vectors)
    dot_error=gamma(2*n+2)
    positive_dot_upper=maximum_fraction(absolute_gram)/(1-dot_error)
    gram_residual=max(abs(float_fraction(gram[i,j])-int(i==j)) for i in range(n) for j in range(n))
    gram_norm_upper=n*(gram_residual+dot_error*positive_dot_upper)
    if gram_norm_upper>F(1,2):raise ArithmeticError('numerical basis cannot be certified invertible')
    reconstruction=(vectors*values)@vectors.T
    absolute_rows=abs(vectors)@abs(vectors).T
    diagonal_norm=maximum_fraction(values)
    product_error=gamma(2*n+4)*diagonal_norm*maximum_fraction(absolute_rows)/(1-dot_error)
    direct_residual=max(abs(float_fraction(matrix[i,j])-float_fraction(reconstruction[i,j]))
                        for i in range(n) for j in range(n))
    # Polar factor U: ||V-U|| <= epsilon, ||V||+1 <3 when epsilon<=1/2.
    spectral_error=n*(F(entry_error)+direct_residual+product_error)+3*diagonal_norm*gram_norm_upper
    eigenvalue=float_fraction(values[-2]);gap=(1-eigenvalue-spectral_error,1-eigenvalue+spectral_error)
    if gap[0]<=0:raise ArithmeticError('positive component gap not separated by certified arithmetic')
    return {'second_eigenvalue_interval':[str(eigenvalue-spectral_error),str(eigenvalue+spectral_error)],
        'gap_interval':[str(x) for x in gap],'gap_lower_numeric':float(gap[0]),'gap_upper_numeric':float(gap[1]),
        'matrix_entry_input_error_bound':str(entry_error),'eigenbasis_Gram_norm_error_bound':str(gram_norm_upper),
        'full_spectrum_matching_error_bound':str(spectral_error),
        'error_method':'exact entry intervals; gamma bounded Gram and reconstruction products; polar orthogonalization; Weyl',
        'floating_domain':'binary64 normal products and sums without overflow or underflow; all bound arithmetic rational'}


def component_matrix(q,component):
    source=bath.source();rows,_,_=bath.instrument(q);index={x:i for i,x in enumerate(component)}
    n=len(component);matrix=np.zeros((n,n));error=F(0)
    for x in component:
        for y,weight in rows[x].items():
            if y not in index:raise ArithmeticError('bath crosses a retained component')
            if x==y:lo=hi=weight
            else:
                exponent=abs(source['energy'][y]-source['energy'][x])//2
                a,b=sqrt_box(F(q)**exponent)
                factor=source['proposal'][(x,y)]/source['uniformization'];lo,hi=factor*a,factor*b
            midpoint=float((lo+hi)/2);matrix[index[x],index[y]]=midpoint
            represented=float_fraction(midpoint);error=max(error,abs(represented-lo),abs(represented-hi))
    if not np.array_equal(matrix,matrix.T):raise ArithmeticError('symmetrized interval midpoint is not symmetric')
    return matrix,error


def integer_log2_upper(value):
    value=F(value)
    if value<=1:return 0
    b=max(0,value.numerator.bit_length()-value.denominator.bit_length())
    while F(1<<b)<value:b+=1
    return b


def mixing_steps(gap,pi_min,epsilon=F(1,1000)):
    # exp(-2 gamma n)<=2^-b since log(2)<1; avoids unproved numeric logarithms.
    b=integer_log2_upper(1/(4*epsilon**2*pi_min))
    x=F(b)/(2*gap)
    return (x.numerator+x.denominator-1)//x.denominator


def lag_interval(probabilities,target,l1_error):
    central=sum((abs(float_fraction(x)-y) for x,y in zip(probabilities,target,strict=True)),F(0))/2
    radius=F(l1_error)/2
    return [str(max(F(0),central-radius)),str(min(F(1),central+radius))]


def certificate():
    parent=authenticate();pins=bindings();source=bath.source();temperatures=[]
    for q in bath.QUANTA:
        _,records,_=bath.instrument(q);components=[]
        for number,component in enumerate(source['components']):
            energies=[source['energy'][i] for i in component];floor=min(energies)
            raw=[q**((e-floor)//2) for e in energies];partition=sum(raw,F(0));minimum=min(raw)/partition
            row={'component':number,'state_indices':component,'states':len(component),
                'preserved_component_weight_exact':str(bath.component_weights()[number]),'minimum_conditional_Gibbs_weight':str(minimum)}
            if len(component)==1:
                row.update({'nonconstant_spectrum_empty':True,'mixing_steps_for_TV_1e_minus3':0})
            else:
                matrix,error=component_matrix(q,component);bound=eigensystem_enclosure(matrix,error)
                row.update(bound);row['mixing_steps_for_TV_1e_minus3']=mixing_steps(F(bound['gap_interval'][0]),minimum)
            components.append(row)
        encoded=[[x,proposal,target,label,str(p)] for x,proposal,target,label,p in records]
        fingerprint=hashlib.sha256(json.dumps(encoded,separators=(',',':')).encode()).hexdigest()
        original=next(r for r in parent['temperature_rows'] if r['q']==str(q))
        if fingerprint!=original['record_sha256']:raise ArithmeticError('original complete bath instrument changed')
        temperatures.append({'q_exp_minus_two_beta':str(q),'components':components,
            'all_original_proposal_rejection_idle_records':len(records),'original_record_sha256':fingerprint,
            'zero_records_retained':sum(p==0 for *_,p in records)})
    schedules=[]
    for experiment in parent['cooling_experiments']:
        stages=[]
        for stage in experiment['stages']:
            q=F(stage['q_exp_minus_two_beta']);target,_=bath.equilibrium(q)
            bounds=lag_interval(stage['complete_state_probabilities'],target,stage['state_probability_L1_roundoff_bound'])
            stages.append({'q':str(q),'attempts':stage['steps'],'TV_to_component_Gibbs_interval':bounds,
                'TV_lower_numeric':float(F(bounds[0])),'TV_upper_numeric':float(F(bounds[1])),
                'original_state_probability_L1_roundoff_bound':stage['state_probability_L1_roundoff_bound'],
                'complete_original_state_probabilities':stage['complete_state_probabilities'],
                'exact_component_Gibbs_probabilities':[str(x) for x in target]})
        schedules.append({'steps_per_setting':experiment['steps_per_temperature'],'stages':stages})
    result={'schema':'oph.native_bath_mixing.v1','source_bindings':pins,'all_checks_pass':True,
        'complete_state_count':1296,'retained_components':26,'temperature_rows':temperatures,
        'certified_cooling_lag':schedules,
        'conditional_TV_bound':'TV(p_C P_C^n, pi_C)<=min(1,sqrt(1/pi_C,min-1)*(1-gap_C)^n/2)',
        'global_bound':'sum_C w_C times the conditional bound, only for the same preserved component weights',
        'scope':{'unchanged_bath_instrument':True,'native_descending_repair_changed':False,
                 'global_unique_equilibrium_claimed':False,'physical_temperature_or_cosmic_time_derived':False,
                 'complete_bath_history_mixes_to_input_independence_claimed':False}}
    if bindings()!=pins:raise ArithmeticError('mixing sources changed during execution')
    return result


def equal(a,b):
    if isinstance(b,dict):return isinstance(a,dict) and set(a)==set(b) and all(equal(a[k],v) for k,v in b.items())
    if isinstance(b,list):return isinstance(a,list) and len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b,strict=True))
    return type(a) is type(b) and a==b


def main():
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=parser.parse_args()
    result=certificate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif not equal(json.loads(path.read_text()),result):raise ArithmeticError('certified mixing receipt differs')
    print('PASS: all component gaps and complete cooling-lag enclosures')
