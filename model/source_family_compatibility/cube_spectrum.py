"""Complete finite-clock spectral records on the centered full-cube ground.

Walsh labels specify prepared coefficients on the open cube. They are not
eigenvalues of a spatial translation generator. Every probability retains
its original addition/removal weight.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import factorial, isqrt

import numpy as np
from scipy import sparse

from . import cube_response as reference
from model.centered_source_ground import particles

BITS = 46
ROOT_BITS = 42
DEGREE = 32
PHASE_BINS = 256
INTMAX = int(np.iinfo(np.int64).max)
SOURCE_UPPER = F(45, 8)
GROUND_ERROR = F(1, 1920)
NUMERICAL_ERROR = F(14, 10**9)


def _norm2(x):
    return sum(int(y)**2 for y in x)


def _dot(x, y):
    return sum(int(a)*int(b) for a,b in zip(x,y,strict=True))


def _maximum(x):
    return max(abs(int(x.min())),abs(int(x.max())))


def _rounded(x, denominator):
    if type(denominator)is not int or denominator<=0:
        raise ValueError('positive integer rounding denominator required')
    if _maximum(x)+denominator//2>INTMAX:
        raise OverflowError('dyadic rounding addition exceeds int64')
    return (x+denominator//2)//denominator


def _action(h, x, row_bound):
    if _maximum(x)*row_bound>INTMAX:
        raise OverflowError('exact clock action exceeds int64')
    return h@x


def phase_step(h, real, imag, row_bound):
    """Round the degree-32 source Taylor polynomial over exact integers."""
    tr,ti=real.copy(),imag.copy();out_r,out_i=tr.copy(),ti.copy()
    for k in range(1,DEGREE+1):
        nr=_rounded(_action(h,ti,row_bound),32*k)
        ni=-_rounded(_action(h,tr,row_bound),32*k)
        if max(_maximum(out_r)+_maximum(nr),_maximum(out_i)+_maximum(ni))>INTMAX:
            raise OverflowError('Taylor sum exceeds int64')
        out_r+=nr;out_i+=ni;tr,ti=nr,ni
    return out_r,out_i


def step_error(dimension):
    """Source norm <=45/8 and duration1/2; every vector has norm <=4."""
    if type(dimension)is not int or dimension<1:
        raise ValueError('positive source dimension required')
    root=isqrt(2*dimension)
    if root*root<2*dimension:root+=1
    rounding=F(root,2**(BITS+1))
    # e^3<21; total rounded-term error <=degree*e^3*rounding.
    return 4*F(3**(DEGREE+1),factorial(DEGREE+1))+DEGREE*21*rounding


def pi_interval():
    """Machin's formula and alternating rational arctangent remainders."""
    def atan_inverse(n):
        total=sum((F((-1)**j,(2*j+1)*n**(2*j+1))for j in range(40)),F())
        return total,total+F(1,81*n**81)
    al,au=atan_inverse(5);bl,bu=atan_inverse(239)
    lo,hi=16*al-4*bu,16*au-4*bl
    scale=2**80
    lower=F((lo*scale).numerator//(lo*scale).denominator,scale)
    upper=F(-((-hi*scale).numerator//(-hi*scale).denominator),scale)
    if not F(3)<lower<upper<F(22,7):raise ArithmeticError('Machin interval failed')
    return lower,upper


@lru_cache(maxsize=4)
def roots(count=PHASE_BINS):
    """Exact dyadic Fourier coefficients, each complex error <2^-41."""
    if type(count)is not int or count<2 or count&(count-1):
        raise ValueError('power-of-two phase register required')
    pl,pu=pi_interval();pm=(pl+pu)/2;scale=2**ROOT_BITS
    rows=[]
    for j in range(count):
        signed=j if j<=count//2 else j-count
        angle=2*pm*F(signed,count)
        cosine=sum(((-1)**k*angle**(2*k)/factorial(2*k)for k in range(31)),F())
        sine=sum(((-1)**k*angle**(2*k+1)/factorial(2*k+1)for k in range(31)),F())
        # Lagrange remainders and the Lipschitz error in the Machin midpoint.
        error=max(abs(angle)**62/F(factorial(62)),abs(angle)**63/F(factorial(63)))+(pu-pl)/2
        integers=tuple((value*scale+F(1,2)).numerator//(value*scale+F(1,2)).denominator
                       for value in(cosine,sine))
        if any(abs(F(integer,scale)-value)+error>F(1,scale)
               for integer,value in zip(integers,(cosine,sine),strict=True)):
            raise ArithmeticError('Fourier coefficient enclosure failed')
        rows.append(integers)
    return tuple(rows)


def source_profiles(data):
    """Actual N9, Sz=1/2 action and all eight normalized real wavepackets."""
    points=tuple(reference.global_ground.cube.GEOMETRY.sites)
    if not F(3,4)<=F(data['norm_squared'],2**46)<=2:
        raise ValueError('the certified 23-bit ground normalization is required')
    labels=tuple(product((0,1),repeat=3))
    masks=np.asarray([m for m in range(1<<16)
                      if m.bit_count()==9 and (m&0x5555).bit_count()==5],dtype=np.int64)
    if len(masks)!=3920:raise ArithmeticError('complete charged source block changed')
    row_of={int(m):i for i,m in enumerate(masks)}
    h=data['centered_numerator'][masks][:,masks].tocsr()
    if h.dtype!=np.dtype(np.int64):
        raise ValueError('exact int64 source clock required')
    if data['centered_numerator'][:,masks].nnz!=h.nnz:
        raise ArithmeticError('charged block is not reducing for the actual source')
    h=h+66*sparse.eye(len(masks),dtype=np.int64,format='csr')
    if (h-h.T).nnz:raise ArithmeticError('source clock is not symmetric')
    row_bound=max(map(int,np.asarray(abs(h).sum(axis=1)).ravel()))
    if row_bound>90:raise ArithmeticError('source clock exceeds 45/8 bound')
    entries=[([],[],[])for _ in labels]
    for column,mask in enumerate(data['balanced_indices']):
        mask=int(mask)
        for x,point in enumerate(points):
            mode=2*x
            if not mask&(1<<mode):
                sign=(-1)**(mask&((1<<mode)-1)).bit_count()
                row=row_of[mask|(1<<mode)]
                for col,k in enumerate(labels):
                    entries[col][0].append(row);entries[col][1].append(column)
                    entries[col][2].append(sign*(-1)**sum(a*b for a,b in zip(k,point,strict=True)))
    creators=tuple(sparse.csr_matrix((np.asarray(values,dtype=np.int64),(rows,cols)),
                       shape=(len(masks),4900))for rows,cols,values in entries)
    vectors=np.column_stack([a@data['trial']for a in creators])
    if any(F(_norm2(vectors[:,j]),2**46)>9 for j in range(8)):
        raise ArithmeticError('unnormalized computation vectors exceed norm three')
    return {'clock':h,'row_bound':row_bound,'vectors':vectors,'labels':labels,
            'addresses':masks,'points':points,'unnormalized_creators':creators}


def symmetry_certificate(data):
    """Bind the charge-conjugate table to the literal full-source involution."""
    complement,phase=particles.particle_hole_data()
    h=data['centered_numerator']
    transformed=h[complement][:,complement].multiply(phase[:,None]).multiply(phase[None,:]).tocsr()
    if (transformed-h).nnz:
        raise ArithmeticError('the literal charge conjugation changes the source clock')
    points=tuple(reference.global_ground.cube.GEOMETRY.sites)
    for mode in range(16):
        lower=(1<<mode)-1
        for mask in range(1<<16):
            if mask&(1<<mode):continue
            comp=int(complement[mask]);expected=(-1)**(sum(points[mode//2])+(mask&lower).bit_count())
            actual=int(phase[mask])*int(phase[comp^(1<<mode)])*(-1)**(comp&lower).bit_count()
            if actual!=expected:raise ArithmeticError('charge-conjugate creator sign changed')
    return {'complete_source_clock_invariance':True,'exact_creator_maps':16,
            'ground_unique_spin_singlet':True,'full_number_and_spin_source_proof_consumed':True,
            'physical_particle_conjugation_identified':False}


def moment_intervals(data,source):
    """Original moments of K+33/8, before any charge conditioning."""
    vectors=source['vectors'];h=source['clock'];bound=source['row_bound']
    first=_action(h,vectors,bound);second=_action(h,first,bound)
    ground_error=data['half_trace_distance_upper'];den=8*data['norm_squared']
    rows=[]
    for j in range(8):
        v,u,w=vectors[:,j],first[:,j],second[:,j]
        numerators=(_norm2(v),_dot(v,u),_norm2(u),_dot(u,w),_norm2(w))
        moments=tuple(F(value,den*16**power)for power,value in enumerate(numerators))
        errors=tuple(2*ground_error*SOURCE_UPPER**power for power in range(5))
        rows.append(tuple((max(F(0),m-e),m+e)for m,e in zip(moments,errors,strict=True)))
    return tuple(rows)


def excitation_moments(data,source):
    """Ground-energy cancellation uses the exact source commutator."""
    one_indices=1<<np.arange(16)
    one=data['centered_numerator'][one_indices][:,one_indices].toarray()+8*np.eye(16,dtype=np.int64)
    if not np.array_equal(one@one,3*np.eye(16,dtype=np.int64)):
        raise ArithmeticError('actual cube one-body kinetic norm changed')
    # ||h f||=sqrt(3)/16 <7/64. The dressed b family has exact CAR.
    commutator_norm=F(39,64)
    indices=data['balanced_indices']
    neutral=data['centered_numerator'][indices][:,indices]+66*sparse.eye(4900,dtype=np.int64,format='csr')
    row_bound=max(map(int,np.asarray(abs(neutral).sum(axis=1)).ravel()))
    hp=_action(neutral,data['trial'],row_bound)
    transported=np.column_stack([a@hp for a in source['unnormalized_creators']])
    acted=_action(source['clock'],source['vectors'],source['row_bound'])
    if _maximum(transported)+_maximum(acted)>INTMAX:
        raise OverflowError('exact source commutator subtraction exceeds int64')
    commutator=acted-transported
    denominator=8*data['norm_squared'];delta=data['half_trace_distance_upper']
    result=[]
    for j in range(8):
        v,w=source['vectors'][:,j],commutator[:,j]
        values=(F(_norm2(v),denominator),F(_dot(v,w),16*denominator),F(_norm2(w),256*denominator))
        errors=(delta,2*delta*commutator_norm,delta*commutator_norm**2)
        result.append(tuple((max(F(0),x-e),x+e)for x,e in zip(values,errors,strict=True)))
    return tuple(result)


def phase_records(data,source,count=PHASE_BINS):
    """Every clock record, certified without identifying a Ritz value as a pole."""
    if type(count)is not int or count>PHASE_BINS:
        raise ValueError('at most256 source clock addresses admitted')
    fourier=roots(count);v=source['vectors'];h=source['clock'];bound=source['row_bound']
    real=v*(1<<(BITS-23));imag=np.zeros_like(real)
    correlations=[]
    for step in range(count):
        correlations.append(tuple((_dot(v[:,j],real[:,j]),_dot(v[:,j],imag[:,j]))for j in range(8)))
        if step+1<count:
            real,imag=phase_step(h,real,imag,bound)
            if any(_norm2(real[:,j])+_norm2(imag[:,j])>16*2**(2*BITS)for j in range(8)):
                raise ArithmeticError('approximate clock vector exceeds its certified norm bound')
    base=8*data['norm_squared']*2**(BITS-23)
    denominator=count*count*base*2**ROOT_BITS
    values=[]
    for col in range(8):
        row=[]
        for record in range(count):
            numerator=count*correlations[0][col][0]*2**ROOT_BITS
            for lag in range(1,count):
                ar,ai=correlations[lag][col];cr,ci=fourier[(lag*record)%count]
                numerator+=2*(count-lag)*(ar*cr-ai*ci)
            row.append(F(numerator,denominator))
        values.append(tuple(row))
    e=step_error(v.shape[0])
    trial_scale=F(data['norm_squared'],2**46)
    if not F(3,4)<=trial_scale<=2:
        raise ValueError('the certified clock correlation normalization is required')
    numerical=e*F(count*count-1,8*count)/trial_scale+F(1,2**(ROOT_BITS-2))
    if numerical>=F(1,10**6):
        raise ArithmeticError('clock precision fails its admitted probability tolerance')
    return {'values':tuple(values),'per_record_numerical_error':numerical,
            'ground_effect_error':data['half_trace_distance_upper'],
            'one_step_vector_error':e,'count':count}


def outward(value,grid=10**6):
    x=F(value)*grid
    return (F(x.numerator//x.denominator,grid),F(-((-x.numerator)//x.denominator),grid))


def interval(lo,hi,grid=10**6):
    if lo>hi:raise ValueError('ordered interval required')
    return (outward(lo,grid)[0],outward(hi,grid)[1])


def window_tail(count,distance):
    """Wrapped source-energy error tail; duration1/2 and pi<22/7."""
    d=F(distance)
    if type(count)is not int or count<2 or count&(count-1) or d<=0:
        raise ValueError('power-of-two clock and positive energy tolerance required')
    q=F(7*count,88)*d
    return min(F(1),F(1,2)/(q-1)) if q>1 else F(1)


def summarize(data,source,readout,moments,local_moments):
    n=readout['count'];numerr=readout['per_record_numerical_error'];gerr=readout['ground_effect_error']
    if n!=PHASE_BINS:raise ValueError('the256-record source partition is required')
    pl,pu=pi_interval();elo,ehi=data['ground_energy_interval'];offset=(elo+F(33,8),ehi+F(33,8))
    groups=((0,3),(3,6),(6,9),(9,12),(12,15),(15,18),(18,24),(24,32),(32,64),(64,128),(128,n))
    rows=[]
    for index,label in enumerate(source['labels']):
        values=readout['values'][index];m=moments[index]
        probabilities=[]
        for left,right in groups:
            if right<=left:continue
            value=sum(values[left:right],F());error=(right-left)*numerr+gerr
            probabilities.append({'phase_records':[left,right],
                'original_probability':interval(max(F(0),value-error),min(F(1),value+error))})
        local=local_moments[index];weight=local[0]
        if weight[0]<=0 or local[1][0]<0:
            raise ArithmeticError('positive original weight and energy moment required')
        mean=interval(local[1][0]/weight[1],local[1][1]/weight[0])
        variance=interval(max(F(0),local[2][0]/weight[1]-(local[1][1]/weight[0])**2),
                          local[2][1]/weight[0]-(local[1][0]/weight[1])**2)
        windows=[]
        for left,right,distance in((6,15,F(1,5)),(3,21,F(1,4)),(0,24,F(1,4))):
            if right>=n:continue
            value=sum(values[left:right+1],F());error=(right-left+1)*numerr+gerr
            tail=window_tail(n,distance)
            lower=max(F(0),value-error-tail*weight[1])
            window=(max(F(1,8),F(4*left,n)*pl-distance-offset[1]),
                    F(4*right,n)*pu+distance-offset[0])
            windows.append({'excitation_energy_window':interval(*window),
                            'original_weight_lower':outward(lower)[0],
                            'phase_tail_upper':tail,'phase_records_inclusive':[left,right]})
        rows.append({'profile_pi_multiple':label,'original_addition_weight':interval(*weight),
                     'conditional_excitation_mean':mean,'conditional_energy_variance':variance,
                     'original_excitation_moment_intervals':tuple(interval(*x)for x in local),
                     'original_shifted_clock_moments_through_four':tuple(interval(*x)for x in m),
                     'all_phase_record_probability_intervals':tuple(
                         interval(max(F(0),v-numerr-gerr),min(F(1),v+numerr+gerr))for v in values),
                     'complete_phase_record_groups':probabilities,'true_spectral_windows':windows})
    return {'cube_sites':8,'c':F(1,8),'U':F(1),'source_full_Fock_dimension':65536,
            'computed_charge_block_dimension':3920,'phase_bins':n,'clock_step':F(1,2),
            'clock_reference':'K+33/8; excitation energy subtracts E0+33/8',
            'ground_energy_interval':interval(elo,ehi,10**8),'ground_half_trace_error':gerr,
            'phase_record_numerical_error_upper':outward(numerr,10**12)[1],
            'profiles':rows,'removal_profile_correspondence':'k maps to k xor (1,1,1)',
            'spin_components_have_equal_ground_response':True,
            'complete_addition_removal_instrument':True,'phase_labels_are_exact_momentum':False,
            'phase_record_is_sharp_spectral_projector':False,'native_large_clock_tape_executed':False,
            'coherent_preparation_and_clock_service_supplied':True,
            'physical_service_failures_retained_with_separate_error_budget':True,
            'physical_mass_assigned':None,'continuum_spectral_atom_inferred':False}


def _admit_profile(actual):
    """Check every fixed reported envelope against the fresh source calculation."""
    proposed=_ADMITTED_PROFILES[actual['profile_pi_multiple'][2]]
    def admit(value,expected,key=''):
        if key=='profile_pi_multiple':return tuple(value)
        if isinstance(value,(list,tuple))and len(value)==2 and all(type(x)is F for x in value):
            low,high=map(F,expected)
            if not low<=value[0]<=value[1]<=high:
                raise ArithmeticError('source interval misses its fixed output envelope: '+key)
            return low,high
        if type(value)is F:
            target=F(expected)
            if key=='original_weight_lower':
                if not 0<=target<=value:raise ArithmeticError('source window misses fixed original weight')
            elif value!=target:raise ArithmeticError('exact source allowance changed')
            return target
        if isinstance(value,dict):
            if list(value)!=list(expected):raise ArithmeticError('source record structure changed')
            return {k:admit(v,expected[k],k)for k,v in value.items()}
        if isinstance(value,(list,tuple)):
            if len(value)!=len(expected):raise ArithmeticError('complete source record count changed')
            values=[admit(v,e,key)for v,e in zip(value,expected,strict=True)]
            return tuple(values)if isinstance(value,tuple)else values
        if value!=expected:raise ArithmeticError('source record label changed')
        return value
    return admit(actual,proposed)


@lru_cache(maxsize=1)
def _certificate():
    data=reference.certified_ground()
    if data['half_trace_distance_upper']>=GROUND_ERROR:
        raise ArithmeticError('reference misses the fixed spectral readout accuracy')
    bounds=(reference.gap.GROUND_ENCLOSURE_LOWER,reference.gap.GROUND_ENCLOSURE_UPPER)
    if not bounds[0]<=data['ground_energy_interval'][0]<=data['ground_energy_interval'][1]<=bounds[1]:
        raise ArithmeticError('source ground interval exceeds the inherited fixed enclosure')
    data['half_trace_distance_upper']=F(GROUND_ERROR.numerator,GROUND_ERROR.denominator)
    data['ground_energy_interval']=bounds
    source=source_profiles(data)
    symmetry=symmetry_certificate(data)
    moments=moment_intervals(data,source);local_moments=excitation_moments(data,source)
    readout=phase_records(data,source)
    if readout['per_record_numerical_error']>=NUMERICAL_ERROR:
        raise ArithmeticError('clock arithmetic misses the fixed output enclosure')
    readout['per_record_numerical_error']=F(NUMERICAL_ERROR.numerator,NUMERICAL_ERROR.denominator)
    result=summarize(data,source,readout,moments,local_moments)
    result['profiles']=[_admit_profile(row)for row in result['profiles']]
    result['charge_and_spin_attachment']=symmetry
    return result


def certificate():
    def detach(x):
        if type(x)is F:return F(x.numerator,x.denominator)
        if isinstance(x,dict):return {k:detach(v)for k,v in x.items()}
        if isinstance(x,tuple):return tuple(detach(v)for v in x)
        if isinstance(x,list):return [detach(v)for v in x]
        return x
    return detach(_certificate())


# Rational output proposals; every entry is independently admitted against
# fresh exact source intervals. They are not trusted spectral input.
_ADMITTED_PROFILES = ({'profile_pi_multiple': (0, 0, 0),
  'original_addition_weight': ('58337/100000', '14611/25000'),
  'conditional_excitation_mean': ('5997/12500', '24141/50000'),
  'conditional_energy_variance': ('493/50000', '139/10000'),
  'original_excitation_moment_intervals': (('58337/100000', '14611/25000'),
                                           ('14019/50000', '3521/12500'),
                                           ('14199/100000', '89/625')),
  'original_shifted_clock_moments_through_four': (('14571/25000', '1828/3125'),
                                                  ('2757/10000', '5749/20000'),
                                                  ('439/4000', '1757/10000'),
                                                  ('0', '26123/100000'),
                                                  ('0', '27123/25000')),
  'all_phase_record_probability_intervals': (('0', '53/50000'),
                                             ('3/20000', '61/50000'),
                                             ('39/100000', '73/50000'),
                                             ('1/1250', '187/100000'),
                                             ('159/100000', '133/50000'),
                                             ('173/50000', '453/100000'),
                                             ('129/12500', '1139/100000'),
                                             ('12633/100000', '637/5000'),
                                             ('3113/50000', '3167/50000'),
                                             ('4083/100000', '419/10000'),
                                             ('1577/20000', '999/12500'),
                                             ('9939/100000', '5023/50000'),
                                             ('3201/100000', '3309/100000'),
                                             ('443/6250', '1439/20000'),
                                             ('1503/50000', '3113/100000'),
                                             ('283/50000', '673/100000'),
                                             ('59/50000', '9/4000'),
                                             ('3/6250', '31/20000'),
                                             ('17/100000', '31/25000'),
                                             ('0', '53/50000'),
                                             ('0', '19/20000'),
                                             ('0', '87/100000'),
                                             ('0', '81/100000'),
                                             ('0', '77/100000'),
                                             ('0', '37/50000'),
                                             ('0', '71/100000'),
                                             ('0', '69/100000'),
                                             ('0', '17/25000'),
                                             ('0', '67/100000'),
                                             ('0', '67/100000'),
                                             ('0', '2/3125'),
                                             ('0', '63/100000'),
                                             ('0', '31/50000'),
                                             ('0', '31/50000'),
                                             ('0', '61/100000'),
                                             ('0', '3/5000'),
                                             ('0', '3/5000'),
                                             ('0', '59/100000'),
                                             ('0', '59/100000'),
                                             ('0', '29/50000'),
                                             ('0', '29/50000'),
                                             ('0', '29/50000'),
                                             ('0', '29/50000'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '29/50000'),
                                             ('0', '29/50000'),
                                             ('0', '29/50000'),
                                             ('0', '29/50000'),
                                             ('0', '59/100000'),
                                             ('0', '59/100000'),
                                             ('0', '3/5000'),
                                             ('0', '3/5000'),
                                             ('0', '61/100000'),
                                             ('0', '61/100000'),
                                             ('0', '31/50000'),
                                             ('0', '63/100000'),
                                             ('0', '2/3125'),
                                             ('0', '13/20000'),
                                             ('0', '33/50000'),
                                             ('0', '17/25000'),
                                             ('0', '69/100000'),
                                             ('0', '9/12500'),
                                             ('0', '37/50000'),
                                             ('0', '39/50000'),
                                             ('0', '41/50000'),
                                             ('0', '11/12500'),
                                             ('0', '3/3125')),
  'complete_phase_record_groups': ({'phase_records': (0, 3), 'original_probability': ('1/625', '67/25000')},
                                   {'phase_records': (3, 6), 'original_probability': ('173/25000', '1/125')},
                                   {'phase_records': (6, 9),
                                    'original_probability': ('19999/100000', '10053/50000')},
                                   {'phase_records': (9, 12),
                                    'original_probability': ('11007/50000', '22121/100000')},
                                   {'phase_records': (12, 15),
                                    'original_probability': ('13403/100000', '1351/10000')},
                                   {'phase_records': (15, 18),
                                    'original_probability': ('839/100000', '473/50000')},
                                   {'phase_records': (18, 24),
                                    'original_probability': ('49/25000', '303/100000')},
                                   {'phase_records': (24, 32),
                                    'original_probability': ('61/100000', '21/12500')},
                                   {'phase_records': (32, 64),
                                    'original_probability': ('31/50000', '169/100000')},
                                   {'phase_records': (64, 128),
                                    'original_probability': ('3/100000', '11/10000')},
                                   {'phase_records': (128, 256),
                                    'original_probability': ('93/25000', '479/100000')}),
  'true_spectral_windows': ({'excitation_energy_window': ('1/8', '93539/100000'),
                             'original_weight_lower': '46633/100000',
                             'phase_tail_upper': '55/338',
                             'phase_records_inclusive': (6, 15)},
                            {'excitation_energy_window': ('1/8', '127991/100000'),
                             'original_weight_lower': '6277/12500',
                             'phase_tail_upper': '11/90',
                             'phase_records_inclusive': (3, 21)},
                            {'excitation_energy_window': ('1/8', '142717/100000'),
                             'original_weight_lower': '25251/50000',
                             'phase_tail_upper': '11/90',
                             'phase_records_inclusive': (0, 24)})},
 {'profile_pi_multiple': (0, 0, 1),
  'original_addition_weight': ('10389/25000', '41663/100000'),
  'conditional_excitation_mean': ('26149/50000', '26369/50000'),
  'conditional_energy_variance': ('83/10000', '91/6250'),
  'original_excitation_moment_intervals': (('10389/25000', '41663/100000'),
                                           ('5447/25000', '10959/50000'),
                                           ('2983/25000', '11973/100000')),
  'original_shifted_clock_moments_through_four': (('1297/3125', '10429/25000'),
                                                  ('2663/12500', '22479/100000'),
                                                  ('1087/12500', '15291/100000'),
                                                  ('0', '317/1250'),
                                                  ('0', '13537/12500')),
  'all_phase_record_probability_intervals': (('0', '79/100000'),
                                             ('0', '17/20000'),
                                             ('0', '19/20000'),
                                             ('1/25000', '111/100000'),
                                             ('33/100000', '141/100000'),
                                             ('97/100000', '41/20000'),
                                             ('31/10000', '417/100000'),
                                             ('3501/100000', '3609/100000'),
                                             ('877/25000', '113/3125'),
                                             ('1557/25000', '1267/20000'),
                                             ('3741/100000', '481/12500'),
                                             ('883/10000', '8937/100000'),
                                             ('489/12500', '4019/100000'),
                                             ('727/20000', '1871/50000'),
                                             ('2577/50000', '5261/100000'),
                                             ('1003/100000', '1111/100000'),
                                             ('3/2500', '57/25000'),
                                             ('17/50000', '71/50000'),
                                             ('3/100000', '11/10000'),
                                             ('0', '47/50000'),
                                             ('0', '21/25000'),
                                             ('0', '39/50000'),
                                             ('0', '73/100000'),
                                             ('0', '7/10000'),
                                             ('0', '67/100000'),
                                             ('0', '33/50000'),
                                             ('0', '2/3125'),
                                             ('0', '63/100000'),
                                             ('0', '63/100000'),
                                             ('0', '61/100000'),
                                             ('0', '61/100000'),
                                             ('0', '3/5000'),
                                             ('0', '3/5000'),
                                             ('0', '59/100000'),
                                             ('0', '29/50000'),
                                             ('0', '29/50000'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '27/50000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '11/20000'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '7/12500'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '57/100000'),
                                             ('0', '29/50000'),
                                             ('0', '29/50000'),
                                             ('0', '29/50000'),
                                             ('0', '59/100000'),
                                             ('0', '3/5000'),
                                             ('0', '3/5000'),
                                             ('0', '61/100000'),
                                             ('0', '31/50000'),
                                             ('0', '63/100000'),
                                             ('0', '2/3125'),
                                             ('0', '33/50000'),
                                             ('0', '17/25000'),
                                             ('0', '7/10000'),
                                             ('0', '37/50000')),
  'complete_phase_record_groups': ({'phase_records': (0, 3), 'original_probability': ('9/20000', '19/12500')},
                                   {'phase_records': (3, 6),
                                    'original_probability': ('121/50000', '349/100000')},
                                   {'phase_records': (6, 9),
                                    'original_probability': ('7427/100000', '3767/50000')},
                                   {'phase_records': (9, 12),
                                    'original_probability': ('3781/20000', '4753/25000')},
                                   {'phase_records': (12, 15),
                                    'original_probability': ('1601/12500', '2583/20000')},
                                   {'phase_records': (15, 18),
                                    'original_probability': ('253/20000', '1373/100000')},
                                   {'phase_records': (18, 24),
                                    'original_probability': ('67/50000', '241/100000')},
                                   {'phase_records': (24, 32),
                                    'original_probability': ('11/50000', '129/100000')},
                                   {'phase_records': (32, 64),
                                    'original_probability': ('19/100000', '127/100000')},
                                   {'phase_records': (64, 128), 'original_probability': ('0', '87/100000')},
                                   {'phase_records': (128, 256),
                                    'original_probability': ('43/25000', '279/100000')}),
  'true_spectral_windows': ({'excitation_energy_window': ('1/8', '93539/100000'),
                             'original_weight_lower': '1341/4000',
                             'phase_tail_upper': '55/338',
                             'phase_records_inclusive': (6, 15)},
                            {'excitation_energy_window': ('1/8', '127991/100000'),
                             'original_weight_lower': '17961/50000',
                             'phase_tail_upper': '11/90',
                             'phase_records_inclusive': (3, 21)},
                            {'excitation_energy_window': ('1/8', '142717/100000'),
                             'original_weight_lower': '3607/10000',
                             'phase_tail_upper': '11/90',
                             'phase_records_inclusive': (0, 24)})})
