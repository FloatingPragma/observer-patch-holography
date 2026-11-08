"""Certified complete finite-source spectra on the joint field/clock parameters."""
from fractions import Fraction as F
from pathlib import Path
from math import ceil
import numpy as np
from model.native_ordered_ground import source
from model.centered_particle_spectrum.certificate import sqrt_bounds,decimal_bounds
from model.spin_dynamical_memory import zero_field
from . import uniform_mode

ROOT=Path(__file__).resolve().parents[2]
WITNESS='model/spin_joint_field_limit/uniform_mode_witness.npz'
BITS=48
SCALES=(2,3,4)


def propose(path=ROOT/WITNESS):
    data={};unit=1<<BITS
    for n in SCALES:
        h,scale=source.matrix(field=F(1,n**5))
        e,v=np.linalg.eigh(h/scale)
        if v[:,0].sum()<0:v[:,0]*=-1
        data[f'e{n}']=np.rint(e*unit).astype(np.int64)
        data[f'v{n}']=np.rint(v*unit).astype(np.int64)
    np.savez_compressed(path,**data)


def verify(n,energies,vectors):
    h,scale=source.matrix(field=F(1,n**5));size=len(h);d=1<<BITS
    if energies.shape!=(size,) or vectors.shape!=(size,size):raise ValueError('complete source eigensystem required')
    if energies.dtype!=np.int64 or vectors.dtype!=np.int64:raise ValueError('literal integer proposals required')
    e=energies.astype(object);v=vectors.astype(object)
    if any(x>y for x,y in zip(e,e[1:])):raise ValueError('ordered source eigenvalues required')
    norm_bound=F(12)+F(4,n**5)
    if any(abs(F(int(x),d))>norm_bound for x in e):raise ValueError('source norm bound exceeded')
    gram=v.T@v-d*d*np.eye(size,dtype=object)
    delta=F(max(sum(abs(int(x))for x in row)for row in gram),d*d)
    if delta>=F(1,10**6):raise ValueError('complete Gram certificate failed')
    residual=h.astype(object)@v*d-scale*v*e[None,:]
    epsilon=sqrt_bounds(F(sum(int(x)**2 for x in residual.flat),scale**2*d**4))[1]
    radius=epsilon+2*norm_bound*delta
    e0,e1=F(int(e[0]),d),F(int(e[1]),d)
    norm=F(sum(int(x)**2 for x in v[:,0]),d*d)
    r0=sqrt_bounds(F(sum(int(x)**2 for x in residual[:,0]),scale**2*d**4)/norm)[1]
    if e1-radius-e0-r0<=0:raise ValueError('source ground not isolated')
    error=2*r0/(e1-radius-e0)
    if error>F(1,10**8):raise ValueError('source ground error too large')
    return {'dimension':size,'vector_bits':BITS,'gram_defect':delta,
            'full_polar_spectrum_radius':radius,'ground_residual':r0,
            'ground_vector_error':error,'ground_energy_interval':[e0-r0,e0+r0],
            'source_gap_interval':[e1-radius-e0-r0,e1+radius-e0+r0],
            'source_operator_norm_upper':norm_bound,
            'every_source_eigenvalue_certified':True,'floating_solver_trusted':False}


def operators(h):
    size=256;indices=np.arange(size);identity=np.eye(size,dtype=np.int64)
    signs=tuple(-1 if x in source.cube()[1]else 1 for x in range(8))
    charge=sum((sign*identity[indices^(1<<x)]for x,sign in enumerate(signs)),np.zeros_like(identity))
    z=np.diag([1-2*(word&1)for word in indices]);xz=z[indices^1,:]
    creators={q:z+q*xz for q in(-1,1)}
    if not np.array_equal(charge@h,h@charge):raise ArithmeticError('actual charge does not commute with source')
    for q,a in creators.items():
        if not np.array_equal(charge@a-a@charge,2*q*a):raise ArithmeticError('actual source creation charge changed')
    if not np.array_equal(sum(a.T@a for a in creators.values()),4*identity):raise ArithmeticError('original creation records incomplete')
    group=(0,3,5,6)
    permutations={a:np.array([sum(((word>>j)&1)<<(j^a)for j in range(8))for word in indices])for a in group}
    projectors=[]
    for k in range(4):
        p=sum(((-1)**((k&a).bit_count())*identity[permutations[a]]for a in group),np.zeros_like(identity))
        if not np.array_equal(p@h,h@p) or not np.array_equal(p@charge,charge@p):raise ArithmeticError('source reflection sector not reducing')
        projectors.append(p)
    if not np.array_equal(sum(projectors),4*identity):raise ArithmeticError('source reflection effects incomplete')
    for i,p in enumerate(projectors):
        for j,q in enumerate(projectors):
            if not np.array_equal(p@q,4*p if i==j else np.zeros_like(p)):raise ArithmeticError('source reflection projections changed')
    m=np.array([sum(signs[x]*(1-2*((word>>x)&1))for x in range(8))for word in indices],dtype=np.int64)
    commutator=charge@np.diag(m)-np.diag(m)@charge
    if not np.array_equal(charge@commutator-commutator@charge,4*np.diag(m)):
        raise ArithmeticError('uniform source preparation has charges beyond plus/minus one')
    if any(not np.array_equal(m[p],m)for p in permutations.values()):
        raise ArithmeticError('uniform source preparation has a nontrivial even-reflection character')
    odd=np.array([sum(((word>>j)&1)<<(j^4)for j in range(8))for word in indices])
    if not np.array_equal(h[np.ix_(odd,odd)],h) or not np.array_equal(charge[np.ix_(odd,odd)],-charge) or not np.array_equal(m[odd],-m):raise ArithmeticError('actual charge-pairing source symmetry changed')
    return creators,projectors,m


def ratio_interval(x,y):
    if y[0]<=0:raise ArithmeticError('positive moment floor missing')
    return min(a/b for a in x for b in y),max(a/b for a in x for b in y)


def _uniform(checked,energies,vectors,m):
    d=1<<BITS;v=vectors.astype(object);g=v[:,0];norm=int(g@g)
    gaps=[F(int(e)-int(energies[0]),d)for e in energies]
    # M/(N/2)=diag(m/8), a contraction with every failure retained.
    image=m.astype(object)*g;coeff=v.T@image
    weights=[F(int(c)**2,64*d*d*norm)for c in coeff]
    eg=checked['ground_vector_error'];delta=checked['gram_defect']
    gap=checked['source_gap_interval'][0];a=min(F(1),gap)/2
    perturb=checked['full_polar_spectrum_radius']+checked['ground_residual']+2*eg+4*delta
    cap=2*checked['source_operator_norm_upper']+1
    shape=2*eg+2*delta+delta*delta
    moments={}
    for power in(1,-1,-3):
        value=sum(w*(e if i else F(1))**power for i,(w,e)in enumerate(zip(weights,gaps)))
        bound=cap if power==1 else a**power
        lipschitz=F(1) if power==1 else abs(power)*a**(power-1)
        error=lipschitz*perturb+shape*bound
        moments[power]=(max(F(),value-error),value+error)
    mean=ratio_interval(moments[-1],moments[-3])
    second=ratio_interval(moments[1],moments[-3])
    variance=(max(F(),second[0]-mean[1]**2),second[1]-mean[0]**2)
    if variance[0]<=F(1,100):raise ArithmeticError('finite source weighted variance not certified positive')
    weight=F(int(image@image),64*norm);wi=(max(F(),weight-2*eg),min(F(1),weight+2*eg))
    return {'preparation':'M3/4, with its complementary outcome',
            'original_preparation_probability_interval':decimal_bounds(*wi,12),
            'unscaled_original_moments':{str(k):decimal_bounds(*x,12)for k,x in moments.items()},
            'unscaled_weighted_energy_square_mean_interval':decimal_bounds(*mean,12),
            'unscaled_weighted_energy_square_variance_interval':decimal_bounds(*variance,12),
            'each_charge_original_probability_interval':decimal_bounds(wi[0]/2,wi[1]/2,12),
            'charge_pair_exactly_equal_spectral_measures':True,
            'uniform_preparation_has_only_charge_plus_minus_one':True,
            'even_cube_reflection_character':0,
            'normalized_preparation_does_not_replace_original_probability':True}


def finite_case(n,energies,vectors):
    checked=verify(n,energies,vectors);h,scale=source.matrix(field=F(1,n**5))
    creators,projectors,m=operators(h);unit=1<<BITS;g=vectors[:,0].astype(object);norm=int(g@g)
    delta=checked['gram_defect'];eg=checked['ground_vector_error']
    eta=checked['full_polar_spectrum_radius']+checked['ground_residual']
    gaps=[F(int(e)-int(energies[0]),unit)for e in energies]
    last=ceil(2*checked['source_operator_norm_upper'])+1
    cuts=[F(1,17)+i for i in range(last)]
    distances=[min(abs(e-cut)for e in gaps)-eta for cut in cuts]
    if min(distances)<=0:raise ArithmeticError('fixed energy bin boundary not certified disjoint')
    # Half-line spectral projections: leakage of 256 approximate eigenvectors.
    cut_errors=[min(F(1),32*eta/distance)for distance in distances]
    intervals=[(F(),cuts[0])]+list(zip(cuts,cuts[1:]))+[(cuts[-1],F(last+1))]
    bands=[]
    for i,(lo,hi)in enumerate(intervals):
        selected=[j for j,e in enumerate(gaps)if (i==0 or e>lo)and e<=hi]
        pe=(cut_errors[i-1]if i else F())+(cut_errors[i]if i<len(cut_errors)else F())
        bands.append((lo,hi,selected,pe))
    channels=[];all_low=F();all_high=F()
    for q,a in sorted(creators.items()):
        for k,p in enumerate(projectors):
            image=p.astype(object)@(a.astype(object)@g)
            coeff=vectors.T.astype(object)@image
            weights=[F(int(c)**2,64*unit**2*norm)for c in coeff]
            records=[]
            for lo,hi,selected,pe in bands:
                value=sum(weights[j]for j in selected)
                error=2*eg+2*delta+delta*delta+pe
                # Empty spectral ranges, or just the charge-zero ground, have exact zero weight.
                wi=(F(),F())if not any(j for j in selected)else(max(F(),value-error),min(F(1),value+error))
                all_low+=wi[0];all_high+=wi[1]
                records.append({'scaled_energy_interval':[str(n*n*lo),str(n*n*hi)],
                                'original_probability_interval':decimal_bounds(*wi,12),
                                'full_source_levels_in_bin':len(selected)})
            channels.append({'source_charge':q,'even_cube_reflection_character':k,'records':records})
    if not all_low<=1<=all_high:raise ArithmeticError('complete source probability enclosures lost normalization')
    uniform=_uniform(checked,energies,vectors,m)
    for key in('mean','variance'):
        values=tuple(F(x)for x in uniform[f'unscaled_weighted_energy_square_{key}_interval'])
        factor=n**(4 if key=='mean'else 8)
        uniform[f'scaled_weighted_energy_square_{key}_interval']=decimal_bounds(values[0]*factor,values[1]*factor,10)
    return {'n':n,'field':str(F(1,n**5)),'p':str(F(1,n*n)),
            'source':'actual 12-edge eight-spin open cube',
            'clock':'(H_B-E_B)/(Jp), J=1','certificate':checked,
            'uniform_response':uniform,'local_charge_reflection_energy_readout':channels,
            'complete_original_records':8*len(bands),'all_zero_slots_retained':True,
            'reflection_characters_are_continuum_momentum':False,
            'every_charge_energy_and_reflection_complement_retained':True,
            'native_spectral_pulse_tape_executed':False,'physical_mass':None}


def certificate():
    with np.load(ROOT/WITNESS,allow_pickle=False)as data:
        if set(data.files)!={f'{kind}{n}'for n in SCALES for kind in('e','v')}:raise ValueError('all source scale witnesses required')
        rows=[finite_case(n,data[f'e{n}'],data[f'v{n}'])for n in SCALES]
    with np.load(zero_field.ROOT/zero_field.WITNESS,allow_pickle=False)as data:
        zero=zero_field.verify_witness({k:data[k]for k in data.files})
    gap=F(zero['gap_interval'][0])
    return {'finite_cases':rows,'same_cube_zero_field_certificate':zero,
            'fixed_box_scaled_gap_lower':{str(n):str(uniform_mode.fixed_box_gap(n,gap))for n in SCALES},
            'fixed_box_table_substituted_for_growing_volume_field':False,
            'uniform_positive_momentum_shell_proved':False,
            'physical_mass_prediction_count':0}
