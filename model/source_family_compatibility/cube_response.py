"""Certified composite response of the complete centered eight-site source.

The reference is the full centered ground at c=1/8 and U=1. Sparse
integer actions and rational expectations retain its entire occupation
space; no aligned-spin reference or fitted dispersion is substituted.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from copy import deepcopy

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

from model.centered_source_ground import gap
from model.interacting_vacuum import global_ground
from model.full_fock_field.field import FullFock


DIMENSION = 1 << 16
MODES = 16
SCALE = 16
INTMAX = int(np.iinfo(np.int64).max)


def _norm2(vector):
    return sum(int(x)**2 for x in vector)


def _dot(left, right):
    """Exact int64 dot after a Python-integer Cauchy prebound."""
    if _norm2(left)*_norm2(right) > INTMAX**2:
        raise OverflowError('integer response dot product exceeds its arithmetic domain')
    return int(left @ right)


def _action(matrix, vector):
    maximum = max((abs(int(x)) for x in vector), default=0)
    rows = np.asarray(abs(matrix).sum(axis=1)).ravel()
    if max((int(x) for x in rows), default=0)*maximum > INTMAX:
        raise OverflowError('integer sparse response action could overflow')
    return matrix @ vector


def upper_sqrt(value, denominator=10**8):
    """Rational outward square-root enclosure, without floating arithmetic."""
    value = F(value)
    if value < 0 or type(denominator) is not int or denominator <= 0:
        raise ValueError('nonnegative rational and positive integer grid required')
    q = value.numerator*denominator**2//value.denominator
    k = isqrt(q)
    if F(k*k, denominator**2) < value:
        k += 1
    return F(k, denominator)


def propose_trial(numerator):
    """Use the inherited deterministic proposal, without a second Gram factor."""
    h = gap._source(numerator)
    if h.shape != (4900, 4900):
        raise ValueError('complete balanced cube source required')
    initial = np.arange(1, 4901, dtype=float)
    initial /= np.linalg.norm(initial)
    _, vectors = eigsh(h.astype(float)/16, k=1, which='SA', tol=1e-11,
                      v0=initial, maxiter=10000)
    rounded = np.rint(vectors[:, 0]*(1 << gap.TRIAL_BITS))
    if not np.isfinite(rounded).all() or np.max(np.abs(rounded)) >= 1 << 30:
        raise ArithmeticError('ground proposal exceeds its exact integer domain')
    return rounded.astype(np.int64)


def verify_trial(numerator, trial):
    """Exact fidelity estimate consuming the inherited full-gap theorem."""
    h = gap._source(numerator)
    p = gap._integer_array(trial, (4900,), 'complete balanced trial')
    if h.shape != (4900, 4900) or (h-h.T).nnz:
        raise ValueError('complete symmetric balanced source required')
    if (h-global_ground.source_block(8)).nnz:
        raise ValueError('the literal original-default cube source is required')
    norm = _norm2(p)
    if not norm:
        raise ValueError('nonzero integer trial required')
    hp = _action(h, p)
    rayleigh = F(_dot(p, hp), 16*norm)-4
    variance = F(_norm2(hp), 256*norm)-(rayleigh+4)**2
    separation = gap.EXCITED_LOWER-rayleigh
    if variance < 0 or rayleigh >= gap.GROUND_UPPER or separation <= 0:
        raise ArithmeticError('trial fails the inherited ground separation')
    infidelity = variance/separation**2
    if infidelity >= gap.INFIDELITY_CEILING:
        raise ArithmeticError('trial fails the inherited fidelity ceiling')
    return {'norm_squared':norm, 'rayleigh':rayleigh,
            'variance':variance, 'infidelity_upper':infidelity,
            'half_trace_distance_upper':upper_sqrt(infidelity),
            'ground_energy_interval':(rayleigh-variance/separation, rayleigh)}


@lru_cache(maxsize=1)
def _ground_data():
    premise = gap.certificate()
    if (premise['full_Fock_dimension'] != DIMENSION or
            premise['full_ground_multiplicity'] != 1 or
            premise['full_centered_gap_strict_lower'] != '1/128'):
        raise ArithmeticError('inherited full centered ground premise changed')
    indices = global_ground.balanced_indices(8)
    h = global_ground.source_block(8)
    trial = propose_trial(h)
    data = verify_trial(h, trial)
    full = np.zeros(DIMENSION, dtype=np.int64)
    full[indices] = trial
    ops = global_ground.cube.operators()
    if ops.scale != SCALE or ops.hubbard_H.imag.nnz:
        raise ArithmeticError('actual source is not the real signed numerator16H')
    original = gap._source(ops.hubbard_H.real)
    number = np.asarray([m.bit_count() for m in range(DIMENSION)], dtype=np.int64)
    centered = original-sparse.diags(8*number, dtype=np.int64, format='csr')
    if (centered[indices][:,indices]-h+64*sparse.eye(4900,dtype=np.int64)).nnz:
        raise ArithmeticError('balanced and full centered clocks differ')
    data.update({'trial':trial, 'full_trial':full, 'balanced_indices':indices,
                 'centered_numerator':centered.tocsr(), 'scale':SCALE,
                 'full_gap_premise':premise})
    return data


def certified_ground():
    """Detached exact reference data for source response and spectral owners."""
    data = _ground_data()
    return deepcopy(data)


@lru_cache(maxsize=1)
def _operators():
    ops = global_ground.cube.operators()
    h = gap._source(ops.hubbard_H.real)
    masks = np.arange(DIMENSION,dtype=np.int64)
    doublon = np.asarray([sum((int(m)>>(2*x))&3 == 3 for x in range(8))
                          for m in masks],dtype=np.int64)
    if not np.array_equal(h.diagonal(),16*doublon):
        raise ArithmeticError('literal source onsite action differs from U=1')
    kinetic = h-sparse.diags(16*doublon,dtype=np.int64,format='csr')
    number=np.asarray([int(m).bit_count()for m in masks],dtype=np.int64)
    centered=h-sparse.diags(8*number,dtype=np.int64,format='csr')
    space = FullFock(MODES)
    a=[];b=[];y=[];norms=[]
    for mode in range(MODES):
        actual = space.annihilator(mode)
        if actual.imag.count_nonzero():
            raise ArithmeticError('real interleaved CAR convention required')
        annihilator = actual.real.astype(np.int64).tocsr()
        annihilator.eliminate_zeros()
        phase = 2*((masks>>(mode^1))&1)-1
        dressed = sparse.diags(phase,format='csr',dtype=np.int64)@annihilator
        force = kinetic@dressed.T-dressed.T@kinetic
        force.eliminate_zeros()
        identity=centered@dressed.T-dressed.T@centered-force-8*annihilator.T
        identity.eliminate_zeros()
        if identity.nnz:
            raise ArithmeticError('dressed onsite derivative differs from the original source')
        # sqrt(||Y||1 ||Y||inf) <= max(||Y||1,||Y||inf).
        row = max(map(int,np.asarray(abs(force).sum(axis=1)).ravel()),default=0)
        col = max(map(int,np.asarray(abs(force).sum(axis=0)).ravel()),default=0)
        norms.append(F(max(row,col),16))
        a.append(annihilator);b.append(dressed.tocsr());y.append(force.tocsr())
    one = kinetic[1<<np.arange(MODES)][:,1<<np.arange(MODES)].toarray()
    return {'kinetic':kinetic.tocsr(),'a':tuple(a),'b':tuple(b),'Y':tuple(y),
            'Y_norm_upper':tuple(norms),'one_numerator':one}


def _bond_numerator(left,right):
    """Actual even numerator 4(S_left.S_right+eta_left.eta_right)."""
    a=_operators()['a'];masks=np.arange(DIMENSION,dtype=np.int64)
    points=global_ground.cube.GEOMETRY.sites
    if sum(abs(x-y)for x,y in zip(points[left],points[right],strict=True))!=1:
        raise ValueError('an actual nearest-neighbour source bond required')
    n={x:tuple((masks>>(2*x+s))&1 for s in(0,1))for x in(left,right)}
    diagonal=(n[left][0]-n[left][1])*(n[right][0]-n[right][1])
    diagonal+=(n[left][0]+n[left][1]-1)*(n[right][0]+n[right][1]-1)
    value=sparse.diags(diagonal,dtype=np.int64,format='csr')
    spin={x:a[2*x].T@a[2*x+1]for x in(left,right)}
    eta={x:(-1)**sum(points[x])*(a[2*x].T@a[2*x+1].T)for x in(left,right)}
    for operators in(spin,eta):
        value+=2*(operators[left]@operators[right].T+operators[left].T@operators[right])
    value.eliminate_zeros()
    return value.tocsr()


def bond_response_on_trial(full_trial,norm):
    p=gap._integer_array(full_trial,(DIMENSION,),'complete full trial')
    if type(norm)is not int or norm<=0 or _norm2(p)!=norm:
        raise ValueError('exact positive trial norm required')
    one=_operators()['one_numerator'];bonds={}
    for x in range(8):
        for y in range(x+1,8):
            if one[2*x,2*y]:
                bonds[x,y]=F(_dot(p,_action(_bond_numerator(x,y),p)),norm)
    C=tuple(tuple(F(int(one[i,j]),16)*bonds.get(tuple(sorted((i//2,j//2))),F())
                  for j in range(MODES))for i in range(MODES))
    return {'four_spin_eta_bond_expectations':bonds,'C':C}


def _square(matrix):
    return tuple(tuple(sum((matrix[i][k]*matrix[k][j]for k in range(len(matrix))),F())
                       for j in range(len(matrix)))for i in range(len(matrix)))


def enclosure(ground_data,response,bonds):
    """Full matrix-norm enclosures from source operator bounds and fidelity."""
    delta=F(1,1920)
    if ground_data['infidelity_upper']>=delta*delta:
        raise ArithmeticError('ground fidelity misses the admitted response bound')
    C=bonds['C'];C2=_square(C);G=response['Y_Gram']
    if any(G[i][j]!=G[j][i] or C[i][j]!=C[j][i]
           for i in range(MODES)for j in range(MODES)):
        raise ArithmeticError('real symmetric response proposals required')
    R=tuple(tuple(G[i][j]-C2[i][j]for j in range(MODES))for i in range(MODES))
    # Each true C entry is h_ij times a bond observable of norm <=3.
    one=_operators()['one_numerator'];ynorm=_operators()['Y_norm_upper']
    C_error=max(sum((6*delta*abs(F(int(one[i,j]),16))for j in range(MODES)),F())
                for i in range(MODES))
    G_error=max(sum((4*delta*ynorm[i]*ynorm[j]for j in range(MODES)),F())
                for i in range(MODES))
    C_bound=F(9,16)  # 3 c degree/2, c=1/8, degree=3.
    residual_error=G_error+2*C_bound*C_error+C_error*C_error
    lower=min(R[i][i]-sum(abs(R[i][j])for j in range(MODES)if j!=i)
              for i in range(MODES))-residual_error
    upper=max(R[i][i]+sum(abs(R[i][j])for j in range(MODES)if j!=i)
              for i in range(MODES))+residual_error
    if lower<=F(3,80) or upper>=F(61,1000):
        raise ArithmeticError('complete cube residual misses its fixed positive interval')
    if any(not F(-157,100)<value-6*delta<value+6*delta<F(-155,100)
           for value in bonds['four_spin_eta_bond_expectations'].values()):
        raise ArithmeticError('actual bond coefficient misses the fixed enclosure')
    # Both even source readouts enclose the same exact-ground C. Their
    # disagreement is tested without imposing exact symmetry on the trial.
    for i in range(MODES):
        for j in range(MODES):
            allowance=4*delta*ynorm[j]+6*delta*abs(F(int(one[i,j]),16))
            if abs(response['C'][i][j]-C[i][j])>allowance:
                raise ArithmeticError('independent force and bond response readouts disagree')
    return {'trial_C':C,'trial_Y_Gram':G,'trial_residual_Gram':R,
            'half_trace_bound':delta,'C_operator_error_upper':C_error,
            'Y_Gram_operator_error_upper':G_error,'residual_operator_error_upper':residual_error,
            'residual_spectrum_interval':(F(3,80),F(61,1000)),
            'four_spin_eta_bond_interval':(F(-157,100),F(-155,100))}


def response_on_trial(full_trial, norm):
    """Rational trial expectations of the actual even response operators."""
    p=gap._integer_array(full_trial,(DIMENSION,),'complete full trial')
    if type(norm) is not int or norm<=0 or _norm2(p)!=norm:
        raise ValueError('exact positive trial norm required')
    ops=_operators()
    bp=tuple(_action(b.T,p)for b in ops['b'])
    bm=tuple(_action(b,p)for b in ops['b'])
    yp=tuple(_action(y,p)for y in ops['Y'])
    ym=tuple(_action(y.T,p)for y in ops['Y'])
    C=tuple(tuple(F(_dot(bp[i],yp[j])+_dot(ym[j],bm[i]),16*norm)
                  for j in range(MODES))for i in range(MODES))
    G=tuple(tuple(F(_dot(yp[i],yp[j])+_dot(ym[j],ym[i]),256*norm)
                  for j in range(MODES))for i in range(MODES))
    return {'C':C,'Y_Gram':G}


@lru_cache(maxsize=1)
def _certified_response():
    ground=certified_ground();p=ground['full_trial'];norm=ground['norm_squared']
    response=response_on_trial(p,norm)
    bonds=bond_response_on_trial(p,norm)
    return enclosure(ground,response,bonds)


def certificate():
    """Fixed admitted rational source results, without solver-noise digits."""
    from .composite import approximation_budget
    data=_certified_response();one=_operators()['one_numerator']
    interval=data['four_spin_eta_bond_interval']
    C=[]
    for row in one:
        C.append(tuple(tuple(sorted(F(int(x),16)*q for q in interval)) if x
                       else (F(),F())for x in row))
    lo,hi=data['residual_spectrum_interval']
    return {'source':'original-default signed Hubbard chart of the complete Weyl/contact cube',
            'sides':[2,2,2],'c':F(1,8),'U':F(1),'centered_clock':'K=H-N/2',
            'full_occupation_dimension':DIMENSION,'balanced_ground_basis_dimension':4900,
            'original_retarded_modes':MODES,'source_sites':8,'source_bonds':12,
            'hidden_spin_and_staggered_eta_invariant_ground':True,
            'full_ground_gap_strict_lower':F(1,128),
            'trial_ground_half_trace_distance_strict_upper':data['half_trace_bound'],
            'four_spin_eta_bond_expectation_interval':interval,
            'C_entry_intervals':tuple(C),
            'C_operator_error_upper':data['C_operator_error_upper'],
            'Y_Gram_operator_error_upper':data['Y_Gram_operator_error_upper'],
            'residual_operator_error_upper':data['residual_operator_error_upper'],
            'residual_squared_operator_interval':(lo,hi),
            'residual_squared_trace_interval':(16*lo,16*hi),
            'fourth_original_moment_defect_operator_interval':(lo/4,hi/4),
            'cyclic_response_dimension_lower':48,
            'each_local_original_addition_removal_weight':F(1,2),
            'normalized_smearing_total_signed_weight':F(1),
            'all_smearings_have_equal_charge_weights':False,
            'original_response_unit_time_error':approximation_budget(hi,F(1),F(1),F(1)),
            'nonzero_residual_in_every_original_mode_direction':True,
            'all_occupation_sectors_in_source_clock':True,
            'complete_local_source_operators_checked':True,
            'numerical_eigenvalues_used_for_admission':False,
            'native_ground_preparation_executed':False,
            'retarded_compression_is_arbitrary_reference_physical_channel':False,
            'extra_independent_CAR_flavours':0,'physical_species_added':0,
            'physical_mass':None,'nonzero_momentum_pole_identified':False,
            'full_ground_continuum_limit_identified':False}
