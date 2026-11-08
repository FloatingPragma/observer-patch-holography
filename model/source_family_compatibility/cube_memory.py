"""Certified frequency-dependent memory of the readable-default full cube.

Sparse linear solves supply proposals. Exact rational residuals and the
same full-ground certificate determine every retained interval.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import lcm
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh,splu
from . import cube_response as reference
from . import cube_spectrum
from . import momentum_response as algebra

P_BITS=31
X_BITS=40
Q_BITS=40
GROUND_HALFTRACE=F(1,400000)
ENERGY_REFERENCE=F(-4124058733378,10**12)
ENERGY_ERROR=F(1,1<<38)
RESIDUAL_CEILING=F(1,1<<27)
SOURCE_SCALE=16


def _norm2(v):return sum(int(x)**2 for x in np.asarray(v).ravel())
def _dot(x,y):return sum(int(a)*int(b)for a,b in zip(x,y,strict=True))
def _sqrt(x):return reference.upper_sqrt(x,10**15)
def _complex(a=0,b=0):return(F(a),F(b))
def _add(a,b):return(a[0]+b[0],a[1]+b[1])
def _mul(a,b):return(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def _neg(a):return(-a[0],-a[1])
def _sub(a,b):return _add(a,_neg(b))
def _abs_upper(a):return abs(a[0])+abs(a[1])


def _matrix_bound(a):
    rows=max(sum((_abs_upper(x)for x in row),F())for row in a)
    cols=max(sum((_abs_upper(a[i][j])for i in range(len(a))),F())for j in range(len(a)))
    return max(rows,cols)


def inverse_proposal(a):
    """Dyadic inverse with an exact residual; no LAPACK accuracy assumption."""
    value=np.asarray([[complex(float(x),float(y))for x,y in row]for row in a])
    inverse=np.linalg.inv(value)
    qr=np.rint(inverse.real*(1<<Q_BITS));qi=np.rint(inverse.imag*(1<<Q_BITS))
    if not np.isfinite(qr).all()or not np.isfinite(qi).all()or max(abs(qr).max(),abs(qi).max())>=2**62:
        raise ArithmeticError('inverse proposal exceeds its integer domain')
    q=tuple(tuple((F(int(x),1<<Q_BITS),F(int(y),1<<Q_BITS))for x,y in zip(r,s,strict=True))for r,s in zip(qr,qi,strict=True))
    residual=[]
    for i in range(len(a)):
        row=[]
        for j in range(len(a)):
            product=sum_complex(_mul(a[i][k],q[k][j])for k in range(len(a)))
            row.append(_sub(_complex(i==j),product))
        residual.append(tuple(row))
    return q,_matrix_bound(q),_matrix_bound(residual)


def sum_complex(items):
    value=_complex()
    for item in items:value=_add(value,item)
    return value


@lru_cache(maxsize=1)
def _source():
    data=reference.certified_ground()
    indices=data['balanced_indices'];neutral=reference.global_ground.source_block(8)
    initial=data['trial'].astype(float);initial/=np.linalg.norm(initial)
    _,vectors=eigsh(neutral.astype(float)/16,k=1,which='SA',tol=1e-13,v0=initial,maxiter=10000)
    rounded=np.rint(vectors[:,0]*(1<<P_BITS))
    if not np.isfinite(rounded).all()or np.max(np.abs(rounded))>=2**32:
        raise ArithmeticError('refined ground proposal exceeds its domain')
    p=rounded.astype(np.int64);cert=reference.verify_trial(neutral,p)
    if cert['half_trace_distance_upper']>GROUND_HALFTRACE:
        raise ArithmeticError('refined ground fails its fixed trace enclosure')
    if max(abs(x-ENERGY_REFERENCE)for x in cert['ground_energy_interval'])>ENERGY_ERROR:
        raise ArithmeticError('actual ground energy misses fixed clock enclosure')
    full=np.zeros(1<<16,dtype=np.int64);full[indices]=p
    masks=np.asarray([m for m in range(1<<16)if m.bit_count()==9 and(m&0x5555).bit_count()==5],dtype=np.int64)
    if len(masks)!=3920:raise ArithmeticError('complete charged source dimension differs')
    rows={int(m):i for i,m in enumerate(masks)};b=np.zeros((3920,16),dtype=np.int64)
    for col,mask in enumerate(indices):
        mask=int(mask)
        for x in range(8):
            mode=2*x
            if mask&(1<<mode):continue
            sign=(-1)**(mask&((1<<mode)-1)).bit_count();row=rows[mask|(1<<mode)]
            b[row,x]+=sign*p[col]
            b[row,8+x]+=sign*p[col]*(2*((mask>>(mode+1))&1)-1)
    centered=data['centered_numerator'];h=centered[masks][:,masks].tocsr()
    if centered[:,masks].nnz!=h.nnz or(h-h.T).nnz:raise ArithmeticError('actual charged clock is not reducing')
    if h.dtype!=np.dtype(np.int64):raise ArithmeticError('integer source clock required')
    # Complete full-source conjugation verifies the reused reversed removal block.
    symmetry=cube_spectrum.symmetry_certificate(data)
    sites=reference.global_ground.cube.GEOMETRY.sites
    eps=tuple((-1)**sum(x)for x in sites);signs=eps+tuple(-x for x in eps)
    one=reference._operators()['one_numerator'][::2,::2]
    norm=cert['norm_squared'];bonds={}
    for x in range(8):
        for y in range(x+1,8):
            if one[x,y]:
                op=reference._bond_numerator(x,y);bonds[x,y]=F(_dot(full,op@full),norm)
    c=tuple(tuple(F(int(one[i,j]),16)*bonds.get(tuple(sorted((i,j))),F())for j in range(8))for i in range(8))
    a=tuple(tuple(_complex(F(int(one[i,j]),16)if i<8 and j<8 else
                                 c[i-8][j-8]if i>=8 and j>=8 else
                                 F(1,2)if abs(i-j)==8 else F())for j in range(16))for i in range(16))
    gram=b.astype(object).T@b.astype(object)
    gram=tuple(tuple(F(int(gram[i,j])*(1+signs[i]*signs[j]),norm)for j in range(16))for i in range(16))
    error=max(sum((abs(gram[i][j]-(i==j))for j in range(16)),F())for i in range(16))
    if error>4*GROUND_HALFTRACE+4*GROUND_HALFTRACE**2:raise ArithmeticError('trial response Gram differs from actual orthogonal embeddings')
    return {'ground':cert,'trial':p,'balanced_indices':indices,'charged_indices':masks,'b':b,'clock_numerator':h,'signs':signs,'A':a,'symmetry':symmetry,
            'norm_squared':norm,'C_error':F(9,8)*GROUND_HALFTRACE}


def _rounded_source_solve(data,eta):
    h=data['clock_numerator'];b=data['b'];n=len(b)
    matrix=(complex(float(ENERGY_REFERENCE),float(eta))*sparse.eye(n)-h.astype(float)/16).tocsc()
    x=splu(matrix).solve(b.astype(float)/(1<<P_BITS))
    rr=np.rint(x.real*(1<<X_BITS));ii=np.rint(x.imag*(1<<X_BITS))
    if not np.isfinite(rr).all()or not np.isfinite(ii).all()or max(abs(rr).max(),abs(ii).max())>=2**50:
        raise ArithmeticError('source solve proposal exceeds integer domain')
    xr,xi=rr.astype(np.int64),ii.astype(np.int64)
    rowbound=max(map(int,np.asarray(abs(h).sum(axis=1)).ravel()))
    if rowbound*max(abs(xr).max(),abs(xi).max())>=2**63:
        raise OverflowError('integer source residual action could overflow')
    hr,hi=h@xr,h@xi
    den=lcm(ENERGY_REFERENCE.denominator,eta.denominator,16)
    ee=int(ENERGY_REFERENCE*den);et=int(eta*den)
    real=ee*xr.astype(object)-(den//16)*hr.astype(object)-et*xi.astype(object)-b.astype(object)*(den*(1<<(X_BITS-P_BITS)))
    imag=ee*xi.astype(object)-(den//16)*hi.astype(object)+et*xr.astype(object)
    squared=F(_norm2(real)+_norm2(imag),(den*(1<<X_BITS))**2)
    residual=_sqrt(squared)
    if residual>RESIDUAL_CEILING:raise ArithmeticError('rounded source solve fails its fixed rational residual')
    if not F(3,4)<=F(data['norm_squared'],1<<(2*P_BITS))<=F(5,4):
        raise ArithmeticError('refined source normalization misses its analytic domain')
    gr=b.astype(object).T@xr.astype(object);gi=b.astype(object).T@xi.astype(object)
    scale=F(1<<P_BITS,(1<<X_BITS)*data['norm_squared'])
    signs=data['signs'];g=[]
    for i in range(16):
        row=[]
        for j in range(16):
            sg=signs[i]*signs[j]
            row.append((scale*int(gr[i,j])*(1-sg),scale*int(gi[i,j])*(1+sg)))
        g.append(tuple(row))
    return tuple(g),residual


@lru_cache(maxsize=8)
def _frequency(eta):
    if type(eta)not in(int,F)or eta<=0:raise ValueError('positive rational source frequency required')
    eta=F(eta);d=_source();g,residual=_rounded_source_solve(d,eta)
    ep=3*GROUND_HALFTRACE
    gerror=(2*ep+ep*ep)/eta+(1+ep)**2*ENERGY_ERROR/eta**2+4*RESIDUAL_CEILING/eta
    q,qbound,qres=inverse_proposal(g);defect=qres+gerror*qbound
    if defect>=1:raise ArithmeticError('source resolvent inverse is not certified')
    inverse_error=qbound*defect/(1-defect)
    memory_error=d['C_error']+inverse_error
    memory=tuple(tuple(_sub(_sub(_complex(0,eta if i==j else 0),d['A'][i][j]),q[i][j])for j in range(16))for i in range(16))
    cheap_matrix=tuple(tuple(_sub(_complex(0,eta if i==j else 0),d['A'][i][j])for j in range(16))for i in range(16))
    cheap,cheap_bound,cheap_res=inverse_proposal(cheap_matrix)
    if cheap_res>=1:raise ArithmeticError('compressed comparison inverse fails')
    cheap_error=cheap_bound*cheap_res/(1-cheap_res)+d['C_error']/eta**2
    original=-eta*g[0][0][1];comparison=-eta*cheap[0][0][1]
    return {'eta':eta,'memory_00_imag_interval':(memory[8][8][1]-memory_error,memory[8][8][1]+memory_error),
            'memory_operator_error':memory_error,'original_filter_interval':(original-eta*gerror,original+eta*gerror),
            'compression_filter_interval':(comparison-eta*cheap_error,comparison+eta*cheap_error),
            'difference_interval':(original-comparison-eta*(gerror+cheap_error),original-comparison+eta*(gerror+cheap_error)),
            'source_linear_solve_residual':residual}


def frequency(eta):
    if type(eta)not in(int,F)or eta<=0:raise ValueError('positive rational source frequency required')
    return reference._detached(_frequency(F(eta)))


def certificate():
    first,second=frequency(F(1,4)),frequency(F(1,8))
    intervals={F(1,4):{'memory':(F(-395,10000),F(-392,10000)),
                      'original':(F(2175,10000),F(2177,10000)),
                      'compression':(F(1943,10000),F(1945,10000)),
                      'difference':(F(23,1000),F(234,10000))},
               F(1,8):{'memory':(F(-239,10000),F(-235,10000)),
                      'original':(F(667,10000),F(669,10000)),
                      'compression':(F(569,10000),F(571,10000)),
                      'difference':(F(97,10000),F(100,10000))}}
    rows=[]
    for row in(first,second):
        fixed=intervals[row['eta']]
        for key,target in [('memory_00_imag_interval','memory'),('original_filter_interval','original'),
                           ('compression_filter_interval','compression'),('difference_interval','difference')]:
            lo,hi=row[key];flo,fhi=fixed[target]
            if not flo<=lo<=hi<=fhi:raise ArithmeticError('actual source response misses fixed '+target+' enclosure')
        rows.append({'eta':row['eta'],'memory_00_imag_interval':fixed['memory'],
                     'original_complete_filter_probability':fixed['original'],
                     'memory_free_filter_probability':fixed['compression'],
                     'original_minus_memory_free_probability':fixed['difference']})
    if first['memory_00_imag_interval'][1]>=second['memory_00_imag_interval'][0]-F(15,1000):
        raise ArithmeticError('actual omitted response does not resolve the frequency dependence')
    return {'source_c':F(1,8),'source_U':F(1),'full_Fock_dimension':65536,
            'complete_ground_block_dimension':4900,'charged_spin_block_dimension':3920,
            'one_spin_retarded_two_block_columns':16,'both_spin_response_equivalence_proved':True,'refined_ground_dyadic_bits':P_BITS,
            'full_ground_halftrace_ceiling':GROUND_HALFTRACE,'ground_energy_reference':ENERGY_REFERENCE,
            'ground_energy_error_ceiling':ENERGY_ERROR,'source_solve_residual_ceiling':RESIDUAL_CEILING,
            'source_clock_and_ground_unchanged':True,'full_PH_and_spin_ground_symmetry_consumed':True,
            'frequencies':rows,'imaginary_memory_change_strict_lower':F(15,1000),
            'local_mode_charge_filter_records':4,'all_failure_complements_retained':True,
            'native_large_filter_tape_executed':False,'arbitrary_reference_channel_for_retarded_columns_claimed':False,
            'thermodynamic_charge_gap_proved':False,'physical_mass_assigned':None}
