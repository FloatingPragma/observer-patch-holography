from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import sys,json
ROOT=Path(__file__).resolve().parents[2]
import numpy as np
from model.protected_particle_readout import circuits,execute,bounds
from model.cube_eta_protection import source,repair
from model.certified_reference import interval as iv
I=iv.Interval;zero=I.point(0);one=I.point(1)
def sqrtup(x):return F(iv.sqrt_fraction(x).hi,iv.SCALE)
def entry(z):return I.point(F(float(complex(z).real))),I.point(F(float(complex(z).imag)))
def matrix_error(actual,target):
 s=F()
 for i in range(actual.shape[0]):
  for j in range(actual.shape[1]):
   r,t=entry(actual[i,j]);a,b=target[i][j]
   s+=(r-a).absolute_upper()**2+(t-b).absolute_upper()**2
 return sqrtup(s)
def blank(d):return [[(zero,zero) for j in range(d)]for i in range(d)]
def ident(d):
 a=blank(d)
 for i in range(d):a[i][i]=(one,zero)
 return a
def field_target(kind,z=0):
 a=ident(4)
 if kind=='fswap':
  a=blank(4);a[0][0]=(one,zero);a[1][2]=a[2][1]=(one,zero);a[3][3]=(-one,zero)
 elif kind=='contact':
  s,c=iv.sine_cosine(-F(1,128));a[3][3]=(c,s)
 else:
  s,c=iv.sine_cosine(F(1,2048));a[1][1]=a[2][2]=(c,zero)
  w=complex(z)/abs(z)
  a[1][2]=(s*int(w.imag),-s*int(w.real));a[2][1]=(-s*int(w.imag),-s*int(w.real))
 return a
from functools import lru_cache

@lru_cache(None)
def certificate():
    root=iv.sqrt_fraction(F(1,2));pointer=blank(8)
    for j in range(8):
     for i in (j,j^4):
      sign=-1 if j&4 and not i&4 else 1
      phase=(1j)**((i&3).bit_count()) if i&4 else 1
      pointer[i][j]=(root*int(sign*phase.real),root*int(sign*phase.imag))
    blocks=circuits.field_binding(F(1,2),32)
    from model.common_source.chiral_candidate import gaussian_source as g
    from model.interacting_field import source as inter
    lookup={execute.key(g.local_gate('fswap').word):('fswap',field_target('fswap'))}
    for group in circuits.GEOMETRY.groups:
     for _,_,z in group:lookup[execute.key(g.local_gate('hopping',z,1/128).word)]=('hopping_'+str(z),field_target('hopping',z))
    _,c=inter.contact(1/128,1.);lookup[execute.key(c['word'])]=('contact',field_target('contact'))
    # The direct interval comparison includes local expm, coefficient, pi-angle,
    # and local composition error. It assumes no accuracy from these algorithms.
    counts=Counter(i.tape for i,_ in blocks);rows=[];betas=[]
    u=F(1,2**53)
    def gamma(k):return k*u/(1-k*u)
    for tape,count in counts.items():
     instrument=next(i for i,_ in blocks if i.tape==tape)
     name,target=lookup[tape];delta=matrix_error(instrument.success,target)
     beta=delta+gamma(8*4)*2*(1+delta)
     betas.append(2*count*beta)
     rows.append({'kind':name,'blocks_both_segments':2*count,'exact_dyadic_Frobenius_error_upper':str(delta),'error_decimal':float(delta),'matrix_product_error_coefficient':str(gamma(32)*2*(1+delta))})
    p=circuits.block(circuits.pointer_prepare(),3,(0,1,16))[0]
    pdelta=matrix_error(p.success,pointer);pbeta=pdelta+gamma(64)*sqrtup(F(8))*(1+pdelta)
    betas.append(pbeta)
    analyzer_errors=[]
    for phase,phase_factor in ((F(),1),(F(1,2),1j)):
        target=blank(8)
        for j in range(8):
            w=phase_factor*((-1j)**((j&3).bit_count())) if j&4 else 1
            for i in (j,j^4):
                sign=-1 if not j&4 and i&4 else 1
                target[i][j]=(root*int(sign*w.real),root*int(sign*w.imag))
        instrument=circuits.block(circuits.pointer_analyzer(phase),3,(0,1,16))[0]
        analyzer_errors.append(str(matrix_error(instrument.success,target)))
    # Directly enclose the actual copied input against exact R/sqrt56.
    copied=source.runs()[1].success[:,0];r=iv.sqrt_fraction(F(1,56));desired={sum(3<<(2*j)for j in range(8)if mask>>j&1):value for mask,value in enumerate(source.reference_values()) if value}
    squared=F()
    for i,z in enumerate(copied):
     a,b=entry(z);target=r*desired.get(i,0)
     squared+=(a-target).absolute_upper()**2+b.absolute_upper()**2
    input_error=sqrtup(squared)
    # Certify numerical V against the exact radical-normalized embedding.
    a,b=repair.exact_embedding();v=repair.embedding().tocoo();s17=iv.sqrt_fraction(F(17));den=4+4*(4-s17)*(4-s17)
    inv=I(iv.ceiling(iv.SCALE*iv.SCALE,den.hi)-1,iv.ceiling(iv.SCALE*iv.SCALE,den.lo))
    # Lower reciprocal floor is used; the -1 is harmless conservative outward rounding.
    err2=F();aa=a.tocsr();bb=b.tocsr()
    assert set(zip(v.row.tolist(),v.col.tolist())) == (set(zip(*a.nonzero()))|set(zip(*b.nonzero())))
    for i,j,z in zip(v.row,v.col,v.data,strict=True):
     av,bv=complex(aa[i,j]),complex(bb[i,j]);target=((int(av.real)+s17*int(bv.real))*inv,(int(av.imag)+s17*int(bv.imag))*inv)
     r0,i0=entry(z);err2+=(r0-target[0]).absolute_upper()**2+(i0-target[1]).absolute_upper()**2
    vdelta=sqrtup(err2);vf=sqrtup(F(24))+vdelta
    k=max(np.diff(repair.embedding().tocsc().indptr));l=max(np.diff(repair.embedding().tocsr().indptr))
    # fl(V* fl(V^*x)): multiply errors plus exact embedding error. All arithmetic
    # bounds use standard no-overflow binary64 complex multiply-add, gamma_(8d).
    cerr=gamma(8*int(k))*vf
    rerr=gamma(8*int(l))*vf*(1+vdelta+cerr)+(1+vdelta)*cerr
    qbeta=vdelta*(2+vdelta)+rerr+gamma(4)*(1+(1+vdelta)**2+rerr)
    betas.append(10*qbeta)
    # Fewer than2^80 primitive real operations occur in the finite state path;
    # even2^-1022 absolute underflow per operation contributes below2^-942.
    # This2^-900 allowance also covers the final sparse projection and norm.
    assert len(blocks)<2**20 and (1<<17)<2**20
    underflow=F(1,2**900)
    S=sum(betas,F())+underflow;assert S<F(1,1000000)
    E=(input_error+S)/(1-S)
    # One compressed native output can only see the propagated error plus
    # the V error and the certified sparse-dot roundoff. Raw normalizations are
    # included conservatively by factor2 after the exact matrix-effect identity.
    compressed=(1+vdelta)*E+vdelta+cerr*(1+E)
    amp=2*(1+vdelta+gamma(8*int(l))*vf)*compressed
    prob=amp**2*(1+gamma(8*131072))+underflow
    result={'schema':'oph.protected_readout_binary64_bound.v1','assumption':'IEEE round-to-nearest ordinary complex dot products; gamma_(8d), no overflow; underflow bounded by an explicit2^-900 total allowance',
     'underflow_absolute_allowance':str(underflow),'local_success_block_bounds':rows,'pointer_success_error_upper':str(pdelta),'analyzer_success_error_uppers':analyzer_errors,'copied_initial_state_error_upper':str(input_error),
     'V_embedding_error_upper':str(vdelta),'V_max_nonzeros_per_column':int(k),'V_max_nonzeros_per_row':int(l),
     'sparse_outside_projector_step_error_coefficient':str(qbeta),
     'global_prefix_state_error_upper':str(E),'global_prefix_state_error_decimal':float(E),
     'per_terminal_raw_probability_numerical_enclosure_upper':str(prob),'per_terminal_raw_probability_numerical_enclosure_upper_decimal':float(prob),
     'bound_is_for_this_executed_backend_and_exact_target_Strang_template':True,'fresh_backend_matrix_bounds_recomputed':True,
     'physical_hardware_error_inferred':False,'replay_comparator_changed':False}
    return result


def primitive_error(kind, angle):
    """Outward bound against the exact six-level proper-code control."""
    from model.observer.source_scalar import compiler
    a=ident(6);s,c=iv.sine_cosine(F(angle));r=iv.sqrt_fraction(F(1,2))
    if kind in ('X','Y','Z','I','H'):
        a[0][0]=a[1][1]=(c,zero)
        if kind=='X':a[0][1]=a[1][0]=(zero,s)
        elif kind=='Y':a[0][1]=(s,zero);a[1][0]=(-s,zero)
        elif kind=='Z':a[0][0]=(c,s);a[1][1]=(c,-s)
        elif kind=='I':a[0][0]=a[1][1]=(c,s)
        elif kind=='H':
            a[0][0]=(c,s*r);a[1][1]=(c,-s*r);a[0][1]=a[1][0]=(zero,s*r)
    elif kind=='J':
        for j in range(3):a[j][j]=(c,s)
    else:raise ValueError('declared primitive required')
    return matrix_error(compiler.native_pulse(kind,angle),a)


@lru_cache(None)
def primitive_certificate():
    from collections import Counter
    words=[source.tree_word(),source.copy_word(),circuits.pointer_prepare(),
           circuits.pointer_analyzer(F()),circuits.pointer_analyzer(F(1,2))]
    counts=Counter((e['kind'],e['angle'])for word in words for e in word)
    for block,n in Counter(i.tape for i,_ in circuits.field_binding(F(1,2),32)).items():
        for kind,angle,_ in block:counts[kind,angle]+=2*n
    return [{'kind':k,'angle_hex':a.hex(),'multiplicity':n,
             'complete_M6_error_upper':str(primitive_error(k,a))}
            for (k,a),n in sorted(counts.items())]


def admitted_budget():
    """Fixed local accuracy caps and an operation-derived full-state bound."""
    u=F(1,2**53)
    gamma=lambda n:n*u/(1-n*u)
    local_cap=F(1,2**46);input_cap=F(1,2**40);v_cap=F(1,2**46)
    b4=local_cap+gamma(32)*2*(1+local_cap)
    b8=local_cap+gamma(64)*sqrtup(F(8))*(1+local_cap)
    vf=sqrtup(F(24))+v_cap;ce=gamma(512)*vf
    re=gamma(8)*vf*(1+v_cap+ce)+(1+v_cap)*ce
    bq=v_cap*(2+v_cap)+re+gamma(4)*(1+(1+v_cap)**2+re)
    underflow=F(1,2**900)
    #26,624 two-bit field blocks, preparation/analyzer pointer blocks and10Q.
    total=26624*b4+2*b8+10*bq+underflow
    assert total<F(1,1000000)
    E=(input_cap+total)/(1-total)
    compressed=(1+v_cap)*E+v_cap+ce*(1+E)
    terminal=(2*(1+v_cap+gamma(8)*vf)*compressed)**2*(1+gamma(8*131072))+underflow
    failures=4*local_cap**2+underflow
    history=sqrtup(E*E+4*681283*failures+120*terminal)
    assert history<F(2,10**8)
    return {'local_matrix_error_cap':str(local_cap),'initial_state_error_cap':str(input_cap),
            'embedding_error_cap':str(v_cap),'complete_state_error_upper':str(E),'complete_orthogonal_history_distance_upper':str(history),
            'one_terminal_probability_upper':str(terminal),'one_export_probability_upper':str(failures),
            'normalization_error_upper':str(2*E+E*E),'underflow_allowance':str(underflow)}
