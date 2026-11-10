"""Complete Gauss and translation-resolved spectra of the periodic source."""
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
import numpy as np
import sympy as s
from .basis import diagonal_basis,physical_basis,coordinate_rows,block_labels,character_multiplicities
from ..two_cell.symmetries import periodic_column,periodic_translation
from ..two_cell.basis import boundary_charges
from ..two_cell.spectrum import exact_eigenbasis_error,exact_dyadic_eigenbasis_error,components
from ..eigenbasis_witness import (SCALE,block_key,load_witness,
                                 numerical_proposal,write_witnesses)

WITNESS_PATH=Path(__file__).with_name('eigenbasis_witnesses.npz')
WITNESS_DOMAIN='periodic_spectrum'


def _kernel(matrix):
    """Exact rational kernel by disjoint row-connected coordinate components."""
    n=matrix.cols;parent=list(range(n));row_columns=defaultdict(list)
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for row,col in matrix.todok():row_columns[row].append(col)
    for cols in row_columns.values():
        for col in cols[1:]:parent[find(col)]=find(cols[0])
    blocks=defaultdict(list);blockrows=defaultdict(list)
    for col in range(n):blocks[find(col)].append(col)
    for row,cols in row_columns.items():blockrows[find(cols[0])].append(row)
    entries={};coordinates=[]
    for cols in sorted(blocks.values(),key=lambda x:x[0]):
        small=s.Matrix(matrix.extract(blockrows[find(cols[0])],cols));rref,pivots=small.rref()
        for free in (j for j in range(len(cols)) if j not in pivots):
            output=len(coordinates);coordinates.append(cols[free]);entries[cols[free],output]=1
            for row,pivot in enumerate(pivots):
                if rref[row,free]:entries[cols[pivot],output]=-rref[row,free]
    result=s.ImmutableSparseMatrix(n,len(coordinates),entries)
    assert matrix*result==s.zeros(matrix.rows,result.cols)
    assert result.extract(coordinates,range(result.cols))==s.eye(result.cols)
    return result,tuple(coordinates)


@lru_cache(maxsize=None)
def physical_block(boundary,L,twoSz):
    keys=diagonal_basis(boundary,L,twoSz);b,g=physical_basis(boundary,L,twoSz)
    positions={key:i for i,key in enumerate(keys)};phases=[s.I**((key[0]&63).bit_count()) for key in keys]
    entries={}
    for j,key in enumerate(keys):
        for target,value in periodic_column(key).items():
            assert target in positions,('omitted full-H output',key,target)
            i=positions[target];value=s.expand(s.conjugate(phases[i])*value*phases[j])
            assert value.is_Rational
            entries[i,j]=value
    h=s.ImmutableSparseMatrix(len(keys),len(keys),entries);assert h==h.T
    hb=h*b;generator=s.ImmutableSparseMatrix(hb.extract(coordinate_rows(boundary,L,twoSz),range(b.cols)))
    assert b*generator==hb
    k=s.ImmutableSparseMatrix(g*generator);assert k==k.T and generator.T*g==g*generator
    embedding=s.ImmutableSparseMatrix(b.rows,b.cols,{(i,j):phases[i]*v for (i,j),v in b.todok().items()})
    return keys,embedding,g,k,generator


@lru_cache(maxsize=None)
def neutral_translation(L,twoSz):
    keys,embedding,g,k,h=physical_block(0,L,twoSz)
    b,_=physical_basis(0,L,twoSz);positions={key:i for i,key in enumerate(keys)}
    phases=[s.I**((key[0]&63).bit_count()) for key in keys];xi=s.I**(L%2)
    entries={}
    for j,key in enumerate(keys):
        for target,value in periodic_translation(key).items():
            i=positions[target];value=s.expand(s.conjugate(phases[i])*value*phases[j]/xi)
            assert value.is_Rational;entries[i,j]=value
    t=s.ImmutableSparseMatrix(len(keys),len(keys),entries)
    tb=t*b;c=s.ImmutableSparseMatrix(tb.extract(coordinate_rows(0,L,twoSz),range(b.cols)))
    assert b*c==tb
    assert c*c==s.eye(c.rows) and c.T*g==g*c and h*c==c*h
    return c,xi


@lru_cache(maxsize=None)
def resolved_block(boundary,L,twoSz,sign=0):
    keys,embedding,g,k,h=physical_block(boundary,L,twoSz)
    if sign==0:return keys,embedding,g,k,h
    if boundary_charges(boundary)!=(0,0) or sign not in (-1,1):
        raise ValueError('bare translation resolves the neutral source only')
    t,xi=neutral_translation(L,twoSz)
    f,coordinates=_kernel(t-sign*s.eye(t.rows))
    assert f.cols==(t.rows+sign*s.trace(t))/2
    target=embedding*f;gram=s.ImmutableSparseMatrix(f.T*g*f)
    hf=h*f;generator=s.ImmutableSparseMatrix(hf.extract(coordinates,range(f.cols)))
    assert f*generator==hf
    energy=s.ImmutableSparseMatrix(gram*generator);assert energy==energy.T
    return keys,s.ImmutableSparseMatrix(target),gram,energy,generator


@lru_cache(maxsize=None)
def rational_block(boundary,L,twoSz,sign=0):
    _,_,g,k,_=resolved_block(boundary,L,twoSz,sign)
    if not g.rows:return np.array([],dtype=object),np.zeros((0,0),dtype=object)
    return load_witness(WITNESS_PATH,WITNESS_DOMAIN,
                        block_key(boundary_charges(boundary),L,twoSz,sign),k,g)


@lru_cache(maxsize=None)
def numerical_block(boundary,L,twoSz,sign=0):
    """Float64 view of the stored dyadic proposal, without an eigensolve."""
    d,z=rational_block(boundary,L,twoSz,sign)
    return np.asarray(d,dtype=float)/SCALE,np.asarray(z,dtype=float)/SCALE


@lru_cache(maxsize=None)
def certified_block(boundary,L,twoSz,sign=0):
    keys,b,g,k,h=resolved_block(boundary,L,twoSz,sign)
    xi=s.I**(L%2)
    result={'boundary':boundary_charges(boundary),'L':L,'twoSz':twoSz,
            'translation_involution_sign':sign,'translation_eigenvalue':str(sign*xi) if sign else 'not resolved',
            'physical_dimension':g.rows,'diagonal_dimension':len(keys),'energy_intervals':[]}
    if not g.rows:return result
    # Connected generalized-energy components reduce the exact complete-
    # eigenbasis certification cost; their direct sum retains every state.
    records=[];upper=s.Integer(0)
    for component in components(k,g):
        kc=k.extract(component,component);gc=g.extract(component,component)
        d,z=load_witness(WITNESS_PATH,WITNESS_DOMAIN,
                        block_key(boundary_charges(boundary),L,twoSz,sign,component),kc,gc)
        error,gram_error=exact_dyadic_eigenbasis_error(kc,gc,d,z)
        upper=max(upper,error)
        fine=1<<44;scale=1<<32;radius=s.Rational(1,1<<26)
        for value in d:
            c=int(round(s.Rational(int(value)*scale,fine)))
            assert abs(s.Rational(int(value),fine)-s.Rational(c,scale))+error<=radius, 'exact eigenbasis error exceeds energy enclosure'
            records.append((int(c),str(s.Rational(int(c),scale)-radius),str(s.Rational(int(c),scale)+radius)))
    records.sort(key=lambda x:x[0]);assert len(records)==g.rows
    # Whole-block numerical vectors used by detector consumers have their
    # separate residual/projector certificates; this receipt certifies all
    # energies via the exact component direct sum.
    result.update({'energy_center_denominator':1<<32,'energy_centers_numerator':[r[0] for r in records],
        'energy_interval_radius':str(s.Rational(1,1<<26)),
        'energy_intervals':[[r[1],r[2]] for r in records],
        'largest_certified_component':max(map(len,components(k,g))),
        'exact_complete_Gauss_and_H_invariance':True})
    return result


def generate_witnesses():
    """Propose whole blocks and their exact connected energy components."""
    def blocks():
        for boundary in (0,-1,1,(-1,0),(1,0)):
            signs=(0,-1,1) if boundary_charges(boundary)==(0,0) else (0,)
            for L,twoSz in block_labels(boundary):
                for sign in signs:
                    _,_,g,k,_=resolved_block(boundary,L,twoSz,sign)
                    if not g.rows:continue
                    d,z=numerical_proposal(k,g)
                    # Full-block consumers execute their individual-column
                    # residual and full-Gram projector checks on replay.
                    error,gram=exact_dyadic_eigenbasis_error(k,g,d,z)
                    assert gram<s.Rational(1,1<<23) and error<s.Rational(1,10**7)
                    yield block_key(boundary_charges(boundary),L,twoSz,sign),k,g,d,z
                    for component in components(k,g):
                        kc=k.extract(component,component);gc=g.extract(component,component)
                        dc,zc=numerical_proposal(kc,gc)
                        error,gram=exact_dyadic_eigenbasis_error(kc,gc,dc,zc)
                        radius=s.Rational(1,1<<26)
                        assert all(abs(s.Rational(int(value),SCALE)-
                                       s.Rational(round(s.Rational(int(value),1<<12)),1<<32))+error<=radius for value in dc)
                        yield block_key(boundary_charges(boundary),L,twoSz,sign,component),kc,gc,dc,zc
    count=write_witnesses(WITNESS_PATH,WITNESS_DOMAIN,blocks())
    rational_block.cache_clear();numerical_block.cache_clear();certified_block.cache_clear()
    return count


def sector_certificate(boundary=0):
    signs=(-1,1) if boundary_charges(boundary)==(0,0) else (0,)
    blocks=[certified_block(boundary,L,spin,sign) for L,spin in block_labels(boundary) for sign in signs]
    records=[(s.Rational(pair[0]),s.Rational(pair[1]),row['L'],row['twoSz'],row['translation_involution_sign'],j)
             for row in blocks for j,pair in enumerate(row['energy_intervals'])]
    lo=min(r[0] for r in records);hi=min(r[1] for r in records)
    ground=[r for r in records if r[0]<=hi];other=[r[0] for r in records if r[0]>hi]
    assert sum(row['physical_dimension'] for row in blocks)==sum(character_multiplicities(boundary).values())
    return {'boundary':boundary_charges(boundary),'physical_dimension':sum(row['physical_dimension'] for row in blocks),
        'blocks':blocks,'global_ground_interval':[str(lo),str(hi)],
        'ground_window':[{'L':r[2],'twoSz':r[3],'translation_sign':r[4],'eigen_index':r[5]} for r in ground],
        'ground_window_rank':len(ground),'ground_certifiably_simple':len(ground)==1,
        'ground_window_to_rest_gap_lower_bound':str(min(other)-hi)}
