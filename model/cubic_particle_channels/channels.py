"""Exact central filters and source-bound classification of every visible line."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import numpy as np
from . import group as g
from model.centered_particle_spectrum import source
from model.interacting_vacuum import ground

ROOT=Path(__file__).resolve().parents[2]

@lru_cache(maxsize=1)
def ground_invariance():
    with np.load(ROOT/'model/centered_particle_spectrum/witness.npz',allow_pickle=False) as w:
        columns=tuple(int(i) for i in w['ground_columns'])
    block=source.reduced('ground',columns)
    full=block['full_basis']
    for axes,shift,_,_ in source.automorphisms():
        target,sign=source.fock_action(4,4,axes,shift)
        image=np.empty_like(full);image[target]=sign[:,None]*full
        if not np.array_equal(image,full):raise ArithmeticError('ground block has nontrivial spatial character')
    for name in ('spin_plus','eta_minus'):
        if np.any(source.ladder(4,4,name)@full):raise ArithmeticError('ground block is not a singlet')
    return {'all_signed_automorphisms':48,'ground_block_dimension':block['dimension'],
      'all_ground_block_columns_fixed_exactly':True,'ground_hidden_spin_zero':True,
      'conclusion':'actual ground spatial irrep A1; ground localization in this block is inherited from FA85/88'}

def commutant_rank(matrices):
    # vec(BM-MB), row-major. A modular nonzero minor proves the rational
    # rank lower bound. Exhibited independent commuting matrices prove equality.
    eye=np.eye(8,dtype=np.int64)
    constraints=np.vstack([np.kron(eye,M.T)-np.kron(M,eye) for M in matrices])
    return source.modular_rank(constraints)

@lru_cache(maxsize=1)
def exact_filters():
    filters={};I=np.eye(16,dtype=np.int64)
    for name,d in zip(g.IRREPS,g.DIMENSIONS,strict=True):
        if name not in ('G1','G2','H'):continue
        real=np.full((16,16),g.Q(),dtype=object);imag=real.copy()
        for q in g.elements():
            S,W=g.hidden(q);a,b=g.spinor_parts(W);chi=g.characters(q)[name]
            # Sum chi(q) D(q); inversion gives the same central projector.
            for x,y in zip(*np.nonzero(S),strict=True):
                for i in range(2):
                    for j in range(2):
                        real[2*x+i,2*y+j]+=d*int(S[x,y])*chi*a[i][j]
                        imag[2*x+i,2*y+j]+=d*int(S[x,y])*chi*b[i][j]
        real=np.array([[v.integer() for v in row] for row in real],dtype=np.int64)
        imag=np.array([[v.integer() for v in row] for row in imag],dtype=np.int64)
        if not np.array_equal(real.T,real) or not np.array_equal(imag.T,-imag):raise ArithmeticError('non-Hermitian central filter')
        if not np.array_equal(real@real-imag@imag,48*real) or not np.array_equal(real@imag+imag@real,48*imag):raise ArithmeticError('central filter is not idempotent')
        filters[name]=(real,imag)
    if not np.array_equal(sum(v[0] for v in filters.values()),48*I) or np.any(sum(v[1] for v in filters.values())):raise ArithmeticError('incomplete charged filters')
    A=np.kron(np.asarray(ground.band_data()[0],dtype=np.int64),np.eye(2,dtype=np.int64))
    for name,(real,imag) in filters.items():
        if not np.array_equal(A@real,real@A) or not np.array_equal(A@imag,imag@A):raise ArithmeticError('filter mixes charged bands')
        if np.trace(A@real)!=0 or np.trace(A@imag)!=0:raise ArithmeticError('unequal band multiplicities')
        for other,(rr,ii) in filters.items():
            if other!=name and (np.any(real@rr-imag@ii) or np.any(real@ii+imag@rr)):raise ArithmeticError('nonorthogonal filters')
    return filters

@lru_cache(maxsize=1)
def structure():
    A=np.asarray(ground.band_data()[0],dtype=np.int64);I=np.eye(8,dtype=np.int64)
    if not np.array_equal(A@A,3*I):raise ArithmeticError('wrong signed hopping bands')
    signed=[];plain=[]
    for axes,shift,p,s in source.automorphisms():
        S=np.zeros((8,8),dtype=np.int64);S[p,np.arange(8)]=s;signed.append(S)
        P=np.zeros((8,8),dtype=np.int64);P[p,np.arange(8)]=1;plain.append(P)
    sr=commutant_rank(signed);pr=commutant_rank(plain)
    if sr!=62 or pr!=60:raise ArithmeticError('commutant dimensions changed')
    Ns=[g.walsh_projector_numerator(k) for k in range(4)]
    for k,N in enumerate(Ns):
        if not np.array_equal(N@N,8*N) or any(np.any(N@M-M@N) for M in plain):raise ArithmeticError('wrong neutral projector')
        if any(np.any(N@N2) for j,N2 in enumerate(Ns) if j!=k):raise ArithmeticError('neutral projectors not orthogonal')
    if not np.array_equal(sum(Ns),8*I):raise ArithmeticError('incomplete neutral filters')
    filters=exact_filters()
    return {'ground_invariance':ground_invariance(),
      'signed_commutant_rank_mod_1009':sr,'ordinary_commutant_rank_mod_1009':pr,
      'signed_commutant_dimension':2,'ordinary_commutant_dimension':4,
      'charged_filter_denominator':48,
      'charged_filters':{k:{'real_numerator':v[0].tolist(),'imaginary_numerator':v[1].tolist(),
            'rank':int(np.trace(v[0]))//48} for k,v in filters.items()},
      'neutral_filter_denominator':8,'neutral_projectors':[N.tolist() for N in Ns]}

def interval(pair):
    if type(pair) is not list or len(pair)!=2 or any(type(v) not in (float,int,str) for v in pair):raise ValueError('invalid retained interval')
    a,b=map(lambda x:F(str(x)),pair)
    if a>b:raise ValueError('reversed retained interval')
    return a,b

def definite_sign(pair):
    a,b=interval(pair)
    if a>0:return 1
    if b<0:return -1
    raise ValueError('interval does not isolate a sign')

def scaled(pair,value):
    a,b=interval(pair);c,d=a*value,b*value
    return min(c,d),max(c,d)

def overlap(a,b):return max(a[0],b[0])<=min(a[1],b[1])

def classify(charged,neutral):
    """Intervals select exact symmetry branches; they are not a new eigensolve."""
    representations=g.certificate();rows=[];bands=Counter()
    for index,row in enumerate(charged['spectrum']['rows']):
        if row['visible_spatial_multiplicity']!=4 or row['strictly_positive_weight_proved'] is not True:raise ValueError('charged rank/positivity premise changed')
        if interval(row['local_branch_weight_interval'])[0]<=0:raise ValueError('charged weight not strictly positive')
        sign=definite_sign(row['normalized_spatial_overlap_intervals'][1]);bands[str(sign)]+=1
        rows.append({'source_index':index,'energy_interval':row['excitation_energy_interval'],
          'local_branch_weight_interval':row['local_branch_weight_interval'],'hopping_band_sign':sign,
          'visible_dimension_each_charge':8,'binary_cubic_irreps':{'G1':1,'G2':1,'H':1},
          'complete_probe_trace_fractions':{'G1':'1/4','G2':'1/4','H':'1/2'}})
    if len(rows)!=228 or dict(bands)!={'-1':114,'1':114}:raise ValueError('charged source coverage changed')
    result={'charged':rows,'charged_band_counts':dict(bands),'neutral':[]}
    for response in neutral['spectrum']['responses']:
        output=[];counts=Counter();dark=0
        # The retained response order is part of the input contract.
        spin=len(result['neutral'])==0
        for row in response['rows']:
            if row['exactly_dark']:
                if row['strictly_positive'] or interval(row['weight_interval'])!=(0,0):raise ValueError('invalid dark branch')
                dark+=1;continue
            if row['strictly_positive'] is not True or interval(row['weight_interval'])[0]<=0:raise ValueError('neutral weight not positive')
            rank=row['visible_spatial_rank'];B=row['spatial_residue_intervals']
            if rank==1:degree=3
            elif rank==3:degree=1 if definite_sign(B[0][1])>0 else 2
            else:raise ValueError('unexpected neutral spatial rank')
            N=g.walsh_projector_numerator(degree)
            # The exact commutant/rank proof fixes B=(w/rank)N; check all
            # retained intervals against it, including off-diagonal zeroes.
            for x in range(8):
                for y in range(8):
                    if not overlap(interval(B[x][y]),scaled(row['weight_interval'],F(int(N[x,y]),rank))):raise ValueError('neutral residue contradicts spatial branch')
            reps=representations['spin_by_walsh_degree' if spin else 'density_by_walsh_degree'][str(degree)]
            dimension=rank*(3 if spin else 1)
            fractions={name:str(F(g.DIMENSIONS[g.IRREPS.index(name)]*mult,dimension)) for name,mult in reps.items()}
            output.append({'source_index':row['index'],'energy_interval':row['energy_interval'],
              'local_component_weight_interval':row['weight_interval'],'walsh_degree':degree,
              'visible_dimension_all_components':dimension,'cubic_irreps':reps,'complete_probe_trace_fractions':fractions})
            counts[str(degree)]+=1
        if len(output)!=135 or dark!=1 or dict(counts)!={'3':24,'2':54,'1':57}:raise ValueError('neutral source coverage changed')
        result['neutral'].append({'channel':'hidden_spin_vector' if spin else 'density_scalar','rows':output,'degree_counts':dict(counts),'exact_dark_candidates':dark})
    if len(result['neutral'])!=2:raise ValueError('wrong neutral response channels')
    return result


def spectrum():
    read=lambda p:json.loads((ROOT/p).read_text())
    return classify(read('model/centered_particle_spectrum/results.json'),read('model/centered_neutral_spectrum/results.json'))
