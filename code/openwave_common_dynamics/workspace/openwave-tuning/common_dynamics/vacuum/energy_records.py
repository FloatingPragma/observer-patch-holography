"""Outward electric-energy records of the complete dressed Slater state."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

import numpy as np

import reference as ref

HERE=Path(__file__).resolve().parent
SCALE=1<<160


def convolve(left,right):
    lo,hi=left;a,b=right
    lower=[sum(lo[j]*a[(k-j)%16] for j in range(16))//SCALE for k in range(16)]
    upper=[(sum(hi[j]*b[(k-j)%16] for j in range(16))+SCALE-1)//SCALE for k in range(16)]
    return lower,upper


def electric_interval(amplitudes):
    geometry=ref.dp.field.current_walk.geometry();C=geometry['C'];P=geometry['P']
    Z=sum(a*a for a in amplitudes)
    lo=[a*a*SCALE//Z for a in amplitudes]
    hi=[(a*a*SCALE+Z-1)//Z for a in amplitudes]
    cache={};total_lo=total_hi=0
    G=ref.reference_charge();offset=-P@G
    for edge,row in enumerate(C):
        key=(int(np.count_nonzero(row==1)),int(np.count_nonzero(row==-1)))
        if key not in cache:
            one=[SCALE]+[0]*15;distribution=(one,one.copy())
            for sign,count in ((1,key[0]),(-1,key[1])):
                kernel=([lo[(sign*j)%16] for j in range(16)],[hi[(sign*j)%16] for j in range(16)])
                for _ in range(count):distribution=convolve(distribution,kernel)
            cache[key]=distribution
        lower,upper=map(list,cache[key])
        # Each negative orbital occupies either endpoint of its disjoint
        # two-mode pair with probability1/2; all16 choices are independent
        # for diagonal occupation/field measurements.
        for a,b,sign in ref.pairs():
            x=int(P[edge,12*(a//4)]);y=int(P[edge,12*(b//4)])
            lower=[(lower[(j-x)%16]+lower[(j-y)%16])//2 for j in range(16)]
            upper=[(upper[(j-x)%16]+upper[(j-y)%16]+1)//2 for j in range(16)]
        squares=[((j+int(offset[edge])+8)%16-8)**2 for j in range(16)]
        total_lo+=sum(p*e for p,e in zip(lower,squares))
        total_hi+=sum(p*e for p,e in zip(upper,squares))
    return F(total_lo,2*SCALE),F(total_hi,2*SCALE)


def run():
    paths=[Path(__file__),Path(ref.__file__),HERE/'REFERENCE.json']
    pins={str(p.relative_to(ref.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    uniform=electric_interval([1]*16)
    assert uniform[0]<=3096<=uniform[1] and uniform[1]-uniform[0]<F(1,10**35)
    rows=[]
    for base in (4,16,64):
        lower,upper=electric_interval(ref.rational_gaussian_amplitudes(base))
        assert 0<lower<=upper<3096 and upper-lower<F(1,10**35)
        rows.append({'base':base,'electric_energy_lower_exact':str(lower),
                     'electric_energy_upper_exact':str(upper),'energy':float((lower+upper)/2),
                     'interval_width':float(upper-lower),'full65536occupation_law_included':True})
    assert pins=={str(p.relative_to(ref.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result={'schema':'oph.dirac.dressed_reference_electric_records.v1',
      'checks':{'all193sharedcycles_marginalized_exactly':True,'complete_Slater_occupation_distribution':True,
               'integer_outward_cyclic_convolutions':True,'uniform_field_energy_control':True,
               'source_bytes_unchanged':True},
      'fixed_point_bits':160,'states':rows,'uniform_energy_exact':'3096',
      'native_scope':'The base16 target is the one compiled by native_preparation.py. Base4/64 are ideal finite preparation comparisons and do not inherit that native certificate.',
      'measurement':'Read all288centeredZ16fields and compute sumF_e²/2. Bounded range0..9216 propagates complete recordTV error delta to an expectation error at most9216delta.',
      'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'ENERGY_RECORDS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'states':rows}))


if __name__=='__main__':run()
