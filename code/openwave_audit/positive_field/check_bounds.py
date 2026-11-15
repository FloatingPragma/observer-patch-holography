"""Exact arithmetic for the positive-field breakup theorem and its margin."""
from fractions import Fraction as F
from pathlib import Path
import json

def compute():
    a0=F(57,16); separated=F(281,48); weak=F(91,16)
    b2=F(3,32); h=F(1,2); variance=F(1,1000)
    delta=3*variance/(4*h)
    assert delta==F(3,2000)
    d0=separated-a0; d1=weak-a0
    assert d0==F(55,24) and d1==F(17,8)
    zero0=4*b2/(d0*d0+4*b2); zero1=4*b2/(d1*d1+4*b2)
    assert zero0==F(216,3241) and zero1==F(24,313)
    q0=F(7,100); q1=F(2,25)
    # Squared, positive quantities certify the radical comparisons.
    assert 4*b2*(1-q0)/q0==F(279,56)<F(9,4)**2
    assert 4*b2*(1-q1)/q1==F(69,16)<F(133,64)**2
    assert d0-F(9,4)==F(1,24)
    assert d1-F(133,64)==F(3,64)
    assert q0/F(24)>delta and 3*q1/64>delta
    assert q1+F(1,1000)<F(1,10)
    # A common-source clock rescaling must also rescale harmonic energy.
    assert 3*F(1,100)/(4*F(1,8))==F(3,50)>F(1,40)
    # The declared squeezed state saturates uncertainty with finite variances.
    assert variance*250==F(1,4)
    # The degree computation is the exact cross-block return count.
    counts=[]
    for L in range(2,6):
        for x in range(L):
            degree=int(x>0)+int(x<L-1)
            gram=F(3*degree,64)
            assert gram<=b2
            counts.append({'L':L,'site':x,'degree':degree,'cross_gram':str(gram)})
    return {'schema':'oph.color_positive_field.exact_bounds.v1',
      'a0':str(a0),'separated_floor':str(separated),'weak_separated_floor':str(weak),
      'cross_norm_squared_upper':str(b2),'h_star':str(h),'initial_field_excess':str(delta),
      'zero_excess_unforced_bound':str(zero0),'zero_excess_weak_bound':str(zero1),
      'finite_energy_unforced_bound_strict':str(q0),'finite_energy_weak_bound_strict':str(q1),
      'unforced_energy_contradiction_margin':str(q0/24-delta),
      'weak_energy_contradiction_margin':str(3*q1/64-delta),
      'native_bad_probability_bound':str(q1+F(1,1000)),
      'cross_gram_by_site':counts,'actual_gaussian_native_program_verified':False,
      'full_row_award':False}

if __name__=='__main__':
    out=compute();Path(__file__).with_name('BOUNDS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='cross_gram_by_site'},sort_keys=True))
