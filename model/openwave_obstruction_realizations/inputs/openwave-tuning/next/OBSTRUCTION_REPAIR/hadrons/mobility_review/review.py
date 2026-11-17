"""Independent raw24-CAR parity and one-hop denominator controls."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import json
HERE=Path(__file__).resolve().parent

def op(w,j,create):
    present=(w>>j)&1
    if present==create:return None,0
    return w^(1<<j),(-1)**((w&((1<<j)-1)).bit_count())

def pair(v,a,b):
    out=defaultdict(int)
    for w,x in v.items():
        z,s=op(w,b,False)
        if z is None:continue
        y,t=op(z,a,True)
        if y is not None:out[y]+=x*s*t
    return {w:x for w,x in out.items() if x}

def total(v,terms):
    out=defaultdict(int)
    for a,b in terms:
        for w,x in pair(v,a,b).items():out[w]+=x
    return {w:x for w,x in out.items() if x}

def hop(v):return total(v,[(a,a+12) for a in range(12)]+[(a+12,a) for a in range(12)])

def cas6(v):
    out=defaultdict(int)
    for a in range(3):
        for b in range(3):
            one=total(v,[(3*x+b,3*x+a) for x in range(4)])
            two=total(one,[(3*x+a,3*x+b) for x in range(4)])
            for w,z in two.items():out[w]+=3*z
    for w,z in v.items():out[w]-=(w&4095).bit_count()**2*z
    return {w:x for w,x in out.items() if x}

def triple_amp(w,modes,create):
    amp=1
    for j in (list(reversed(modes)) if create else modes):
        w,s=op(w,j,create)
        if w is None:return 0,None
        amp*=s
    return amp,w

def main():
    rows=[]
    for spectator in range(4):
        left=(1<<(3*spectator))-1;triple=list(range(3*spectator,3*spectator+3));right=sum(1<<j for j in triple)
        initial=left|(right<<12);final=left|right
        one=hop({initial:1});assert cas6(one)=={w:8*x for w,x in one.items()}
        value=hop(hop(one)).get(final,0)
        bl,fl=triple_amp(left,triple,True);br,fr=triple_amp(right,triple,False);assert fl==final and fr==0
        withparity=6*bl*((-1)**left.bit_count())*br;without=6*bl*br
        assert value==withparity
        rows.append(dict(left_quark_number=left.bit_count(),right_quark_number=3,raw24CAR_Kcubed_matrix_element=value,
            parity_formula=withparity,omitted_parity_formula=without,parity_omission_changes_sign=without!=value,
            one_hop_Casimir6_eigenvalue=8))
    coeff=6*F(1,8)**3/F(4,3)**2;assert coeff==F(27,4096)
    out=dict(schema='oph.hadrons.mobility_independent_review.v1',producer_code_imported=False,
        complete8464_matrix_reconstruction_claim=False,raw24CAR_controls=rows,
        leading_baryon_transfer_amplitude='27/(4096*kappa^2)',
        outside_resolvent_argument='P K Q H0^-1=Delta^-1 P K Q and H0^-1 Q K P=Delta^-1 Q K P; no assumption on all two-hop intermediates',
        normalization_terms='At order3, two hopping factors plus PVP mass have net local transfer0, not3; mixed auxiliary factors cannot replace three hops',
        physical_halo_or_full_Heff3_claim=False,all_checks_passed=True)
    (HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
