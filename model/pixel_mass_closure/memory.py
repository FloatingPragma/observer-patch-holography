"""Certify the complete native-memory inverse throughout each root enclosure."""
from fractions import Fraction as F

import numpy as np

from model.native_flavor_orbit.intervals import I, sqrt_rational, rounded_out
from model.mcpr_native_memory_bridge import construction as old


def sign(x):
    if x.lo>0:return 1
    if x.hi<0:return -1
    return 0


def bracket(function,lo,hi,tolerance):
    lo,hi,tolerance=F(lo),F(hi),F(tolerance)
    left,right=sign(function(I(lo))),sign(function(I(hi)))
    if left*right != -1:
        raise ArithmeticError('uniform endpoint signs required')
    while hi-lo>tolerance:
        middle=(lo+hi)/2
        sm=sign(function(I(middle)))
        if not sm:
            # Return the proved enclosing bracket, not an uncertain point.
            break
        if sm==left:lo=middle
        else:hi=middle
    return I(lo,hi)


def certify(root):
    masses=[I(row['lo'],row['hi']) for row in root['implicit_and_mass_enclosures']['masses_over_v']]
    total=sum(masses,I(0))
    if not 0<total.lo<total.hi<1:
        raise ArithmeticError('positive subunit total mass required')
    lam=[x/total for x in masses]
    e2=sum((lam[i]*lam[j] for i in range(3) for j in range(i)),I(0))
    e3=lam[0]*lam[1]*lam[2]
    determinant=I(F(2,5))-2*sqrt_rational(5)/25
    poly=lambda z:z**3-z*z+F(5,4)*e2*z-e3/determinant
    w=bracket(poly,'.0009','.0011','1e-10')
    u=bracket(poly,'.069','.072','1e-10')
    v=bracket(poly,'.92','.94','1e-10')
    A=v*(v+w); B=v*(1-v+w)+(v+w)*w
    quad=lambda p:-A*p*p+B*p-w
    inverses=[]
    for lo,hi in (('.018','.021'),('.058','.062')):
        p=bracket(quad,lo,hi,'1e-7')
        q=1-v/(1-p*(v+w))
        gain=total/p
        if not 0<p.lo<p.hi<F(1,12) or not 0<q.lo<q.hi<F(1,6):
            raise ArithmeticError('native p,q domain not certified')
        inverses.append({'p':rounded_out(p,12).strings(),'q':rounded_out(q,12).strings(),
                         'gain_over_v':rounded_out(gain,12).strings()})
    return {'root_domain':'entire certified implicit-coordinate enclosure',
            'normalized_mass_intervals':[rounded_out(x,18).strings() for x in lam],
            'projector_weight_intervals':[rounded_out(x,14).strings() for x in (u,v,w)],
            'inverse_branches':inverses,'cubic_root_count':3,'quadratic_root_count':2,
            'total_mass_over_v':rounded_out(total,18).strings(),
            'method':'rational endpoint signs, degree exhaustion, and exact generic native-memory identities',
            'all_complete_record_probabilities_retained':True}


def point_controls(root):
    """Check the literal 512-record transition for a center in each root box."""
    masses=np.asarray([(float(b['lo'])+float(b['hi']))/2 for b in root['implicit_and_mass_enclosures']['masses_over_v']])
    total=masses.sum();lam=masses/total
    det=2/5-2*np.sqrt(5)/25
    weights=np.sort(np.roots([1,-1,1.25*sum(lam[i]*lam[j] for i in range(3) for j in range(i)),-np.prod(lam)/det]))
    w,u,v=weights
    A=v*(v+w);B=v*(1-v+w)+(v+w)*w
    pvals=np.sort(np.roots([A,-B,w]))
    controls=[]
    for p in pvals:
        q=1-v/(1-p*(v+w));gain=total/p
        law=old.direct_memory({'p':p,'q':q})
        actual=gain*law['coefficient_matrix']
        error=float(np.max(np.abs(np.linalg.eigvalsh(actual)-masses)))
        if not error<2e-12:
            raise ArithmeticError('complete-record coefficient does not realize the enclosed source law')
        controls.append({'p':float(p),'q':float(q),'gain_over_v':float(gain),
                         'retained_records_per_input':law['retained_records_per_input'],
                         'retained_old_inputs':law['retained_old_inputs'],
                         'spectral_error_below_2e_minus_12':True})
    return controls
