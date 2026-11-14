"""Exact finite-time cyclic-to-integer field cutoff certificate.

Run with Python's standard library. The hypotheses and operator proof are
in FIELD_REGISTER_PROOF.md. This module does not construct a Maxwell state.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent


def exponential_tail_upper(magnitude, first_order):
    """A rational upper bound on sum_{k>=first_order} magnitude**k/k!."""
    magnitude = F(magnitude)
    if magnitude < 0 or type(first_order) is not int or first_order < 0:
        raise ValueError('nonnegative magnitude and integral order required')
    start = max(first_order, math.floor(magnitude)+1)
    prefix = sum((magnitude**k/math.factorial(k)
                  for k in range(first_order,start)),F())
    return prefix + (magnitude**start/math.factorial(start))/(1-magnitude/(start+1))


def state_error_bound(modulus, initial_cap, matter_motion, magnetic_coefficients):
    """Return a norm bound for a complete finite observation interval.

    Supply every face coefficient at every update, including repeats.
    matter_motion bounds total absolute per-edge transport at all prefixes.
    Coefficients are exact rationals or exact rationalized source floats.
    """
    if type(modulus) is not int or modulus < 2 or modulus%2:
        raise ValueError('an even modulus at least two is required')
    if any(type(v) is not int or v<0 for v in (initial_cap,matter_motion)):
        raise ValueError('nonnegative integer support bounds required')
    magnitude = sum((abs(F(beta)) for beta in magnetic_coefficients),F())
    distance = modulus//2-initial_cap-matter_motion
    if distance <= 0:
        return F(2)
    return min(F(2),2*exponential_tail_upper(magnitude,distance))


def run():
    # Positive exact-series comparisons verify the tail certificate without
    # using the same geometric-remainder implementation as the checker.
    cases=0
    for magnitude in (F(0),F(1,10),F(123,100),F(15,2)):
        for order in (0,1,5,13):
            upper=exponential_tail_upper(magnitude,order)
            partial=sum((magnitude**k/math.factorial(k)
                         for k in range(order,order+160)),F())
            assert upper>=partial
            cases+=1
    beta=F.from_float(.41)
    one_face=[beta]*3
    homogeneous=[beta]*(3*208)
    rows=[]
    for modulus in (16,32,64,128):
        bound=state_error_bound(modulus,0,3,one_face)
        rows.append({'modulus':modulus,'distance_to_first_possible_wrap':modulus//2-3,
                     'state_norm_bound':float(bound),'terminal_total_variation_bound':float(min(F(1),bound)),
                     'bound_exact':str(bound)})
    assert state_error_bound(32,0,3,one_face)<F(52,10**10)
    assert state_error_bound(32,0,3,homogeneous)==2
    assert state_error_bound(16,0,8,[])==2
    assert state_error_bound(16,0,3,[])==0
    assert state_error_bound(32,0,3,[-beta]*3)==state_error_bound(32,0,3,one_face)
    rejected=0
    for args in ((15,0,0,[]),(16,-1,0,[]),(16,0,True,[]),(True,0,0,[])):
        try: state_error_bound(*args)
        except ValueError: rejected+=1
        else: raise AssertionError('invalid field contract admitted')
    result={'schema':'oph.field_register_transfer.certificate.v1',
            'checks':{'exact_positive_series_controls':cases,
                      'invalid_domains_rejected':rejected,
                      'all_faces_and_all_times_required':True,
                      'matter_motion_and_support_retained':True,
                      'zero_magnetic_coefficient_control':True,
                      'coefficient_sign_uses_absolute_magnitude':True},
            'one_face_three_step_controls':rows,
            'homogeneous208_face_modulus32_bound':2,
            'coefficient_ancestry':'beta is the exact rational represented by source binary64 .41; three updates, oneface or208faces at each update',
            'theorem_scope':'Cyclic field to integer compact-rotor evolution on an admitted initial support, arbitrary quantum references, common complete terminal instrument. Native compiler/reset/readout errors are additional.',
            'not_claimed':['Maxwell continuum','prepared interacting vacuum','physical unit calibration','native circuit resource equality between moduli','official platform validation'],
            'source_sha256':{Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}}
    (HERE/'field_register_transfer_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'Z32_oneface_bound':rows[1]['state_norm_bound']}))


if __name__=='__main__':
    run()
