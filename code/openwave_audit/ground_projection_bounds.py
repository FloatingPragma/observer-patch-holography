"""Exact arithmetic for the finite-energy projection theorem, no eigensolve."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import sys

sys.set_int_max_str_digits(20000)
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())


def main():
    manifests={
      'openwave-tuning/common_dynamics/autonomous_car_vacuum/MANIFEST.json':'df34a04cc324363c56f7cc110d0dc240a91fcd707ba3584fac461889c9529368',
      'openwave-tuning/common_dynamics/autonomous_car_response/G1_MANIFEST.json':'281723896dd3e09289b57b7f08543c9d310db261f6b9ffd3346d3bdcf1b99211',
      'openwave-tuning/common_dynamics/autonomous_external_screening/MANIFEST.json':'81fd403ce22e6472795cca74bbbf02b9dc09ee74baaaf9dc41e08cb81b48aae4',
      'openwave-tuning/common_dynamics/autonomous_car_native/MANIFEST.json':'1f510cbd194afa43f4bb37fc1cc45792d582aff67bee201b06504141a6f892d7',
    }
    pins={}
    for name,digest in manifests.items():
        p=ROOT/name;assert sha(p)==digest;pins[name]=digest
        m=read(p)
        for section in ('files','source_pins'):
            for key,value in m.get(section,{}).items():
                assert key not in pins or pins[key]==value,key
                assert sha(ROOT/key)==value,key;pins[key]=value
    fock=ROOT/'openwave-tuning/next/CHAIN3/native_field_connection/inputs/openwave-tuning/common_dynamics/quadratic_field/FOCK.json'
    spec=read(fock)['spectral_certificate']
    omega_squared_lower=96*(F(spec['minimum_eigenvalue'])-F(spec['actual_curvature_operator_error_upper']))
    assert omega_squared_lower>1
    epsilon=F(1,10**7);A=F(358);k=F(384)
    cutoff=A*(1+(k/epsilon)**2)
    delta=A/cutoff
    assert k*k*delta/(1-delta)==epsilon*epsilon
    B=cutoff.numerator//cutoff.denominator
    q2=6*cutoff;Q=math.isqrt(q2.numerator//q2.denominator)+1
    assert Q*Q>q2
    dimension=math.comb(32,16)*(B+1)**190*(2*Q+1)**3
    screening=F(1252815,135895533121)
    # Two further epsilon-sized allowances: finite-vector preparation and
    # the combined complete energy-reading interval. None is called executed.
    remaining=screening-3*epsilon
    assert remaining>0
    result={
      'schema':'oph.N16.ground_projection.v1','physical_source_unchanged':True,
      'uniform_probe_interval':['-4','4'],'normal_ordered_free_field_mean_upper_exact':str(A),
      'kinetic_norm_upper_exact':str(k),'physical_min_frequency_squared_lower_exact':str(omega_squared_lower),
      'requested_Ritz_energy_error_exact':str(epsilon),'field_energy_cutoff_exact':str(cutoff),
      'tail_probability_upper_exact':str(delta),'proved_Ritz_error_squared_exact':str(k*k*delta/(1-delta)),
      'basis_dimension_upper_decimal_digits':len(str(dimension)),
      'basis_dimension_upper_bit_length':dimension.bit_length(),
      'basis_dimension_upper_sha256':hashlib.sha256(str(dimension).encode()).hexdigest(),
      'per_mode_occupation_upper':B,'harmonic_coordinate_upper':Q,
      'screening_margin_exact':str(screening),'unallocated_margin_after_three_epsilons_exact':str(remaining),
      'compression_negative_control':{'physical_ground_energy':0,'exact_Ritz_energy':0,'wrong_exponentiated_compression_ground_energy':-1},
      'ground_gap_assumed':False,'unique_ground_assumed':False,
      'finite_matrix_generated':False,'native_ground_preparation_generated':False,
      'native_ground_preparation_executed':False,'full_criterion_met':False,
      'producer_self_approval':False,'source_pins':pins,
    }
    (HERE/'ground_projection_bounds.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('field_energy_cutoff_exact','basis_dimension_upper_decimal_digits','unallocated_margin_after_three_epsilons_exact','native_ground_preparation_generated')},indent=2))


if __name__=='__main__':main()
