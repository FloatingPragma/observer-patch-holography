"""A fixed dyadic geometry prepared by the declared native-charge extension."""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import numpy as np
import sympy as s
from model.openwave_native_force.force import green_solve, finite_energy_certificate, laplacian
from model.finite_area_response.response import Interval, DENOMINATOR
from math import isqrt


def root(value):
    value = F(value)
    if value < 0: raise ValueError('nonnegative square root required')
    scale = DENOMINATOR
    a = isqrt(value.numerator*scale*scale//value.denominator)
    return Interval(F(a, scale), F(a+1, scale))


@lru_cache(maxsize=1)
def prepared():
    n = 31; denominator = 2**48
    source = np.zeros((n,n,n)); source[15,15,15] = float(F(9,140000))
    integers = np.rint(green_solve(source)*denominator).astype(np.int64)
    assert integers.max() < 0
    radii = tuple(1+F(int(integers[15,j,15]), denominator) for j in range(n))
    assert all(F(1)-F(1,50000) < r < 1 for r in radii)
    cert = finite_energy_certificate(n, ((15,15,15),))
    return integers, radii, cert


def metric_certificate():
    integers, radii, energy = prepared()
    near = radii[15]; far = radii[23]
    travel = Interval.point(0)
    for r in radii:
        travel += Interval.point(1)/(Interval.point(r)*root(r))
    delay = travel-len(radii)
    assert delay.lo > 0
    matrices = []
    for r in (near, far):
        sr = s.Rational(r)
        density = s.diag(-1/sr, 1/sr, sr**2, sr**2)
        metric = s.diag(-sr**2, sr**2, 1/sr, 1/sr)
        assert s.simplify(s.sqrt(-density.det())*density.inv()-metric) == s.zeros(4)
        matrices.append({'r': str(r), 'C': [str(x) for x in density.diagonal()],
                         'g': [str(x) for x in metric.diagonal()]})
    return {'geometry_preparation': 'dyadic trial of the same 31-cube native-charge elastic energy',
            'native_geometry_charge': '9/140000', 'dyadic_bits': 48,
            'integer_geometry_sha256': hashlib.sha256(integers.astype('<i8').tobytes()).hexdigest(),
            'native_force_energy_certificate': energy,
            'full_nonlinear_geometry_minimizer_asserted': False,
            'transverse_slice_r': list(map(str, radii)), 'local_metric_rows': matrices,
            'radial_null_coordinate_speed_squared': '1',
            'transverse_null_coordinate_speed_squared': 'r^3',
            'lapse_ratio_near_to_far': str(near/far),
            'piecewise_metric_transverse_ray_time_interval': travel.encoded(),
            'piecewise_metric_transverse_delay_interval': delay.encoded(),
            'piecewise_metric_delay_decimal': [float(delay.lo), float(delay.hi)],
            'ray_is_exact_finite_lattice_wavefront_claimed': False,
            'geometry_is_static_during_probe': True}


def einstein_test():
    # Independent linearized Ricci construction, without imposing Einstein balance.
    t,x,y,z = s.symbols('t x y z', real=True); coordinates = (t,x,y,z)
    q = s.Function('q')(x,y,z); eta = s.diag(-1,1,1,1)
    h = s.diag(-2*q,2*q,-q,-q); trace = s.trace(eta*h)
    ricci = s.zeros(4)
    for a in range(4):
        for b in range(4):
            ricci[a,b] = s.simplify(sum(eta[c,c]*(
                s.diff(h[c,b],coordinates[c],coordinates[a])
                +s.diff(h[c,a],coordinates[c],coordinates[b])
                -s.diff(h[a,b],coordinates[c],coordinates[c])) for c in range(4))/2
                -s.diff(trace,coordinates[a],coordinates[b])/2)
    scalar = s.trace(eta*ricci); tensor = s.simplify(ricci-eta*scalar/2)
    assert s.simplify(tensor[0,0]-s.diff(q,x,2)+(s.diff(q,y,2)+s.diff(q,z,2))/2) == 0
    # Even a smooth harmonic inverse-radius tail has a nonzero vacuum residual.
    radius = s.sqrt(x*x+y*y+z*z); harmonic = -1/radius
    evaluated = tensor[0,0].subs(q,harmonic).doit().subs({x:3,y:0,z:0})
    assert evaluated == -s.Rational(1,9)
    integers, _, _ = prepared(); denominator = 2**48
    point = (21,15,15)
    derivatives = []
    for axis in range(3):
        left=list(point); right=list(point); left[axis]-=1; right[axis]+=1
        derivatives.append(F(int(integers[tuple(left)])-2*int(integers[point])+int(integers[tuple(right)]), denominator))
    g00 = derivatives[0]-(derivatives[1]+derivatives[2])/2
    residual = sum(derivatives)
    assert abs(residual) <= F(6, denominator) and abs(g00) > F(1,10**8)
    return {'linearized_Einstein_00': str(tensor[0,0]),
            'smooth_harmonic_q_minus_inverse_radius_G00_at_3_0_0': str(evaluated),
            'finite_source_sample': list(point), 'discrete_laplacian_q': str(residual),
            'discrete_linearized_G00': str(g00),
            'outside_pinned_source_matter_density': '0',
            'stress_comparison_domain': 'linear weak field against the pinned pressureless source; no complete covariant stress dictionary',
            'same_source_vacuum_Einstein_premise_passed': False,
            'stress_defined_by_Einstein_tensor_to_force_balance': False,
            'regional_Einstein_decoder_used_without_its_premises': False,
            'native_universal_scalar_mass_coupling_obstruction_preserved': True}
