"""Exact cube residues of the actual native-spin Fourier contraction."""
from functools import lru_cache
import sympy as s
from ..atomic_spin_geometry.cube import magnetization_block, full_operators, pauli
from ..atomic_spin_geometry.momentum import translation
from ..family_cell.quotient import evaluate, divide, outward


@lru_cache(maxsize=1)
def exact_data():
    basis, h = magnetization_block(4)
    x = s.Symbol('x'); polynomial = s.Poly(x**3+20*x**2-16*x-576, x)
    w = s.Matrix.hstack(*(h**3+20*h**2-16*h-576*s.eye(70)).nullspace())
    gram = w.T*w; reduced = gram.inv()*w.T*h*w
    assert h*w == w*reduced and reduced.charpoly(x).as_expr() == polynomial.as_expr()
    triplet = (h+16*s.eye(70)).nullspace(); assert len(triplet) == 1
    t = triplet[0]; tnorm = (t.T*t)[0]
    neel = s.diag(*(sum((-1)**site.bit_count()*(1-2*((mask >> site)&1))
        for site in range(8))/s.Integer(2) for mask in basis))
    u = gram.inv()*w.T*neel*t
    response = u*u.T*gram/tnorm
    total = gram.inv()*w.T*neel**2*w
    projector = (reduced**2+(20+x)*reduced+(x*x+20*x-16)*s.eye(3))/polynomial.diff().as_expr()
    line = s.factor(s.trace(projector*response))
    success = s.factor(s.trace(projector*total))
    assert s.cancel(line-(23*x*x+24*x-720)/(3*(3*x*x+40*x-16))) == 0
    assert s.cancel(success-32*(x*x+6*x-24)/(3*(3*x*x+40*x-16))) == 0
    interval = tuple(map(s.Rational, polynomial.intervals(eps=s.Rational(1,10**70))[0][0]))
    result = {'ground_polynomial': polynomial, 'ground_interval': interval,
        'line': line/3600, 'success': success/3600, 'outside': s.factor((success-line)/3600),
        'line_moments': [s.trace(reduced**j*response) for j in range(3)],
        'total_moments': [s.trace(reduced**j*total) for j in range(3)]}
    for name in ('line','success','outside'):
        numerator, denominator = s.fraction(s.factor(result[name]))
        bounds = divide(evaluate(s.Poly(numerator,x),interval),evaluate(s.Poly(denominator,x),interval))
        assert 0 < bounds[0] <= bounds[1] < 1
        result[name+'_interval'] = tuple(map(s.Rational,outward(bounds,64)))
    assert result['outside_interval'][1] < result['line_interval'][0]
    return result


def coordinate_join():
    plain, twisted, unitary = full_operators()
    checks = phases = 0
    for axis, bit in enumerate((4,2,1)):
        native = sum(((-1)**int(bool(site & bit))*pauli(site,axis)/2
            for site in range(8)),s.SparseMatrix(256,256,{}))
        neel = sum(((-1)**site.bit_count()*pauli(site,axis)/2
            for site in range(8)),s.SparseMatrix(256,256,{}))
        assert unitary.H*native*unitary == neel
        for direction in range(3):
            t=translation(direction)
            assert t*native*t.H == (-1)**int(axis==direction)*native
            phases += 1
        checks += 1
    return {'complete_256_state_native_to_Neel_intertwiners':checks,
        'actual_source_Fourier_phase_checks':phases,
        'Klein_coordinate_change_redefines_source_momentum':False,
        'first_line_selection': 'the actual native X_a spin Fourier component becomes a Neel component in the cube coordinates'}


def certificate():
    data=exact_data()
    return {'ground_4H_polynomial':str(data['ground_polynomial'].as_expr()),
        'ground_root_selection':'least real root',
        'ground_4H_isolating_interval':[str(v) for v in data['ground_interval']],
        'complete_ground_cubic_subspace_dimension':3,
        'complete_magnetization_zero_block_dimension':70,
        'first_triplet_magnetic_zero_rank':1,
        'raw_line_moments':[str(v) for v in data['line_moments']],
        'raw_total_moments':[str(v) for v in data['total_moments']],
        'exact_normalized_line_weight':str(data['line']),
        'exact_normalized_success_weight':str(data['success']),
        'exact_normalized_outside_weight':str(data['outside']),
        'line_probability_interval':[str(v) for v in data['line_interval']],
        'success_probability_interval':[str(v) for v in data['success_interval']],
        'outside_probability_interval':[str(v) for v in data['outside_interval']],
        'reference_vector_channels':3,
        'each_channel_is_one_third_of_total_weight':True,
        'normalized_covariant_pump_divisor_from_raw_single_component':'3600',
        'coordinate_join':coordinate_join()}
