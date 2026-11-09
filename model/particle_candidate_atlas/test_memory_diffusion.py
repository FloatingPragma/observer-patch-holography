"""Independent three-ray, full Nambu and CP-variance controls."""
from itertools import combinations
import json

import pytest
import sympy as s

from . import memory_diffusion as m


def clean(a):
    return s.Matrix(a).applyfunc(lambda x: s.factor(s.radsimp(s.expand(x))))


def comm(a, b):
    return a*b-b*a


def independent(duty):
    # Literal source frame from slow_geometry.raw_frame, with independent
    # halfturn construction and the proved stationary three-ray lumping.
    phi = (1+s.sqrt(5))/2
    vectors = [s.Matrix(v) for v in ((-1,-phi,0),(1,-phi,0),(-phi,0,-1))]
    p = [clean(v*v.T/(v.T*v)[0]) for v in vectors]
    def halfturn(v): return clean(2*v*v.T/(v.T*v)[0]-s.eye(3))
    r = clean(halfturn(vectors[0]+vectors[2])*halfturn(vectors[0]+vectors[1]))
    a, b, q = s.Rational(5,72), s.Rational(1,72), s.Rational(1,6)
    t = s.Matrix([[1-a,1-b,1-a],[a,0,0],[0,b,a]])
    pi = s.Matrix([1-a,a*(1-a),a*b])/(1-a*a+a*b)
    g = (s.eye(3)-t+pi*s.ones(1,3)).inv()-pi*s.ones(1,3)
    e = {'always':t, 'C':s.Matrix([[b,a,b],[a,0,0],[0,b,a]]),
         'E':s.Matrix([[q,q-b,q],[0,0,0],[0,b,0]]),
         'CE':s.Matrix([[b,0,b],[0,0,0],[0,b,0]])}[duty]
    w = e*pi; k = s.diag(*w)/2-w*w.T/2+e*g*s.diag(*w); d = k+k.T
    def action(x):
        return clean(-sum((d[i,j]*comm(p[i],comm(p[j],x))
                          for i in range(3) for j in range(3)),s.zeros(3))/2)
    return r,p,d,action


@pytest.mark.parametrize('duty', m.POLICIES)
def test_independent_three_ray_source_matches_all_twelve_inputs(duty):
    r,p,d,action = independent(duty)
    row = m.exact_response(duty)
    assert row['R'] == r
    assert row['raw_records_over_all_old_inputs'] == 1536
    assert all(sum(row['T'][:,j]) == 1 for j in range(12))
    for x in (r,r+r*r,s.I*(r-r*r)):
        assert action(x) == m.diffusion(x,duty)
    for n in (1,2,3):
        assert all(s.factor(d.extract(ix,ix).det()) >= 0 for ix in combinations(range(3),n))


@pytest.mark.parametrize('duty', m.POLICIES)
def test_full_six_dimensional_nambu_double_commutator(duty):
    r,p,d,action = independent(duty)
    x = s.Matrix([[0,1],[1,0]])
    h = [s.kronecker_product(x,pi) for pi in p]
    b = r+r*r
    for spin in (s.eye(2),x):
        observable = s.kronecker_product(spin,b)
        direct = clean(-sum((d[i,j]*comm(h[i],comm(h[j],observable))
                             for i in range(3) for j in range(3)),s.zeros(6))/2)
        assert clean(direct-s.kronecker_product(spin,action(b))) == s.zeros(6)


def test_exact_always_and_C_duty_witnesses():
    r = m.exact_response()['R']; h = r+r*r; j = s.I*(r-r*r)
    assert s.expand(clean(comm(m.diffusion(h),r))[0,1]+33292836*(s.sqrt(5)-1)/2151685171) == 0
    assert clean(comm(m.diffusion(j),r))[0,0] == -s.I*s.Rational(52427031,34426962736)
    assert s.expand(clean(comm(m.diffusion(h,'C'),r))[0,0]+8729149*(s.sqrt(5)-1)/4800041280) == 0


@pytest.mark.parametrize('duty', m.POLICIES)
def test_real_and_imaginary_leaks_cannot_cancel_in_root_or_actual_effect(duty):
    r = m.exact_response(duty)['R']; h = r+r*r; j = s.I*(r-r*r)
    lh = clean(comm(m.diffusion(h,duty),r)); lj = clean(comm(m.diffusion(j,duty),r))
    assert lh != s.zeros(3) and lj != s.zeros(3)
    assert clean(lh.conjugate()-lh) == s.zeros(3)
    assert clean(lj.conjugate()+lj) == s.zeros(3)
    # Exact nonzero phase inside the positive chamber: tan(theta/2)=1/10.
    z = s.Rational(99,101)+s.I*s.Rational(20,101)
    a,b = s.sqrt(10/s.Integer(72)),s.sqrt(5/s.Integer(72))
    root = a*s.eye(3)+b*(z*r+s.conjugate(z)*r*r)
    for observable in (root,clean(root*root/3)):
        assert clean(comm(m.diffusion(observable,duty),r)) != s.zeros(3)


@pytest.mark.parametrize('duty', m.POLICIES)
def test_root_diffusion_and_squared_root_diffusion_are_distinct(duty):
    r,p,d,action = independent(duty)
    h = r+r*r
    variance = clean(action(h*h)-action(h)*h-h*action(h))
    gram = clean(sum((d[i,j]*comm(p[i],h).H*comm(p[j],h)
                      for i in range(3) for j in range(3)),s.zeros(3)))
    assert variance == gram
    assert s.simplify(s.trace(variance)) > 0
    for n in (1,2,3):
        assert all(s.simplify(variance.extract(ix,ix).det()) >= 0
                   for ix in combinations(range(3),n))


@pytest.mark.parametrize('matrix,duty',[(s.eye(2),'always'),(s.eye(3),'conditioned_success'),
                                      (s.eye(3)*s.Float('.5'),'C')])
def test_invalid_domain_is_rejected(matrix,duty):
    with pytest.raises(ValueError): m.diffusion(matrix,duty)


def test_certificate_is_json_safe_detached_and_scoped():
    result = m.certificate(); expected = json.dumps(result,sort_keys=True)
    assert len(result['duties']) == 4
    assert result['chi_response_rate_derived'] is False
    assert result['isospectral_match_after_family_rotation_excluded'] is False
    assert result['physical_mass_assigned'] is None
    result['duties'][0]['diffusion_H'][0][0] = 'forged'
    assert json.dumps(m.certificate(),sort_keys=True) == expected
