"""Exact finite spin rotations and the source force outside linear spins."""
from functools import lru_cache
from fractions import Fraction as F
from itertools import product
import sympy as s


def comm(left, right):
    return left*right-right*left


def spins(sites):
    pauli = (s.Matrix([[0,1],[1,0]]), s.Matrix([[0,-s.I],[s.I,0]]), s.diag(1,-1))
    return [[s.kronecker_product(*[pauli[a]/2 if x == y else s.eye(2)
                                  for y in range(sites)]) for a in range(3)]
            for x in range(sites)]


@lru_cache(maxsize=1)
def source_frame_certificate():
    from model.native_ordered_ground import source
    rows = []
    for sites, edges in ((2,((0,1),)),(3,((0,1),(1,2)))):
        spin = spins(sites)
        size = 2**sites
        h = sum((sum((spin[x][a]*spin[y][a] for a in range(3)), s.zeros(size))
                 -s.eye(size)/4 for x,y in edges), s.zeros(size))
        h -= sum((s.Rational(1,4)*(-1)**x*spin[x][0]
                  for x in range(sites)), s.zeros(size))
        marshall = s.eye(size)
        for x in range(sites):
            if x % 2:
                marshall = marshall*(2*spin[x][2])
        h = marshall*h*marshall
        order = [sum(((word>>x)&1)<<(sites-1-x) for x in range(sites))
                 for word in range(size)]
        actual, scale = source.matrix(sites,edges,F(1),F(1,4))
        assert h.extract(order,order)*scale == s.Matrix(actual)
        rows.append({'sites':sites,'edges':[list(edge) for edge in edges],
                     'actual_source_integer_scale':scale,
                     'full_selected_field_matrix_matches':True})
    return rows


@lru_cache(maxsize=1)
def local_certificate():
    spin = spins(2)
    h = sum((spin[0][a]*spin[1][a] for a in range(3)), s.zeros(4))-s.eye(4)/4
    left, right = s.symbols('theta_x theta_y', real=True)
    delta = left-right
    rows = []
    for a in range(3):
        generator = left*spin[0][a]+right*spin[1][a]
        transverse = sum((spin[0][b]*spin[1][b] for b in range(3) if b != a), s.zeros(4))
        current = sum((s.LeviCivita(a,b,c)*spin[0][b]*spin[1][c]
                       for b,c in product(range(3), repeat=2)), s.zeros(4))
        first = s.I*comm(h,generator)
        second = comm(generator,comm(h,generator))
        third = comm(generator,comm(generator,comm(generator,h)))
        assert (first+delta*current).applyfunc(s.expand) == s.zeros(4)
        assert (second+delta**2*transverse).applyfunc(s.expand) == s.zeros(4)
        assert (third+s.I*delta**3*current).applyfunc(s.expand) == s.zeros(4)
        assert sorted(transverse.eigenvals(), key=s.default_sort_key) == [-s.Rational(1,2),0,s.Rational(1,2)]
        assert sorted(current.eigenvals(), key=s.default_sort_key) == [-s.Rational(1,2),0,s.Rational(1,2)]
        # The bond rotates exactly within the transverse/current span.
        assert (comm(generator,transverse)+s.I*delta*current).applyfunc(s.expand) == s.zeros(4)
        assert (comm(generator,current)-s.I*delta*transverse).applyfunc(s.expand) == s.zeros(4)
        rows.append({'axis': a+1, 'transverse_norm': '1/2', 'current_norm': '1/2',
                     'first_term': '-delta current', 'second_term': '-delta^2 transverse',
                     'third_nested_norm': '|delta|^3/2'})
    return {'exact_two_site_axis_checks': rows,
            'source_bond': 'S_x dot S_y - 1/4',
            'scaled_third_order_remainder': 'sum_edges |delta|^3/(12p)'}


@lru_cache(maxsize=1)
def source_force_certificate():
    spin = spins(3)
    h = sum((sum((spin[x][a]*spin[y][a] for a in range(3)), s.zeros(8))-s.eye(8)/4
             for x,y in ((0,1),(1,2))), s.zeros(8))
    rows = []
    for axis in range(3):
        a = spin[0][axis]
        force = comm(h,comm(h,a))-(spin[0][axis]-spin[1][axis])/2
        assert all(s.trace(force*local) == 0 for site in spin for local in site)
        norm = s.simplify(s.trace(force.conjugate().T*force)/8)
        assert norm == s.Rational(1,16)
        rows.append({'axis':axis+1, 'force_trace_norm_squared_over_J4':'1/16',
                     'orthogonal_to_every_one_site_spin':True})
    return {'actual_three_site_path':[[0,1],[1,2]], 'endpoint_test':[1,0,0],
            'source_force_rows':rows,
            'finite_source_linear_spin_equation_is_exact':False,
            'finite_force_proves_nonzero_continuum_force':False}


def hidden_force_control():
    # Exact Hessian of two coupled canonical oscillators. Both observed
    # pure-direction third derivatives vanish, with scalar second derivatives.
    q,x,p,y,t = s.symbols('q x p y t', real=True)
    coupling = s.Rational(1,2)
    h = (p*p+y*y+q*q+x*x+2*coupling*q*x)/2
    shift_q = s.expand(h.subs(q,q+t)-h)
    shift_p = s.expand(h.subs(p,p+t)-h)
    assert s.expand(shift_q-t*(q+x/2)-t*t/2) == 0
    assert s.expand(shift_p-t*p-t*t/2) == 0
    stiffness = s.Matrix([[1,coupling],[coupling,1]])
    assert stiffness.eigenvals() == {s.Rational(1,2):1,s.Rational(3,2):1}
    return {'model_is_the_ordered_spin_source':False,
            'positive_coupled_oscillator_stiffness':[['1','1/2'],['1/2','1']],
            'observed_pure_direction_second_BCH_coefficients':['1','1'],
            'all_pure_direction_BCH_terms_above_second_vanish':True,
            'observed_second_time_derivative':'-q-x/2',
            'observed_response_frequency_squares':['1/2','3/2'],
            'both_observed_response_residues_are_nonzero':True,
            'scalar_quadratic_twists_imply_one_observed_frequency':False}


def reset_energy_certificate():
    """Exact complete reset of one actual bond inside a three-site source."""
    spin = spins(3)
    inside = sum((spin[0][a]*spin[1][a] for a in range(3)),s.zeros(8))-s.eye(8)/4
    crossing = sum((spin[1][a]*spin[2][a] for a in range(3)),s.zeros(8))-s.eye(8)/4
    h = inside+crossing
    singlet = s.Matrix([0,1,-1,0])/s.sqrt(2)
    kraus = [s.kronecker_product(singlet*s.eye(4)[:,i].T,s.eye(2)) for i in range(4)]
    assert sum((a.conjugate().T*a for a in kraus),s.zeros(8)) == s.eye(8)
    image = lambda matrix:sum((a.conjugate().T*matrix*a for a in kraus),s.zeros(8))
    assert image(inside) == -s.eye(8)
    delta = image(crossing)-crossing
    stability = sum((a.conjugate().T*comm(h,a) for a in kraus),s.zeros(8))
    assert stability == -s.eye(8)-inside+delta
    assert all(-1 <= value <= 1 for value in delta.eigenvals())
    assert all(value >= 0 for value in (inside+s.eye(8)).eigenvals())
    return {'actual_path_edges':[[0,1],[1,2]],'reset_region':[0,1],
            'complete_reset_Kraus_count':4,
            'finite_region_ground_energy_over_J':'-1',
            'boundary_channel_change_spectrum_over_J':sorted(map(str,delta.eigenvals())),
            'ground_stability_identity_exact':True,
            'negative_projector_bonds_reduce_reset_cost_to_J_per_boundary_edge':True}
