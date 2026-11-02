"""Exact finite controls of continuity, torque and current spectral weights."""
from functools import lru_cache
import copy
import sympy as s


def _spins(n):
    one = (s.Matrix([[0, 1], [1, 0]])/2,
           s.Matrix([[0, -s.I], [s.I, 0]])/2, s.diag(1, -1)/2)
    return [[s.kronecker_product(*[one[a] if j == x else s.eye(2)
                                  for j in range(n)]) for a in range(3)] for x in range(n)]


def _comm(a, b):
    return a*b-b*a


def _zero(a):
    return all(s.simplify(x) == 0 for x in a)


@lru_cache(maxsize=1)
def _checks():
    spins = _spins(3)
    dim = 8
    J, B = s.symbols('J B', real=True)
    f = s.symbols('f0:3', real=True)
    edges = [(0, 1), (1, 2)]
    h = J*sum((sum((spins[x][a]*spins[y][a] for a in range(3)), s.zeros(dim))
               -s.eye(dim)/4 for x, y in edges), s.zeros(dim))
    hb = h-B*sum(((-1)**x*spins[x][0] for x in range(3)), s.zeros(dim))
    continuity, torque = [], []
    for a in range(3):
        b, c = (a+1)%3, (a+2)%3
        observable = sum((f[x]*spins[x][a] for x in range(3)), s.zeros(dim))
        current = s.I*J*sum(((f[x]-f[y])*(spins[x][b]*spins[y][c]
                                -spins[x][c]*spins[y][b]) for x, y in edges), s.zeros(dim))
        field = -s.I*B*sum(((-1)**x*f[x]*sum((s.LeviCivita(0,a,b)*spins[x][b]
                        for b in range(3)), s.zeros(dim)) for x in range(3)), s.zeros(dim))
        continuity.append(_zero(_comm(h, observable)-current))
        torque.append(_zero(_comm(hb, observable)-current-field))
    # Periodic four-site singlet, J=1. This is a finite algebra control,
    # not a three-dimensional ordered-ground or small-momentum experiment.
    spins = _spins(4); dim = 16
    edges = [(x, (x+1)%4) for x in range(4)]
    h = sum((sum((spins[x][a]*spins[y][a] for a in range(3)), s.zeros(dim))
              -s.eye(dim)/4 for x,y in edges), s.zeros(dim))
    null = (h+3*s.eye(dim)).nullspace()
    assert len(null) == 1
    g = null[0]; g = g/s.sqrt((g.conjugate().T*g)[0])
    kgen = h+3*s.eye(dim)
    levels = (0, 1, 2, 3)
    assert _zero(kgen*(kgen-s.eye(dim))*(kgen-2*s.eye(dim))*(kgen-3*s.eye(dim)))
    weights = []
    for label, phase, midpoint, factor in [('pi/2', -s.I, (1-s.I)/s.sqrt(2), s.sqrt(2)),
                                           ('pi', -1, -s.I, s.Integer(2))]:
        a = sum((phase**x*spins[x][2]/2 for x in range(4)), s.zeros(dim))
        j = sum((midpoint*phase**x*(spins[x][0]*spins[y][1]-spins[x][1]*spins[y][0])/2
                 for x,y in edges), s.zeros(dim))
        assert _zero(_comm(h,a)+factor*j)
        sv, jv = a*g, j*g
        row = {'momentum':label,'lattice_factor':str(factor**2),'lines':[]}
        for e in levels:
            projector = s.eye(dim)
            for other in levels:
                if other != e:projector = projector*(kgen-other*s.eye(dim))/(e-other)
            sw = s.simplify((sv.conjugate().T*projector*sv)[0])
            jw = s.simplify((jv.conjugate().T*projector*jv)[0])
            assert s.simplify(e*e*sw-factor**2*jw) == 0
            row['lines'].append({'energy':str(e),'spin_weight':str(sw),'current_weight':str(jw)})
        row['current_inverse_moment'] = str(s.simplify(sum(s.Rational(x['current_weight'])/s.Rational(x['energy'])
                                                        for x in row['lines'] if x['energy'] != '0')))
        assert row['current_inverse_moment'] == '1/6'
        weights.append(row)
    # One open bond is a pure-gauge twist: its exact singlet energy is constant.
    sp = _spins(2); dim=4
    h = sum((sp[0][a]*sp[1][a] for a in range(3)), s.zeros(dim))-s.eye(dim)/4
    g = s.Matrix([0,1,-1,0])/s.sqrt(2)
    hp = s.I*_comm(sp[0][2],h)
    hpp = -_comm(sp[0][2],_comm(sp[0][2],h))
    diamagnetic = s.simplify((g.T*hpp*g)[0])
    chi = s.simplify((g.T*hp*(h+s.eye(dim))*hp*g)[0])
    assert diamagnetic == s.Rational(1,2) and chi == s.Rational(1,4)
    return {'all_axes_continuity':all(continuity),'all_axes_staggered_torque':all(torque),
            'four_site_ring_ground_energy':'-3','four_site_ring_kappa':'1/3',
            'four_site_ring_spectral_lines':weights,
            'open_bond_twist':{'bare_second_derivative':str(diamagnetic),
                               'twice_inverse_current_moment':str(2*chi),
                               'ground_energy_curvature':str(diamagnetic-2*chi)},
            'finite_controls_are_ordered_ground_scaling_data':False}


def exact_checks():
    value=copy.deepcopy(_checks())
    if value['all_axes_continuity'] is not True or value['all_axes_staggered_torque'] is not True:
        raise ArithmeticError('exact current identity failed')
    return value
