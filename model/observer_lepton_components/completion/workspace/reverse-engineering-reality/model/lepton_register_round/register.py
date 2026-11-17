"""The two-object generation register: the 432-state tick, its branch weights and its compression (section 1).

A seam object has a load in {tail, head, idle}, a record r in {0, 1} and an address o in {+, -}; o = + is the screen
orientation tail -> head of the seam about the face f, o = - its reverse. Basis index 4 l + 2 r + o. The object's
event is `Observer.repair` of the public program at the two addresses of its seam, read off basis state by basis
state. The neutral register F holds |x_0>, |x_1>, |x_2> with R|x_j> = |x_(j+1)>.
"""
from fractions import Fraction as F
from functools import lru_cache

import mpmath as mp
import numpy as np
import sympy as s

from model.coherent_seam_repair.program import Observer
from . import geometry

TAIL, HEAD, IDLE = 0, 1, 2
ACTIVE = (TAIL, HEAD)


def idx(l, r, o):
    return 4 * l + 2 * r + o


def screen_object(seam):
    """(tail, head) of a seam in the screen orientation about the face f."""
    a, b = seam
    return (a, b) if geometry.screen_sign(geometry.FACE, a, b) > 0 else (b, a)


@lru_cache(maxsize=None)
def event(seam):
    """The object's event read from `Observer.repair` at the two addresses of the seam.

    Returns (Q, off_seam_identity): Q as a 12 x 12 tuple of fractions, and whether a load on any other port is left
    unchanged at both addresses.
    """
    obs = Observer()
    tail, head = screen_object(seam)
    port = {TAIL: tail, HEAD: head}
    address = {0: obs.addresses.index((tail, head)), 1: obs.addresses.index((head, tail))}
    q = [[F(0)] * 12 for _ in range(12)]
    supported = True
    for l in range(3):
        for r in range(2):
            for o in range(2):
                col = idx(l, r, o)
                if l == IDLE:
                    q[col][col] = F(1)
                    continue
                local = {obs.index(port[l2], r2, address[o]): idx(l2, r2, o) for l2 in ACTIVE for r2 in range(2)}
                out = obs.repair(obs.basis_state(port[l], r, address[o])).amplitudes
                for i, amp in enumerate(out):
                    if amp == 0:
                        continue
                    if i not in local:
                        supported = False
                        continue
                    q[local[i]][col] = F(amp)
    idle_fixed = True
    for p in range(12):
        if p in (tail, head):
            continue
        for r in range(2):
            for o in range(2):
                state = obs.basis_state(p, r, address[o])
                idle_fixed &= obs.repair(state).amplitudes == state.amplitudes
    return tuple(tuple(row) for row in q), supported and idle_fixed


def _mat(rows):
    if any(F(x).denominator != 1 for row in rows for x in row):
        raise ArithmeticError('integer matrix expected')
    return np.array([[int(x) for x in row] for row in rows], dtype=np.int64)


def object_operators(seam):
    """Integer matrices 2Q, R2 (orientation exchange), 2W = 2 Q R2 and the oriented record V."""
    q, _ = event(seam)
    two_q = _mat([[2 * x for x in row] for row in q])
    r2 = np.zeros((12, 12), dtype=np.int64)
    v = np.zeros((12, 12), dtype=np.int64)
    for l in range(3):
        for r in range(2):
            for o in range(2):
                r2[idx(l, r, 1 - o), idx(l, r, o)] = 1
                v[idx(l, r, o), idx(l, r, o)] = -1 if (r == 1 and o == 1) else 1
    return two_q, r2, two_q @ r2, v


def step_of(i):
    """The register shift applied by T: R when C repairs at +, R^-1 at -, nothing otherwise."""
    l, r, o = i // 4, (i // 2) % 2, i % 2
    return (1 if o == 0 else -1) if r == 1 else 0


def full_tick(oriented=True):
    """4 U0 on C (x) E (x) F as an integer 432 x 432 matrix, U0 = T V (W_C (x) W_E); index 36 iC + 3 iE + x."""
    _, _, w_c, v_c = object_operators(geometry.C_SEAM)
    _, _, w_e, v_e = object_operators(geometry.E_SEAM)
    if not oriented:
        v_c = v_e = np.eye(12, dtype=np.int64)
    t = np.zeros((432, 432), dtype=np.int64)
    for ic in range(12):
        for ie in range(12):
            for x in range(3):
                t[36 * ic + 3 * ie + (x + step_of(ic)) % 3, 36 * ic + 3 * ie + x] = 1
    i3 = np.eye(3, dtype=np.int64)
    return t @ np.kron(np.kron(v_c, v_e), i3) @ np.kron(np.kron(w_c, w_e), i3)


def preparation():
    """sqrt6 |in_obj> = t + sqrt2 i, with t and i the tail and idle indicators at record 0 and both addresses."""
    t = np.zeros(12, dtype=np.int64)
    i = np.zeros(12, dtype=np.int64)
    for o in range(2):
        t[idx(TAIL, 0, o)] = 1
        i[idx(IDLE, 0, o)] = 1
    return t, i


def analyzer_parts():
    """4 sqrt3 Phi = sqrt2 h + 2 (e_+ + e_-) on C (x) E: h = C idle (x) E agrees, e_+- = C repairs +- (x) E idle."""
    h = np.zeros(144, dtype=np.int64)
    e_plus = np.zeros(144, dtype=np.int64)
    e_minus = np.zeros(144, dtype=np.int64)
    for oc in range(2):
        for le in ACTIVE:
            for oe in range(2):
                h[12 * idx(IDLE, 0, oc) + idx(le, 0, oe)] = 1
    for lc in ACTIVE:
        for oe in range(2):
            e_plus[12 * idx(lc, 1, 0) + idx(IDLE, 0, oe)] = 1
            e_minus[12 * idx(lc, 1, 1) + idx(IDLE, 0, oe)] = 1
    return h, e_plus, e_minus


def compression(oriented=True):
    """M(theta) = M_0 + z M_+ + z^-1 M_-, each M_c = (A_c + sqrt2 B_c)/(96 sqrt3) with integer A_c, B_c."""
    n = full_tick(oriented)
    t, i = preparation()
    r0 = np.kron(t, t) + 2 * np.kron(i, i)         # 6 |in_C (x) in_E> = r0 + sqrt2 r1
    r1 = np.kron(t, i) + np.kron(i, t)
    h, e_plus, e_minus = analyzer_parts()
    parts = {c: (np.zeros((3, 3), dtype=np.int64), np.zeros((3, 3), dtype=np.int64)) for c in ('0', '+', '-')}
    for x in range(3):
        ex = np.zeros(3, dtype=np.int64)
        ex[x] = 1
        out0, out1 = n @ np.kron(r0, ex), n @ np.kron(r1, ex)
        for xp in range(3):
            sel = np.zeros(3, dtype=np.int64)
            sel[xp] = 1
            hp, ep, em = np.kron(h, sel), np.kron(e_plus, sel), np.kron(e_minus, sel)
            parts['0'][0][xp, x] = 2 * hp @ out1
            parts['0'][1][xp, x] = hp @ out0
            parts['+'][0][xp, x] = 2 * ep @ out0
            parts['+'][1][xp, x] = 2 * ep @ out1
            parts['-'][0][xp, x] = 2 * em @ out0
            parts['-'][1][xp, x] = 2 * em @ out1
    return parts


def shift(power):
    m = np.zeros((3, 3), dtype=np.int64)
    for j in range(3):
        m[(j + power) % 3, j] = 1
    return m


def object_weights(seam):
    """Weights of idle, agree (record 0 on the seam), repair + and repair - after V W on the prepared object."""
    _, _, w, v = object_operators(seam)
    t, i = preparation()
    a, b = v @ w @ t, v @ w @ i                         # 2 sqrt6 psi = a + sqrt2 b
    disjoint = bool(np.all(a * b == 0))
    weight = [F(int(x * x + 2 * y * y), 24) for x, y in zip(a, b)]
    out = {'idle': F(0), 'agree': F(0), 'repair +': F(0), 'repair -': F(0)}
    for j, wj in enumerate(weight):
        l, r, o = j // 4, (j // 2) % 2, j % 2
        key = 'idle' if l == IDLE else 'agree' if r == 0 else 'repair +' if o == 0 else 'repair -'
        out[key] += wj
    return out, disjoint


def joint_weights():
    """P(C repairs +), P(E agrees) and P(both) from the full 432-state vector U0 |in, x_0>."""
    n = full_tick(True)
    t, i = preparation()
    e0 = np.array([1, 0, 0], dtype=np.int64)
    a = n @ np.kron(np.kron(t, t) + 2 * np.kron(i, i), e0)
    b = n @ np.kron(np.kron(t, i) + np.kron(i, t), e0)   # 24 U0 |in, x_0> = a + sqrt2 b
    disjoint = bool(np.all(a * b == 0))
    p_rp = p_h = p_both = total = F(0)
    for j in range(432):
        w = F(int(a[j] * a[j] + 2 * b[j] * b[j]), 576)
        if not w:
            continue
        ic, ie = j // 36, (j // 3) % 12
        rp = ic // 4 in ACTIVE and (ic // 2) % 2 == 1 and ic % 2 == 0
        hold = ie // 4 in ACTIVE and (ie // 2) % 2 == 0
        total += w
        p_rp += w * rp
        p_h += w * hold
        p_both += w * (rp and hold)
    return dict(total=total, repair_plus=p_rp, hold=p_h, both=p_both), disjoint


def eigen_checks():
    """Characters of R, the eigenvalues of M(theta) on them, the Koide value and the positive chamber (sympy)."""
    th = s.Symbol('theta', real=True)
    z = s.Symbol('z')                                       # z = exp(i theta)
    w = (-1 + s.sqrt(3) * s.I) / 2                          # exp(2 pi i/3)
    R = s.Matrix(shift(1).tolist())
    a, b = s.sqrt(s.Rational(1, 6)), s.sqrt(s.Rational(1, 12))
    M = s.sqrt(2) / 3 * (a * s.eye(3) + b * (z * R + R.inv() / z))
    ok_vec = True
    for k in range(3):
        v = s.Matrix([w ** ((-j * k) % 3) for j in range(3)]) / s.sqrt(3)
        ok_vec &= (R * v - w ** k * v).applyfunc(lambda e: s.expand(e)) == s.zeros(3, 1)
        mu = s.sqrt(2) / 3 * (a + b * (z * w ** k + w ** ((-k) % 3) / z))
        ok_vec &= (M * v - mu * v).applyfunc(lambda e: s.simplify(s.expand(e))) == s.zeros(3, 1)
        # z w^k + conj(z w^k) = 2 cos(theta + 2 pi k/3) on z = exp(i theta)
        ok_vec &= s.simplify(s.expand_complex((z * w ** k + w ** ((-k) % 3) / z).subs(z, s.exp(s.I * th)))
                             - 2 * s.cos(th + 2 * s.pi * k / 3)) == 0
    lam = [s.sqrt(2) / 3 * (1 / s.sqrt(6) + s.cos(th + 2 * s.pi * k / 3) / s.sqrt(3)) for k in range(3)]
    c_lam = [x * 3 / s.sqrt(2) for x in lam]
    sum_sq = s.simplify(s.expand_trig(sum(x ** 2 for x in c_lam)))
    total = s.simplify(s.expand_trig(sum(c_lam)))
    edge = [s.simplify(x.subs(th, s.pi / 12)) for x in c_lam]
    return dict(
        eigenvectors_are_characters=bool(ok_vec),
        sum_of_squared_eigenvalues_of_C=str(sum_sq), sum_of_eigenvalues_of_C=str(total),
        koide_value=str(s.nsimplify(sum_sq / total ** 2)),
        eigenvalues_of_C_at_pi_over_12=[str(x) for x in edge])


def register_outputs(theta=F(2, 9), digits=40):
    """Squared eigenvalue ratios of C(theta): electron, muon, tau are the characters w^1, w^2, w^0."""
    with mp.workdps(digits):
        th = mp.mpf(theta.numerator) / theta.denominator
        lam = [1 / mp.sqrt(6) + mp.cos(th + 2 * mp.pi * k / 3) / mp.sqrt(3) for k in range(3)]
        e, mu, tau = lam[1] ** 2, lam[2] ** 2, lam[0] ** 2
        return dict(theta=str(theta), eigenvalues_of_C=[mp.nstr(x, 12) for x in lam],
                    all_positive=all(x > 0 for x in lam), below_pi_over_12=bool(th < mp.pi / 12),
                    mu_over_e=mp.nstr(mu / e, 15), tau_over_e=mp.nstr(tau / e, 15),
                    _mu=float(mu / e), _tau=float(tau / e))


def read_effect_trace(parts):
    """tr(M^dag M) for |z| = 1 when the three parts have disjoint supports and A_c B_c = 0 entrywise."""
    supports = [((p[0] != 0) | (p[1] != 0)) for p in parts.values()]
    disjoint = not any(np.any(x & y) for i, x in enumerate(supports) for y in supports[i + 1:])
    no_cross = all(np.all(p[0] * p[1] == 0) for p in parts.values())
    total = sum(int(np.sum(p[0] * p[0] + 2 * p[1] * p[1])) for p in parts.values())
    return F(total, 96 * 96 * 3), disjoint and no_cross


def result():
    geo = geometry.summary()
    q_c, ok_c = event(geometry.C_SEAM)
    q_e, ok_e = event(geometry.E_SEAM)
    two_q, r2, two_w, v = object_operators(geometry.C_SEAM)
    n = full_tick(True)
    n_plain = full_tick(False)
    parts = compression(True)
    plain = compression(False)
    R, Ri, I3 = shift(1), shift(-1), np.eye(3, dtype=np.int64)
    zero = np.zeros((3, 3), dtype=np.int64)
    w_c, dis_c = object_weights(geometry.C_SEAM)
    w_e, dis_e = object_weights(geometry.E_SEAM)
    joint, dis_j = joint_weights()
    p_rp, p_h, p_both = joint['repair_plus'], joint['hold'], joint['both']
    count = p_rp + p_h - 2 * p_both
    eig = eigen_checks()
    out = register_outputs()
    h, e_plus, e_minus = analyzer_parts()
    exact_form = (np.array_equal(parts['0'][0], 32 * I3) and np.array_equal(parts['0'][1], zero)
                  and np.array_equal(parts['+'][0], zero) and np.array_equal(parts['+'][1], 16 * R)
                  and np.array_equal(parts['-'][0], zero) and np.array_equal(parts['-'][1], 16 * Ri))
    # a = A_0/(96 sqrt3) = 32/(96 sqrt3) and |b| = 16 sqrt2/(96 sqrt3), so (b/a)^2 = 2 * 16^2/32^2
    b_over_a_sq = F(2 * int(parts['+'][1][1, 0]) ** 2, int(parts['0'][0][0, 0]) ** 2)
    trace, trace_exact = read_effect_trace(parts)
    checks = dict(
        face_is_fixed_by_the_two_seam_half_turns=(geo['half_turn_product_order'] == 3 and
                                                  geo['faces_fixed_by_half_turn_product'] == [[4, 5, 11], [6, 7, 8]]),
        hexagon_and_screen_orientations=(geo['hexagon_in_screen_sense'] == [2, 3, 9, 1, 0, 10]
                                         and geo['c_seam_screen_orientation'] == [2, 3]
                                         and geo['e_seam_screen_orientation'] == [10, 2]
                                         and geo['third_turn_advances_hexagon_by'] == [2]),
        object_event_is_observer_repair_on_the_seam=ok_c and ok_e and q_c == q_e,
        object_event_is_an_involution=bool(np.array_equal(two_q @ two_q, 4 * np.eye(12, dtype=np.int64))),
        object_factors_unitary=bool(np.array_equal(two_w.T @ two_w, 4 * np.eye(12, dtype=np.int64))
                                    and np.array_equal(r2 @ r2, np.eye(12, dtype=np.int64))
                                    and np.array_equal(v @ v, np.eye(12, dtype=np.int64))),
        full_432_state_tick_unitary=bool(np.array_equal(n.T @ n, 16 * np.eye(432, dtype=np.int64))),
        object_weights_2_3_1_6_1_12_1_12=(all(w == {'idle': F(2, 3), 'agree': F(1, 6), 'repair +': F(1, 12),
                                                    'repair -': F(1, 12)} for w in (w_c, w_e))
                                          and dis_c and dis_e),
        analyzer_flags_orthonormal=(int(h @ h) == 8 and int(e_plus @ e_plus) == 4 and int(e_minus @ e_minus) == 4
                                    and int(h @ e_plus) == int(h @ e_minus) == int(e_plus @ e_minus) == 0),
        compression_is_sqrt2_over_3_C_theta=exact_form,
        compression_commutes_with_R=all(np.array_equal(m @ R, R @ m) for p in parts.values() for m in p),
        eigenvalues_on_the_characters=eig['eigenvectors_are_characters'],
        balance_b_over_a_is_one_over_sqrt2=exact_form and b_over_a_sq == F(1, 2),
        koide_two_thirds_in_the_positive_chamber=(eig['sum_of_squared_eigenvalues_of_C'] == '1'
                                                  and eig['koide_value'] == '2/3'
                                                  and eig['eigenvalues_of_C_at_pi_over_12'][1] == '0'),
        plain_record_gives_quadrature_form=(np.array_equal(plain['0'][0], 32 * I3)
                                            and np.array_equal(plain['+'][1], 16 * R)
                                            and np.array_equal(plain['-'][1], -16 * Ri)),
        repair_plus_1_12_hold_1_6_joint_1_72=(p_rp, p_h, p_both, joint['total']) == (F(1, 12), F(1, 6), F(1, 72), 1)
                                             and dis_j,
        exclusive_read_count_is_2_9=count == F(2, 9) == F(1, 12) + F(1, 6) - 2 * F(1, 72),
        trace_of_the_read_effect_equals_the_count=trace_exact and trace == count,
        register_outputs_at_theta_2_9=(out['all_positive'] and out['below_pi_over_12']
                                       and abs(out['_mu'] - 206.770316) < 5e-7 and abs(out['_tau'] - 3477.4728) < 5e-5))
    data = dict(
        geometry=geo,
        object_event_C_seam={'basis': 'index 4 l + 2 r + o; l tail, head, idle; o + screen, - reversed',
                             'two_Q': [[int(x) for x in row] for row in two_q]},
        object_weights={'C': {k: str(x) for k, x in w_c.items()}, 'E': {k: str(x) for k, x in w_e.items()}},
        compression={c: {'A': p[0].tolist(), 'B_sqrt2': p[1].tolist()} for c, p in parts.items()},
        compression_reading='M = sum_c z^c (A_c + sqrt2 B_c)/(96 sqrt3); A_0 = 32 I, B_+ = 16 R, B_- = 16 R^-1 gives '
                            '(sqrt2/3) C(theta), C = sqrt(1/6) I + sqrt(1/12)(z R + z^-1 R^-1), z = exp(i theta)',
        plain_record_compression={c: {'A': p[0].tolist(), 'B_sqrt2': p[1].tolist()} for c, p in plain.items()},
        b_over_a_squared=str(b_over_a_sq),
        eigenvalues=eig,
        read_counts={'P(R+)': str(p_rp), 'P(H)': str(p_h), 'P(R+ and H)': str(p_both), 'P(R+ xor H)': str(count)},
        trace_of_M_dag_M=str(trace),
        register_outputs={k: v for k, v in out.items() if not k.startswith('_')})
    return data, checks
