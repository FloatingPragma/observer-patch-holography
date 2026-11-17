"""The orientation sign of the public repair law, the kick as the commutator of the oriented toggle with the pair's
tick, and the two-tick crossing of a double load on its own seam (section 5).
"""
from fractions import Fraction as F

from model.coherent_seam_repair.program import Observer
from . import geometry, hexagon as hx, register as rg, walk as wk


def _matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def _lin(*terms):
    """sum of coefficient * matrix."""
    n = len(terms[0][1])
    return [[sum(c * m[i][j] for c, m in terms) for j in range(n)] for i in range(n)]


def seam_exchange():
    """K = |p,1><m,0| + |m,0><p,1| on the basis (a,0), (b,0), (a,1), (b,1) of one seam; m = (a - b)/sqrt2, p = (a + b)/sqrt2."""
    m0 = [F(1), F(-1), F(0), F(0)]
    p1 = [F(0), F(0), F(1), F(1)]
    return [[(p1[i] * m0[j] + m0[i] * p1[j]) / 2 for j in range(4)] for i in range(4)]


def orientation_sign():
    """On every seam: Observer.repair at (a, b) is 1 - K^2 + K and at (b, a) it is 1 - K^2 - K."""
    obs = Observer()
    K = seam_exchange()
    K2 = _matmul(K, K)
    I4 = [[F(int(i == j)) for j in range(4)] for i in range(4)]
    plus, minus = _lin((1, I4), (-1, K2), (1, K)), _lin((1, I4), (-1, K2), (-1, K))
    cubic = _matmul(K2, K) == K
    ok = True
    for a, b in obs.seams:
        basis = [(a, 0), (b, 0), (a, 1), (b, 1)]
        for address, expected in (((a, b), plus), ((b, a), minus)):
            k = obs.addresses.index(address)
            block = [[F(0)] * 4 for _ in range(4)]
            for col, (port, par) in enumerate(basis):
                out = obs.repair(obs.basis_state(port, par, k)).amplitudes
                support = {i for i, x in enumerate(out) if x != 0}
                inside = {obs.index(p2, r2, k) for p2, r2 in basis}
                ok &= support <= inside
                for row, (p2, r2) in enumerate(basis):
                    block[row][col] = F(out[obs.index(p2, r2, k)])
            ok &= block == expected
    return dict(all_seams=ok, seams=len(obs.seams), exchange_cubed_is_exchange=cubic,
                K=[[str(x) for x in row] for row in K])


def clocked_register():
    """The C object: J = |1><1|_r sigma_z(o) = (V_- - V_+)/2, Z_r = V_+ V_-, and the oriented law on blank records."""
    diag = {}
    for i in range(12):
        l, r, o = i // 4, (i // 2) % 2, i % 2
        diag[i] = dict(J=r * (1 - 2 * o), Vm=-1 if (r == 1 and o == 1) else 1, Vp=-1 if (r == 1 and o == 0) else 1,
                       Zr=1 - 2 * r)
    j_form = all(2 * d['J'] == d['Vm'] - d['Vp'] for d in diag.values())
    zr_form = all(d['Zr'] == d['Vp'] * d['Vm'] for d in diag.values())
    flag_values = {'C idle': {diag[rg.idx(rg.IDLE, 0, o)]['J'] for o in range(2)},
                   'C repairs +': {diag[rg.idx(l, 1, 0)]['J'] for l in rg.ACTIVE},
                   'C repairs -': {diag[rg.idx(l, 1, 1)]['J'] for l in rg.ACTIVE}}
    # the object's event from Observer.repair: Q = 1 - K^2 + sigma_z(o) K; V Q = 1 - K^2 + K on blank records
    q, _ = rg.event(geometry.C_SEAM)
    K = [[F(0)] * 12 for _ in range(12)]
    ks = seam_exchange()
    for o in range(2):
        cells = [rg.idx(rg.TAIL, 0, o), rg.idx(rg.HEAD, 0, o), rg.idx(rg.TAIL, 1, o), rg.idx(rg.HEAD, 1, o)]
        for i, ri in enumerate(cells):
            for j, cj in enumerate(cells):
                K[ri][cj] = ks[i][j]
    K2 = _matmul(K, K)
    I12 = [[F(int(i == j)) for j in range(12)] for i in range(12)]
    sz = [[F((1 - 2 * (i % 2)) * (i == j)) for j in range(12)] for i in range(12)]
    plain = _lin((1, I12), (-1, K2), (1, _matmul(sz, K)))
    oriented = _lin((1, I12), (-1, K2), (1, K))
    V = [[F(diag[i]['Vm'] * (i == j)) for j in range(12)] for i in range(12)]
    vq = _matmul(V, [list(row) for row in q])
    blank = [j for j in range(12) if (j // 2) % 2 == 0]
    return dict(j_form=j_form, zr_form=zr_form, flag_values={k: sorted(v) for k, v in flag_values.items()},
                plain_law=[list(r) for r in q] == plain,
                oriented_on_blank=all(vq[i][j] == oriented[i][j] for i in range(12) for j in blank),
                oriented_differs_on_written=any(vq[i][j] != oriented[i][j] for i in range(12) for j in range(12)
                                                if j not in blank))


def crossing():
    """Two ticks of a seam's own two-address walk on a blank double, on every seam from both ends."""
    seams, _, _, _ = geometry.carrier()
    plain_cross = fixed_return = four_ticks = 0
    for a, b in seams:
        plain = wk.Walk(((a, b), (b, a)))
        fixed = wk.Walk(((a, b), (b, a)), ((a, b), (a, b)))
        for x, y in ((a, b), (b, a)):
            for o in range(2):
                start = {(wk.config(x, x), 0, o): F(1)}
                two = plain.step(plain.step(start))
                plain_cross += two == {(wk.config(y, y), 0, o): F(1)}
                four_ticks += plain.step(plain.step(two)) == start
                fixed_return += fixed.step(fixed.step(start)) == start
    return dict(starts=4 * len(seams), plain_law_crosses=plain_cross, plain_law_returns_in_four=four_ticks,
                fixed_order_returns_in_two=fixed_return)


def result():
    sign = orientation_sign()
    reg = clocked_register()
    pair = {f: hx.round_numbers(f) for f in hx.LAWS}
    at2 = hx.round_numbers('plain', (hx.D2,), (hx.D2,))
    at10 = hx.round_numbers('plain', (hx.D10,), (hx.D2,))
    cross = crossing()
    c2, c10 = at2['overlap'], at10['overlap']
    checks = dict(
        repair_is_1_minus_K2_plus_sigma_z_K_on_every_seam=sign['all_seams'] and sign['exchange_cubed_is_exchange'],
        object_event_is_1_minus_K2_plus_sigma_z_K=reg['plain_law'],
        oriented_record_gives_1_minus_K2_plus_K_on_blank_records=reg['oriented_on_blank'] and reg['oriented_differs_on_written'],
        current_is_half_difference_of_oriented_records=reg['j_form'] and reg['zr_form'],
        flags_carry_current_0_plus_1_minus_1=reg['flag_values'] == {'C idle': [0], 'C repairs +': [1], 'C repairs -': [-1]},
        pair_toggle_identity_SS_minus_2X=all(p['toggle_identity'] for p in pair.values()),
        kick_is_the_toggle_commutator_with_h_2_eta=all(p['unflagged_start'] and p['toggle_is_kick'] for p in pair.values()),
        toggle_round_amplitude_1_minus_P_plus_P_exp_minus_2i_eta_j=all(p['overlap'] == F(2, 9) for p in pair.values()),
        double_crosses_its_seam_in_two_plain_ticks=(cross['plain_law_crosses'] == cross['starts']
                                                    and cross['plain_law_returns_in_four'] == cross['starts']),
        fixed_order_law_returns_the_double=cross['fixed_order_returns_in_two'] == cross['starts'],
        one_port_round_coefficients_5_18_and_1_6_mean_2_9=(c2, c10, (c2 + c10) / 2) == (F(5, 18), F(1, 6), F(2, 9)),
        one_port_non_returns_stay_at_the_port=all(
            r[o]['outside_the_blank_doubles'] == 0 and r[o]['uniform_address_part'] == 0
            for r in (at2, at10) for o in ('A', 'B')),
        one_port_no_weights_65_324_and_5_36=(at2['norm'] - c2 ** 2, at10['norm'] - c10 ** 2) == (F(65, 324), F(5, 36)))
    data = dict(
        seam_exchange_K=sign['K'], seams_checked=sign['seams'],
        current_on_flags=reg['flag_values'],
        toggle_round={f: {'<start|X>': str(p['overlap'])} for f, p in pair.items()},
        crossing=cross,
        one_port_round={'A at port 2, B at port 2': {'coefficient': str(c2), 'no weight over abs(1 - z)^2': str(at2['norm'] - c2 ** 2)},
                        'A at port 10, B at port 2': {'coefficient': str(c10), 'no weight over abs(1 - z)^2': str(at10['norm'] - c10 ** 2)},
                        'mean coefficient': str((c2 + c10) / 2)})
    return data, checks
