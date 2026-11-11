"""Double loads on the hexagon of the face f, and the weak round with the inverse tick (sections 2 and 3).

The hexagon object has the six hexagon ports, two loads, the four-valued record of the public completion and an
address register over the twelve ordered hexagon seams, reflected about its own uniform state. The plain law applies
the completion in the order of the address. The oriented law applies, at both addresses of a seam, the completion
in one fixed order whose first port lies in the orbit {2, 9, 0} of the face third-turn.
"""
from fractions import Fraction as F

import mpmath as mp
import sympy as s

from model.coherent_seam_repair.program import Observer
from . import geometry, walk as wk

RING = geometry.hexagon_seams()
ADDRESSES = tuple(p for (x, y) in RING for p in ((x, y), (y, x)))
_RHO = geometry.third_turn(geometry.FACE)
ORBIT = frozenset({2, _RHO[2], _RHO[_RHO[2]]})
D2, D10 = wk.config(2, 2), wk.config(10, 10)
PAIR_IN, PAIR_FWD, PAIR_BWD = wk.config(2, 10), wk.config(2, 3), wk.config(10, 0)
LAWS = ('oriented', 'plain')


def hexagon_walk(frame):
    if frame == 'plain':
        return wk.Walk(ADDRESSES)
    pairs = []
    for x, y in ADDRESSES:
        pairs.append((x, y) if x in ORBIT else (y, x))
    return wk.Walk(ADDRESSES, pairs)


def flag(key):
    cfg, rec, _ = key
    if rec in (2, 3):
        return {PAIR_FWD: 'forward', PAIR_IN: 'in-dimer', PAIR_BWD: 'backward'}.get(cfg, 'other')
    if rec == 0 and cfg in (D2, D10):
        return 'double'
    return 'other'


def standing(walk_, cfgs=(D2, D10)):
    """The blank double, uniform over the given ports and over the addresses; unit amplitudes and its norm squared."""
    state = walk_.uniform({c: F(1) for c in cfgs})
    return state, F(len(state))


def tick_weights(walk_, cfgs):
    state, n2 = standing(walk_, cfgs)
    return wk.weights(walk_.step(state), n2, flag)


def section_double_loads():
    u = wk.completion()
    involution = all(sum(u[i][k] * u[k][j] for k in range(24)) == (1 if i == j else 0)
                     for i in range(24) for j in range(24))
    symmetric = all(u[i][j] == u[j][i] for i in range(24) for j in range(24))
    rows, ok_s, ok_one = {}, True, True
    expected = {'double': F(2, 3), 'in-dimer': F(1, 6), 'forward': F(1, 12), 'backward': F(1, 12)}
    for frame in LAWS:
        w = hexagon_walk(frame)
        ws, w2, w10 = tick_weights(w, (D2, D10)), tick_weights(w, (D2,)), tick_weights(w, (D10,))
        p, q = ws['forward'], ws['in-dimer']
        rows[frame] = {'symmetric double': {k: str(v) for k, v in sorted(ws.items())},
                       'double at port 2': {k: str(v) for k, v in sorted(w2.items())},
                       'double at port 10': {k: str(v) for k, v in sorted(w10.items())},
                       'P(R+ xor H) of two objects': str(p + q - 2 * p * q)}
        ok_s &= ws == expected and p + q - 2 * p * q == F(2, 9)
        ok_one &= (w2 == {'double': F(2, 3), 'forward': F(1, 6), 'in-dimer': F(1, 6)}
                   and w10 == {'double': F(2, 3), 'backward': F(1, 6), 'in-dimer': F(1, 6)})
    checks = dict(
        completion_is_a_real_symmetric_involution=involution and symmetric,
        symmetric_double_tick_weights_2_3_1_6_1_12_1_12_and_count_2_9=ok_s,
        one_port_doubles_emit_1_6_across_each_of_their_two_seams=ok_one)
    data = dict(hexagon_ports=list(geometry.hexagon()), screen_seams=[list(x) for x in RING],
                third_turn_orbit=sorted(ORBIT), tick_one=rows)
    return data, checks, expected


# ---------------------------------------------------------------------------------------------------------- weak round
def flags_of(key_pair):
    """(f, i): object A read for the forward flag, object B read for the in-dimer flag."""
    return int(flag(key_pair[0]) == 'forward'), int(flag(key_pair[1]) == 'in-dimer')


def round_numbers(frame, a_cfgs=(D2, D10), b_cfgs=(D2, D10)):
    """Exact numbers of one weak round on two independent objects: A read for R+ (forward), B read for H (in-dimer).

    X = (W^-1 x W^-1) Qx (W x W)|start>, Qx the projector on (R+ of A) xor (H of B). The round with the kick
    exp(-i h j Qx) maps the start to |start> + (z - 1)|X>, z = exp(-i h j).
    """
    w = hexagon_walk(frame)
    SA, nA = standing(w, a_cfgs)
    SB, nB = standing(w, b_cfgs)
    WA, WB = w.step(SA), w.step(SB)
    back = w.inverse_step(WA) == SA and w.inverse_step(WB) == SB
    a = w.inverse_step(wk.restrict(WA, lambda k: flag(k) == 'forward'))      # W^-1 Pi_F W S_A
    b = w.inverse_step(wk.restrict(WB, lambda k: flag(k) == 'in-dimer'))     # W^-1 Pi_I W S_B
    s_a, s_b = wk.dot(SA, a) / nA, wk.dot(SB, b) / nB
    start, pn = wk.tensor(SA, SB), nA * nB
    ticked = wk.tensor(WA, WB)
    flagged = {k: v for k, v in ticked.items() if flags_of(k)[0] != flags_of(k)[1]}
    X = wk.on_both(flagged, w.inverse_step)
    overlap, norm = wk.dot(start, X) / pn, wk.norm2(X) / pn
    # (W^-1 x W^-1)(Z_F x Z_I)(W x W)|start> = |start> - 2|X>
    toggled = {k: v * (1 - 2 * flags_of(k)[0]) * (1 - 2 * flags_of(k)[1]) for k, v in ticked.items()}
    toggle_identity = wk.on_both(toggled, w.inverse_step) == wk.combine(start, X, -2)
    # The toggle round W^-1 T(eta) W T(eta)^-1 as a polynomial in u = exp(-i eta j): T(eta)^-1 multiplies the
    # unflagged start by u, T(eta) multiplies a ticked branch with flags (f, i) by u^(-(1 - 2f)(1 - 2i)).
    unflagged_start = all(flags_of(k) == (0, 0) for k in start)
    by_power = {}
    for k, v in ticked.items():
        f, i = flags_of(k)
        by_power.setdefault(1 - (1 - 2 * f) * (1 - 2 * i), {})[k] = v
    toggle_round = {m: wk.on_both(part, w.inverse_step) for m, part in by_power.items()}
    toggle_is_kick = toggle_round == {0: wk.combine(start, X, -1), 2: X}

    def content(v, s_, S, n2, cfgs, size=w.size):
        perp = wk.combine(v, S, -s_)
        total = wk.norm2(perp) / n2
        off = sum((x * x for (cfg, rec, k), x in perp.items() if rec != 0 or cfg not in cfgs), F(0)) / n2
        uniform = F(0)
        for cfg in cfgs:
            tot = sum((perp.get((cfg, 0, k), F(0)) for k in range(size)), F(0))
            uniform += tot * tot / size / n2
        return dict(norm2=total, outside_the_blank_doubles=off, uniform_address_part=uniform,
                    non_uniform_address_part=total - off - uniform)
    return dict(inverse_restores=back, s_A=s_a, s_B=s_b, P=s_a + s_b - 2 * s_a * s_b, overlap=overlap, norm=norm,
                toggle_identity=toggle_identity, unflagged_start=unflagged_start, toggle_is_kick=toggle_is_kick,
                A=content(a, s_a, SA, nA, a_cfgs), B=content(b, s_b, SB, nB, b_cfgs))


def phase_table(P, Ns=(1, 2, 10, 100, 1000)):
    """Phase per unit clock time and kept modulus of (1 - P + P exp(-i/N))^N, and the distance to exp(-i P)."""
    rows = {}
    with mp.workdps(40):
        pm = mp.mpf(P.numerator) / P.denominator
        for N in Ns:
            z = mp.exp(-1j / mp.mpf(N))
            amp = (1 - pm + pm * z) ** N
            rows[N] = {'phase per unit clock time': mp.nstr(-mp.arg(amp), 15), 'modulus kept': mp.nstr(abs(amp), 12),
                       'distance to exp(-i P)': mp.nstr(abs(amp - mp.exp(-1j * pm)), 6),
                       'expected resets per unit clock time':
                           mp.nstr(N * 4 * mp.sin(1 / (2 * mp.mpf(N))) ** 2 * pm * (1 - pm), 8),
                       '_phase': float(-mp.arg(amp)), '_distance': float(abs(amp - mp.exp(-1j * pm)))}
    return rows


def series(P):
    """theta(x) = atan(P sin x/(1 - P + P cos x))/x and log|1 - P + P exp(-ix)|/x^2 as series in x = 1/N."""
    x = s.Symbol('x')
    p = s.Rational(P.numerator, P.denominator)
    theta = s.series(s.atan(p * s.sin(x) / (1 - p + p * s.cos(x))) / x, x, 0, 4).removeO()
    modulus = s.series(s.log((1 - p + p * s.cos(x)) ** 2 + (p * s.sin(x)) ** 2) / (2 * x ** 2), x, 0, 2).removeO()
    return s.expand(theta), s.expand(modulus)


def forward_only(frame='oriented', ticks=6):
    w = hexagon_walk(frame)
    state, n2 = standing(w)
    start = dict(state)
    rows = []
    for _ in range(ticks):
        state = w.step(state)
        wt = wk.weights(state, n2, flag)
        p, q = wt.get('forward', F(0)), wt.get('in-dimer', F(0))
        rows.append(dict(weights={k: str(v) for k, v in sorted(wt.items())}, count=p + q - 2 * p * q))
    return rows, state == start


def public_round(port=0):
    """The strong round on the public one-load program with its own methods only."""
    obs = Observer()
    D = obs.address_count
    s0 = obs.embedded_port(port)
    s1 = obs.step(s0)
    flagged = [i for i, x in enumerate(s1.amplitudes) if x != 0 and (i // D) % 2 == 1]
    p = sum((s1.amplitudes[i] ** 2 for i in flagged), F(0)) / s1.norm_denominator
    s2 = s1
    for i in flagged:
        s2, _ = obs.quarter_turn(s2, i)
    s3 = obs.inverse_step(s2)
    amp = s.nsimplify(obs.embedded_return_amplitude(s3, port))
    branches = {rec['outcome']: s.nsimplify(st.norm_squared()) for rec, st in obs.embedded_return_measurement(s3, port)}
    return dict(inverse_of_step=obs.inverse_step(s1).amplitudes == s0.amplitudes, quarter_turns=len(flagged),
                p=p, amplitude=amp, yes=branches[True], no=branches[False])


def section_weak_round():
    nums = {f: round_numbers(f) for f in LAWS}
    P = nums['oriented']['P']
    table = phase_table(P)
    theta, modulus = series(P)
    x = s.Symbol('x')
    fwd, period = forward_only('oriented')
    counts = [r['count'] for r in fwd]
    pub = public_round()
    stated = {1: 0.205334030970278, 10: 0.222062099670503, 1000: 0.222222206218563}
    same = all(nums['plain'][k] == nums['oriented'][k] for k in ('s_A', 's_B', 'P', 'overlap', 'norm', 'A', 'B'))
    checks = {}
    for frame, n in nums.items():
        checks[f'{frame}_inverse_tick_restores_the_pair'] = n['inverse_restores']
    o = nums['oriented']
    checks.update(
        every_round_reads_the_tick_one_weights=o['s_A'] == F(1, 12) and o['s_B'] == F(1, 6) and P == F(2, 9),
        zero_return_amplitude_is_1_minus_P_plus_P_z=o['overlap'] == F(2, 9) and o['norm'] == F(2, 9),
        non_return_weight_is_14_81_per_unit_abs_1_minus_z_squared=o['norm'] - o['overlap'] ** 2 == F(14, 81),
        non_return_norms_11_144_and_5_36=o['A']['norm2'] == F(11, 144) and o['B']['norm2'] == F(5, 36),
        non_return_vectors_are_blank_doubles=(o['A']['outside_the_blank_doubles'] == 0
                                              and o['B']['outside_the_blank_doubles'] == 0
                                              and o['A']['uniform_address_part'] == F(1, 144)
                                              and o['A']['non_uniform_address_part'] == F(10, 144)
                                              and o['B']['uniform_address_part'] == 0),
        plain_and_oriented_rounds_agree=same,
        phase_per_unit_clock_time_at_N_1_10_1000=all(abs(table[N]['_phase'] - v) < 1e-12 for N, v in stated.items()),
        phase_series_2_9_minus_35_over_2187_N2=s.expand(theta - (s.Rational(2, 9) - s.Rational(35, 2187) * x ** 2)) == 0,
        modulus_series_minus_7_over_81_N=s.expand(modulus + s.Rational(7, 81)) == 0,
        distance_to_the_limit_at_most_1_over_N=all(r['_distance'] <= 1 / N for N, r in table.items()),
        forward_only_counts_2_9_35_81_38_81_35_81_2_9_0=counts == [F(2, 9), F(35, 81), F(38, 81), F(35, 81), F(2, 9), F(0)],
        forward_only_mean_8_27_and_period_six=sum(counts, F(0)) / 6 == F(8, 27) and period,
        public_one_load_round_returns_11_12_plus_i_12=(pub['inverse_of_step'] and pub['p'] == F(1, 12)
                                                       and s.simplify(pub['amplitude'] - (F(11, 12) + s.I / 12)) == 0
                                                       and pub['no'] == s.Rational(11, 72) and pub['yes'] == s.Rational(61, 72)))
    data = dict(
        round={f: {'<S|W^-1 Pi_F W|S>': str(n['s_A']), '<S|W^-1 Pi_I W|S>': str(n['s_B']), 'P': str(n['P']),
                   '<SS|X>': str(n['overlap']), '<X|X>': str(n['norm']),
                   'non_return_A': {k: str(v) for k, v in n['A'].items()},
                   'non_return_B': {k: str(v) for k, v in n['B'].items()}} for f, n in nums.items()},
        zero_return_amplitude='<SS|Psi(z)> = 1 + (z - 1) <SS|X> = 1 - P + P z, z = exp(-i h j)',
        rounds={str(N): {k: v for k, v in r.items() if not k.startswith('_')} for N, r in table.items()},
        phase_series=str(theta), log_modulus_times_N=str(modulus),
        forward_only_oriented=[{'tick': t + 1, 'weights': r['weights'], 'count': str(r['count'])}
                               for t, r in enumerate(fwd)],
        forward_only_mean=str(sum(counts, F(0)) / 6),
        public_one_load_round={'quarter_turns': pub['quarter_turns'], 'P(record 1)': str(pub['p']),
                               'return amplitude': str(pub['amplitude']), 'yes weight': str(pub['yes']),
                               'no weight': str(pub['no'])})
    return data, checks, nums
