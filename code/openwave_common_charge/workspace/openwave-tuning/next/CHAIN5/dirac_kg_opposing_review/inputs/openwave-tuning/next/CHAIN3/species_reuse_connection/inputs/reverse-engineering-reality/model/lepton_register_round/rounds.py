"""Many rounds of the weak round on the sixty-address carrier (PROOF.md section 8).

The complete record of the two return tests, the unconditional phase of the clocked system, partners that share the
kick, the round without the reset, the forward-only control, the register tick with record flags, and an execution
of ten thousand rounds. Rational arithmetic unless a tolerance is given. Every walk operation is the public one:
`Walk.step`, `Walk.inverse_step` and `images` of `walk.py` on the sixty ordered seams of `Observer.addresses`.
"""
from fractions import Fraction as F
from itertools import product
import cmath
import math

import mpmath as mp
import numpy as np
import sympy as s

try:                                    # inside model/lepton_register_round
    from . import carrier, walk as wk
except ImportError:                     # run from a working copy outside the package
    from model.lepton_register_round import carrier, walk as wk

PORT = 2
P = F(2, 9)


# ------------------------------------------------------------------------------------------------ cells
def standing(w, port=PORT):
    """The blank double at the port with amplitude one on every address (norm squared 60)."""
    return w.uniform({wk.config(port, port): F(1)})


def records_by_address(w, port=PORT):
    """The record that the event writes at each address from the blank double; the event maps |D, 0, k> to one
    basis state and back."""
    cfg0 = wk.config(port, port)
    recs, exact = [], True
    for a, b in w.addresses:
        im = wk.images(cfg0, 0, a, b)
        exact &= len(im) == 1 and im[0][2] == 1
        back = wk.images(im[0][0], im[0][1], a, b)
        exact &= len(back) == 1 and back[0] == (cfg0, 0, F(1))
        recs.append(im[0][1])
    return tuple(recs), exact


def pair_array(pair, port=PORT, size=60):
    """A pair dictionary as a size x size table over the addresses of blank doubles, and the weight elsewhere."""
    cfg0 = wk.config(port, port)
    table = [[F(0)] * size for _ in range(size)]
    outside = F(0)
    for ((c1, r1, k1), (c2, r2, k2)), v in pair.items():
        if c1 == cfg0 and c2 == cfg0 and r1 == 0 and r2 == 0:
            table[k1][k2] = v
        else:
            outside += v * v
    return table, outside


# ------------------------------------------------------------------------------------------------ complete record
def complete_record(w, port=PORT):
    """|start> -> |start> + (z - 1)|X>. Weights of the four outcomes of the two return tests in |X>, each test the
    projector on the uniform address state of its cell."""
    start = standing(w, port)
    ticked = w.step(start)
    flagged = {k: v for k, v in wk.tensor(ticked, ticked).items() if (k[0][1] == 2) != (k[1][1] in (2, 3))}
    X = wk.on_both(flagged, w.inverse_step)
    n = w.size
    table, outside = pair_array(X, port, n)
    n2 = F(n * n)
    total = sum((v * v for row in table for v in row), F(0)) / n2
    overlap = sum((v for row in table for v in row), F(0)) / n2
    col = [sum(row, F(0)) / n for row in table]                       # Y on its uniform state
    row_ = [sum((table[i][j] for i in range(n)), F(0)) / n for j in range(n)]
    mc, mr = sum(col, F(0)) / n, sum(row_, F(0)) / n
    x_no = sum(((c - mc) ** 2 for c in col), F(0)) * n / n2            # X no, Y yes
    y_no = sum(((r - mr) ** 2 for r in row_), F(0)) * n / n2           # X yes, Y no
    return dict(outside=outside, overlap=overlap, norm=total, yes_yes=overlap * overlap, x_no=x_no, y_no=y_no,
                both_no=total - overlap * overlap - x_no - y_no)


def multiplier_identity():
    """(1 - P + P z)(1 - P + P/w) + P (1 - P)(z - 1)(1/w - 1) = 1 - P + P z/w for unimodular z, w."""
    z, w_ = s.symbols('z w')
    p = s.Rational(2, 9)
    lhs = (1 - p + p * z) * (1 - p + p / w_) + p * (1 - p) * (z - 1) * (1 / w_ - 1)
    return s.simplify(lhs - (1 - p + p * z / w_)) == 0


# ------------------------------------------------------------------------------------------------ no reset
def without_reset(w, port=PORT):
    """Read weights of a cell whose return test said no, and the count of the pair in the four test states."""
    start = standing(w, port)
    ticked = w.step(start)
    n2 = F(w.size)
    out = {}
    for name, pred in (('X', lambda k: k[1] == 2), ('Y', lambda k: k[1] in (2, 3))):
        back = w.inverse_step(wk.restrict(ticked, pred))
        sv = wk.dot(start, back) / n2
        perp = wk.combine(back, start, -sv)                           # the non-return vector
        again = w.step(perp)
        norm = wk.norm2(again)
        out[name] = (sv, sum((v * v for k, v in again.items() if pred(k)), F(0)) / norm)
    (sx, tx), (sy, ty) = out['X'], out['Y']

    def count(a, b):
        return a + b - 2 * a * b
    return dict(returned=(sx, sy), not_returned=(tx, ty),
                counts=(count(sx, sy), count(tx, sy), count(sx, ty), count(tx, ty)))


# ------------------------------------------------------------------------------------------------ partners
def shared_kick_polynomial(recs, n):
    """Coefficients c_m of the yes amplitude sum_m c_m u^m of one round with n partners and the generator
    sum_l Qx_l (u = exp(-i h j/n) for the shared kick), by enumeration of the flags of all cells."""
    size = len(recs)
    fx = [F(sum(1 for r in recs if (r == 2) == bool(f)), size) for f in (0, 1)]
    hy = [F(sum(1 for r in recs if (r in (2, 3)) == bool(h)), size) for h in (0, 1)]
    coeff = [F(0)] * (n + 1)
    for f in (0, 1):
        for hs in product((0, 1), repeat=n):
            weight = fx[f]
            for h in hs:
                weight *= hy[h]
            coeff[sum(f ^ h for h in hs)] += weight
    return coeff


def enumerate_two_partners(recs):
    """The same coefficients for n = 2 by a sum over all 60^3 address triples."""
    f = np.array([r == 2 for r in recs], dtype=np.int64)
    h = np.array([r in (2, 3) for r in recs], dtype=np.int64)
    m = (f[:, None, None] ^ h[None, :, None]) + (f[:, None, None] ^ h[None, None, :])
    return [F(int((m == v).sum()), len(recs) ** 3) for v in range(3)]


def at_once_count(n):
    return F(1, 12) + F(5, 6) * (1 - F(5, 6) ** n)


def shared_phase(coeff, n, N, digits=30):
    with mp.workdps(digits):
        u = mp.exp(-1j / (mp.mpf(N) * n))
        a = sum((mp.mpf(c.numerator) / c.denominator) * u ** m for m, c in enumerate(coeff))
        return float(-N * mp.arg(a))


# ------------------------------------------------------------------------------------------------ forward only
def forward_only(w, ticks=4, port=PORT):
    state = standing(w, port)
    n2 = F(w.size)
    counts = []
    for _ in range(ticks):
        state = w.step(state)
        p = sum((v * v for k, v in state.items() if k[1] == 2), F(0)) / n2
        q = sum((v * v for k, v in state.items() if k[1] in (2, 3)), F(0)) / n2
        counts.append(p + q - 2 * p * q)
    return counts


# ------------------------------------------------------------------------------------------------ register
def register_branches(w, port=PORT):
    """Double loads with record flags: weights of the analyzer branches and equality of the amplitudes inside a
    branch. M = a I + b (z R + z^-1 R^-1) with 3 a^2 = P(X blank) P(Y written), 3 b^2 = P(X record 2) P(Y blank)."""
    ticked = w.step(standing(w, port))
    n2 = F(w.size)
    groups = {'blank': [], 'record 2': [], 'record 3': []}
    for (cfg, rec, k), v in ticked.items():
        groups['blank' if rec == 0 else f'record {rec}'].append(v)
    weight = {g: sum((v * v for v in vs), F(0)) / n2 for g, vs in groups.items()}
    equal_positive = all(len(set(vs)) == 1 and vs[0] > 0 for vs in groups.values())
    written = weight['record 2'] + weight['record 3']
    a2 = weight['blank'] * written / 3
    b2 = weight['record 2'] * weight['blank'] / 3
    b2m = weight['record 3'] * weight['blank'] / 3
    return dict(weight=weight, equal_positive=equal_positive, a2=a2, b2=b2, b2_minus=b2m,
                trace=3 * (a2 + b2 + b2m))


def native_pair_branches(w, ports):
    """Native pairs with the oriented record: 3 a^2 = P(X untouched) P(Y agrees), 3 b^2 = P(X repairs +) P(Y untouched)."""
    wt, signs = carrier.native_pair(w, ports)
    a2 = wt['untouched'] * wt['agree'] / 3
    b2 = wt['repair, load first'] * wt['untouched'] / 3
    b2m = wt['repair, load second'] * wt['untouched'] / 3
    return dict(a2=a2, b2=b2, b2_minus=b2m, trace=3 * (a2 + b2 + b2m), signs=signs)


def return_operator(recs, theta, h):
    """Float check. Register F explicit, arrays over (position, address of X, address of Y). Returns the largest
    entry of G - (1 + (exp(-i h) - 1) M^dag M), the off-diagonal part of G on the characters and M."""
    n = len(recs)
    rp = np.array([r == 2 for r in recs])
    rm = np.array([r == 3 for r in recs])
    hy = rp | rm
    kp = cmath.exp(1j * theta)

    def reflect(x):
        x = 2 * x.mean(axis=-2, keepdims=True) - x
        return 2 * x.mean(axis=-1, keepdims=True) - x

    def tick(x):
        x = reflect(x)
        out = x.copy()
        out[:, rp, :] = kp * np.roll(x[:, rp, :], 1, axis=0)
        out[:, rm, :] = kp.conjugate() * np.roll(x[:, rm, :], -1, axis=0)
        return out

    def untick(x):
        out = x.copy()
        out[:, rp, :] = kp.conjugate() * np.roll(x[:, rp, :], -1, axis=0)
        out[:, rm, :] = kp * np.roll(x[:, rm, :], 1, axis=0)
        return reflect(out)
    blank = ~(rp | rm)
    phi = (np.outer(blank / math.sqrt(blank.sum()), hy / math.sqrt(hy.sum()))
           + np.outer(rp / math.sqrt(rp.sum()), ~hy / math.sqrt((~hy).sum()))
           + np.outer(rm / math.sqrt(rm.sum()), ~hy / math.sqrt((~hy).sum()))) / math.sqrt(3)
    uu = np.full((n, n), 1.0 / n, dtype=np.complex128)
    M = np.zeros((3, 3), dtype=np.complex128)
    G = np.zeros((3, 3), dtype=np.complex128)
    for f in range(3):
        e = np.zeros(3, dtype=np.complex128)
        e[f] = 1
        t = tick(e[:, None, None] * uu)
        c = np.tensordot(t, phi, axes=([1, 2], [0, 1]))
        M[:, f] = c
        back = untick(t + (cmath.exp(-1j * h) - 1) * c[:, None, None] * phi)
        G[:, f] = back.sum(axis=(1, 2)) / n
    w3 = cmath.exp(2j * math.pi / 3)
    V = np.array([[w3 ** (-(j * k)) for k in range(3)] for j in range(3)]) / math.sqrt(3)
    Gc = V.conj().T @ G @ V
    return (float(np.abs(G - (np.eye(3) + (cmath.exp(-1j * h) - 1) * M.conj().T @ M)).max()),
            float(np.abs(Gc - np.diag(np.diag(Gc))).max()), M)


# ------------------------------------------------------------------------------------------------ execution
def executed_rounds(w, recs, N, rounds, check=5, port=PORT):
    """The round on complex128 arrays over the 60 x 60 addresses of the pair, with the projected state carried from
    round to round; the first `check` rounds are also executed through the public walk on dictionaries."""
    n = len(recs)
    xor = np.array([r == 2 for r in recs])[:, None] ^ np.array([r in (2, 3) for r in recs])[None, :]
    kick = np.where(xor, cmath.exp(-1j / N), 1.0)
    uu = np.full((n, n), 1.0 / n, dtype=np.complex128)

    def reflect(x):
        x = 2 * x.mean(axis=0, keepdims=True) - x
        return 2 * x.mean(axis=1, keepdims=True) - x
    amp = 1.0 / math.sqrt(n)
    cell = {(wk.config(port, port), 0, k): amp for k in range(n)}
    start = wk.tensor(cell, cell)
    pd = dict(start)
    psi = uu.copy()
    g = 1.0 + 0j
    phase = 0.0
    worst = count_dev = 0.0
    z = cmath.exp(-1j / N)
    for r in range(rounds):
        t = reflect(psi)
        if r < 100 or r % 1000 == 0:
            count_dev = max(count_dev, abs(float((abs(t) ** 2 * xor).sum() / (abs(t) ** 2).sum()) - 2 / 9))
        out = reflect(t * kick)
        a = out.sum() / n
        if r < check:
            ticked = wk.on_both(pd, w.step)
            kicked = {k: (v * z if (k[0][1] == 2) != (k[1][1] in (2, 3)) else v) for k, v in ticked.items()}
            back = wk.on_both(kicked, w.inverse_step)
            a_pub = sum(complex(start[k]).conjugate() * v for k, v in back.items() if k in start)
            worst = max(worst, abs(a_pub - a))
            pd = {k: (a_pub / abs(a_pub)) * v for k, v in start.items()}
        rel = a * g.conjugate()
        phase += -math.atan2(rel.imag, rel.real)
        g = a / abs(a)
        psi = g * uu
    return dict(phase_per_unit_clock_time=phase / (rounds / N), public_walk_difference=worst,
                count_deviation=count_dev)


def closed_phase(N, digits=30):
    with mp.workdps(digits):
        p = mp.mpf(2) / 9
        return float(-N * mp.arg(1 - p + p * mp.exp(-1j / mp.mpf(N))))


# ------------------------------------------------------------------------------------------------ result
def result():
    w = carrier.carrier_walk()
    recs, exact_event = records_by_address(w)
    rec = {p: complete_record(w, p) for p in (2, 10)}
    r2 = rec[2]
    nr = without_reset(w)
    polys = {n: shared_kick_polynomial(recs, n) for n in (1, 2, 3, 6)}
    means = {n: sum((m * c for m, c in enumerate(polys[n])), F(0)) / n for n in polys}
    closed = {n: [F(0)] * (n + 1) for n in polys}
    for n in polys:                                   # (11/12)(5/6 + u/6)^n + (1/12)(1/6 + 5u/6)^n
        for m in range(n + 1):
            closed[n][m] = (F(11, 12) * math.comb(n, m) * F(1, 6) ** m * F(5, 6) ** (n - m)
                            + F(1, 12) * math.comb(n, m) * F(5, 6) ** m * F(1, 6) ** (n - m))
    fwd = forward_only(w)
    reg = register_branches(w)
    nb = carrier.neighbours()
    other = min(p for p in range(12) if p != PORT and p not in nb[PORT])
    nat = native_pair_branches(w, (PORT, other))
    g_dev, g_off, M = return_operator(recs, 2 / 9, 0.01)
    run = executed_rounds(w, recs, 100, 10000)
    phases = {(n, N): shared_phase(polys[n], n, N) for n, N in ((1, 1), (1, 10), (2, 1), (3, 1), (2, 10))}
    checks = dict(
        event_maps_each_blank_double_address_to_one_basis_state=exact_event and len(recs) == 60,
        yes_amplitude_overlap_and_norm_are_2_9=all(r['overlap'] == r['norm'] == P and r['outside'] == 0
                                                   for r in rec.values()),
        test_outcome_weights_4_81_11_324_125_1296_55_1296=all(
            (r['yes_yes'], r['x_no'], r['y_no'], r['both_no']) == (F(4, 81), F(11, 324), F(125, 1296), F(55, 1296))
            for r in rec.values()),
        no_weights_sum_to_14_81=r2['x_no'] + r2['y_no'] + r2['both_no'] == F(14, 81) == P * (1 - P),
        complete_record_multiplier_is_1_minus_P_plus_P_z_over_w=multiplier_identity(),
        not_returned_cells_read_11_12_and_5_6=nr['returned'] == (F(1, 12), F(1, 6))
                                              and nr['not_returned'] == (F(11, 12), F(5, 6)),
        counts_without_reset_are_2_9_7_9_7_9_2_9=nr['counts'] == (F(2, 9), F(7, 9), F(7, 9), F(2, 9)),
        shared_kick_polynomial_is_the_product_form=polys == closed,
        two_partner_polynomial_by_enumeration_of_216000_triples=enumerate_two_partners(recs) == polys[2],
        shared_kick_generator_has_mean_2_9_for_1_2_3_6_partners=all(v == P for v in means.values()),
        partners_read_at_once_count_73_216_and_563_1296=(at_once_count(1), at_once_count(2), at_once_count(3))
                                                        == (F(2, 9), F(73, 216), F(563, 1296)),
        shared_kick_phases_at_N_1_and_10=(abs(phases[(1, 1)] - 0.205334030970278) < 1e-12
                                          and abs(phases[(1, 10)] - 0.222062099670503) < 1e-12
                                          and abs(phases[(2, 1)] - 0.215732955808) < 1e-11
                                          and abs(phases[(3, 1)] - 0.217602641007) < 1e-11
                                          and abs(phases[(2, 10)] - 0.222158623247) < 1e-11),
        forward_only_counts_on_sixty_addresses=fwd == [F(2, 9), F(3544109, 7290000), F(21968948741, 41006250000),
                                                       F(1396558119030851, 2989355625000000)],
        double_load_register_has_a2_5_108_b2_5_216_trace_5_18=(reg['equal_positive'] and reg['a2'] == F(5, 108)
                                                               and reg['b2'] == reg['b2_minus'] == F(5, 216)
                                                               and reg['trace'] == F(5, 18)),
        native_pair_register_has_a2_1_27_b2_1_54_trace_2_9=(nat['a2'] == F(1, 27) and nat['b2'] == nat['b2_minus'] == F(1, 54)
                                                            and nat['trace'] == F(2, 9)),
        both_registers_are_balanced=reg['b2'] / reg['a2'] == nat['b2'] / nat['a2'] == F(1, 2),
        mass_round_return_operator_is_1_plus_e_minus_1_M_dag_M=g_dev < 1e-13 and g_off < 1e-13,
        compression_on_sixty_addresses_is_sqrt_5_18_C=float(np.abs(
            M - math.sqrt(5 / 18) * (math.sqrt(1 / 6) * np.eye(3) + math.sqrt(1 / 12) * (
                cmath.exp(2j / 9) * np.roll(np.eye(3), 1, axis=0)
                + cmath.exp(-2j / 9) * np.roll(np.eye(3), -1, axis=0)))).max()) < 1e-14,
        ten_thousand_rounds_keep_the_phase_and_the_count=(
            abs(run['phase_per_unit_clock_time'] - closed_phase(100)) < 1e-11
            and run['public_walk_difference'] < 1e-12 and run['count_deviation'] < 1e-14))
    data = dict(
        test_outcome_weights={k: str(r2[k]) for k in ('yes_yes', 'x_no', 'y_no', 'both_no')},
        read_weights_of_a_cell_after_a_no={'X': str(nr['not_returned'][0]), 'Y': str(nr['not_returned'][1])},
        counts_in_the_four_test_states=[str(c) for c in nr['counts']],
        shared_kick_polynomials={str(n): [str(c) for c in polys[n]] for n in polys},
        shared_kick_generator_means={str(n): str(v) for n, v in means.items()},
        at_once_counts={str(n): str(at_once_count(n)) for n in (1, 2, 3, 6)},
        shared_kick_phases={f'n={n}, N={N}': mp.nstr(mp.mpf(v), 13) for (n, N), v in phases.items()},
        forward_only_counts=[str(c) for c in fwd],
        double_load_register={'a^2': str(reg['a2']), 'b^2': str(reg['b2']), 'trace': str(reg['trace'])},
        native_pair_register={'a^2': str(nat['a2']), 'b^2': str(nat['b2']), 'trace': str(nat['trace'])},
        executed_rounds={'N': 100, 'rounds': 10000,
                         'phase_per_unit_clock_time': mp.nstr(mp.mpf(run['phase_per_unit_clock_time']), 12),
                         'closed_form': repr(closed_phase(100))})
    return data, checks
