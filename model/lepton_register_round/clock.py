"""The unit of the declared clock (section 4): the edge anchor on the public root-SWAP and the attempt count of the
public seam program.

The declared clock is H = sum_x (I - Re U_x) + sum_edges (pi/4) dGamma(I - S_xy) with evolution exp(-i t H)
(`model/coherent_defect_wave/energy/PROOF.md`). The tick's edge gate is E = a I + b S with a = (1 - i)/2 and
b = (1 + i)/2 (`model/coherent_defect_wave/spatial/family.py`).
"""
from fractions import Fraction as F
import random

import sympy as s

from model.coherent_defect_wave.model import JointObserver
from model.coherent_defect_wave.spatial.family import A as GATE_A, B as GATE_B, STAR, PairState, dressed_swap, exchange
from model.coherent_seam_repair.program import Observer, State

TAUS = (s.Integer(1), s.Integer(5), s.Integer(-3), s.Integer(9), s.Integer(-7), s.Integer(-1), s.Integer(3),
        s.Integer(2), s.Integer(0), s.Integer(4), s.Rational(1, 2), s.Rational(4, 5), s.pi / 2, 2 / s.pi)


def clock_coefficients(tau):
    """exp(-i tau (pi/4)(I - S)) = c0 I + c1 S for a Hermitian involution S."""
    phase = s.exp(-s.I * s.pi * tau / 2)
    return s.expand(s.expand_complex((1 + phase) / 2)), s.expand(s.expand_complex((1 - phase) / 2))


def pair_state():
    """A two-cell state on actual link words, with frame signs on both cells."""
    words = JointObserver().link_words
    upper = [g for g in words if g > (g ^ STAR)]
    lower = [g for g in words if g < (g ^ STAR)]
    return PairState({((upper[0], 5, 3, 0, 0), (lower[0], 6, 0, 1, 8)): 1,
                      ((upper[1], 5, 4, 1, 59), (lower[1], 5, 3, 0, 4)): s.I,
                      ((lower[2], 0, 0, 0, 7), (upper[2], 5, 4, 1, 18)): 2}, 6)


def clock_edge(state, tau):
    c0, c1 = clock_coefficients(tau)
    src, swp = state.as_dict(), dressed_swap(state).as_dict()
    return PairState({k: s.expand(c0 * src.get(k, 0) + c1 * swp.get(k, 0)) for k in src.keys() | swp.keys()},
                     state.norm_denominator)


def edge_anchor():
    pair = pair_state()
    gate, gate_inverse, swapped = exchange(pair), exchange(pair, inverse=True), dressed_swap(pair)
    rows = []
    for tau in TAUS:
        st = clock_edge(pair, tau)
        integer = tau.is_Integer
        rows.append(dict(tau=str(tau), gate=st == gate, inverse_gate=st == gate_inverse, swap=st == swapped,
                         identity=st == pair, in_1_plus_4Z=bool(integer and int(tau) % 4 == 1)))
    n = s.Symbol('n', integer=True)
    family = s.simplify(s.exp(-s.I * s.pi * (1 + 4 * n) / 2) + s.I) == 0
    principal = [t for t in range(-23, 26, 4) if abs(t) < 2]
    return dict(rows=rows, family=family, principal=principal, swap_nontrivial=swapped != pair,
                involution=dressed_swap(swapped) == pair, gate_squares_to_swap=exchange(gate) == swapped,
                gate_coefficients=[str(GATE_A), str(GATE_B)], clock_at_one=[str(x) for x in clock_coefficients(1)])


def laplacian(obs):
    lap = [[0] * 12 for _ in range(12)]
    for a, b in obs.seams:
        lap[a][a] += 1
        lap[b][b] += 1
        lap[a][b] -= 1
        lap[b][a] -= 1
    return lap


def attempt_count(samples=6, seed=11):
    obs = Observer()
    D = obs.address_count
    lap = laplacian(obs)

    def inner(x, y):
        return sum(a * b for a, b in zip(x.amplitudes, y.amplitudes)) / x.norm_denominator

    def record_one(y):
        return sum(y.amplitudes[obs.index(p, 1, k)] ** 2 for p in range(12) for k in range(D)) / y.norm_denominator

    def consensus(vec):
        amps = [F(0)] * obs.dimension
        for port, c in enumerate(vec):
            if c:
                for i, x in enumerate(obs.embedded_port(port).amplitudes):
                    amps[i] += F(c) * x
        return State(tuple(amps), D)

    ports = [obs.embedded_port(p) for p in range(12)]
    steps = [obs.step(x) for x in ports]
    identity = {(1 - inner(x, y) / inner(x, x), record_one(y) / inner(x, x)) for x, y in zip(ports, steps)}
    block = [[inner(ports[p], steps[q]) for q in range(12)] for p in range(12)]
    block_ok = all(block[p][q] == (1 if p == q else 0) - F(lap[p][q], 60) for p in range(12) for q in range(12))
    invariant = True
    for p, x in enumerate(ports):
        y, yi = steps[p], obs.inverse_step(x)
        lhs = [a - (b + c) / 2 for a, b, c in zip(x.amplitudes, y.amplitudes, yi.amplitudes)]
        rhs = consensus([F(lap[q][p], 60) for q in range(12)])
        invariant &= lhs == list(rhs.amplitudes) and x.norm_denominator == rhs.norm_denominator
    rng = random.Random(seed)
    random_ok = True
    for _ in range(samples):
        x = consensus([F(rng.randint(-9, 9), rng.randint(1, 7)) for _ in range(12)])
        y = obs.step(x)
        random_ok &= (1 - inner(x, y) / inner(x, x)) == record_one(y) / inner(x, x)
    # off the consensus sector: a blank record with a non-uniform address
    amps = [F(0)] * obs.dimension
    for p in range(12):
        for k in range(D):
            amps[obs.index(p, 0, k)] = F(rng.randint(-3, 3))
    x = State(tuple(amps))
    y = obs.step(x)
    off = (1 - inner(x, y) / inner(x, x), record_one(y) / inner(x, x))
    eig = s.Matrix(lap).eigenvals()
    return dict(identity=identity, block_ok=block_ok, invariant=invariant, random_ok=random_ok, off=off,
                spectrum={str(k): int(v) for k, v in eig.items()})


def result():
    edge = edge_anchor()
    att = attempt_count()
    rows = edge['rows']
    checks = dict(
        swap_is_a_nontrivial_involution_and_gate_squares_to_it=(edge['swap_nontrivial'] and edge['involution']
                                                                and edge['gate_squares_to_swap']),
        clock_at_unit_time_has_the_gate_coefficients=edge['clock_at_one'] == edge['gate_coefficients'],
        clock_equals_the_edge_gate_exactly_on_1_plus_4Z=all(r['gate'] == r['in_1_plus_4Z'] for r in rows),
        clock_gives_inverse_gate_and_swap_on_minus_1_and_2_plus_4Z=all(
            r['inverse_gate'] == (r['tau'] in ('-1', '3')) and r['swap'] == (r['tau'] in ('2',)) for r in rows),
        exp_minus_i_pi_1_plus_4n_over_2_is_minus_i=edge['family'],
        principal_branch_member_is_tau_1=edge['principal'] == [1],
        one_minus_overlap_is_record_one_weight_1_12_on_every_port=att['identity'] == {(F(1, 12), F(1, 12))},
        identity_on_random_consensus_inputs=att['random_ok'],
        identity_needs_the_consensus_address=att['off'][0] != att['off'][1],
        consensus_block_is_I_minus_L_over_60=att['block_ok'],
        onsite_term_keeps_the_consensus_sector=att['invariant'])
    data = dict(
        gate_coefficients=edge['gate_coefficients'], clock_coefficients_at_tau_1=edge['clock_at_one'],
        edge_rows=rows, principal_branch_members_of_1_plus_4Z=edge['principal'],
        port_laplacian_spectrum=att['spectrum'],
        one_minus_overlap_and_record_one_weight=[[str(a), str(b)] for a, b in sorted(att['identity'])],
        blank_record_with_non_uniform_address=[str(x) for x in att['off']])
    return data, checks
