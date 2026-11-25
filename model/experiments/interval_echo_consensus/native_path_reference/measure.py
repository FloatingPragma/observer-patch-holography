"""Execute native multiplicity references, their classical projection and quantum boundary."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from math import log
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import sympy as sp

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.clock import Registers
from model.candidate_laws.interval_echo.projection import Branch
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
PORTS = (0, 1, 5)
DETECTORS = ((0, 10), (1, 8), (5, 9))
EDGES = ((0, 1), (1, 5))


class SourceTape:
    """One interior uniform variate on the requested side of either native threshold."""

    def __init__(self, bits):
        self.bits = iter(bits)

    def random(self):
        return 0.125 if next(self.bits) else 0.875


def check(condition, message):
    if not condition:
        raise ArithmeticError(message)


def push(table, function):
    result = Counter()
    for key, p in table.items():
        result[function(key)] += p
    return result


def conditional(table, event):
    total = sum(p for h, p in table.items() if event(h))
    check(total > 0, 'empty admissible event')
    return Counter({h: p / total for h, p in table.items() if event(h)}), total


def kl(p, q):
    check(all(q.get(h, 0) > 0 for h, v in p.items() if v), 'relative entropy outside support')
    return sum(float(v) * log(float(v / q[h])) for h, v in p.items() if v)


def prepare(board):
    r = Registers(board, seed=0)
    r.loads[0] = r.loads[2] = 2
    check(all(abs(p - .5) < 2e-13 for p in r.rates.values()), 'native source is not fair on both sharp inputs')
    return r


def read_marker(r):
    flags = [r.flag(a, b, {0}) for a, b in DETECTORS]
    check(sum(flags) == 1, 'one-hot physical marker record failed')
    return flags.index(1), flags


def write_marker(r, label):
    current = PORTS[read_marker(r)[0]]
    target = PORTS[label]
    if current != target:
        for a, b in r.route(current, target, {0}):
            r.swap(a, b, {0})


def execute(board, bits, depth, menu='reference'):
    check(menu in ('reference', 'biased_source', 'always_stay', 'mean_standard'), 'unknown native menu')
    r = prepare(board)
    r.rng = SourceTape(bits)
    records = []
    for stage in range(depth):
        draws = []
        for edge in EDGES:
            if menu == 'biased_source':
                # The sharp source starts each pi/3 pulse on its first rail.
                a, b, vacancy = 2, 3, 4
                r.native(a, b, {0})
                check((r.loads[a], r.loads[b], r.loads[vacancy]) == (2, 0, 0), 'unprepared biased source')
                outcome = r.rng.random() < .25
                r.loads[a], r.loads[b] = (0, 2) if outcome else (2, 0)
                r.gates['source_pi_over_three'] += 1
                bit = int(r.flag(b, vacancy, {0}))
                if bit:
                    r.swap(a, b, {0})
                    r.gates['source_resets'] += 1
            else:
                bit = r.fair(0, {0})
            draws.append(bit)
            if menu in ('reference', 'biased_source') or (menu == 'always_stay' and stage == 0):
                if bit:
                    r.swap(*edge, {0})
        if menu == 'mean_standard':
            # Two native record bits define the integer observable before readout.
            write_marker(r, sum(draws))
        observed, flags = read_marker(r)
        records.append(dict(source_bits=draws, position=observed, complete_flags=flags))
    check(r.gates['half_contacts'] + r.gates['source_pi_over_three'] == 2 * depth, 'a nominal draw was discarded')
    return tuple(q['position'] for q in records), dict(gates=dict(r.gates), records=records)


def all_native_records(board, depth, menu='reference'):
    counts = Counter()
    raw = []
    for bits in product((0, 1), repeat=2 * depth):
        record, audit = execute(board, bits, depth, menu)
        probability = (Fraction(1,4)**sum(bits) * Fraction(3,4)**(len(bits)-sum(bits))
                       if menu == 'biased_source' else Fraction(1, 4**depth))
        counts[record] += probability
        raw.append(dict(tape=list(bits), output=list(record), tape_probability=str(probability), **audit))
    return counts, raw


def cover_weights(board, depth):
    counts = Counter()
    for bits in product((0, 1), repeat=max(0, depth - 1)):
        r = prepare(board)
        r.rng = Branch(bits)
        chosen = depth
        for j in range(1, depth):
            if r.fair(0, {0}):
                chosen = j
                break
        counts[chosen] += 1
    return [Fraction(counts[j], 2**max(0, depth - 1)) for j in range(1, depth + 1)]


def weighted_cost(p, reference, weights):
    return sum(float(w) * kl(push(p, lambda h: h[:j]), push(reference, lambda h: h[:j]))
               for j, w in enumerate(weights, 1))


def cover_rank(histories, mode):
    rows = []
    depth = len(histories[0])
    projections = [(lambda h, j=j: h[:j]) for j in range(1, depth + 1)] if mode == 'prefix' else [
        (lambda h, j=j: h[j:j+1]) for j in range(depth)]
    for projection in projections:
        for label in sorted({projection(h) for h in histories}):
            rows.append([int(projection(h) == label) for h in histories])
    matrix = sp.Matrix(rows)
    return dict(rank=int(matrix.rank()), reachable_joints=len(histories),
                invisible_classical_dimension=len(histories) - int(matrix.rank()))


def symbolic_mean_projection():
    root = sp.sqrt(2)
    p = sp.Matrix([(4-root)/7, (2*root-1)/7, (4-root)/7])
    tau = sp.Matrix([sp.Rational(1, 2), sp.Rational(1, 4), sp.Rational(1, 4)])
    check(sp.simplify(sum(p)-1) == 0 and sp.simplify(p[1]+2*p[2]-1) == 0, 'projection not feasible')
    check(all(v.is_positive for v in p), 'projection not faithful')
    # Its exact ratio is c*(sqrt2)^E: this is the finite KL Pythagorean certificate.
    ratios = [sp.simplify(p[i]/tau[i]) for i in range(3)]
    check(all(sp.simplify(ratios[i] - ratios[0]*root**i) == 0 for i in range(3)), 'tilt ratio not affine in the observable')
    a = sp.Symbol('a')
    check(sp.simplify((7*a*a - 8*a + 2).subs(a, p[0])) == 0, 'optimum polynomial differs')
    rational_standard = [sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(1, 4)]
    excess = sum(q*sp.log(q/t) for q, t in zip(rational_standard, tau, strict=True)) - sum(q*sp.log(q/t) for q, t in zip(p, tau, strict=True))
    return dict(reference=[str(v) for v in tau], observable=[0, 1, 2], native_mean_standard=[str(v) for v in rational_standard],
                optimizer=[str(v) for v in p], optimizer_numerical=[float(v) for v in p],
                ratio_certificate=[str(v) for v in ratios], tilt_base='sqrt(2)',
                optimum_zero_minimal_polynomial='7*a^2-8*a+2', finite_denominator_error_bound='1/(8*k^2), with accepted tapes k<=2^N',
                native_feasible_standard_KL_excess_nats=float(excess.evalf()),
                limit='Every bounded fair-read sampler and every finite retained-event conditioning has rational masses; exact optimizer does not')


def interval_decision(address, denominator):
    """Classify a dyadic cell using the decreasing polynomial on [0, 1/2]."""
    left = 2 * address < denominator
    lower = address if left else denominator - 1 - address

    def polynomial(n):
        return 7 * n * n - 8 * n * denominator + 2 * denominator * denominator

    if polynomial(lower + 1) > 0:
        return 0 if left else 2
    if polynomial(lower) < 0:
        return 1
    return None


def sample_mean(registers, cap=None):
    """Execute the polynomial stopping rule on the public source and terminal marker."""
    if cap is not None and cap < 1:
        raise ValueError('the cutoff requires at least one source read')
    bits = []
    address = 0
    label = None
    while label is None and (cap is None or len(bits) < cap):
        bit = registers.fair(0, {0})
        bits.append(bit)
        address = 2 * address + bit
        label = interval_decision(address, 2**len(bits))
    # Both unresolved cells publish the middle label at a finite cutoff.
    write_marker(registers, 1 if label is None else label)
    observed, flags = read_marker(registers)
    check(observed == (1 if label is None else label), 'native terminal marker differs')
    return dict(tape=bits, cutoff=label is None, terminal=observed, complete_flags=flags)


def adaptive_sampler(board, cap):
    """Enumerate the live stopping prefixes and physically replay every terminal leaf."""
    if cap < 1:
        raise ValueError('the cutoff requires at least one source read')
    pending = [()]
    leaves = []
    pending_by_depth = []
    capped_law = [Fraction(0) for _ in PORTS]
    expected_reads = Fraction(0)
    for depth in range(1, cap + 1):
        denominator = 2**depth
        expected_reads += Fraction(len(pending), 2**(depth - 1))
        next_pending = []
        for parent in pending:
            for bit in (0, 1):
                tape = (*parent, bit)
                r = prepare(board)
                r.rng = Branch(tape)
                record = sample_mean(r, cap=depth)
                check(tuple(record['tape']) == tape, 'adaptive replay passed an earlier stop')
                if record['cutoff']:
                    next_pending.append(tape)
                    if depth < cap:
                        continue
                capped_law[record['terminal']] += Fraction(1, denominator)
                leaves.append(dict(**record, probability=str(Fraction(1, denominator))))
        pending = next_pending
        check(len(pending) == 2, 'two reflected threshold cells must remain unresolved')
        pending_by_depth.append(str(Fraction(len(pending), denominator)))
    check(sum(capped_law) == 1 and capped_law[0] == capped_law[2], 'capped law lost mass or the native mean')
    return dict(cutoff_reads=cap, stopping_leaves=leaves, unresolved_probability_by_depth=pending_by_depth,
                finite_output_law=[str(p) for p in capped_law], finite_expected_reads=str(expected_reads),
                finite_total_variation_upper_bound=pending_by_depth[-1], unbounded_expected_reads=3,
                controller='Public fair reads, integer prefix and polynomial comparison; unresolved cutoffs publish label one')


def quantum_audit():
    # Independently use actual 99-state generator columns, not an installed qubit.
    law = candidate.native_law
    index = [law.SRC.index((2, 0)), law.SRC.index((0, 2))]
    columns = np.eye(law.DIM)[:, index]
    h0, h1 = law.local_term_parts()
    h = h0 + candidate.KAPPA * h1
    contrast_shift = float(np.trace(columns.T @ h1 @ columns).real / 2)
    contrast_residual = float(np.linalg.norm(h1 @ columns - contrast_shift * columns))
    check(contrast_residual < 1e-12, 'contrast distinguishes the two source rails')
    plus = columns @ (np.array([1., 1j]) / np.sqrt(2))
    minus = columns @ (np.array([1., -1j]) / np.sqrt(2))
    half = expm(-.5j * np.pi * h)
    biased = expm(-1j * np.pi * h / 3)
    biased_columns = np.array([[abs((biased @ v)[i])**2 for i in index] for v in columns.T])
    source_leakage = float(np.linalg.norm(np.delete(biased @ columns, index, axis=0)))
    check(source_leakage < 1e-11, 'biased source leaked from the two sharp rails')
    check(np.max(abs(biased_columns - [[.75,.25],[.25,.75]])) < 1e-11, 'native biased source reduction differs')
    p = np.array([[abs((half @ v)[i])**2 for i in index] for v in (plus, minus)])
    computational = np.array([[abs(v[i])**2 for i in index] for v in (plus, minus)])
    h_zero = h0
    changed = expm(-.5j * np.pi * h_zero)
    control = np.array([[abs((changed @ v)[i])**2 for i in index] for v in (plus, minus)])
    check(np.max(abs(computational - .5)) < 1e-12, 'diagonal quantum witness not indistinguishable')
    check(np.max(abs(p - [[0, 1], [1, 0]])) < 1e-11, 'native coherent cover witness not distinguished')
    check(np.max(abs(control - p)) < 1e-11, 'contrast-zero control should preserve this block read')
    return dict(actual_native_pi_over_three_source_probabilities=biased_columns.tolist(),
                biased_source_off_block_norm=source_leakage, contrast_common_shift=contrast_shift,
                contrast_nonscalar_residual=contrast_residual,
                native_sharp_diagonal_probabilities=computational.tolist(),
                actual_native_half_contact_read=p.tolist(), contrast_zero_same_half_contact_read=control.tolist(),
                arbitrary_d_diagonal_trace_zero_kernel_dimension='d*(d-1)',
                statement='All sharp permutation/diagonal-read path effects are diagonal; native half-contact on the data is outside that grammar and distinguishes these states',
                changed_law_scope='Contrast removed, alphabet/routing/pulse durations unchanged: the fair-read reference does not select kappa8')


def fractions(table):
    return [dict(record=list(h), probability=str(p)) for h, p in sorted(table.items())]


def marker_departed(history):
    return history[0] != 0


def run():
    panels = []
    prior = None
    for level, depth in ((0, 1), (1, 2), (2, 3)):
        board = federation(level)
        ref, raw = all_native_records(board, depth)
        changed, raw_changed = all_native_records(board, depth, 'always_stay')
        biased, raw_biased = all_native_records(board, depth, 'biased_source')
        event = marker_departed
        optimum, admission = conditional(ref, event)
        rival, changed_admission = conditional(changed, event)
        biased_optimum, biased_admission = conditional(biased, event)
        weights = cover_weights(board, depth)
        check(sum(ref.values()) == sum(changed.values()) == sum(biased.values()) == 1, 'a source law lost mass')
        if depth == 1:
            check([ref[(i,)] for i in range(3)] == [Fraction(1, 2), Fraction(1, 4), Fraction(1, 4)],
                  'native counting reference differs from the exact mean projection input')
        check(admission == changed_admission == Fraction(1, 2), 'native acceptance differs')
        check(sum(weights) == 1 and all(w > 0 for w in weights), 'observer cover weights not positive')
        check(all(p > 0 for p in ref.values()), 'reachable reference not faithful')
        if prior:
            check(push(ref, lambda h: h[:-1]) == prior['ref'], 'native reference refinement fails')
            check(push(optimum, lambda h: h[:-1]) == prior['optimum'], 'optimizer refinement fails')
            check(weights[:-2] + [sum(weights[-2:])] == prior['weights'], 'reader weights do not coarsen')
        full_reference = Counter({tuple(r['tape']): Fraction(r['tape_probability']) for r in raw})
        accepted_tapes = {tuple(r['tape']) for r in raw if event(tuple(r['output']))}
        full_optimum, _ = conditional(full_reference, lambda t: t in accepted_tapes)
        full_biased = Counter({tuple(r['tape']): Fraction(r['tape_probability']) for r in raw_biased})
        full_control, _ = conditional(full_biased, lambda t: t in accepted_tapes)
        raw_weights = [Fraction(0) if j % 2 else weights[j//2-1] for j in range(1, 2*depth+1)]
        full_gap = weighted_cost(full_control, full_reference, raw_weights) - weighted_cost(full_optimum, full_reference, raw_weights)
        check(full_gap > .1, 'same-alphabet biased source did not reject native optimizer')
        lift_error = None
        if prior:
            coarse = Counter(prior['optimum'])
            names = list(coarse)
            amount = coarse[names[0]]/2
            coarse[names[0]] -= amount
            coarse[names[1]] += amount
            lifted = Counter({h: coarse[h[:-1]] * value / prior['ref'][h[:-1]]
                              for h, value in ref.items() if h[:-1] in coarse})
            check(push(lifted, lambda h: h[:-1]) == coarse, 'arbitrary feasible reference lift not surjective')
            check(all(event(h) for h, value in lifted.items() if value), 'lift loses visible constraint')
            lift_error = str(sum(lifted.values())-1)
        # General feasible lift checked on an asymmetric, non-Markov joint competitor.
        q = Counter(optimum)
        keys = list(q)
        if len(keys) >= 2:
            mass = q[keys[0]] / 2
            q[keys[0]] -= mass
            q[keys[1]] += mass
        identity_error = weighted_cost(q, ref, weights) - weighted_cost(optimum, ref, weights) - weighted_cost(q, optimum, weights)
        check(abs(identity_error) < 3e-13, 'weighted exact Pythagorean check fails')
        gap = weighted_cost(rival, ref, weights) - weighted_cost(optimum, ref, weights)
        check(depth == 1 or gap > .1, 'feasible changed controller did not expose selection')
        panels.append(dict(board_cells=board.cells, depth=depth, complete_tapes=len(raw),
                           source_records=raw, feasible_changed_controller_records=raw_changed,
                           biased_native_source_records=raw_biased, biased_source_admission=str(biased_admission),
                           full_source_conditional_reference=fractions(full_optimum),
                           complete_source_cover=dict(rank=len(raw), reachable_joints=len(raw), invisible_classical_dimension=0),
                           path_only_cover_on_full_source=dict(rank=len(ref), reachable_joints=len(raw), invisible_classical_dimension=len(raw)-len(ref)),
                           complete_source_diagonal_quantum_kernel_dimension=len(raw)*(len(raw)-1),
                           full_source_biased_pulse_excess_nats=full_gap,
                           quotient_biased_pulse_excess_nats=weighted_cost(biased_optimum, ref, weights)-weighted_cost(optimum, ref, weights),
                           asymmetric_feasible_lift_normalization_error=lift_error,
                           reference=fractions(ref), conditional_optimizer=fractions(optimum),
                           retained_rejection_probability=str(1-admission), observer_weights=[str(w) for w in weights],
                           weighted_reference_cost_nats=weighted_cost(optimum, ref, weights),
                           feasible_controller_excess_nats=gap,
                           asymmetric_competitor_Pythagorean_residual=identity_error,
                           prefix_cover=cover_rank(list(ref), 'prefix'),
                           individual_time_cover=cover_rank(list(ref), 'single'),
                           diagonal_quantum_tangent_kernel_dimension=len(ref)*(len(ref)-1)))
        prior = dict(ref=ref, optimum=optimum, weights=weights)
    # Independent physical source generates the rational mean-one constraint datum.
    standard, raw_standard = all_native_records(federation(0), 1, 'mean_standard')
    check(sum(h[0] * p for h, p in standard.items()) == 1, 'native mean target not one')
    check([standard[h] for h in ((0,), (1,), (2,))] == [Fraction(1,4), Fraction(1,2), Fraction(1,4)], 'mean standard map differs')
    return dict(source='Interval Echo kappa8, native sharp fair source plus pi contact permutations and echo/read/undo',
                levels=panels, native_constraint_standard=raw_standard,
                exact_mean_projection=symbolic_mean_projection(), adaptive_sampler=adaptive_sampler(federation(0), 12),
                quantum_boundary=quantum_audit(),
                scope='General native bounded-fair grammar has a multiplicity reference and conditional-support projection; universal exact mean projection and quantum joint separation do not follow')


def compare(a, b, path='native_path_reference'):
    if isinstance(b, dict):
        check(isinstance(a, dict) and a.keys() == b.keys(), path)
        for k in b:
            compare(a[k], b[k], path+'.'+k)
    elif isinstance(b, list):
        check(isinstance(a, list) and len(a) == len(b), path)
        for i, (c, d) in enumerate(zip(a, b, strict=True)):
            compare(c, d, path+f'[{i}]')
    elif isinstance(b, float):
        check(np.isfinite(a) and np.isclose(a, b, rtol=1e-10, atol=3e-12), path)
    else:
        check(type(a) is type(b) and a == b, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    encoded = json.dumps(run(), indent=2, allow_nan=False)+'\n'
    if args.output:
        with args.output.open('x') as stream:
            stream.write(encoded)
    else:
        compare(json.loads(encoded), json.loads((HERE/'results.json').read_text()))
        print('Verified native path multiplicities, classical projection and exact sampling/cover boundaries.')


if __name__ == '__main__':
    main()
