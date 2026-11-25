"""Send an independent load-one tag across the clock ball while both residents move."""
from fractions import Fraction
from itertools import combinations
import json
import math

import numpy as np
from scipy.linalg import eigh, expm
from scipy.sparse import csr_matrix
from scipy.stats import beta

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.clock import Registers, _dimer_pulse
from model.executable_model import federation
from .. import clock_ball, measure as regional

SHOTS = 2**26
OPTICAL_SHOTS = 2**34
ZERO = ((0, 1),)
ONE = ((0, 2),)


def certify_acquisition(basis, rim, outside):
    """Check native ternary pointer flags and all sharp endpoint transpositions of the shuffle."""
    law = candidate.native_law
    echo = law.echo()
    pointer = np.arange(law.DIM) % len(law.SRC) == law.SRC.index((1, 1))
    h0, h1 = law.local_term_parts()
    swap = expm(-1j*np.pi*(h0+candidate.KAPPA*h1))
    identity = np.eye(law.DIM)
    error = 0.
    for i, (terminal, standard) in enumerate(law.SRC):
        column = echo[:, i]
        error = max(error, np.linalg.norm(column[pointer != (terminal + standard == 2)]),
                    np.linalg.norm(echo @ column - identity[:, i]),
                    np.linalg.norm(swap[:, i] - identity[:, law.SRC.index((standard, terminal))]))
    assert error < 1e-11
    reg = Registers(federation(0), seed=8477)
    ports = sorted({tag for tag, _ in basis})
    inner = sorted(set(ports) - set(rim) - {outside})
    _, _, border, _ = regional.apparatus(reg.car, [0])
    gateway = min(border)
    for i, j in reg.route(outside, gateway, {0}):
        reg.native(i, j, {0})
    reg.native(gateway, border[gateway], {0, border[gateway]//12})
    count = 0
    for region in (inner, rim):
        for a, b in combinations(region, 2):
            route = reg.route(a, b, {0})
            check = np.arange(len(ports))
            for i, j in route:
                reg.native(i, j, {0})
                check[i], check[j] = check[j], check[i]
            target = np.arange(len(ports))
            target[a], target[b] = target[b], target[a]
            np.testing.assert_array_equal(check, target)
            count += 1
    return dict(verified_routed_transpositions=count, ternary_pointer_error=float(error),
                sharp_swap_error=reg.error)


def material_hand(time, rng):
    """Two analyzer settings on independent copies of the registered, isolated dimer."""
    probabilities = [abs(_dimer_pulse(time+delay)[0, 0])**2 for delay in (0, 1.5*np.pi)]
    counts = rng.binomial(SHOTS, probabilities)
    readings = 2*counts/SHOTS-1
    return dict(analyzer_counts=counts.tolist(), phase_radians=float(np.arctan2(readings[1], readings[0])),
                theoretical_pointer_rates=probabilities)


def determinant(x, y):
    return (x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0])


def calibrate_from_counts(matrix, reverse, arms, rng):
    """Simultaneous bounds using independently counted class standards and held banks."""
    train_counts = rng.binomial(OPTICAL_SHOTS, np.vstack([matrix[1:], reverse]))
    rates = train_counts/OPTICAL_SHOTS
    measured = np.vstack([np.ones(3), rates[:2]])
    charge = np.linalg.solve(measured.T, np.log([1., 5., 10.]))
    outward = np.linalg.solve(measured.T, rates[2])
    # Nine pure-class standards, four held inward banks and two held outward checks.
    eps = math.sqrt(math.log(2*15/.01)/(2*OPTICAL_SHOTS))
    estimates = {}
    for name, arm in arms.items():
        observed = np.array(arm["exterior_hits"])/OPTICAL_SHOTS
        read = np.r_[1., observed[:2]]
        estimates[name] = dict(charge_nats=float(charge @ read),
            predicted_outward_rate=float(outward @ read), observed_outward_rate=float(observed[2]))
    loss = estimates["spread"]["charge_nats"] - estimates["split"]["charge_nats"]
    radius = 4*eps*sum(abs(charge[1:]))
    return dict(standard_counts=train_counts.tolist(), shots_per_bank=OPTICAL_SHOTS,
        arms=estimates, charge_change_nats=loss, simultaneous_frequency_error=eps,
        simultaneous_99_charge_interval=[loss-radius, loss+radius],
        simultaneous_99_outward_discrepancy_bound=float(2*eps*(1+sum(abs(outward[1:])))),
        least_singular_value=float(np.linalg.svd(measured, compute_uv=False)[-1]),
        matrix_error_bound=math.sqrt(6)*eps)


def optical_readout(tables, leading, order4, rim, rng):
    """Infer the three unknown class weights from two exterior arrival frequencies."""
    groups = [[w for w in leading if len(set(w).intersection(rim)) == k] for k in range(3)]
    def coefficients(values, scale):
        return [Fraction(sum(values[w][0] for w in group), scale*len(group)) for group in groups]
    a = coefficients(leading, 4)
    b = coefficients(order4, 192)
    def responses(table, channel, scale):
        return np.array([sum(table[w][channel] for w in group)/(scale*len(group)) for group in groups])
    # Time rows are fixed by the regional exposure before a stock is measured.
    mid, late = [responses(x, 0, 1) for x in tables[1:]]
    reverse = responses(tables[0], 1, 9)
    matrix = np.vstack([np.ones(3), mid, late])
    costs = np.log([1., 5., 10.])
    decoder = np.linalg.solve(matrix.T, costs)
    reverse_decoder = np.linalg.solve(matrix.T, reverse)
    eps = math.sqrt(math.log(2*6/.01)/(2*OPTICAL_SHOTS))
    arms = {}
    for name, stock in (("split", clock_ball.SPLIT), ("spread", clock_ball.SPREAD)):
        actual = np.array([sum(w in group for w in stock)/len(stock) for group in groups])
        rates = np.array([mid @ actual, late @ actual, reverse @ actual])
        counts = rng.binomial(OPTICAL_SHOTS, rates)
        read = np.r_[1., counts[:2]/OPTICAL_SHOTS]
        arms[name] = dict(exterior_hits=counts.tolist(), expected_rates=rates.tolist(),
            decoded_distribution=np.linalg.solve(matrix, read).tolist(),
            decoded_charge_nats=float(decoder @ read),
            predicted_outward_rate=float(reverse_decoder @ read),
            observed_outward_rate=float(counts[2]/OPTICAL_SHOTS))
    loss = arms["spread"]["decoded_charge_nats"]-arms["split"]["decoded_charge_nats"]
    radius = 2*sum(abs(decoder[1:]))*eps
    return dict(counted_calibration=calibrate_from_counts(matrix, reverse, arms, rng),
        class_sizes=[len(group) for group in groups],
        times=[regional.EXPOSURE, 2*regional.EXPOSURE], held_outward_time=regional.EXPOSURE/2,
        classes=[0, 1, 2], class_calibration_matrix=matrix.tolist(),
        responses_outward=reverse.tolist(), class_t2=[str(i) for i in a], class_t4=[str(i) for i in b],
        t6_determinant=str(12*determinant(a,b)), condition_number=float(np.linalg.cond(matrix)),
        nats_per_record=decoder.tolist(), shots_per_bank=OPTICAL_SHOTS, arms=arms,
        charge_change_nats=loss, simultaneous_99_charge_interval=[loss-radius,loss+radius],
        simultaneous_99_outward_discrepancy_bound=float((1+sum(abs(reverse_decoder[1:])))*eps))


def source(contrast):
    """Twice the exact blank-memory native generator on one tag and two double loads."""
    car = federation(0)
    ports, contacts, _, words = regional.apparatus(car, [0])
    chains = clock_ball.shortest_chains(ports, contacts, clock_ball.CENTER)
    rim = tuple(p for p in ports if len(chains[p]) == clock_ball.RADIUS)
    outside, = (p for p in ports if len(chains[p]) > clock_ball.RADIUS)
    basis = [(tag, pair) for pair in words for tag in ports if tag not in pair]
    index = {event: i for i, event in enumerate(basis)}
    local = candidate.native_law.seam_term(contrast)
    symbols = candidate.native_law.SRC
    transitions = {}
    for i, (a, b) in enumerate(symbols):
        exact = np.zeros(len(symbols))
        exact[i] = 2 * contrast * (a-b)**2
        if a != b:
            exact[i] += 1
            exact[symbols.index((b, a))] -= 1
        np.testing.assert_array_equal(2*local[:, i], exact)
        transitions[a, b] = int(exact[i])
    twice = np.zeros((len(basis), len(basis)), dtype=np.int64)
    for j, (tag, pair) in enumerate(basis):
        values = dict.fromkeys(ports, 0)
        values[tag] = 1
        for p in pair:
            values[p] = 2
        twice[j, j] = 2*contrast*sum(values[p]**2 * (6-int(car.degree[p])) for p in ports)
        for a, b in contacts:
            left, right = values[a], values[b]
            twice[j, j] += transitions[left, right]
            if left != right:
                exchange = {a: b, b: a}
                moved = exchange.get(tag, tag), tuple(sorted(exchange.get(p, p) for p in pair))
                twice[index[moved], j] -= 1
    np.testing.assert_array_equal(twice, twice.T)
    return twice, basis, rim, outside


def meters(kernel, basis, outside):
    """Unconditional joint hit flags. The reverse launch samples nine vacant sites."""
    exterior = np.array([tag == outside for tag, _ in basis])
    eligible = np.array([outside not in pair for _, pair in basis])
    inside = eligible & ~exterior
    rows = {}
    for w in sorted({w for _, w in basis if outside not in w}):
        stock = np.array([pair == w for _, pair in basis])
        src, = np.flatnonzero(exterior & stock)
        targets = np.flatnonzero(inside & stock)
        assert len(targets) == 9
        # For integer Taylor coefficients summing outgoing flags avoids any rounding.
        rows[w] = (kernel[inside, src].sum().item(),
                   kernel[np.ix_(exterior, targets)].sum().item(),
                   kernel[~eligible, src].sum().item())
    return rows


def taylor(twice):
    k = csr_matrix(twice)
    k2 = k @ k
    k3 = k2 @ k
    # Only off-diagonal entries feed the three disjoint outcome comparisons.
    return (k.multiply(k).toarray(),
            (3*k2.multiply(k2) - 4*k.multiply(k3)).toarray())


def comparison(scores):
    def mean(stock, i):
        return sum(scores[w][i] for w in stock)/len(stock)
    return [[mean(b,i)-mean(a,i) for i in range(2)]
            for a,b in ((ZERO,ONE),(clock_ball.SPLIT,clock_ball.SPREAD))]


def prediction_intervals(banks):
    """Two-sided binomial intervals, simultaneous over four arms and both directions."""
    alpha = .01/8
    rows = {}
    for name, bank in banks.items():
        pairs = []
        for k in (bank["inward_counts"][0], bank["outward_count"]):
            lo = beta.ppf(alpha/2, k, SHOTS-k+1) if k else 0.
            hi = beta.ppf(1-alpha/2, k+1, SHOTS-k) if k < SHOTS else 1.
            pairs.append((lo, hi))
        rows[name] = pairs
    def contrast(after, before, i):
        a, b = rows[after][i], rows[before][i]
        return [a[0]-b[1], a[1]-b[0]]
    intervals = []
    for i in range(2):
        train = contrast("rim_standard", "interior_standard", i)
        held = contrast("spread", "split", i)
        if train[0]*train[1] <= 0:
            intervals.append(None)
        else:
            endpoints = [math.log(5)*x/y for x in held for y in train]
            intervals.append([min(endpoints), max(endpoints)])
    return intervals


def record(contrast, seed):
    twice, basis, rim, outside = source(contrast)
    vals, vecs = eigh(twice / 2)
    first, fourth = taylor(twice)
    leading = meters(first, basis, outside)
    order4 = meters(fourth, basis, outside)
    exact = [[str(Fraction(int(sum(order4[w][i] for w in b)), len(b)*192) -
                  Fraction(int(sum(order4[w][i] for w in a)), len(a)*192)) for i in range(2)]
             for a,b in ((ZERO, ONE), (clock_ball.SPLIT, clock_ball.SPREAD))]
    rng = np.random.default_rng(seed)
    exposures = []
    tables = []
    for time in (regional.EXPOSURE/2, regional.EXPOSURE, regional.EXPOSURE*2):
        u = (vecs * np.exp(-1j*time*vals)) @ vecs.T
        probabilities = abs(u)**2
        np.testing.assert_allclose(probabilities.sum(axis=0), 1, atol=2e-12, rtol=0)
        scores = meters(probabilities, basis, outside)
        tables.append(scores)
        banks = {}
        stocks = {"interior_standard":ZERO, "rim_standard":ONE,
                  "split":clock_ball.SPLIT, "spread":clock_ball.SPREAD}
        for name, stock in stocks.items():
            f, rev_scaled, displaced = np.mean([scores[w] for w in stock], axis=0)
            # Includes all starts. The tag outside implies both residents are inside.
            law = [float(f), float(displaced), float(1-f-displaced)]
            reverse = rev_scaled/9
            banks[name] = dict(words=stock, inward_categories=law, outward_hit=float(reverse),
                inward_counts=rng.multinomial(SHOTS, law).tolist(),
                outward_count=int(rng.binomial(SHOTS, reverse)))
        calibration, held = comparison(scores)
        predictions = [math.log(5)*v/u if abs(u)>1e-12 else None
                       for u,v in zip(calibration, held, strict=True)]
        conservation = np.mean(list(scores.values()), axis=0)
        np.testing.assert_allclose(conservation[0], conservation[1], rtol=0, atol=2e-12)
        exposures.append(dict(exposure=time, banks=banks,
            sampled_99_held_nats=prediction_intervals(banks),
            inward_and_scaled_outward_at_uniform_reference=conservation[:2].tolist(),
            calibration_change=calibration, held_change=held, held_predicted_nats=predictions,
            two_rim_words=[dict(word=w, inward=scores[w][0], outward=scores[w][1]/9)
                           for w in combinations(rim, 2)]))
    optical = optical_readout(tables, leading, order4, rim, rng) if contrast == candidate.KAPPA else None
    return dict(contrast=contrast, acquisition=certify_acquisition(basis, rim, outside),
        randomized_transmission=optical, joint_records=len(basis), shell=rim, exterior=outside,
        inward_leading_coefficients=sorted({str(Fraction(s[0],4)) for s in leading.values()}),
        scaled_outward_leading_coefficients=sorted({str(Fraction(s[1],4)) for s in leading.values()}),
        exact_t4_calibration_and_held=exact, exposures=exposures)


def measure():
    rng = np.random.default_rng(20716)
    return dict(shots_per_bank=SHOTS,
                independent_hand=[material_hand(t, rng) for t in
                                  (regional.EXPOSURE/2, regional.EXPOSURE, 2*regional.EXPOSURE)],
                true_held_completion_nats=math.log(2/5)/2,
                active=record(candidate.KAPPA, 36997), exchange_control=record(0, 41219))


if __name__ == "__main__":
    print(json.dumps(measure(), indent=2, allow_nan=False))
