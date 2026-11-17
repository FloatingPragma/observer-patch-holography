"""Exact finite predictions distinguishing current and contrast instruments.

Both candidates implement the existing native endpoint repair on 78 load
states. They differ only in the coherent grouping of their recorded Kraus
branches. Preparation and the final source interference measurement are
supplied operations; this is a conditional apparatus prediction, not an
empirical verdict or a derivation of those services.

The protocol orders one repair and one binary readout. It assigns no physical
duration to either operation. All records and both readout outcomes survive.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache

from model.observer import native

CANDIDATES = ("C", "E")
SETTINGS = ("phase_0", "phase_pi", "phase_pi_over_2", "diagonal_control")
STATES = native.states()
INDEX = {state: i for i, state in enumerate(STATES)}
ZERO = (Q(0), Q(0))


def _occupation(*ports):
    state = [0] * 12
    for port in ports:
        state[port] += 1
    return tuple(state)


X = INDEX[_occupation(2, 3)]
Z = INDEX[_occupation(1, 2)]


def _candidate(candidate):
    if candidate not in CANDIDATES:
        raise ValueError("Candidate must be C or E")


def _add(target, key, value, scale=Q(1)):
    previous = target.get(key, ZERO)
    updated = tuple(previous[i] + scale * value[i] for i in (0, 1))
    if updated == ZERO:
        target.pop(key, None)
    else:
        target[key] = updated


def prepared_density(setting):
    """Sparse density in exact Gaussian rationals (real, imaginary).

    Pure settings prepare (|X> + exp(i theta)|Z>)/sqrt(2). The diagonal
    control is their incoherent equal mixture. These share the same load
    probabilities, independently of the selected channel.
    """
    if setting not in SETTINGS:
        raise ValueError("Unknown fixed preparation")
    result = {(X, X): (Q(1, 2), Q(0)), (Z, Z): (Q(1, 2), Q(0))}
    if setting == "diagonal_control":
        return result
    phase = {"phase_0": (Q(1, 2), Q(0)),
             "phase_pi": (Q(-1, 2), Q(0)),
             "phase_pi_over_2": (Q(0), Q(-1, 2))}[setting]
    result[X, Z] = phase
    result[Z, X] = (phase[0], -phase[1])
    return result


@lru_cache(maxsize=2)
def _images(candidate):
    _candidate(candidate)
    result = []
    for state in STATES:
        images = []
        for coin in (0, 1):
            output, current_record = native.native(state, native.E, coin)
            record = current_record if candidate == "C" else (
                coin, state[native.E[0]] - state[native.E[1]])
            images.append((record, INDEX[output]))
        result.append(tuple(images))
    return tuple(result)


def kraus_support(candidate):
    """Every nonzero Kraus entry as (input, output), of amplitude 1/sqrt(2)."""
    grouped = defaultdict(list)
    for source, images in enumerate(_images(candidate)):
        for record, target in images:
            grouped[record].append((source, target))
    return {record: tuple(entries) for record, entries in sorted(grouped.items())}


def cptp_certificate(candidate):
    """Check the full Gram identity, including arbitrary input coherences.

    Each input column has two entries of squared amplitude 1/2. Distinct
    columns have disjoint (record, output) support. Thus W*W=I exactly;
    recording W's label defines a completely positive, trace-preserving
    instrument on every density operator, including entangled inputs.
    """
    images = _images(candidate)
    seen = set()
    for column in images:
        if len(column) != 2 or len(set(column)) != 2:
            raise AssertionError("Invalid Kraus column norm")
        if seen.intersection(column):
            raise AssertionError("Distinct Kraus columns overlap")
        seen.update(column)
    return {"candidate": candidate, "input_dimension": len(STATES),
            "record_count": len(kraus_support(candidate)),
            "nonzero_kraus_entries": len(seen), "gram_identity_exact": True}


def instrument_density(candidate, density):
    """Apply the recorded instrument to a sparse operator exactly.

    The caller supplies Gaussian-rational entries on the 78-state carrier.
    This is a linear operation, so it also accepts off-diagonal matrix units
    used to check a channel. Each returned record block is unnormalized.
    """
    images = _images(candidate)
    blocks = {record: {} for record in kraus_support(candidate)}
    for (left, right), value in density.items():
        if (type(left) is not int or type(right) is not int
                or not 0 <= left < len(STATES) or not 0 <= right < len(STATES)
                or not isinstance(value, tuple) or len(value) != 2
                or any(not isinstance(part, Q) for part in value)):
            raise ValueError("Expected carrier indices and Gaussian-rational entries")
        for coin in (0, 1):
            left_record, a = images[left][coin]
            right_record, b = images[right][coin]
            if left_record == right_record:
                _add(blocks[left_record], (a, b), value, Q(1, 2))
    return blocks


def trace(density):
    real = sum((value[0] for (a, b), value in density.items() if a == b), Q())
    imag = sum((value[1] for (a, b), value in density.items() if a == b), Q())
    if imag:
        raise ValueError("Trace is not real")
    return real


def readout_density(density, outcome):
    """P rho P for P=|X+Z><X+Z|/2 or its orthogonal complement.

    The projection acts only on the working source, separately inside each
    old-record block. It neither mixes nor erases committed record values.
    """
    if type(outcome) is not int or outcome not in (-1, 1):
        raise ValueError("Readout outcome must be -1 or +1")

    def column(source):
        plus = {X: Q(1, 2), Z: Q(1, 2)} if source in (X, Z) else {}
        if outcome == 1:
            return plus
        result = {target: -coefficient for target, coefficient in plus.items()}
        result[source] = result.get(source, Q()) + 1
        return result

    output = {}
    for (left, right), value in density.items():
        for a, coefficient_a in column(left).items():
            for b, coefficient_b in column(right).items():
                _add(output, (a, b), value, coefficient_a * coefficient_b)
    return output


def _record(candidate, record):
    coin, value = record
    if candidate == "C":
        return {"coin": coin, "current": value}
    current = (value - (2 * coin - 1) * (abs(value) % 2)) // 2
    return {"coin": coin, "contrast": value, "current": current}


def event_outputs(candidate, setting):
    """All old-record/readout events, including zero-probability events.

    Post-readout matrices are unnormalized: their traces equal event
    probabilities. Fractions serialize as strings without roundoff.
    E's contrast label stays present even though both models expose current.
    """
    events = []
    for record, block in instrument_density(candidate, prepared_density(setting)).items():
        for outcome in (1, -1):
            output = readout_density(block, outcome)
            probability = trace(output)
            if probability < 0:
                raise AssertionError("Negative event probability")
            events.append({"record": _record(candidate, record), "outcome": outcome,
                           "probability": str(probability),
                           "source_density": [
                               {"row": a, "column": b, "real": str(value[0]),
                                "imaginary": str(value[1])}
                               for (a, b), value in sorted(output.items())]})
    if sum((Q(event["probability"]) for event in events), Q()) != 1:
        raise AssertionError("Incomplete readout instrument")
    return events


def prediction(candidate, setting):
    return sum((Q(event["probability"]) for event in event_outputs(candidate, setting)
                if event["outcome"] == 1), Q())


def public_distribution(candidate, setting):
    """Common (coin, current, readout) probabilities, without postselection.

    E's additional recorded contrast is summed for this reporting interface;
    no coherent amplitudes are combined and its decoherence is retained.
    Every common bin is present, including zero-probability bins.
    """
    result = {(coin, current, outcome): Q()
              for coin in (0, 1) for current in (-1, 0, 1) for outcome in (-1, 1)}
    for event in event_outputs(candidate, setting):
        record = event["record"]
        result[record["coin"], record["current"], event["outcome"]] += Q(event["probability"])
    return result


def prediction_table():
    return {candidate: {setting: str(prediction(candidate, setting)) for setting in SETTINGS}
            for candidate in CANDIDATES}
