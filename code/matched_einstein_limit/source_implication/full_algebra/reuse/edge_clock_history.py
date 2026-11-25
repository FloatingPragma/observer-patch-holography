"""Reconstruct a whole head from one directed radar leg pair and one probe charge."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np

from . import clock_history as clock, edge_clock as edge, measure

EVENTS = tuple((a, t) for t in range(2, 13) for a in range(1, t))
ARCHIVED_ROUTE = (2, 1)


def block_norm(matrix):
    return edge.upper_sqrt(Fraction(int((matrix*matrix).sum()), clock.F**2))


def calibrated_events():
    board = clock.station.iso.pairing.formation.federation(1)
    ports = clock.station.iso.pairing.formation.region(board)
    addresses = np.array([4*p+i for p in ports for i in range(4)])
    weights = clock.native.prior.weights()
    rows = np.zeros((2, 12, len(EVENTS)+11, 31), dtype=object)
    offsets = np.zeros(rows.shape[:3], dtype=object)
    errors = np.full(rows.shape[:3], Fraction(0), dtype=object)
    rounding = Fraction(8, clock.F)
    sample = Fraction(1, 10**19) + rounding
    first_errors = [Fraction(0)]
    for t in range(1, 12):
        first_errors.append(sample + sum(sample + (1+sample)*first_errors[a] + rounding
                                         for a in range(1, t)))
    for qi, charge in enumerate((-1, 1)):
        for parity in (0, 1):
            h = clock.native.generator(board, ports, charge, bool(parity))
            u = clock.native.kernel_panels(h, addresses, weights)
            hits = {r: clock.first_hits(u, *r) for r in measure.ROUTES}
            norms = {r: [block_norm(x) for x in v] for r, v in hits.items()}
            for i, route in enumerate(measure.ROUTES):
                back = route[::-1]
                arrivals = []
                for j, (a, t) in enumerate(EVENTS):
                    b = t-a
                    amp = clock.native.prior.multiply(hits[back][b], hits[route][a])
                    ap_error = (norms[back][b]*first_errors[a] + norms[route][a]*first_errors[b]
                                + first_errors[a]*first_errors[b] + rounding)
                    arrivals.append((amp, ap_error))
                arrivals.extend((hits[route][a], first_errors[a]) for a in range(1, 12))
                for j, (amp, ap_error) in enumerate(arrivals):
                    error = ap_error*(2*block_norm(amp)+ap_error) + rounding
                    errors[qi, i, j] = max(errors[qi, i, j], error)
                    effect = clock.native.prior.multiply(clock.native.prior.adjoint(amp), amp)
                    moments, trace = edge.coordinates(effect)
                    rows[qi, i, j, 15*parity:15*(parity+1)] = 2*moments
                    rows[qi, i, j, -1] += (2 if parity else -2)*trace
                    offsets[qi, i, j] += trace
    return rows, offsets, errors


def certificate(rows, errors, scale_rows=False):
    """Exact residual of a candidate left inverse including event-specific error."""
    parameters = rows.shape[1]
    numerical = np.asarray(rows, dtype=float) / edge.DEN
    scales = 1/np.asarray(errors, dtype=float) if scale_rows else np.ones(len(rows))
    pinv = np.linalg.pinv(numerical*scales[:, None], rcond=1e-17)*scales[None, :]
    unit = 2**40
    proposal = np.array([[round(float(x)*unit) for x in row] for row in pinv], dtype=object)
    norm = edge.upper_sqrt(Fraction(int((proposal*proposal).sum()), unit**2))
    residual = (edge.DEN*unit)*np.eye(parameters, dtype=object) - proposal@rows
    center_error = edge.upper_sqrt(Fraction(int((residual*residual).sum()), (edge.DEN*unit)**2))
    # Each coefficient of the true row differs by at most twice the effect error.
    # Bound the product of the reconstruction with that error one event at a time.
    columns = proposal.T
    transfer = 2*edge.upper_sqrt(Fraction(parameters))*sum(
        err*edge.upper_sqrt(Fraction(int(col@col), unit**2))
        for err, col in zip(errors, columns, strict=True))
    defect = center_error + transfer
    if defect >= 1 and not scale_rows:
        return certificate(rows, errors, scale_rows=True)
    out = dict(left_inverse_defect_upper=edge.upper(defect))
    if defect < 1:
        out.update(rank=parameters, smallest_singular_value_lower=edge.lower((1-defect)/norm),
                   coordinate_inverse_gain_upper=edge.upper(norm/(1-defect)))
    return out


def run():
    rows, offsets, errors = calibrated_events()
    certs = []
    archived = []
    for qi, q in enumerate((-1, 1)):
        for i, route in enumerate(measure.ROUTES):
            row = dict(charge=q, route=list(route))
            for end in (9, 11, 12):
                count = sum(t <= end for _, t in EVENTS)
                row[f"through_tick_{end}"] = certificate(rows[qi, i, :count], errors[qi, i, :count])
            compact = [j for j, (_, t) in enumerate(EVENTS) if t<=8] + list(range(66, 74))
            row["through_tick_8_with_relay"] = certificate(rows[qi, i, compact], errors[qi, i, compact])
            certs.append(row)
            if route == ARCHIVED_ROUTE:
                assert row["through_tick_9"].get("rank") == 31
                archived.append(dict(charge=q, route=list(route),
                                     matrix=[[str(x) for x in r] for r in rows[qi, i, :66]],
                                     offsets=[str(x) for x in offsets[qi, i, :66]],
                                     relay_matrix=[[str(x) for x in r] for r in rows[qi, i, 66:]],
                                     relay_offsets=[str(x) for x in offsets[qi, i, 66:]],
                                     event_effect_error_upper=[edge.upper(x) for x in errors[qi,i]]))
    hidden = []
    for i, route in enumerate(measure.ROUTES):
        delta = offsets[1, i] - offsets[0, i]
        matrix = np.column_stack((rows[0, i], rows[1, i], delta))
        max_errors = [max(a, b) for a, b in zip(errors[0, i], errors[1, i], strict=True)]
        early = [j for j, (_, t) in enumerate(EVENTS) if t<=11] + list(range(66, 77))
        hidden.append(dict(route=list(route), **certificate(matrix[:66], max_errors[:66]),
                           through_tick_11_with_relay=certificate(matrix[early], [max_errors[j] for j in early])))
    assert hidden[measure.ROUTES.index(ARCHIVED_ROUTE)].get("rank") == 63
    return dict(event_labels=[[a,t] for a,t in EVENTS],
                affine_denominator=str(edge.DEN), origin="I8/8",
                normalized_head_coordinate_order="full_clock.json",
                reported_routes=certs, native_affine_clock=archived,
                unknown_charge_joint_state=hidden)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    text = json.dumps(run(), indent=2, allow_nan=False)+"\n"
    if args.write:
        Path(__file__).with_suffix(".json").write_text(text)
    else:
        print(text, end="")
