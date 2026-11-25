"""Finite clock and two-path load records on one perturbed thermal cell."""

import json
import numpy as np

from ..measure import nested_regions, reduced_gibbs_states, spectrum
from .short_flight import PARITY_WORD, parity_terminals, paths, shell_averages, short_coefficients


def flight_records(blocks, kappa, times, terminals=None):
    """Exact-ensemble Born probabilities and their response on the two signed tangents."""
    a, b, sender, receiver = paths()
    if terminals is not None:
        sender, receiver = terminals
    bit = 1 << sender
    # At half filling on the standalone 5-regular cell mu=kappa.
    weights = [np.exp(-(energies - kappa * n)/16)
               for n, (energies, _, _, _) in enumerate(blocks)]
    partition = sum(np.sum(w) for w in weights)
    thermal = [(vec * (w / partition)) @ vec.T
               for (_, vec, _, _), w in zip(blocks, weights, strict=True)]
    floor = min(w.min() for w in weights) / partition
    # The unreset full tangent has norm 1/1024. Both signs are faithful at this amplitude.
    amplitude = min(1/32, 512 * floor)
    # Axes: times, launch, baseline and its t-derivative and three state derivatives.
    readings = np.zeros((len(times), 2, 5))
    for n, (energies, vec, masks, index) in enumerate(blocks):
        clicked = (masks & (1 << receiver)) != 0
        detector = vec[clicked].T @ vec[clicked]
        differences = energies[:, None] - energies[None, :]
        for occupied in (0, 1):
            keep = np.flatnonzero(((masks & bit) != 0) == bool(occupied))
            if not len(keep):
                continue
            bath = masks[keep] & ~bit
            # Reset the sender by a full swap with a fresh sharp rail. The old
            # sender is traced out including its correlations with the cell.
            left = [blocks[n - occupied][3][int(m)] for m in bath]
            right = [blocks[n - occupied + 1][3][int(m | bit)] for m in bath]
            base = thermal[n - occupied][np.ix_(left, left)]
            base = base + thermal[n - occupied + 1][np.ix_(right, right)]
            v = vec[keep]
            baseline = (v.T @ (base @ v)) * detector
            single = ((bath >> a) ^ (bath >> b)) & 1
            diag = .5 * (2 * single - 1)
            odd = np.flatnonzero(single)
            partner = [index[int(m ^ (1 << a) ^ (1 << b))] for m in masks[keep[odd]]]
            # After reset there are nine maximally mixed spectator ports. These
            # matrices use the same endpoint tangents as the regional log.
            flip_entries = v[odd].T @ vec[partner]
            exchange = (-flip_entries / 512) * detector
            aligned = ((v.T @ (diag[:, None] * v) + .5 * flip_entries) / 512) * detector
            parity = np.zeros_like(detector)
            if sender not in PARITY_WORD:
                z = np.ones(len(bath), dtype=int)
                for port in PARITY_WORD:
                    z *= 1 - 2 * ((bath >> port) & 1)
                parity = ((v.T @ (z[:, None] * v)) / 2048) * detector
            for i, t in enumerate(times):
                cos = np.cos(differences * t)
                deriv = -differences * np.sin(differences * t)
                readings[i, occupied] += [np.vdot(baseline, cos), np.vdot(baseline, deriv),
                                          np.vdot(exchange, cos), np.vdot(aligned, cos), np.vdot(parity, cos)]

    rows = []
    for t, records in zip(times, readings, strict=True):
        pulse, slope, exchange, aligned, parity = records[1] - records[0]
        # The material clock is a separate isolated blank one-load dimer.
        rows.append(dict(exposure=t, clock_transfer=float(np.sin(t / 2)**2),
                         receiver_empty_and_full=records[:, 0].tolist(),
                         receiver_exchange_derivatives=records[:, 2].tolist(),
                         receiver_aligned_derivatives=records[:, 3].tolist(),
                         receiver_parity_derivatives=records[:, 4].tolist(), parity_response=float(parity),
                         threshold_time_parity=float(-parity/slope),
                         signal=float(pulse), signal_time_derivative=float(slope),
                         exchange_response=float(exchange), aligned_response=float(aligned),
                         threshold_time_exchange=float(-exchange/slope),
                         threshold_time_aligned=float(-aligned/slope)))
    return dict(thermal_floor=float(floor), valid_signed_amplitude=float(amplitude), rows=rows)


def measure():
    times = (1/32, 1/16, 1/8, 1/4)
    reports = {}
    for kappa in (8, 0):
        blocks = spectrum(kappa)
        report = flight_records(blocks, kappa, times)
        terminals = parity_terminals()
        secondary = flight_records(blocks, kappa, times, terminals)
        report["parity_terminals"] = terminals
        report["parity_route"] = secondary["rows"]
        report["parity_route_even_taylor"] = {
            key: [str(x) for x in term] for key, term in short_coefficients(kappa, 4, terminals).items()}
        if kappa:
            report["regional_modular"] = nested_regions(blocks, 1/16, kappa)[-1]
            ports = report["regional_modular"]["ports"]
            state, = reduced_gibbs_states(blocks, 1/16, [ports], kappa)
            values, vec = np.linalg.eigh(state)
            modular_diagonal = (vec**2) @ -np.log(values)
            mask = sum(1 << ports.index(site) for site in PARITY_WORD)
            signs = np.array([(-1)**((i & mask).bit_count()) for i in range(len(state))])
            report["regional_modular"]["parity_response"] = float(np.mean(signs * modular_diagonal))
        report["even_taylor_coefficients"] = {
            key: [str(x) for x in term] for key, term in short_coefficients(kappa).items()}
        report["quartic_shell_average"] = [
            dict(distance=distance, ordered_pairs=count, exchange=str(e), aligned=str(c))
            for distance, (count, e, c) in shell_averages(kappa).items()]
        reports[str(kappa)] = report
    return dict(source_and_terminals=paths(), by_contrast=reports)


if __name__ == "__main__":
    print(json.dumps(measure(), indent=2, allow_nan=False))
