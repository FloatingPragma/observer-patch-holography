"""Complete monitored round trips with one retained quantum collar.

The source has two reducing sectors, blank and symmetrically written. A click
projects the scalar onto one station but leaves the collar's conditional
quantum state intact. Every terminal history is read at the same declared
deadline. The ideal detectors, capture swap, phase analyzers and sampling
schedule are supplied apparatus; this module supplies no autonomous clock.

The initial collar is (blank + written)/sqrt(2). Its written free energy is
kappa. The common energies of the six station collars cancel from these
records. A successful return swaps the scalar into an isolated capture mode
with equal energy in the two collar sectors. Timeout branches evolve freely
on the native scalar source to the deadline, without further detection.
"""

import numpy as np
from scipy.linalg import eigh

from ..carrier import require


def labels(ticks):
    """(0,0): outbound timeout; (j,0): return timeout; (j,k): two hits."""
    return np.array([(0, 0)] + [(j, k) for j in range(1, ticks + 1)
                               for k in range(ticks + 1)], dtype=np.int64)


def first_hit(unitary, sender, receiver, ticks):
    """Subnormalized first-hit amplitudes and the complete no-click vector.

    The hit instrument is P U (Q U)^(j-1), and the timeout operator is
    (Q U)^ticks. These alternatives are orthogonal classical records.
    """
    state = np.eye(len(unitary), dtype=complex)[:, sender].copy()
    amplitudes = np.zeros(ticks, dtype=complex)
    for j in range(ticks):
        state = unitary @ state
        amplitudes[j] = state[receiver]
        state[receiver] = 0
    require(abs(np.vdot(state, state).real + np.vdot(amplitudes, amplitudes).real - 1) < 2e-11,
            "first-hit instrument lost probability")
    return amplitudes, state


def _quadratures(weights, overlap):
    """X and Y plus/minus probabilities, without conditioning on the history.

    overlap is <blank scalar|written scalar>, including the collar phase.
    Y+ is (|blank> + i |written>)/sqrt(2), so its signed mean is Im overlap.
    The factor 1/2 from the equal initial collar is included here.
    """
    total = weights.sum(axis=0) / 2
    signed = np.stack([overlap.real, overlap.imag])
    result = np.stack([total[None, :] + signed, total[None, :] - signed], axis=-1) / 2
    require(result.min() >= -2e-11 and np.max(abs(result.sum(axis=(1, 2)) - 1)) < 2e-11,
            "incomplete terminal collar instrument")
    return result


def pair_records(unitaries, sender, receiver, ticks, dt, kappa):
    """Exact records for a directed pair and the independent classical reset.

    The reset control starts from an equal classical collar mixture. After
    each outbound click it discards that collar and prepares an independent
    equal classical mixture at the same seam. This is a specified local
    intervention, not an alternative evaluation of the retained-collar law.
    Its X/Y means vanish, and its two-leg history probabilities factor.
    """
    require(len(unitaries) == 2 and unitaries[0].shape == unitaries[1].shape,
            "two equal scalar sectors required")
    require(type(ticks) is int and ticks > 0 and dt > 0 and sender != receiver,
            "invalid finite radar menu")
    dimension = len(unitaries[0])
    require(0 <= sender < dimension and 0 <= receiver < dimension,
            "radar station outside scalar source")
    outcome_labels = labels(ticks)
    forward = [first_hit(u, sender, receiver, ticks) for u in unitaries]
    reverse = [first_hit(u, receiver, sender, ticks) for u in unitaries]

    # Return-timeout residuals at different j evolve for N-j further steps.
    # Keep the whole scalar vector; multiplying its norms loses the overlap.
    forward_tail = []
    reverse_tail = []
    for sector, u in enumerate(unitaries):
        state = forward[sector][1].copy()
        for _ in range(ticks):
            state = u @ state
        forward_tail.append(state)
        states = [reverse[sector][1].copy()]
        for _ in range(ticks - 1):
            states.append(u @ states[-1])
        reverse_tail.append(states)

    weights = np.zeros((2, len(outcome_labels)))
    overlap = np.zeros(len(outcome_labels), dtype=complex)
    independent = np.zeros(len(outcome_labels))
    phase = np.exp(-1j * kappa * (2 * ticks * dt))
    forward_weights = np.array([abs(row[0]) ** 2 for row in forward])
    reverse_weights = np.array([abs(row[0]) ** 2 for row in reverse])
    reverse_misses = np.array([np.vdot(row[1], row[1]).real for row in reverse])

    for index, (j, k) in enumerate(outcome_labels):
        if j == 0:
            states = forward_tail
            independent[index] = sum(np.vdot(v, v).real for v in states) / 2
        elif k == 0:
            states = [forward[r][0][j - 1] * reverse_tail[r][ticks - j] for r in range(2)]
            independent[index] = forward_weights[:, j - 1].mean() * reverse_misses.mean()
        else:
            # The return click and a local capture swap leave a common isolated
            # scalar mode. No source potential acts on that captured quantum.
            states = [np.array([forward[r][0][j - 1] * reverse[r][0][k - 1]]) for r in range(2)]
            independent[index] = forward_weights[:, j - 1].mean() * reverse_weights[:, k - 1].mean()
        weights[:, index] = [np.vdot(v, v).real for v in states]
        overlap[index] = phase * np.vdot(states[0], states[1])

    require(np.max(abs(weights.sum(axis=1) - 1)) < 2e-11 and abs(independent.sum() - 1) < 2e-11,
            "round-trip outcomes lost probability")
    coherent = _quadratures(weights, overlap)
    dephased = _quadratures(weights, np.zeros_like(overlap))
    reset = np.repeat(np.repeat((independent / 2)[None, :, None], 2, axis=0), 2, axis=2)
    return dict(coherent=coherent, dephased=dephased, independent_reprepared=reset,
                sector_weights=weights, overlap=overlap)


def run(spec):
    """Return a summary and compact exact arrays; this performs no sampling.

    The complete instrument is stable under a source-generator perturbation
    bounded by epsilon: half-diamond distance is at most min(1, T epsilon),
    T=2N dt. Insert identical ideal instruments in the Duhamel/telescoping
    bound and retain their records. Capture can only shorten source action.
    For changing this single contact by delta_g, epsilon=abs(delta_g), since
    the written-arm projector times one scalar occupation has norm one.
    The bound also holds with an untouched quantum reference.
    """
    from .joint_source import radar_source

    endpoints = spec["endpoints"]
    ticks, dt, kappa = spec["ticks"], spec["dt"], spec["kappa"]
    require(len(endpoints) == len(set(endpoints)) and spec["mark"] not in endpoints,
            "clock stations and coherent mark must be distinct")
    unitaries = []
    for marked in (0, 1):
        h = radar_source(spec, marked)
        values, modes = eigh(h)
        unitaries.append((modes * np.exp(-1j * dt * values)) @ modes.conj().T)
    pairs = [(a, b) for a in endpoints for b in endpoints if a != b]
    records = [pair_records(unitaries, a, b, ticks, dt, kappa) for a, b in pairs]
    laws = {key: np.stack([row[key] for row in records]) for key in records[0]}
    coherent = laws["coherent"]
    dephased_tv = abs(coherent - laws["dephased"]).sum(axis=(2, 3)) / 2
    reset_tv = abs(coherent - laws["independent_reprepared"]).sum(axis=(2, 3)) / 2
    histories = coherent.sum(axis=-1)
    reset_histories = laws["independent_reprepared"].sum(axis=-1)
    history_tv = abs(histories - reset_histories).sum(axis=-1) / 2
    require(np.max(abs(histories - laws["dephased"].sum(axis=-1))) < 2e-12,
            "mark dephasing changed the scalar history marginal")
    require(dephased_tv.max() <= .5 + 2e-11,
            "terminal instrument increased initial coherent/dephased trace distance")
    outcome_labels = labels(ticks)
    successful = (outcome_labels[:, 0] > 0) & (outcome_labels[:, 1] > 0)
    summary = {
        "law": "fixed +1 contact; blank/written coherent mark retained across turnaround",
        "preparation": "equal coherent mark, one scalar quantum at sender, six written clock stations",
        "apparatus": "ideal station projections, immediate local capture after return, fixed finite query schedule",
        "terminal_mark_read": "X and Y, each with plus/minus and every first-hit or timeout history",
        "deadline_in_exchange_units": 2 * ticks * dt,
        "directed_pairs": len(pairs), "histories_per_pair": len(outcome_labels),
        "quadratures": ["X", "Y"],
        "independent_reprepared_control": "classical equal mixture before each leg; local reset at outbound click",
        "maximum_coherent_vs_dephased_TV": float(dephased_tv.max()),
        "maximum_coherent_vs_independent_reprepared_TV": float(reset_tv.max()),
        "maximum_history_TV_from_independent_repreparation": float(history_tv.max()),
        "roundtrip_probability_by_pair": histories[:, 0, successful].sum(axis=-1).tolist(),
        "outbound_timeout_probability_by_pair": histories[:, 0, 0].tolist(),
        "return_timeout_probability_by_pair": histories[:, 0, (outcome_labels[:, 0] > 0) &
                                                                 (outcome_labels[:, 1] == 0)].sum(axis=-1).tolist(),
        "source_perturbation_bound": "half-diamond <= min(1, 2*N*dt*abs(delta_g)); identical apparatus and input",
        "statistics": "exact finite Born probabilities; no sampled trials or statistical discovery claim",
        "physical_metric_or_autonomous_schedule": False,
    }
    arrays = {"coherent_radar_pairs": np.asarray(pairs, dtype=np.int64),
              "coherent_radar_labels": outcome_labels}
    arrays.update({"coherent_radar_" + key: value for key, value in laws.items()})
    return summary, arrays
