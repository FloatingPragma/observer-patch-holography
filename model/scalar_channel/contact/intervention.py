"""Same-law, local scalar phase pulse versus equal-duration sham with two joint pointers."""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from ..carrier import require
from .joint_source import arm_sources


def step(h, time, frame):
    return expm_multiply(-1j * time * csr_matrix(h), frame, traceA=-1j * time * np.trace(h))


def laws(spec):
    settings = spec["phase_kick"]
    seeds = [*spec["endpoints"], spec["mark"]]
    # Two mark sectors, two pulse choices, seven seeds, two waits, two axes,
    # two local collar signs, six station pointers and a complete 'elsewhere' outcome.
    p = np.empty((2, 2, 7, 2, 2, 2, 7))
    for contact in (0, 1):
        arms = arm_sources(spec, contact, 1)
        before = [step(h, settings["pre_wait"], np.eye(len(h))[:, seeds]) for h in arms]
        for control, strength in enumerate(settings["strengths"]):
            during = []
            for h, state in zip(arms, before):
                driven = h.copy()
                driven[spec["mark"], spec["mark"]] += strength
                during.append(step(driven, settings["duration"], state))
            for j, wait in enumerate(settings["post_waits"]):
                u, v = [step(h, wait, w) for h, w in zip(arms, during)]
                for axis, phase in enumerate((1, 1j)):
                    # Initial equal arms and terminal projective read each contribute 1/sqrt(2).
                    full = np.array([abs((u + sign * phase * v)/2)**2 for sign in (1, -1)])
                    mask = np.ones(full.shape[1], dtype=bool)
                    mask[spec["endpoints"]] = False
                    reduced = np.concatenate((full[:, spec["endpoints"]], full[:, mask].sum(axis=1, keepdims=True)), axis=1)
                    p[contact, control, :, j, axis] = reduced.transpose(2, 0, 1)
    require(np.min(p) >= 0 and np.max(abs(p.sum(axis=(-1, -2)) - 1)) < 1e-11,
            "incomplete phase intervention instrument")
    return p


def _intervals(counts, shots, tail):
    low = np.zeros(counts.shape)
    high = np.ones(counts.shape)
    pos, below = counts > 0, counts < shots
    low[pos] = beta.ppf(tail, counts[pos], shots-counts[pos]+1)
    high[below] = beta.isf(tail, counts[below]+1, shots-counts[below])
    # Marginals contain at most 14 probabilities on the acquisition lattice.
    return low - 24e-12, high + 24e-12


def read(spec, counts):
    settings = spec["phase_kick"]
    shots = settings["shots"]
    require(counts.shape == (2, 2, 7, 2, 2, 2, 7) and counts.dtype == np.int64 and
            np.all((0 <= counts) & (counts <= shots)) and np.all(counts.sum(axis=(-1,-2)) == shots),
            "invalid phase intervention receipt")
    # Read material '+' and all six station click frequencies from the same joint counts.
    observables = np.concatenate([counts[..., 0, :].sum(axis=-1, keepdims=True),
                                  counts[..., :6].sum(axis=-2)], axis=-1)
    lower, upper = _intervals(observables, shots, settings["alpha"] / (2 * observables.size))
    delta = (observables[:, 1] - observables[:, 0]) / shots
    low, high = lower[:, 1] - upper[:, 0], upper[:, 1] - lower[:, 0]
    # The difference of differences removes scalar scattering that would occur
    # with a spectator clock. Selection uses simultaneous bounds over all settings.
    extra = delta[1] - delta[0]
    def certified(lo, hi):
        return np.maximum(0, np.maximum(lo, -hi))
    absolute = certified(low[1], high[1])
    interaction = certified(low[1]-high[0], high[1]-low[0])
    floor = np.minimum(absolute, interaction)
    simultaneous = np.minimum(floor[..., :1], floor[..., 1:])
    for emitter in range(6):
        simultaneous[emitter, ..., emitter] = 0
    qualified = int(np.count_nonzero(simultaneous))
    result = {"menus_with_remote_and_material_response": qualified,
              "strongest_joint_lower_bound": float(simultaneous.max()),
              "alpha": settings["alpha"], "shots_per_joint_setting": shots,
              "settings": int(np.prod(counts.shape[:-2]))}
    if not qualified:
        # With no certified remote response, even the argmax can be the emitter.
        return dict(result, seed_seam=None, receiver_seam=None, post_wait=None, analyzer=None,
                    affirmative_clock_and_station_rate_changes=None, contact_specific_changes=None,
                    simultaneous_absolute_lower_bounds=None)
    seed, wait, axis, receiver = map(int, np.unravel_index(simultaneous.argmax(), simultaneous.shape))
    pair = [0, receiver+1]
    menu = (seed, wait, axis)
    return dict(result, seed_seam=[*spec["endpoints"], spec["mark"]][seed],
                receiver_seam=spec["endpoints"][receiver], post_wait=settings["post_waits"][wait],
                analyzer=["X", "Y"][axis], affirmative_clock_and_station_rate_changes=delta[1][menu][pair].tolist(),
                contact_specific_changes=extra[menu][pair].tolist(),
                simultaneous_absolute_lower_bounds=floor[menu][pair].tolist())


def run(spec):
    p = laws(spec)
    counts = np.zeros(p.shape, dtype=np.int64)
    rng = np.random.Generator(np.random.PCG64([spec["seed"], 5880]))
    for idx in np.ndindex(p.shape[:-2]):
        raw = p[idx].ravel()
        units = np.rint(raw / raw.sum() * 10**12).astype(np.int64)
        units[np.argmax(raw)] += 10**12 - units.sum()
        counts[idx] = rng.multinomial(spec["phase_kick"]["shots"], units/10**12).reshape(2, 7)
    return read(spec, counts), counts
