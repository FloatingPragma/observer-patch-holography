"""A prepared pair that persists in free evolution.

The protected state of the torus law is a signed circulation of a bound pair around one
plaquette: the four edges of a square carry the amplitudes +1, -1, +1, -1, so that at each
corner the two occupied slots cancel and the overlap with the equal-slot vector S vanishes.
On the complement of S the contact coin acts as the number -exp(i phi / 2), and the routing
swaps the two ends of each edge, which leaves the circulation unchanged. The circulation is
therefore an eigenvector of every update of the law, at every age of the pair, and its
capture projection, the overlap a returning marker would need to release it, is zero at every
address. Under independent transport of the two constituents the same state keeps only a
quarter of its matched coincident weight after one update.
"""
from __future__ import annotations

import time

import numpy as np

from .law import EQUAL, LAWS, SLOTS, step_contact, step_maturing, step_pair

# Recorded by the observers in the pair-persistence experiment on the torus law.
RECORDED = {
    "cohorts": 5,
    "trials_per_cohort": 4194304,
    "separations": 0,
    "control_separated": [3144922, 4194304],
    "ordinary_return_separation": 0.0047253,
    "route_cut_return_separation": 0.0116733,
}


def circulation(size: int = 8, corner=(3, 3)) -> np.ndarray:
    """The signed plaquette circulation of a bound pair, shape (N, N, 1, 1, 8), unit norm."""
    x, y = corner
    chi = np.zeros((size, size, 1, 1, SLOTS), complex)
    edges = [((x, y), 0, (x + 1, y), 1, 1.0),        # x edge at the bottom
             ((x + 1, y), 2, (x + 1, y + 1), 3, -1.0),  # y edge on the right
             ((x + 1, y + 1), 1, (x, y + 1), 0, 1.0),   # x edge at the top
             ((x, y + 1), 3, (x, y), 2, -1.0)]          # y edge on the left
    for (a, sa, b, sb, c) in edges:
        chi[a[0], a[1], 0, 0, sa] = c
        chi[b[0], b[1], 0, 0, sb] = c
    return chi / np.linalg.norm(chi)


def capture_projection(chi: np.ndarray) -> float:
    """Largest overlap with the equal-slot vector over all addresses."""
    return float(np.max(np.abs(chi @ EQUAL.conj())))


def embed(chi: np.ndarray) -> np.ndarray:
    """The matched coincident two-constituent state whose bound-sector amplitude is chi."""
    carrier = chi.shape[:4]
    volume = int(np.prod(carrier))
    psi = np.zeros((volume, SLOTS, volume, SLOTS), complex)
    flat = chi.reshape(volume, SLOTS)
    for v in range(volume):
        psi[v, np.arange(SLOTS), v, np.arange(SLOTS)] = flat[v]
    return psi.reshape(carrier + (SLOTS,) + carrier + (SLOTS,))


def matched_weight(psi: np.ndarray) -> float:
    carrier = psi.shape[:4]
    volume = int(np.prod(carrier))
    flat = psi.reshape(volume, SLOTS, volume, SLOTS)
    here = np.arange(volume)
    return float(np.sum(np.abs(flat[here[:, None], np.arange(SLOTS)[None, :], here[:, None], np.arange(SLOTS)[None, :]]) ** 2))


def pair_persistence_run(size: int = 8, updates: int = 24, ages=(0, 64, 127, 128, 8191)):
    """Persistence of the circulation under the torus law, and the independent-transport control."""
    started = time.time()
    law = LAWS["torus"]
    phi = law.rail_contact_phases[0]
    start = circulation(size)
    # The bound sector at several ages, including the opening and the mature pair.
    worst_profile, worst_projection, eigen = 0.0, 0.0, []
    for age in ages:
        chi = start.copy()
        for n in range(updates):
            nxt = step_maturing(chi, law, age + n)
            ratio = np.vdot(chi.ravel(), nxt.ravel())
            eigen.append(complex(ratio))
            worst_profile = max(worst_profile, float(np.max(np.abs(np.abs(nxt) - np.abs(start)))))
            worst_projection = max(worst_projection, capture_projection(nxt))
            chi = nxt
    # Two constituents: the contact law keeps the matched pair; independent transport loses it.
    pair = embed(start)
    contact, control = pair.copy(), pair.copy()
    contact_weights, control_weights, control_density, pair_density = [1.0], [1.0], [], []
    for n in range(updates + 1):
        control_density.append(np.sum(np.abs(control) ** 2, axis=(2, 3, 4, 5, 6, 7, 8, 9)))
        pair_density.append(np.sum(np.abs(contact) ** 2, axis=(2, 3, 4, 5, 6, 7, 8, 9)))
        if n == updates:
            break
        contact = step_contact(contact, phi)
        control = step_pair(control, law)
        contact_weights.append(matched_weight(contact))
        control_weights.append(matched_weight(control))
    expected = -np.exp(0.5j * phi)
    return {"law": law, "phi": phi, "size": size, "updates": updates, "ages": list(ages),
            "profile_change": worst_profile, "capture_projection": worst_projection,
            "eigenvalue": [float(expected.real), float(expected.imag)],
            "eigen_error": float(max(abs(e - expected) for e in eigen)),
            "contact_matched": contact_weights, "control_matched": control_weights,
            "pair_density": np.array(pair_density), "control_density": np.array(control_density),
            "recorded": RECORDED, "seconds": time.time() - started}
