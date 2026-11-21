"""Stationary local readout and accessible modes of a marked, matched torus pair.

For a Fourier routing involution F and equal-slot ray S, the span of S
and FS is invariant under both the contact coin D and every partial hop
(a linear polynomial in F). A newborn contributes a scalar times S.
The two F parity spaces split orthogonally into their bright ray
(I +/- F)S when nonzero, and their dark complement. A nonzero first hop
generates FS from S, so the plane is the minimal *linear invariant span*.
This need not give independent control over its two mature pole weights.
A marker-only read or a sum of natal histories preserves the span.

On that plane the mature eigenvalues are exp(+/- iE) with
cos E = cos(mu) sum_j cos(k_j)/4 and mu = phi/2. At an endpoint FS = +/-S
the accessible pole is just +/- exp(-i mu). On the two dark spaces
the mature poles are -sigma exp(i mu), sigma the routing parity.
A local rail quadrature compares amplitudes by arg(h0 conj(h1)).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .. import formation, law, pair_persistence


PROTOCOL = Path(__file__).with_name("clock_modes_protocol.json")


def contact_fibre(angles):
    """Return routing, bright projector and its two dark parity complements."""
    k = np.asarray(angles, dtype=float)
    if k.shape != (4,) or not np.isfinite(k).all():
        raise ValueError("four finite Fourier angles required")
    angle = np.fmod(k[law.SLOT_AXIS] * law.SLOT_SIGN, 2 * np.pi)
    eye = np.eye(law.SLOTS, dtype=complex)
    routing = np.exp(1j * angle)[:, None] * eye[law.PARTNER]
    # These half-angle forms retain small differences between nearby rays.
    cos_half = np.where(np.abs(angle) == np.pi, 0, np.cos(angle / 2))
    sin_half = np.sin(angle / 2)
    phase = np.exp(0.5j * angle)
    rays = (phase * cos_half, -1j * phase * sin_half)
    bright = np.zeros_like(routing)
    dark = []
    for parity, ray in zip((1, -1), rays):
        project = (eye + parity * routing) / 2
        scale = np.max(np.abs(ray))
        if scale:
            unit = ray / scale
            unit /= np.linalg.norm(unit)
            one = np.outer(unit, unit.conj())
            bright += one
            project -= one
        dark.append(project)
    return routing, bright, tuple(dark)


def mature_born_poles(angles, rail):
    """Signed mature pole projectors on the birth span, including both endpoints."""
    routing, bright, _ = contact_fibre(angles)
    phi = law.LAWS["torus"].rail_contact_phases[rail]
    mu = phi / 2
    a = np.mean(np.cos(angles))
    energy = np.arccos(a * np.cos(mu))
    step = routing @ law.contact_coin(np.eye(law.SLOTS, dtype=complex), phi)
    # The two formal roots are distinct even at a collapsed one-ray endpoint.
    # Multiplying by bright removes the unused root and all dark complements.
    negative = (np.exp(1j * energy) * bright - step @ bright) / (2j * np.sin(energy))
    return energy, negative, bright - negative


def _strip_amplitudes(length, updates, rail):
    """Equal-slot Fourier amplitude from a delta birth through the frozen opening."""
    member = law.LAWS["torus"]
    mu = member.rail_contact_phases[rail] / 2
    opening = member.maturation[2]
    k = 2 * np.pi * np.fft.fftfreq(length)
    a = (3 + np.cos(k)) / 4
    # In the orthonormal frame S,(FS-aS)/b, F has rows (a,b),(b,-a).
    b = np.abs(np.sin(k / 2)) * np.sqrt((1 + a) / 2)
    along = np.ones(length, dtype=complex)
    other = np.zeros_like(along)
    bright_coin = np.exp(-1j * mu)
    dark_coin = -np.exp(1j * mu)
    for age in range(min(updates, opening)):
        theta = law.hop_angle(age, member.maturation)
        leave = np.exp(-1j * theta) * np.cos(theta)
        hop = 1j * np.exp(-1j * theta) * np.sin(theta)
        x, y = bright_coin * along, dark_coin * other
        along, other = leave * x + hop * (a * x + b * y), leave * y + hop * (b * x - a * y)
    n = updates - opening
    if n > 0:
        energy = np.arccos(a * np.cos(mu))
        first = a * bright_coin * along + b * dark_coin * other
        along = (np.cos(n * energy) * along +
                 np.sin(n * energy) * (first - np.cos(energy) * along) / np.sin(energy))
    return k, along


def born_aperture(length, updates, rail):
    """Amplitude conditional on birth and rail at the normalized two-neighbour S port."""
    if not isinstance(length, int) or length < 2 or not isinstance(updates, int) or not 0 <= updates < 8192:
        raise ValueError("distinct receiver sites and a pre-wrap integer age required")
    if rail not in (0, 1):
        raise ValueError("the captured pair has two rails")
    k, along = _strip_amplitudes(length, updates, rail)
    return complex(np.mean((1 + np.exp(1j * k)) * along) / np.sqrt(2))


def leading_aperture(n, rail):
    """Leading unwrapped mature-power amplitude at the registered detector."""
    member = law.LAWS["torus"]
    mu = member.rail_contact_phases[rail] / 2
    curvature = 1 / (4 * np.tan(mu))
    phase = -(member.maturation[2] + n) * mu - np.pi / 4
    return np.exp(1j * phase) / np.sqrt(np.pi * n * curvature)


def natal_acceptance():
    """Probability of the enrolled first marker in a sharp literal collision."""
    size = 3
    carrier = (1, 1, size, 1, law.SLOTS)
    free = np.zeros(carrier + carrier, dtype=complex)
    incoming, opposed = formation.COLLISION
    free[0, 0, 0, 0, incoming, 0, 0, 0, 0, opposed] = 1
    refused, branches = formation.capture(free)
    (site, born), = branches
    expected = np.zeros_like(born)
    expected[0, 0, 0, 0] = np.sin(formation.ANGLE) * law.EQUAL
    birth = float(np.vdot(born, born).real)
    if (site != 0 or np.max(np.abs(born - expected)) > 1e-14 or
            abs(np.vdot(refused, refused).real + birth - 1) > 1e-14):
        raise ArithmeticError("unexpected natal slot, marker or refusal")
    return birth


def flat_control():
    """Material rail analyzer at one corner of a supplied persistent circulation."""
    from ..spectral import rail_phase

    member = law.LAWS["torus"]
    initial = pair_persistence.circulation()
    site = (3, 3, 0, 0)
    records = []
    for age in (0, 63, 127, 128):
        rails = np.stack([law.step_maturing(initial, member, age, rail)[site]
                          for rail in (0, 1)]) / np.sqrt(2)
        p = np.array([np.sum(np.abs(rails[0] + a * rails[1]) ** 2) / 4
                      for a in (1, 1j, -1, -1j)])
        records.append({"starting_age": age, "rail_ports": p.tolist(), "miss": float(1 - p.sum()),
                        "hand_increment": rail_phase(rails)})
    return {"capture_overlap": pair_persistence.capture_projection(initial), "one_update": records}



def station_records():
    """Return every original-start terminal outcome, comparison and held decision."""
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    opening = law.LAWS["torus"].maturation[2]
    birth = natal_acceptance()
    offset = protocol["offset"]
    sections = []
    for wait, tolerance in zip(protocol["mature_waits"], protocol["complex_relative_error_gates"], strict=True):
        appointments = []
        for n in (wait, wait + offset):
            updates = opening + n
            size = 2 * updates + 3
            h = np.array([born_aperture(size, updates, r) for r in (0, 1)])
            leading = np.array([leading_aperture(n, r) for r in (0, 1)])
            error = np.abs(h / leading - 1)
            # The two input rails have amplitudes 1/sqrt(2) relative to the birth.
            rails = np.sqrt(birth / 2) * h
            ports = np.abs(rails[0] + np.array([1, 1j, -1, -1j]) * rails[1]) ** 2 / 4
            admission = float(ports.sum())
            word = (ports[0] - ports[2]) + 1j * (ports[1] - ports[3])
            appointments.append({
                "mature_updates": n, "age": updates, "size": size,
                "birth_refusal": 1 - birth, "born_port_miss": birth - admission,
                "rail_ports": ports.tolist(), "admission": admission,
                "conditional_rail_amplitudes": [[z.real, z.imag] for z in h],
                "leading_amplitudes": [[z.real, z.imag] for z in leading],
                "relative_errors": error.tolist(), "complex_gate": bool(max(error) <= tolerance),
                "hand": float(np.angle(word)), "visibility": float(2 * abs(word) / admission),
            })
        change = float(np.angle(np.exp(1j * (appointments[1]["hand"] - appointments[0]["hand"]))))
        sections.append({"appointments": appointments, "tolerance": tolerance, "advance": change})
    forecast = offset * (law.LAWS["torus"].rail_contact_phases[1] - law.LAWS["torus"].rail_contact_phases[0]) / 2
    final = sections[-1]["advance"]
    return {"prepared_flat_control": flat_control(), "pairs": sections, "advance_forecast": forecast, "final_advance_error": abs(final - forecast),
            "final_phase_gate": bool(final > 0 and abs(final - forecast) <= 0.03)}


if __name__ == "__main__":
    print(json.dumps(station_records(), indent=2, allow_nan=False))
