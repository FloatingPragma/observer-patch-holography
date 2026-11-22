"""Refine the range and duration of a native density contact on two propagating chart registers."""

import argparse
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo.candidate import KAPPA
from .clocked_field import stopping_law
from .range_contact import native_check

HERE = Path(__file__).resolve().parent
G = (2 * KAPPA + 1) / 8
ARMS = (
    "single_site",
    "square_root_width",
    "half_width",
    "whole_chart",
    "weak_half",
    "sham",
)
LABELS = (
    "both_zero",
    "B_plus_x_A_minus_x",
    "B_plus_y_A_minus_y",
    "B_plus_z_A_minus_z",
    "other_opposite_momenta",
    "timeout",
)


def aperture(n, arm):
    widths = dict(
        single_site=1,
        square_root_width=isqrt(n),
        half_width=n // 2,
        whole_chart=n,
        weak_half=n // 2,
        sham=n // 2,
    )
    width = widths[arm]
    members = (np.arange(n) + width // 2) % n < width
    window = members[:, None, None] & members[None, :, None] & members[None, None, :]
    strength = G / n if arm == "weak_half" else (0 if arm == "sham" else G)
    return width, window, strength


def free_step(state):
    """One walk on each register, at zero total momentum and private order zero."""
    # state[r,i,j], r=x_B-x_A. Both coin reflections precede the two signed shifts.
    coin = 2 * state.mean(axis=-1, keepdims=True) - state
    coin = 2 * coin.mean(axis=-2, keepdims=True) - coin
    out = np.empty_like(state)
    for i in range(8):
        for j in range(8):
            shift = [0, 0, 0]
            if i < 6:
                shift[i // 2] -= 1 if i % 2 == 0 else -1
            if j < 6:
                shift[j // 2] += 1 if j % 2 == 0 else -1
            out[..., i ^ 1, j ^ 1] = np.roll(coin[..., i, j], shift, axis=(0, 1, 2))
    return out


def read(state, phase):
    n = len(state)
    v = n**3
    zero = state.sum(axis=(0, 1, 2)) / np.sqrt(v)
    probabilities = [float(np.vdot(zero, zero).real)]
    for axis in range(3):
        reduced = state.sum(axis=tuple(i for i in range(3) if i != axis))
        mode = np.einsum("x,xij->ij", phase, reduced) / np.sqrt(v)
        probabilities.append(float(np.vdot(mode, mode).real))
    norm = float(np.vdot(state.ravel(), state.ravel()).real)
    if abs(norm - 1) > 1e-10 or sum(probabilities) > 1 + 1e-10:
        raise ArithmeticError(
            "The fully retained two-register state lost normalization"
        )
    return np.array([*probabilities, max(0.0, 1 - sum(probabilities)), 0.0])


def held_law(window, strength, time):
    """No transport. This independent Fourier read includes all coin outcomes."""
    phase = np.exp(1j * strength * time * window)
    amplitudes = np.fft.fftn(phase) / window.size
    p = [
        abs(amplitudes[index]) ** 2
        for index in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
    ]
    return np.array([*p, max(0.0, 1 - sum(p)), 0.0])


def simulate(n, arm):
    if type(n) is not int or n < 4 or n % 4:
        raise ValueError("Chart sides must be positive multiples of four")
    width, window, strength = aperture(n, arm)
    dt = 2 / n
    twist = np.exp(1j * strength * dt * window)[..., None, None]
    state = np.full((n, n, n, 8, 8), 1 / (8 * np.sqrt(n**3)), dtype=complex)
    phase = np.exp(-2j * np.pi * np.arange(n) / n)
    fixed = {}
    clock, late = stopping_law(n // 4, n)
    stopped = np.zeros(len(LABELS))
    for epoch in range(n + 1):
        # A constant joint phase preserves the uniform input exactly.
        constant = strength == 0 or width == n
        record = (
            np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]) if constant else read(state, phase)
        )
        stopped += clock[epoch] * record
        if epoch in (n // 4, n // 2):
            fixed[str(epoch * dt)] = record.tolist()
        if epoch < n and not constant:
            state = free_step(free_step(twist * state))
    stopped[-1] = late
    return dict(
        n=n,
        arm=arm,
        width=width,
        aperture_fraction=float(window.mean()),
        strength=strength,
        paired_contacts_per_epoch=64 * n**3 * width**3,
        fixed=fixed,
        clocked_at_one=stopped.tolist(),
        held_at_one=held_law(window, strength, 1.0).tolist(),
    )


def dicke_control(n, filling):
    sites = 8 * n**3
    occupied = int(sites * filling)
    both = occupied * (occupied - 1) / (sites * (sites - 1))
    return dict(
        n=n,
        slots=sites,
        filling=filling,
        distinct_slot_joint=both,
        covariance=-filling * (1 - filling) / (sites - 1),
    )


def run():
    sides = (4, 8, 16, 32)
    return dict(
        protocol_sha256=sha256(
            (HERE / "range_scaling_protocol.json").read_bytes()
        ).hexdigest(),
        source="conditional uniform equal-coin packets in distinct registers",
        bins=LABELS,
        native_contact=[native_check(1 / 16, arm) for arm in (True, False)],
        refinements=[simulate(n, arm) for n in sides for arm in ARMS],
        symmetric_density=[
            dicke_control(n, p) for n in sides for p in (0.25, 0.5, 0.75)
        ],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.parse_args().output.write_text(
        json.dumps(run(), indent=2, allow_nan=False) + "\n"
    )
