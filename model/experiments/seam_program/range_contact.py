"""Compile a pure occupation interaction out of the addressed Interval Echo blank seam."""

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1.sectors import (
    PairSector,
    SeamEvents,
)
from model.candidate_laws.interval_echo.candidate import KAPPA
from .clocked_apparatus import patches
from .compiler import Compiler
from .preparation import raw_event, rotation


def rail_sign(comp, rail, helpers):
    """Z on one binary rail, with 31 initially vacant work rails returned to vacancy."""
    bank = [rail, *helpers]
    if KAPPA != 8 or len(bank) != 32 or len(set(bank)) != 32:
        raise ValueError(
            "The sign word requires the contrast-eight bank of 32 distinct rails"
        )
    for mask in range(1, 32):
        for i, port in enumerate(bank):
            if i < (i ^ mask):
                yield from rotation(comp, port, bank[i ^ mask], np.pi / 8)


def phase_word(comp, a, b, blanks, helpers, tau, coupled=True):
    """Refocus the exchange, then undo the two singleton phases on separate vacant rails."""
    unique = [a, b, *blanks, *helpers]
    if len(blanks) != 2 or len(set(unique)) != len(unique) or comp.marker in unique:
        raise ValueError("Separate data, vacancy and work rails are required")
    if not np.isfinite(tau):
        raise ValueError("Exposure must be finite")
    signed_words = {
        a: list(rail_sign(comp, a, helpers)),
        b: list(rail_sign(comp, b, helpers)),
    }
    # In the sham the forward and backward services both use separate vacancies.
    forward = [(a, b, tau)] if coupled else [(a, blanks[0], tau), (b, blanks[1], tau)]
    for x, y, time in [*forward, (a, blanks[0], -tau), (b, blanks[1], -tau)]:
        half = (time / 2) % (2 * np.pi)
        for _ in range(2):
            yield from rotation(comp, x, y, half)
            yield from signed_words[x]


def native_check(tau, coupled=True):
    _, data, _, _, _ = patches(0)
    comp = Compiler(data, 0, [1, 1, 1, 1])
    sector = PairSector(data)
    events = SeamEvents(sector, KAPPA)
    a, b = 1, 13
    helpers = list(range(30, 61))
    # Two spectator occupations pad the four input columns to the common binary pair sector.
    indices = [
        sector.index(20, 21),
        sector.index(20, b),
        sector.index(a, 21),
        sector.index(a, b),
    ]
    state = np.zeros((sector.P, 4), dtype=complex)
    state[indices, np.arange(4)] = 1
    word = list(phase_word(comp, a, b, [3, 14], helpers, tau, coupled))
    for edge, time in word:
        raw_event(events, state, edge, time)
    phases = state[indices, np.arange(4)].copy()
    expected = np.ones(4, complex)
    if coupled:
        expected[-1] = np.exp(1j * (2 * KAPPA + 1) * tau)
    state[indices, np.arange(4)] -= expected
    return dict(
        physical_seam_pulses=len(word),
        coupled=coupled,
        tau=tau,
        maximum_column_error=float(np.linalg.norm(state, axis=0).max()),
        diagonal_real=phases.real.tolist(),
        diagonal_imag=phases.imag.tolist(),
    )
