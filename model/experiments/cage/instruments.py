"""Instruments on a profile, and what an observer reads from their records.

Three experiments of the construction manual are executed here on one cell:

    the swing of a cage (section 7.4)    preparation, the observer's own reference and the seam record in one trial
    the fuel dial (section 6.6)          the two levels of the pair, read from a fuel whose level is set by a dial
    a cage kept (section 7.5)            three preparations, five waits, four analyzers, one capture

A function whose name ends in `_law` computes the probabilities of records from the simulator's state. A function
whose name starts with `read_` takes counts and settings and nothing of the simulator. States are states of
independent orbitals, held as occupation matrices `S[i, j] = <c_j^+ c_i>`, and an evolution `U` acts as
`U S U^+`. The record laws are exact for such states: the probability of a later reading of the fuel is linear
in the occupation matrix that an earlier reading leaves.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.linalg import block_diag, expm
from scipy.stats import chi2 as chi_square_law

from model.candidate_laws.scalar_channel.carrier import require
from model.candidate_laws.scalar_channel.instruments import prepare_source
from .construction import ENERGY, REST, SIGMA, STEP, STRENGTH, T_COHERENT, T_SHARP

QUERIES = 256                                # waits of the clock experiment, in steps
SETTINGS = 4                                 # quarter turns of the fuel's phase
DIAL_STRENGTH = 0.1                          # the converter of the fuel dial
DIALS = np.linspace(-3.5, 3.5, 141)          # levels of the fuel under the dial
DIAL_WAITS = (128, 256)                      # steps
REFERENCE_QUERIES = 513                      # the observer's reference is read at 0 to 512 steps
WAITS = (0, 8, 32, 128, 512)                 # steps of the persistence experiment
PREPARATIONS = (('coherent', -2 * ENERGY, T_COHERENT),
                ('sharp', -2 * ENERGY ** 2 / REST, T_SHARP),
                ('none', 0.0, T_COHERENT))
ANALYZERS = (('none', 0.0, 0.0),
             ('swap', -2 * REST, math.pi / (2 * SIGMA)),
             ('half swap', -2 * REST, math.pi / (4 * SIGMA)),
             ('free quarter', 0.0, math.pi / (4 * ENERGY)))
SHOTS = dict(clock=10 ** 8, reference=10 ** 9, dial=4 * 10 ** 6, persistence=10 ** 8)
PHASE_STEP_BOUNDS = (1 / 16, 1 / 4)          # what the observer takes as given about its reference, per step
LINE_DOMAIN = (0.25, 12.0)                   # where the reader looks for a line, in units of the reference
SIGMA_DOMAIN = (0.0, 4.0)                    # where the reader looks for the coupling


# ------------------------------------------------------------------ the simulator's side

def schedules():
    """The waits of the clock experiment in the simulator's parameter: steps of equal and of unequal length."""
    equal = np.arange(QUERIES) * STEP
    unequal = np.r_[0.0, np.cumsum(STEP * (1 + 0.2 * np.sin(np.arange(1, QUERIES) * 0.317)))]
    return np.stack([equal, unequal])


def reference(steps=256):
    """The prepared reference of section 5.7 as an occupation matrix on the 150 modes of one cell."""
    frame, info = prepare_source(steps)
    return frame @ frame.conj().T, float(info['distance_from_exact_sea'])


def rotation(rest=REST, sigma=SIGMA):
    """Rows: the lower and the upper mode of a pair in the basis (e, e')."""
    energy = math.hypot(rest, sigma)
    co, si = math.sqrt((energy + rest) / (2 * energy)), math.sqrt((energy - rest) / (2 * energy))
    return np.array([[-si, co], [co, si]])


def sectors(G, rest=REST, sigma=SIGMA):
    """The four number sectors of a pair in a state of independent orbitals, from its occupation matrix in (e, e')."""
    R = rotation(rest, sigma)
    Gc = R @ G @ R.T
    two = float(np.linalg.det(Gc).real)
    return dict(empty=float(np.linalg.det(np.eye(2) - Gc).real), reference=float(Gc[0, 0].real - two),
                excitation=float(Gc[1, 1].real - two), two=two)


def pair_in_the_reference(Phi, occupation):
    """The occupation matrix of the pair in the prepared reference and its distance from the ideal lower mode."""
    G = Phi.T @ occupation @ Phi
    lower = rotation()[0]
    return G, float(np.linalg.norm(G - np.outer(lower, lower), 2))


def after_the_coherent_pulse(G):
    """The pair after the pulse of section 6.6: occupation and rate of e, the excitation, the work of the pulse."""
    h = np.array([[REST, SIGMA], [SIGMA, -REST]])
    U = expm(-1j * T_COHERENT * (h - 2 * ENERGY * np.diag([1.0, 0.0])))
    after = U @ G @ U.conj().T
    return dict(occupation_of_e=float(after[0, 0].real), rate_of_e=float((-1j * (h @ after - after @ h))[0, 0].real),
                excitation=sectors(after)['excitation'], work_of_the_pulse=float(np.trace(h @ (after - G)).real))


def swing(times, rest=REST, sigma=SIGMA):
    """The occupation of e for the object placed on e: computed from the block h, and the closed form."""
    energy = math.hypot(rest, sigma)
    h = np.array([[rest, sigma], [sigma, -rest]])
    computed = np.array([abs(expm(-1j * t * h)[0, 0]) ** 2 for t in times])
    return computed, 1 - (sigma / energy) ** 2 * np.sin(energy * np.asarray(times)) ** 2


def pulses():
    """What the two pulses leave outside their target, for the ideal lower mode: the coherent and the sharp pulse."""
    h = np.array([[REST, SIGMA], [SIGMA, -REST]])
    lower, upper = rotation()
    coherent = expm(-1j * T_COHERENT * (h - 2 * ENERGY * np.diag([1.0, 0.0]))) @ lower
    sharp = expm(-1j * T_SHARP * (h - (2 * ENERGY ** 2 / REST) * np.diag([1.0, 0.0]))) @ lower
    return float(abs(coherent[1])), float(np.linalg.norm(sharp - (upper @ sharp) * upper))


def fock_operators(n):
    """Operators that remove an object from each of n modes, on the 2^n states of those modes."""
    out = []
    for j in range(n):
        a = np.zeros((2 ** n, 2 ** n), complex)
        for x in range(2 ** n):
            if (x >> j) & 1:
                a[x ^ (1 << j), x] = (-1) ** bin(x & ((1 << j) - 1)).count('1')
        out.append(a)
    return out


def lift(h, ops):
    """The law of many objects that exclude each other, from the law h of one object."""
    return sum(h[i, j] * ops[i].conj().T @ ops[j] for i in range(len(ops)) for j in range(len(ops)) if h[i, j] != 0)


def excitation_identities():
    """[H, X^+] = 2E X^+ and [H, P_X] = 0 on all states of the pair, the reference circulation and the fuel."""
    ops = fock_operators(4)
    h = block_diag(np.array([[REST, SIGMA], [SIGMA, -REST]]), -REST, 0.0)
    H = lift(h, ops)
    lower, upper = rotation()
    a_minus = lower[0] * ops[0] + lower[1] * ops[1]
    a_plus = upper[0] * ops[0] + upper[1] * ops[1]
    X = a_plus.conj().T @ a_minus
    n_plus, n_minus = a_plus.conj().T @ a_plus, a_minus.conj().T @ a_minus
    P = n_plus @ (np.eye(16) - n_minus)
    return float(np.abs(H @ X - X @ H - 2 * ENERGY * X).max()), float(np.abs(H @ P - P @ H).max())


def _free(H):
    values, vectors = np.linalg.eigh(H)
    return lambda t: (vectors * np.exp(-1j * values * t)) @ vectors.conj().T


def joint_law(H, e, r, f, start, waits, potential, duration, wait_law=None):
    """The law of the two readings (c, z) of the trial of section 7.4, for every wait and every phase setting.

    `H` is the law of one object on the modes in use, with the fuel at index `f`. `e` is the profile and `r`
    the observer's reference. `start` is the occupation matrix before the trial. The preparation is the
    potential `potential` on `e` for the time `duration`. `wait_law` replaces the law during the wait, for the
    control with the two modes of the cage stopped. Returns `p[wait, setting, c, z]`.

    The trial: preparation; half pulse between `r` and the fuel; wait; phase setting on the fuel; one step;
    second half pulse; first reading `c` of the fuel; one step; capture of `e` into the fuel as found; second
    reading `z`. The rows of the fuel are carried backwards through the program.
    """
    n = len(H)
    fuel = np.zeros(n)
    fuel[f] = 1
    free = _free(H)
    prepare = expm(-1j * duration * (H + potential * np.outer(e, e)))
    half = expm(-1j * (math.pi / (4 * STRENGTH)) * (H + STRENGTH * (np.outer(r, fuel) + np.outer(fuel, r))))
    capture = expm(-1j * (math.pi / (2 * STRENGTH)) * (H + STRENGTH * (np.outer(e, fuel) + np.outer(fuel, e)))) @ free(STEP)
    first = half @ prepare
    S = first @ start @ first.conj().T
    row = capture[f].copy()                    # the fuel after the capture, as a row on the state after the first reading
    kept = abs(row[f]) ** 2
    row[f] = 0
    values, vectors = np.linalg.eigh(H if wait_law is None else wait_law)
    turn = np.exp(-1j * np.outer(np.asarray(waits, float), values))
    out = np.zeros((len(waits), SETTINGS, 2, 2))
    for quarter in range(SETTINGS):
        phase = np.ones(n, complex)
        phase[f] = np.exp(-1j * quarter * math.pi / 2)
        tail = half @ free(STEP) @ np.diag(phase)
        X = (((row @ tail) @ vectors) * turn) @ vectors.conj().T
        Y = ((tail[f] @ vectors) * turn) @ vectors.conj().T
        alpha = np.einsum('ti,ij,tj->t', X, S, X.conj()).real
        beta = np.einsum('ti,ij,tj->t', X, S, Y.conj())
        full = np.einsum('ti,ij,tj->t', Y, S, Y.conj()).real          # the first reading finds the fuel full
        z0 = (1 - full) * alpha + abs(beta) ** 2
        z1 = full * alpha - abs(beta) ** 2 + kept * full
        out[:, quarter, 0, 1], out[:, quarter, 0, 0] = z0, (1 - full) - z0
        out[:, quarter, 1, 1], out[:, quarter, 1, 0] = z1, full - z1
    return out


def fock_joint_law(start, waits, potential, duration):
    """The same law on all sixteen states of the four modes e, e', the reference and the fuel, by projection."""
    ops = fock_operators(4)
    h = block_diag(np.array([[REST, SIGMA], [SIGMA, -REST]]), -REST, 0.0)
    H = lift(h, ops)
    free = _free(H)
    number = ops[3].conj().T @ ops[3]
    P = [np.eye(16) - number, number]
    hp, hh, hc = h.copy(), h.copy(), h.copy()
    hp[0, 0] += potential
    hh[2, 3] = hh[3, 2] = STRENGTH
    hc[0, 3] = hc[3, 0] = STRENGTH
    first = expm(-1j * (math.pi / (4 * STRENGTH)) * lift(hh, ops)) @ expm(-1j * duration * lift(hp, ops))
    half = expm(-1j * (math.pi / (4 * STRENGTH)) * lift(hh, ops))
    capture = expm(-1j * (math.pi / (2 * STRENGTH)) * lift(hc, ops)) @ free(STEP)
    occupation, orbitals = np.linalg.eigh(start)
    rho = np.eye(16, dtype=complex)
    for k in range(4):
        a = sum(np.conj(orbitals[i, k]) * ops[i] for i in range(4))
        nk = a.conj().T @ a
        rho = rho @ ((1 - occupation[k]) * (np.eye(16) - nk) + occupation[k] * nk)
    rho = first @ rho @ first.conj().T
    out = np.zeros((len(waits), SETTINGS, 2, 2))
    for i, t in enumerate(waits):
        waited = free(t) @ rho @ free(t).conj().T
        for quarter in range(SETTINGS):
            tail = half @ free(STEP) @ expm(-1j * quarter * (math.pi / 2) * number)
            read = tail @ waited @ tail.conj().T
            for c in (0, 1):
                branch = P[c] @ read @ P[c]
                after = capture @ branch @ capture.conj().T
                p = float(np.trace(number @ after).real)
                out[i, quarter, c] = [float(np.trace(branch).real) - p, p]
    return out


def reference_law(H, r, f, start, waits):
    """The trials that read the observer's own reference (section 6.2): `p[preparation, setting, wait]`.

    Two pulses of area pi/4 between `r` and the fuel with a wait between them, four phase settings, and the fuel
    prepared empty (0) or full (1). The entry is the probability of finding the fuel full.
    """
    n = len(H)
    fuel = np.zeros(n)
    fuel[f] = 1
    pulse = expm(-1j * (math.pi / (4 * STRENGTH)) * (H + STRENGTH * (np.outer(r, fuel) + np.outer(fuel, r))))
    values, vectors = np.linalg.eigh(H)
    turn = np.exp(-1j * np.outer(np.asarray(waits, float), values))
    empty = pulse @ start @ pulse.conj().T
    states = (empty, empty + np.outer(pulse[:, f], pulse[:, f].conj()))
    out = np.zeros((2, SETTINGS, len(waits)))
    for quarter in range(SETTINGS):
        phase = np.ones(n, complex)
        phase[f] = np.exp(-1j * quarter * math.pi / 2)
        rows = (((pulse[f] * phase) @ vectors) * turn) @ vectors.conj().T
        for k, S in enumerate(states):
            out[k, quarter] = np.einsum('ti,ij,tj->t', rows, S, rows.conj()).real
    return out


def dial_law(G, dials=DIALS, waits=DIAL_WAITS, strength=DIAL_STRENGTH):
    """The fuel dial on the pair and the fuel: `p[wait, fuel prepared empty or full, dial]` of a full fuel."""
    out = np.zeros((len(waits), 2, len(dials)))
    for j, level in enumerate(dials):
        H = np.array([[REST, SIGMA, strength], [SIGMA, -REST, 0.0], [strength, 0.0, level]])
        values, vectors = np.linalg.eigh(H)
        for i, q in enumerate(waits):
            row = (vectors[2] * np.exp(-1j * values * q * STEP)) @ vectors.T
            for k in (0, 1):
                out[i, k, j] = float((row @ block_diag(G, float(k)) @ row.conj()).real)
    return out


def persistence_law(D, e, Phi, start, term=None, converter=True):
    """The experiment of section 7.5 on the modes of one cell: `p[preparation, wait, analyzer]` of a full fuel.

    `term` is added to the law in every step of the run, for the changed law. With `converter` false the capture
    is switched off. Also returns the probability of the excitation in the state before the analyzer, for every
    preparation and wait.
    """
    law = D if term is None else D + term
    n = len(law)
    free = _free(law)
    number = np.outer(e, e)
    coupled = np.zeros((n + 1, n + 1))
    coupled[:n, :n] = law
    if converter:
        coupled[:n, n] = coupled[n, :n] = STRENGTH * e
    capture = expm(-1j * (math.pi / (2 * STRENGTH)) * coupled)[n, :n] @ free(STEP)
    rows = [capture @ expm(-1j * time * (law + potential * number)) for _, potential, time in ANALYZERS]
    out = np.zeros((len(PREPARATIONS), len(WAITS), len(ANALYZERS)))
    excitation = np.zeros((len(PREPARATIONS), len(WAITS)))
    for ip, (_, potential, time) in enumerate(PREPARATIONS):
        U = free(STEP) @ expm(-1j * time * (law + potential * number)) @ free(STEP)
        S = U @ start @ U.conj().T
        for iq, q in enumerate(WAITS):
            W = free(q * STEP)
            waited = W @ S @ W.conj().T
            excitation[ip, iq] = sectors(Phi.T @ waited @ Phi)['excitation']
            for ia, row in enumerate(rows):
                out[ip, iq, ia] = float((row @ waited @ row.conj()).real)
    return out, excitation


def weight_left(D, Phi, start, times, term=None):
    """The weight that the vector `start` keeps in the span of the columns of `Phi`, under the law and a term."""
    values, vectors = np.linalg.eigh(D if term is None else D + term)
    amplitudes = vectors.T @ start
    psi = vectors @ (np.exp(-1j * np.outer(values, np.asarray(times, float))) * amplitudes[:, None])
    return np.sum(np.abs(Phi.T @ psi) ** 2, axis=0)


def pointers(full):
    """Two pointers copy the fuel after the coupling is off: the records 00, 01, 10, 11 from the probability of a full fuel."""
    full = np.asarray(full, float)
    return np.stack([1 - full, np.zeros_like(full), np.zeros_like(full), full], axis=-1)


def draw(rng, law, shots):
    """Counts of `shots` trials for every setting, from a law whose last axis lists the outcomes."""
    law = np.asarray(law, float)
    flat = law.reshape(-1, law.shape[-1])
    require(flat.min() > -1e-11 and float(np.abs(flat.sum(axis=1) - 1).max()) < 1e-9, 'a record law is no probability law')
    flat = np.clip(flat, 0.0, None)
    flat = flat / flat.sum(axis=1, keepdims=True)
    return rng.multinomial(shots, flat).reshape(law.shape).astype(np.int64)


def counts_against_a_law(counts, law, shots):
    """How counts sit in a record law: chi-square, its degrees of freedom, its two tails, and the largest cell.

    Outcomes with an expectation below 50 counts are left out of the sum. The last axis lists the outcomes.
    """
    counts = np.asarray(counts, float).reshape(-1, np.shape(law)[-1])
    law = np.asarray(law, float).reshape(counts.shape)
    expected = shots * law
    use = expected > 50
    value = float((((counts - expected) ** 2) / np.where(use, expected, 1.0))[use].sum())
    freedom = int(use.sum()) - int(use.any(axis=1).sum())
    variance = expected * (1 - law)
    cells = use & (variance > 0)
    largest = float((np.abs(counts - expected)[cells] / np.sqrt(variance[cells])).max()) if cells.any() else 0.0
    return dict(chi_square=value, degrees_of_freedom=freedom, upper_tail=float(chi_square_law.sf(value, freedom)),
                lower_tail=float(chi_square_law.cdf(value, freedom)), largest_cell_in_standard_errors=largest,
                counts_on_outcomes_of_probability_zero=int(np.asarray(counts)[law <= 0].sum()))


# ------------------------------------------------------------------ the observer's side

def _golden(loss, a, b, rounds):
    ratio = (math.sqrt(5) - 1) / 2
    c, d = b - ratio * (b - a), a + ratio * (b - a)
    for _ in range(rounds):
        if loss(c) < loss(d):
            b = d
        else:
            a = c
        c, d = b - ratio * (b - a), a + ratio * (b - a)
    return (a + b) / 2


def read_clock(counts, shots):
    """The observer's phase and the seam record from the counts `counts[wait, setting, c, z]` of one schedule.

    The phase of the reference at every wait comes from the first reading and the four settings. The seam record
    is the mean of the second reading over the four settings. The steps of the phase must lie inside the bounds
    that the observer takes as given, since counts cannot show a missed full turn.
    """
    first = counts[:, :, 1, :].sum(axis=-1) / shots
    phasor = (first[:, 0] - first[:, 2] - 1j * (first[:, 1] - first[:, 3])) / 4
    steps = np.angle(phasor[1:] * phasor[:-1].conj())
    require(bool(np.all((steps > PHASE_STEP_BOUNDS[0]) & (steps < PHASE_STEP_BOUNDS[1]))), 'a step of the reference lies outside the declared bounds')
    return np.r_[0.0, np.cumsum(steps)], counts[:, :, :, 1].sum(axis=(-1, -2)) / (SETTINGS * shots)


def read_line(phase, record, domain=LINE_DOMAIN):
    """One line fitted to a record against the observer's phase: frequency, standard error, amplitude, share of the variance."""
    phase, record = np.asarray(phase, float), np.asarray(record, float)

    def fit(w):
        X = np.column_stack([np.ones(len(phase)), np.cos(w * phase), np.sin(w * phase)])
        beta = np.linalg.lstsq(X, record, rcond=None)[0]
        return float(((record - X @ beta) ** 2).sum()), beta
    grid = np.arange(domain[0], domain[1], 0.004)
    k = int(np.argmin([fit(w)[0] for w in grid]))
    w = _golden(lambda x: fit(x)[0], grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)], 70)
    rest, beta = fit(w)
    total = float(((record - record.mean()) ** 2).sum())
    jacobian = np.column_stack([np.ones(len(phase)), np.cos(w * phase), np.sin(w * phase),
                                phase * (beta[2] * np.cos(w * phase) - beta[1] * np.sin(w * phase))])
    try:
        error = math.sqrt(rest / (len(phase) - 4) * float(np.linalg.inv(jacobian.T @ jacobian)[3, 3]))
    except np.linalg.LinAlgError:
        error = float('inf')
    return dict(frequency=float(w), standard_error=error, amplitude=float(math.hypot(beta[1], beta[2])),
                share_of_the_variance=(1 - rest / total) if total > 0 else 0.0)


def read_reference(counts, shots):
    """The phase of the observer's own reference at every count of steps, from `counts[preparation, setting, wait, outcome]`."""
    full = counts[..., 1] / shots
    difference = full[1] - full[0]
    phasor = (difference[0] - difference[2] - 1j * (difference[1] - difference[3])) / 4
    steps = np.angle(phasor[1:] * phasor[:-1].conj())
    require(bool(np.all((steps > PHASE_STEP_BOUNDS[0]) & (steps < PHASE_STEP_BOUNDS[1]))), 'a step of the reference lies outside the declared bounds')
    return np.r_[0.0, np.cumsum(steps)]


def read_dial(counts, shots, phases, dials=DIALS, strength=DIAL_STRENGTH, rest=REST, domain=SIGMA_DOMAIN):
    """The coupling of the pair and its level difference from the fuel dial.

    `counts[wait, fuel prepared empty or full, dial, record]`; `phases` are the readings of the observer's own
    reference at the two waits. The difference between the two preparations is the survival of one object in
    the fuel, whatever the state of the pair is. The mass term and the strength of the converter are known in
    units of the reference.
    """
    observed = (counts[:, 1, :, 3] - counts[:, 0, :, 3]) / shots

    def survival(s):
        out = np.zeros((len(phases), len(dials)))
        for j, level in enumerate(dials):
            values, vectors = np.linalg.eigh(np.array([[rest, s, strength], [s, -rest, 0.0], [strength, 0.0, level]]))
            for i, t in enumerate(phases):
                out[i, j] = abs((vectors[2] * np.exp(-1j * values * t / rest)) @ vectors[2]) ** 2
        return out
    loss = lambda s: float(((survival(s) - observed) ** 2).sum())
    grid = np.linspace(domain[0], domain[1], 201)
    k = int(np.argmin([loss(s) for s in grid]))
    s = _golden(loss, grid[max(k - 1, 0)], grid[min(k + 1, 200)], 50)
    slope = (survival(s + 1e-5) - survival(s - 1e-5)) / 2e-5
    error = math.sqrt(loss(s) / (observed.size - 1) / float((slope ** 2).sum()))
    level = 2 * math.hypot(rest, s)
    return dict(coupling=float(s), coupling_standard_error=error, level_difference=level,
                level_difference_standard_error=2 * s / math.hypot(rest, s) * error,
                largest_misfit=float(np.abs(survival(s) - observed).max()))


def effects(sigma, rest=REST):
    """The four analyzers as effects on the occupation matrix of the pair, from the observer's settings and coupling.

    Row j gives the probability of a full fuel after analyzer j as a combination of
    (G_ee, G_e'e', Re G_ee', Im G_ee').
    """
    h = np.array([[rest, sigma], [sigma, -rest]])
    coupled = block_diag(h, 0.0)
    coupled[0, 2] = coupled[2, 0] = STRENGTH
    capture = expm(-1j * (math.pi / (2 * STRENGTH)) * coupled)[2, :2] @ expm(-1j * STEP * h)
    out = []
    for _, potential, time in ANALYZERS:
        row = capture @ expm(-1j * time * (h + potential * np.diag([1.0, 0.0])))
        M = np.outer(row.conj(), row)
        out.append([M[0, 0].real, M[1, 1].real, 2 * M[0, 1].real, 2 * M[0, 1].imag])
    return np.array(out)


def read_sectors(counts, shots, sigma, rest=REST):
    """The four sectors of the pair for every preparation and wait, from `counts[preparation, wait, analyzer, record]`.

    The four counts of a preparation and a wait give the occupation matrix of the pair through the inverse of
    the four effects. An estimate is linear in the counts and may leave the range of a probability by its
    sampling error.
    """
    A = effects(sigma, rest)
    estimate = np.einsum('ij,pqj->pqi', np.linalg.inv(A), counts[..., 3] / shots)
    out = []
    for per_preparation in estimate:
        out.append([sectors(np.array([[a, x + 1j * y], [x - 1j * y, b]]), rest, sigma) for a, b, x, y in per_preparation])
    return out, float(np.linalg.cond(A))
