"""Prepared native seams, coarse signed read and phase analyser for PR-177."""

from itertools import product

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import CELL_SEAMS

LAW = candidate.native_law
PORTS = (0, 1, 5)
FINE = (0, 1)
DUAL = ((0, 7), (1, 8))
# These three independent dimers can live on the internal seams of another cell.
FLAGS = ((0, 1), (2, 3))
HAND = (5, 11)
for edge in (*DUAL, *FLAGS, HAND, (0, 5), (1, 5)):
    assert edge in CELL_SEAMS


def repair_blocks(coarse=False):
    """The 99-state echoes before and after a physical collar projection."""
    unitary = LAW.echo()
    blank = unitary[:, : LAW.N_SRC]
    labels = np.array([d for _c, d in LAW.LABELS]).repeat(LAW.N_SRC)
    if coarse:
        labels = abs(labels)
    blocks = {}
    for label in sorted(set(labels)):
        mask = np.r_[np.zeros(LAW.N_SRC, dtype=bool), labels == label]
        recovered = unitary @ (mask[:, None] * blank)
        if np.max(abs(recovered[LAW.N_SRC :])) > 2e-14:
            raise ArithmeticError(
                "Inverse write left a material arm outside the blank collar"
            )
        blocks[int(label)] = recovered[: LAW.N_SRC]
    return blocks


def triangle_blocks(seam, coarse=False):
    """Restrict the recovered projections to the three address states.

    The spectator address is outside the exposed pair and its occupations never
    move. Check the complete local outputs before using this invariant sector.
    """
    effects = {}
    for label, block in repair_blocks(coarse).items():
        if np.max(abs(block - np.diag(np.diag(block)))) > 2e-14:
            raise ArithmeticError("An address changed during publication")
        keys = [(int(port == seam[0]), int(port == seam[1])) for port in PORTS]
        diagonal = [block[LAW.SRC.index(key), LAW.SRC.index(key)] for key in keys]
        if max(abs(x) for x in diagonal) > 1e-14:
            effects[label] = np.diag(diagonal)
    return effects


def dimer(duration, strength=1.0):
    """Exact one-load two-port hold; only the exchange coefficient changes."""
    swap, contrast = LAW._load_ops()
    h = strength * (np.eye(LAW.N_SRC) - swap) / 2 + candidate.KAPPA * contrast
    ids = [LAW.SRC.index(k) for k in ((1, 0), (0, 1))]
    return expm(-1j * duration * h[np.ix_(ids, ids)])


def pair_word(i, j, phase):
    u = np.eye(3, dtype=complex)
    split = dimer(np.pi / 2)
    u[np.ix_([i, j], [i, j])] = split
    u[j, :] *= np.exp(1j * (phase - np.pi / 2))
    return u


def hand_probabilities(n):
    hand = np.linalg.matrix_power(dimer(np.pi / 16), n)[:, 0]
    sine_hand = dimer(np.pi / 2).conj().T @ hand
    return abs(hand) ** 2, abs(sine_hand) ** 2


def flag_probability(duration, strength):
    return float(abs(dimer(duration, strength)[1, 0]) ** 2)


def deposit(bank, record, rho):
    if record in bank:
        bank[record] += rho
    else:
        bank[record] = rho.copy()


def appointment(bank, effects, p, index):
    """Append the latest public label on a positive flag, keep all old labels otherwise."""
    out = {}
    for record, rho in bank.items():
        deposit(out, record, (1 - p) * rho)
        for label, block in effects.items():
            branch = p * block @ rho @ block.conj().T
            if np.linalg.norm(branch) < 1e-14:
                continue
            latest = list(record)
            latest[index] = label
            deposit(out, tuple(latest), branch)
    return out


def history(rho, effects, probabilities, n, latched):
    initial = {("blank",) * len(effects): rho}
    if not latched:
        bank = initial
        for _ in range(n):
            for i, (blocks, p) in enumerate(zip(effects, probabilities, strict=True)):
                bank = appointment(bank, blocks, p, i)
        return bank
    result = {}
    for flags in product((0, 1), repeat=len(effects)):
        weight = np.prod(
            [p if f else 1 - p for f, p in zip(flags, probabilities, strict=True)]
        )
        bank = initial
        for _ in range(n):
            for i, (blocks, f) in enumerate(zip(effects, flags, strict=True)):
                if f:
                    bank = appointment(bank, blocks, 1, i)
        for key, branch in bank.items():
            deposit(result, key, weight * branch)
    return result


def raw_total_two():
    """Address laws after an uncorrected write with all collar outputs retained in the trace."""
    ket = LAW.echo()[:, : LAW.N_SRC].reshape(-1, LAW.N_SRC, LAW.N_SRC)
    ids = [LAW.SRC.index(k) for k in ((0, 2), (1, 1), (2, 0))]
    return (abs(ket[:, ids][:, :, ids]) ** 2).sum(axis=0)
