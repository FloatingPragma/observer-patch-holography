"""The seam reader: two port observers, a courier and an agreed public record, run with the law acting.

The reader block of one seam has eleven bits: observer A (value, valid), observer B (value, valid), the courier
(value, valid), two decision bits and three public record bits. Each observer copies the predicate
`B_e = 1{n_a = n_b = 1}` into its private pointer; the courier carries A's copy to B; the comparison writes the
public triplet `7 B_e` only when both copies are valid and agree. Decisions: 1 agreement, 2 mismatch,
3 unavailable. Every stage is a reversible permutation on all 2048 words, dirty words included.

Finite readout: one cell, the full sector of two loads (78 states, doubled loads included) and the fuel state,
tensored with the 2048 words of the reader of seam (0, 1): 161,792 states. The input is the emitted state at
`kappa = 128`, `g = 64`, after the law acts for 4 units of the evolution parameter. Each of the four stages lasts
`tau_read/4` with Hamiltonian `pi (I - G_j)/(2 tau_j)` added to the law, which acts throughout. The comparison
applies the exact reader at the start of the window with the diagonal part of the law; Duhamel's formula bounds the
difference by `tau_read ||K||`.
"""
from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from . import carrier as cr
from .common import clean
from .sectors import G_EMIT, MultisetSector

BITS = 11
APP = 1 << BITS
CONTROLS = {'nominal': {}, 'mismatch': {'flip_message': True}, 'unavailable': {'read_b': False},
            'wrong_agreement': {'common_flip': True}}


def gate(flag, word, read_a=True, read_b=True, communicate=True, flip_message=False, common_flip=False,
         inverse=False):
    """The complete reversible reader on one 11-bit word; `flag` is the source predicate `B_e`."""
    def final(w):
        b, c = (w >> 2) & 3, (w >> 4) & 3
        valid = (b & 2) != 0 and (c & 2) != 0
        equal = (b & 1) == (c & 1)
        decision = 1 if valid and equal else (2 if valid else 3)
        public = 7 * (b & 1) if valid and equal else 0
        return w ^ (decision << 6) ^ (public << 8)

    w = word
    if inverse:
        w = final(w)
        if flip_message:
            w ^= 1 << 4
        if communicate:
            w ^= (w & 3) << 4
        if read_b:
            w ^= (2 | (flag ^ common_flip)) << 2
        if read_a:
            w ^= 2 | (flag ^ common_flip)
        return w
    if read_a:
        w ^= 2 | (flag ^ common_flip)
    if read_b:
        w ^= (2 | (flag ^ common_flip)) << 2
    if communicate:
        w ^= (w & 3) << 4
    if flip_message:
        w ^= 1 << 4
    return final(w)


def stages(flag, read_a=True, read_b=True, communicate=True, flip_message=False, common_flip=False):
    """The four stages as permutations of the 2048 words; each is an involution."""
    words = np.arange(APP, dtype=np.int64)
    A = words ^ (2 | (flag ^ common_flip)) if read_a else words.copy()
    B = words ^ ((2 | (flag ^ common_flip)) << 2) if read_b else words.copy()
    C = words ^ ((words & 3) << 4) if communicate else words.copy()
    if flip_message:
        C ^= 1 << 4
    b, c = (words >> 2) & 3, (words >> 4) & 3
    valid = ((b & 2) != 0) & ((c & 2) != 0)
    equal = (b & 1) == (c & 1)
    decision = np.where(valid & equal, 1, np.where(valid, 2, 3))
    public = np.where(valid & equal, 7 * (b & 1), 0)
    D = words ^ (decision << 6) ^ (public << 8)
    out = [A, B, C, D]
    for p in out:
        if not np.array_equal(p[p], words):
            raise ValueError('reader stage is not involutive')
    return out


def audit():
    """Every word, every flag, every control: the stages compose to the gate, and the gate has an inverse."""
    lines, columns = [], 0
    words = np.arange(APP)
    for label, kw in CONTROLS.items():
        for flag in (0, 1):
            p = words.copy()
            for st in stages(flag, **kw):
                p = st[p]
                columns += APP
            exact = np.array([gate(flag, int(w), **kw) for w in words])
            inverse = all(gate(flag, int(exact[w]), inverse=True, **kw) == w for w in words)
            out = int(p[0])
            lines.append({'control': label, 'flag': flag, 'stages_compose_to_gate': bool(np.array_equal(p, exact)),
                         'gate_invertible_on_all_words': inverse,
                         'decision': (out >> 6) & 3, 'public_triplet': (out >> 8) & 7})
    return {'bits': BITS, 'words': APP, 'stage_columns_checked': columns, 'fresh_word_readings': lines}


def finite_readout(kappa=128.0, wait=4.0, tau_read=1 / 16384, seam=(0, 1)):
    car = cr.cell()
    sec = MultisetSector(car, 2)
    H = sec.hamiltonian(kappa)
    dim = len(sec.states)
    occ = [0] * 12
    occ[seam[0]] = occ[seam[1]] = 1
    p = sec.index[tuple(occ)]
    Ef = float(H[p, p])
    joined = np.zeros((dim + 1, dim + 1))
    joined[:dim, :dim] = H
    joined[-1, -1] = Ef
    joined[p, -1] = joined[-1, p] = G_EMIT
    w, v = np.linalg.eigh(joined)
    emitted = v @ (np.exp(-1j * math.pi / (2 * G_EMIT) * w) * v[-1, :].conj())
    wh, vh = np.linalg.eigh(H)
    psi = np.zeros(dim + 1, complex)
    psi[:dim] = vh @ (np.exp(-1j * wait * wh) * (vh.conj().T @ emitted[:dim]))
    psi[-1] = np.exp(-1j * wait * Ef) * emitted[-1]
    K = sec.K
    exchange_norm = float(np.linalg.eigvalsh(K)[-1])
    Hs = sparse.block_diag((sparse.csr_matrix(H), sparse.csr_matrix([[Ef]])), format='csr')
    V = np.r_[kappa * sec.C, Ef]
    flag = np.zeros(dim + 1, dtype=int)
    flag[p] = 1
    st0, st1 = stages(0), stages(1)
    N = (dim + 1) * APP
    joint = sparse.kron(Hs, sparse.identity(APP, format='csr'), format='csr')
    ident = sparse.identity(N, format='csr')
    columns = np.arange(N)
    actual = np.zeros(N, complex)
    actual[::APP] = psi
    ideal = actual.copy()
    stage_time = tau_read / 4
    norms = []
    for j in range(4):
        perm = (np.arange(dim + 1)[:, None] * APP + np.where(flag[:, None], st1[j][None, :], st0[j][None, :])).ravel()
        if not np.array_equal(perm[perm], columns):
            raise ValueError('source-controlled reader stage is not involutive')
        Gj = sparse.csr_matrix((np.ones(N), (perm, columns)), shape=(N, N))
        gen = -1j * (stage_time * joint + np.pi / 2 * (ident - Gj))
        actual = expm_multiply(gen, actual, traceA=complex(gen.diagonal().sum()))
        ideal = ideal[perm]
        norms.append(float(np.vdot(actual, actual).real))
    ideal = ideal * np.repeat(np.exp(-1j * tau_read * V), APP)

    def public_law(vec):
        weights = np.abs(vec.reshape(dim + 1, APP)) ** 2
        private = weights.sum(axis=0)
        words = np.arange(APP)
        return np.bincount(((words >> 6) & 3) * 8 + ((words >> 8) & 7), weights=private, minlength=32)

    pa, pi = public_law(actual), public_law(ideal)
    tv = float(np.abs(pa - pi).sum() / 2)
    err = float(np.linalg.norm(actual - ideal))
    emit = Fraction(64009, 25600000)
    read = Fraction(12) * Fraction(1, 16384)
    after = emit + Fraction(40, 250 ** 2) + read
    return {'seam': list(seam), 'kappa': kappa, 'wait': wait, 'read_duration': tau_read,
            'source_states_with_fuel': dim + 1, 'joint_dimension': N, 'largest_norm_drift': max(abs(x - 1) for x in norms),
            'exchange_norm': exchange_norm, 'joint_vector_error': err, 'public_total_variation_error': tv,
            'vector_error_bound_on_the_cell': exchange_norm * tau_read, 'family_bound_12_tau': float(read),
            'agreed_true_probability_actual': float(pa[1 * 8 + 7]), 'agreed_true_probability_ideal': float(pi[1 * 8 + 7]),
            'public_failure_bound_family': f'{(emit + read).numerator}/{(emit + read).denominator}',
            'public_failure_bound_family_value': float(emit + read),
            'after_one_read_round_bound_value': float(after),
            'within_bounds': err <= exchange_norm * tau_read + 1e-9 and tv <= err + 1e-9}


def run(full=False):
    return clean({'audit': audit(), 'finite_readout': finite_readout()})


if __name__ == '__main__':
    import json
    print(json.dumps(run(), indent=1))
