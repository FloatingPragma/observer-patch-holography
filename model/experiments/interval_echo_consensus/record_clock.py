"""Native parity references and a time-independent lift of the record-change circuit."""
import json
from itertools import product
from math import sqrt

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import Carrier
from .record_action import addresses, changes, signed_table


def parity_clock():
    # Both stored signs have the same absolute contrast. Encode in the physical echo basis.
    ids = [(1+law.LABELS.index((0, d)))*law.N_SRC + law.SRC.index((0, 0)) for d in (1, -1)]
    columns = law.echo()[:, ids]
    flip = np.eye(2)[::-1]
    parity = (np.eye(2)-flip)/2
    h0, h1 = law.local_term_parts()
    np.testing.assert_allclose(h0 @ columns, columns @ parity, atol=3e-15, rtol=0)
    np.testing.assert_allclose(h1 @ columns, columns, atol=3e-15, rtol=0)
    plus, minus = np.array([1., 1.])/sqrt(2), np.array([1., -1.])/sqrt(2)
    ready = np.kron(plus, minus)
    swap = np.eye(4)[:, [0, 2, 1, 3]]
    worst = 0.
    for donor, t in product((0, -8, .5, 8, 16), (1/3, sqrt(2)/3)):
        pulse = expm(-1j*t*(h0 + donor*h1))
        native = columns.T @ pulse @ columns
        np.testing.assert_allclose(pulse @ columns, columns @ native, atol=6e-14, rtol=0)
        active = np.kron(native, np.eye(2))
        # A classical Fredkin permutation routes one of the stored parity states to the exposed rail.
        for bit, route in enumerate((np.eye(4), swap)):
            final = route.T @ active @ route @ ready
            worst = max(worst, float(np.linalg.norm(final - np.exp(-1j*t*(donor+bit))*ready)))
    if worst >= 1e-13:
        raise ArithmeticError("The donor left a relative phase or failed to restore the rails")
    return worst


def prepare_marker():
    """Write a sharp sign, recover the temporary loads, and prepare the clock character."""
    echo = law.echo()
    initial = np.eye(law.DIM)[:, law.SRC.index((1, 0))]
    written = echo @ initial
    unbalance = np.arange(law.DIM)
    for ci, (coin, _) in enumerate(law.LABELS, start=1):
        if coin == 0:
            for n, (a, b) in enumerate(law.SRC):
                unbalance[ci*law.N_SRC+n] = ci*law.N_SRC+law.SRC.index((b, a))
    separate = written[unbalance].reshape(1+law.N_LAB, law.N_SRC)
    sharp = np.zeros(1+law.N_LAB)
    for coin in (0, 1):
        sharp[1+law.LABELS.index((coin, 1))] = 1/np.sqrt(2)
    np.testing.assert_allclose(separate, np.outer(sharp, initial[:law.N_SRC]), atol=3e-16, rtol=0)
    # The temporary sharp load has been recovered and is parked. Only the sign memory is occupied.
    vacuum = np.eye(law.N_SRC)[:, law.SRC.index((0, 0))]
    stored = np.kron(sharp, vacuum)
    sign = np.array([0, *(d for _, d in law.LABELS)])
    exchange, contrast = law.local_term_parts()
    error = 0.
    for donor in (0, -8, .5, 8, 16):
        half = expm(-.5j*np.pi*(exchange+donor*contrast)) @ stored
        for odd, tag in ((True, 1j), (False, -1j)):
            out = np.repeat(np.where(sign < 0, tag, 1), law.N_SRC)*half
            target = np.zeros_like(sharp)
            for i, (_, d) in enumerate(law.LABELS, start=1):
                if abs(d) == 1:
                    target[i] = (-1 if odd and d < 0 else 1)/2
            ray = np.kron(target, vacuum)
            error = max(error, float(np.linalg.norm(out-ray*np.vdot(ray, out))))
    if error >= 1e-13:
        raise ArithmeticError("A donor or its recovered preparation load affected the marker character")
    return error


def literal_lift():
    """All 768 microscopic states of a three-port path with four active/parked clock pairs."""
    board = Carrier("open path", 3, [(0, 1), (1, 2)])
    pairs = [(0, 1), (0, 2), (1, 2)]
    reference = np.max(board.degree)-board.degree
    changed = []
    words = signed_table()
    for pair in pairs:
        loads = np.array([int(v in pair) for v in range(board.n)])
        bits = []
        for edge in board.seams:
            other = [f for f in addresses(board, edge) if f != edge]
            bits.extend(x != 0 for x in changes(other, edge, loads, words))
        bits.extend(bool(loads[i]) for i in np.flatnonzero(reference))
        changed.append(bits)
    bits_per_state = len(changed[0])
    bank = 4**bits_per_state
    dim = bank*len(pairs)
    raw = np.arange(bank)
    identity = sparse.eye(dim, format="csr", dtype=np.int64)

    def laplacian(permutation):
        if not np.array_equal(permutation[permutation], np.arange(dim)):
            raise ValueError("Each microscopic update must be a reversible involution")
        return identity - sparse.csr_matrix((np.ones(dim, dtype=np.int64),
                                              (permutation, np.arange(dim))), shape=(dim, dim))

    kinetic = np.zeros((len(pairs), len(pairs)), dtype=int)
    moved = sparse.csr_matrix(identity.shape, dtype=np.int64)
    undressed = moved.copy()
    for a, b in board.seams:
        destinations, plain = [], []
        for j, pair in enumerate(pairs):
            exchange = {a: b, b: a}
            k = pairs.index(tuple(sorted(exchange.get(v, v) for v in pair)))
            kinetic[j, j] += 1
            kinetic[k, j] -= 1
            routed = raw.copy()
            for i, (old, new) in enumerate(zip(changed[j], changed[k], strict=True)):
                if old != new:
                    unequal = ((raw >> (2*i)) ^ (raw >> (2*i+1))) & 1
                    routed ^= unequal * (3 << (2*i))
            destinations.append(k*bank+routed)
            plain.append(k*bank+raw)
        moved += laplacian(np.concatenate(destinations))
        undressed += laplacian(np.concatenate(plain))
    clocks = sparse.csr_matrix(identity.shape, dtype=np.int64)
    for i in range(bits_per_state):
        clocks += laplacian(np.arange(dim) ^ (1 << (2*i)))
    # Unnormalized Fourier columns have integer entries; each has squared norm 'bank'.
    embedding = np.zeros((dim, len(pairs)), dtype=np.int64)
    for j, bits in enumerate(changed):
        odd_mask = sum(1 << (2*i+int(not bit)) for i, bit in enumerate(bits))
        embedding[j*bank:(j+1)*bank, j] = [(-1)**(int(n) & odd_mask).bit_count() for n in raw]
    effective = kinetic + 2*np.diag(np.sum(changed, axis=1))
    if not np.array_equal(embedding.T @ embedding, bank*np.eye(len(pairs))):
        raise ArithmeticError("Invalid material code")
    native = moved+clocks                         # Twice the microscopic Hamiltonian.
    if not np.array_equal(native @ embedding, embedding @ effective):
        raise ArithmeticError("The checker clocks and material exchange do not compose")
    control = (undressed+clocks) @ embedding
    if not np.any(bank*control - embedding @ (embedding.T @ control)):
        raise ArithmeticError("Omitting the clock updates no longer exposes the stale record")
    error = 0.
    for t in (1/3, sqrt(2)/3):
        error = max(error, float(np.linalg.norm(expm_multiply((-1j*t/2)*native, embedding/sqrt(bank))
                                                - embedding @ expm((-1j*t/2)*effective)/sqrt(bank), ord=2)))
    if error > 2e-13:
        raise ArithmeticError("Literal pulses disagree on the material code")
    return dict(microscopic_states=dim, material_states=len(pairs), flags=bits_per_state,
                full_column_error=error)


if __name__ == "__main__":
    print(json.dumps(dict(donor_error=parity_clock(), preparation_error=prepare_marker(),
                          static_lift=literal_lift()), allow_nan=False))
