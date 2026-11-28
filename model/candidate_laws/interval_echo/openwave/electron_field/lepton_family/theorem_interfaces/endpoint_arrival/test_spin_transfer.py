"""Cross-check the two resident spin readers with interleaved spectator matter."""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply


def annihilator(mode):
    one = sparse.eye(2, format="csr")
    parity = sparse.diags([1.0, -1.0], dtype=float, format="csr")
    lower = sparse.csr_matrix([[0, 1], [0, 0]])
    result = sparse.eye(1, format="csr")
    for slot in reversed(range(7)):
        row = one if slot > mode else lower if slot == mode else parity
        result = sparse.kron(result, row, format="csr")
    return result


def local_spin(first, second):
    up, down = first.T @ first, second.T @ second
    return (up + down - 2 * up @ down, first.T @ second + second.T @ first, up - down)


def average(vector, observable):
    return np.vdot(vector, observable @ vector).real


def test_observers_read_the_transported_spin():
    # Insert an occupied or empty spectator physically between the holder's
    # spin modes. The Jordan-Wigner strings therefore cross a genuine mode.
    matter = [annihilator(i) for i in range(7)]
    holder, launch, receiver = [matter[i] for i in (0, 2)], matter[3:5], matter[5:7]
    ph, xh, zh = local_spin(*holder)
    pr, xr, zr = local_spin(*receiver)
    ready = ph @ pr
    mark = sparse.eye(128, format="csr") - 2 * launch[1].T @ launch[1]
    write_link = sparse.csr_matrix([[0, 0], [1, 0]])
    excess = sum(x.T @ y for x, y in zip(launch, holder, strict=True))
    forward = sum(z.T @ x for z, x in zip(receiver, launch, strict=True))
    for sign in (1, -1):
        create_pair, move = (excess, forward) if sign == 1 else (excess.T, forward.T)
        forward_with_link = sparse.kron(write_link, move, format="csr")
        generator = (-sign / 2) * (forward_with_link + forward_with_link.T)
        enrolled_bits = (0, 2) if sign == 1 else (3, 4, 5, 6)
        for spectator in (0, 1):
            stock = np.zeros(128)
            stock[sum(1 << j for j in enrolled_bits) + 2 * spectator] = 1
            launched = create_pair @ stock / np.sqrt(2)
            # With the bond closed the fresh occupation-one spin meter refuses.
            np.testing.assert_allclose(average(launched, pr), 0, atol=1e-12)
            for marked in (False, True):
                pulse_input = mark @ launched if marked else launched
                start = np.concatenate([pulse_input, np.zeros(128)])
                flight = expm_multiply(-1j * np.pi * generator, start)
                constant_word, written_word = np.split(flight, 2)
                np.testing.assert_allclose(np.linalg.norm(constant_word), 0, atol=1e-12)
                np.testing.assert_allclose(average(written_word, ready), 1, atol=1e-12)
                np.testing.assert_allclose(
                    average(written_word, zh @ zr), -1, atol=1e-12
                )
                np.testing.assert_allclose(
                    average(written_word, xh @ xr), 1 if marked else -1, atol=1e-12
                )
