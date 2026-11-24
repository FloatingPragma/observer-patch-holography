"""Audit which written-collar preparations realize the conditional wall coupling."""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate


def audit():
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    term = exchange + candidate.KAPPA * contrast
    words = np.eye(law.N_SRC)
    memory = np.eye(1 + law.N_LAB)
    binary = list(product((0, 1), repeat=2))

    def sharp(coin, contrast):
        return memory[1 + law.LABELS.index((coin, contrast))]

    def parity(coin_parity, sign_parity):
        return sum(coin_parity**c * sign_parity**int(d < 0) * sharp(c, d) / 2
                   for c, d in product((0, 1), (-2, 2)))

    table = {}
    largest_residual = 0.
    for coin, sign in product((1, -1), repeat=2):
        energies = []
        for ab in binary:
            col = np.kron(parity(coin, sign), words[law.SRC.index(ab)])
            observed = float(col @ term @ col)
            energy = round(observed)
            expected = 0 if coin == 1 and ab == (1, 1) else 32 + (1 - sign) // 2
            assert energy == expected
            largest_residual = max(largest_residual, np.linalg.norm(term @ col - expected * col))
            energies.append(energy)
        table[f'{coin:+d},{sign:+d}'] = energies
    idle_residual = max(np.linalg.norm(term @ np.kron(sharp(coin, 0), word))
                        for coin in (0, 1) for word in words)
    assert max(largest_residual, idle_residual) < 1e-12

    def coin_and_sign(c, d):
        return c, d

    partitions = {
        'unread': lambda c, d: (),
        'read_sign': lambda c, d: d,
        'read_coin': lambda c, d: c,
        'read_both': coin_and_sign,
    }
    left = law.SRC.index((0, 1))
    right = law.SRC.index((1, 1))
    # First endpoint is one end of a loaded dimer. Its other end keeps
    # the complementary load. Only this connector is addressed until X readout.
    probabilities = {}
    for time in (np.pi / 32, np.pi):
        evolve = expm(-1j * time * term)
        readings = {}
        for name, pointer in partitions.items():
            branches = {}
            for c, d in product((0, 1), (-2, 2)):
                key = pointer(c, d)
                branches[key] = branches.get(key, 0) + sharp(c, d) / 2
            p_plus = 0.
            for state in branches.values():
                a = (evolve @ np.kron(state, words[left])).reshape(-1, law.N_SRC)
                b = (evolve @ np.kron(state, words[right])).reshape(-1, law.N_SRC)
                np.testing.assert_allclose(np.delete(a, left, axis=1), 0., atol=1e-13)
                np.testing.assert_allclose(np.delete(b, right, axis=1), 0., atol=1e-13)
                p_plus += np.linalg.norm(a[:, left] + b[:, right])**2 / 4
            readings[name] = round(float(p_plus), 12)
        probabilities['pi/32' if time < 1 else 'pi'] = readings
    assert abs(probabilities['pi/32']['unread']) < 1e-12
    assert probabilities['pi/32']['read_coin'] == .5
    return dict(binary_words=[''.join(map(str, ab)) for ab in binary],
                energies_by_coin_and_sign_parity=table,
                eigenvector_and_idle_residuals_below=1e-12,
                detector_plus_probability=probabilities)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='write a new receipt')
    path = parser.parse_args().output
    result = audit()
    if path is None:
        assert result == json.loads(Path(__file__).with_suffix('.json').read_text())
    else:
        with path.open('x') as stream:
            stream.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))
