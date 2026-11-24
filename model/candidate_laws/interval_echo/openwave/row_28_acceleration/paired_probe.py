"""Decoding a bit from two local port readings while preserving the written seam."""

import argparse
import json

import numpy as np
from scipy.linalg import expm

from .cut_signal import LAW, deposited_collar, receive
from .experiment import HERE


def sequential(sign, coin, bit, direction, half, contacts):
    first = receive(deposited_collar(sign, coin), bit, direction, contacts[0], half)
    joint = np.zeros((2, 2))
    for pointer in (0, 1):
        # Store the first pointer, reload the receiver from its private rails
        # and reuse the unnormalised conditional memory left on the same seam.
        branch = first[:, bit, pointer, 1 - pointer]
        second = receive(branch, bit, direction, contacts[1], half)
        joint[pointer] = (abs(second[:, bit, (0, 1), (1, 0)]) ** 2).sum(axis=0)
    assert abs(joint.sum() - 1) < 1e-12
    return joint


def pointer_law(bit, coin, kappa, scale):
    quiet = np.array([[1.0, 0.0], [0.0, 0.0]])
    if bit == 0:
        return quiet
    answer = np.zeros((2, 2))
    for parity in (0, 1):
        angles = (4 * kappa + parity) * scale * np.pi / (2 * np.array([32, 33]))
        p = np.sin(angles) ** 2
        answer += 0.5 * np.outer([1 - p[0], p[0]], [1 - p[1], p[1]])
    return answer if coin is None else (quiet + answer) / 2


def run():
    exchange, contrast = LAW.local_term_parts()
    half = expm(-0.5j * np.pi * LAW.seam_term(8))
    records = []
    rng = np.random.default_rng(281121)
    for kappa, scale in ((8, 1), (8, 1.01), (0, 1)):
        h = exchange + kappa * contrast
        contacts = [expm(-1j * np.pi * scale * h / x) for x in (32, 33)]
        for sign in (-2, 2):
            for coin in (None, 0, 1):
                for bit in (0, 1):
                    for direction in ("left", "right"):
                        joint = sequential(sign, coin, bit, direction, half, contacts)
                        probabilities = pointer_law(bit, coin, kappa, scale)
                        np.testing.assert_allclose(
                            joint, probabilities, atol=2e-13, rtol=0
                        )
                        records.append(
                            dict(
                                kappa=kappa,
                                duration_scale=scale,
                                sign=sign,
                                coin=coin,
                                bit=bit,
                                sender=direction,
                                p_quiet=round(float(joint[0, 0]), 12),
                                joint_counts=rng.multinomial(
                                    2**20, probabilities.ravel()
                                ).tolist(),
                            )
                        )
    return dict(shots=2**20, order=["00", "01", "10", "11"], records=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    path = HERE / "paired_probe.json"
    res = run()
    if parser.parse_args().record:
        with path.open("x") as f:
            json.dump(res, f, separators=(",", ":"))
            f.write("\n")
    else:
        assert res == json.loads(path.read_text())
        print(
            "Verified all paired local pointer records under coherent, dephased and detuned controls."
        )


if __name__ == "__main__":
    main()
