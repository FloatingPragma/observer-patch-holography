"""Interferometric readiness for private sign records, before any pair is prepared."""
import json

import numpy as np


def parity_instrument():
    # Prepare a balanced route, flip the sign on route one and recombine.
    eye = np.eye(2)
    flip = eye[::-1]
    h = np.array([[1, 1], [1, -1]])/np.sqrt(2)
    bs = np.kron(h, eye)
    route = np.block([[eye, np.zeros((2,2))], [np.zeros((2,2)), flip]])
    bright, dark = np.split((bs @ route @ bs)[:, :2], 2)
    np.testing.assert_allclose(bright+dark, eye, atol=5e-16, rtol=0)
    np.testing.assert_allclose(bright.T@bright + dark.T@dark, eye, atol=6e-16, rtol=0)
    return bright, dark


def prepare_responsive():
    """Keep the outcome bit and flip a bright sign's phase; both branches exit dark."""
    bright, dark = parity_instrument()
    tag = np.diag([1., -1.])
    branches = (tag @ bright, dark)
    isometry = np.vstack(branches)
    np.testing.assert_allclose(isometry.T @ isometry, np.eye(2), atol=6e-16, rtol=0)
    # Each column including its coherence with every other column is retained.
    np.testing.assert_allclose(np.kron(np.eye(2), dark) @ isometry, isometry, atol=6e-16, rtol=0)
    return isometry


def check():
    bright, dark = parity_instrument()
    even = np.array([1., 1.])/np.sqrt(2.)
    odd = np.array([1., -1.])/np.sqrt(2.)
    prepared = prepare_responsive()
    # A sign entangled with a spectator reference exercises both kept outcomes.
    bell = np.eye(2).ravel()/np.sqrt(2.)
    records = (np.kron(prepared, np.eye(2)) @ bell).reshape(2, 2, 2)
    output = {"feedback_log": [float(np.linalg.norm(r)**2) for r in records]}
    np.testing.assert_allclose(output["feedback_log"], [.5, .5], atol=8e-16, rtol=0)
    for name, state in (("even", even), ("odd", odd), ("sharp", np.array([1., 0.]))):
        output[name] = [float(np.linalg.norm(k @ state)**2) for k in (bright, dark)]
    np.testing.assert_allclose([output[k] for k in ("even", "odd", "sharp")],
                               [[1,0],[0,1],[.5,.5]], rtol=0, atol=7e-16)
    # Deliberately correlated failures on a coherent bank. Check the whole effect
    # inequality, which also covers arbitrary mixtures and quantum references.
    dim = 8
    all_dark = np.kron(dark, np.kron(dark, dark))
    sum_faults = sum(np.kron(np.eye(2**j), np.kron(bright, np.eye(2**(2-j)))) for j in range(3))
    if (np.linalg.matrix_rank(all_dark, tol=1e-12) != 1 or
            np.linalg.eigvalsh(sum_faults + all_dark - np.eye(dim)).min() < -1e-14):
        raise ArithmeticError("Local parity checks do not lower-bound joint readiness")
    return output


if __name__ == "__main__":
    print(json.dumps(check(), allow_nan=False))
