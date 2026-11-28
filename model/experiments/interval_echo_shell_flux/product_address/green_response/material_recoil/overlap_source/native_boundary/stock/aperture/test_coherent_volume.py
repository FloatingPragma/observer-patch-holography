import numpy as np
from scipy.linalg import expm
from scipy.special import entr

from . import coherent_control as control
from . import regional_tomography as meter


def test_conditional_entropy_slope_with_repeated_membership():
    w = np.array([0.15, 0.35, 0.1, 0.4])
    dw = np.array([-0.04, 0.02, 0.08, -0.06])
    membership = np.array([0, 1, 1, 0])
    material = np.array([[0.4, 0.45], [0.2, 0.38], [0.03, 0.06], [0.64, 0.66]])
    scores, radii = meter.reference_scores(
        w, np.zeros(4), np.array([3, 5, 5, 3]), membership, material, np.zeros(4)
    )
    assert not radii.any()

    def direct(weights):
        return (
            sum(entr(weights))
            - sum(entr(np.bincount(membership, weights)))
            + weights @ material
        )

    t = 1e-5
    np.testing.assert_allclose(
        dw @ scores[:, 1:],
        (direct(w + t * dw) - direct(w - t * dw)) / (2 * t),
        atol=1e-10,
        rtol=0,
    )


def test_pinching_propagation_against_direct_reduced_entropy_and_cut():
    h = np.array(
        [[1, -0.5, 0, 0], [-0.5, 2, -0.5, 0], [0, -0.5, 1, -0.5], [0, 0, -0.5, 2]]
    )
    r = expm(-0.4 * h)
    r /= r.trace()
    e, v = np.linalg.eigh(h)
    mask = np.array([True, True, False, False])
    cut = (np.array([1]), np.array([2]))
    for x in (0, 3):
        pinched = r.copy()
        pinched[x] = 0
        pinched[:, x] = 0
        pinched[x, x] = r[x, x]
        d = pinched - r
        work, width = control.work(r, np.zeros_like(r), h, x)
        np.testing.assert_allclose(work, np.trace(h @ d))
        assert width == 0
        slopes, fringe, _, bound = control.readouts(
            mask, cut, h, r, x, control.propagation(e, v, r, x), np.zeros_like(r), r
        )
        assert (bound > 0) == (x == 0)
        for i, t in enumerate(control.AGES):
            u = expm(-1j * t * h)
            delta = u @ d @ u.conj().T
            step = 1e-5
            plus, minus = [
                (r + sign * step * delta)[np.ix_(mask, mask)] for sign in (1, -1)
            ]
            derivative = (
                meter.entropy(plus, False)[0] - meter.entropy(minus, False)[0]
            ) / (2 * step)
            np.testing.assert_allclose(slopes[i], derivative, atol=5e-10, rtol=0)
            np.testing.assert_allclose(
                fringe[i], 2 * delta[cut][0].real, atol=5e-14, rtol=0
            )
