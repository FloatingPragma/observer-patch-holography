"""Transfer a measured ordinary Poisson calibration to a coherently pulsed native write."""

import argparse
from fractions import Fraction as Q
import math
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import eigh, schur
from scipy import sparse

from . import rest_limit as certificate

from . import poisson_response as model

lab = model.lab
HERE = lab.HERE
READS = 2**48
DUTY = 1e-4
STEP = 1 / 8
CYCLES = 2**21
FILTER_LENGTH = 2**21


def clock_average(x, length):
    """Mean of exp(i j x) for j=0,...,length-1, including narrow spectral gaps."""
    x = (x + np.pi) % (2 * np.pi) - np.pi
    return (
        np.exp(1j * (length - 1) * x / 2)
        * np.sinc(length * x / (2 * np.pi))
        / np.sinc(x / (2 * np.pi))
    )


def unitary(h, time):
    e, v = eigh(h)
    return (v * np.exp(-1j * e * time)) @ v.T


def filtered(e, v):
    uniform = np.ones(len(e)) / np.sqrt(len(e))
    target = round(float(e[0]), 8)
    response = clock_average(-STEP * (e - target), FILTER_LENGTH)
    coeff = (v.T @ uniform) * response**2
    success = float(np.vdot(coeff, coeff).real)
    pure = v @ coeff / np.sqrt(success)
    # All initially excited amplitudes have a phase separated from the window center.
    # For any trial of the rectangular filter |f| <= 1/(M |sin(phase/2)|).
    factor = min(
        1.0, 1 / (FILTER_LENGTH * min(abs(np.sin(STEP * (e[1:] - target) / 2))))
    )
    bound = factor**2 / np.sqrt(success)
    return pure, success, target, bound


def preparation_certificate(h, target):
    centers, matrix, gram = certificate.spectral_disks(sparse.csr_matrix(h))
    scale, radius = certificate.SCALE, certificate.RADIUS
    reference = Q(str(target))
    rest = Q(centers[0], 2 * scale)
    separation = Q(centers[1] - centers[0], scale)
    assert separation > 2 * radius
    detuning = abs(rest - reference) + radius / 2
    excited = Q(centers[1], 2 * scale) - radius / 2 - reference
    phase = excited * Q(1, 16)
    assert phase > 0 and Q(centers[-1], 2 * scale) + radius / 2 - reference < 8
    eta = radius / (separation - radius)
    overlap = Q(sum(map(int, matrix[:, 0])) ** 2, len(centers) * int(gram[0, 0])) - eta
    window = 1 - (FILTER_LENGTH * detuning * Q(1, 8)) ** 2 / 2
    assert window > 0
    pass_floor = overlap * window**4
    factor = 1 / (FILTER_LENGTH * (phase - phase**3 / 6))
    # sqrt(pass_floor) >= pass_floor for a probability floor in (0,1).
    assert 0 < pass_floor <= 1 and factor < 1
    excited_bound = factor**2 / pass_floor
    subdivisions = 2**54
    contraction_error = Q(9 * (3 * FILTER_LENGTH + 1), 128 * subdivisions)
    return dict(
        excited_amplitude_bound=math.ceil(excited_bound * 10**10) / 10**10,
        pass_floor=math.floor(pass_floor * 10**8) / 10**8,
        source_monitors_per_fragment=subdivisions,
        monitored_contraction_error=math.ceil(contraction_error * 10**13) / 10**13,
    )


def time_density(cycle, psi):
    tri, modes = schur(cycle, output="complex")
    diag = tri.diagonal()
    np.testing.assert_allclose(tri, np.diag(diag), atol=4e-13)
    phase = np.angle(diag)
    c = modes.conj().T @ psi
    rho = np.outer(c, c.conj()) * clock_average(phase[:, None] - phase[None, :], CYCLES)
    return modes @ rho @ modes.conj().T


def time_bank(cycle, psi):
    p = np.diag(time_density(cycle, psi)).real
    np.testing.assert_allclose(p.sum(), 1, atol=5e-12)
    return p / p.sum()


def fixture(tape):
    board = lab.lab.cold.federation(0)
    marks, _ = lab.lab.cold.read_map(board, tape["spectral_tape"])
    shells = lab.lab.cold.distances(board, marks[0])
    y = int(np.argmax(shells))
    keep = np.arange(board.n) != y
    h, lap, _ = lab.lab.cold.motion(board, [], 1, model.K)
    h0 = h[keep][:, keep].toarray()
    hw, live, _ = lab.lab.cold.motion(board, marks, 1, model.K)
    # The pinned standard uses the checked native joint block in ordinary_geometry.
    near = board.adjacency()[y, :].toarray().ravel()[keep]
    ho = h0 - np.diag((4 * model.K + 0.5) * near)
    return (
        board,
        keep,
        marks,
        y,
        h0,
        ho,
        hw[keep][:, keep].toarray(),
        (lap - live)[keep][:, keep].toarray(),
        near,
    )


def calibrated(counts, gain, prep_error, fixture, name):
    board, keep, marks, y, h0, _, _, cut, near = fixture
    p = np.array(counts) / READS
    assert sum(counts) == READS and min(counts) >= 0
    # One background event per site; trace distance of the selected blank
    # from its ideal ground bounds the probability error for any event.
    width = np.sqrt(lab.lab.RISK / (2 * READS)) + prep_error
    low, high = p - width, p + width
    assert min(low) > 0
    phi = np.sqrt(p)
    m = np.isin(np.arange(board.n), marks[0]).astype(float)[keep]
    i, j = np.where(np.triu(h0, 1) < 0)
    incidence = np.zeros((len(i), len(p)))
    incidence[np.arange(len(i)), i] = 1
    incidence[np.arange(len(i)), j] = -1
    c = -h0[i, j] * phi[i] * phi[j]
    dc = np.maximum(
        c - (-h0[i, j]) * np.sqrt(low[i] * low[j]),
        (-h0[i, j]) * np.sqrt(high[i] * high[j]) - c,
    )

    def source(pp, g):
        v = np.sqrt(pp)
        return (
            (g + 0.5) * pp * near
            if name == "ordinary"
            else (g * m * pp / 4 + v * (cut @ v) / 2)
        )

    middle = source(p, gain.mean())
    if name == "ordinary":
        sl, sh = source(low, gain[0]), source(high, gain[1])
    else:
        # For the cut component phi_x*(phi_x-phi_y), retain both endpoint intervals.
        endpoints = np.flatnonzero(m)
        u, v = endpoints
        product_lo = np.sqrt(low[u] * low[v])
        product_hi = np.sqrt(high[u] * high[v])
        sl = gain[0] * m * low / 4
        sh = gain[1] * m * high / 4
        sl[endpoints] += (low[endpoints] - product_hi) / 2
        sh[endpoints] += (high[endpoints] - product_lo) / 2
    assert sl.sum() > 0
    rhs = middle - middle.sum() * p
    bl = sl - sh.sum() * high
    bh = sh - sl.sum() * low
    drhs = np.maximum(rhs - bl, bh - rhs)
    D = incidence[:, 1:]
    green = np.linalg.inv((D.T * c) @ D)
    u = np.r_[0.0, green @ rhs[1:]]
    # Grounded interval perturbation, retaining the difference across each bond.
    propagate = abs(green @ D.T) * dc
    matrix = propagate @ abs(D)
    margin = 1 - np.linalg.norm(matrix, np.inf)
    assert margin > 0
    residual = abs(green) @ drhs[1:] + propagate @ abs(incidence @ u)
    eu = np.max(residual) / margin
    # Center u with its measured occupation. The exact derivative is 2 p (u - <u>).
    centered = u - p @ u
    mean_error = 2 * eu + width * sum(abs(centered))
    derivative = 2 * p * centered
    bound = 2 * (width * abs(centered) + high * mean_error)
    rad = lab.lab.cold.distances(board, [y] if name == "ordinary" else marks[0])[keep]
    return derivative, bound, rad


def decode(data, tape):
    fixture_data = fixture(tape)
    gain = lab.ordinary_gain(tape)
    outputs = {}
    error = float(data["filter"]["excited_amplitude_bound"])
    for name in ("ordinary", "written"):
        derivative, bound, rad = calibrated(
            data["reference"], gain, error, fixture_data, name
        )
        blank = np.array(data["arms"]["blank"])
        obs = np.array(data["arms"][name])
        assert blank.sum() == obs.sum() == READS and min(blank) >= 0 and min(obs) >= 0
        delta = (obs - blank) / READS
        width = 2 * np.sqrt(lab.lab.RISK / (2 * READS))
        regional = []
        for radius in (0, 1, 3, 6):
            inside = rad <= radius
            if not np.any(inside):
                continue
            point = DUTY * derivative[inside].sum()
            # Zero total derivative lets a large region use the complementary bound.
            uncertainty = DUTY * min(bound[inside].sum(), bound[~inside].sum())
            observed = delta[inside].sum()
            regional.append(
                dict(
                    radius=radius,
                    prediction=float(point),
                    linear_input_band=[
                        float(point - uncertainty),
                        float(point + uncertainty),
                    ],
                    measured=float(observed),
                    measured_interval=[
                        float(observed - width),
                        float(observed + width),
                    ],
                )
            )
        outputs[name] = dict(
            regions=regional,
            relative_l2_difference=float(
                np.linalg.norm(delta - DUTY * derivative)
                / np.linalg.norm(DUTY * derivative)
            ),
            predicted_change=(DUTY * derivative).tolist(),
        )
    return outputs


def run():
    reference = HERE / "mediated_release.json"
    tape = json.loads(reference.read_bytes())
    board, keep, marks, y, h0, ho, hw, _, _ = fixture(tape)
    e, v = eigh(h0)
    psi, success, target, estimated_error = filtered(e, v)
    certified = preparation_certificate(h0, target)
    assert estimated_error <= certified["excited_amplitude_bound"]
    rng = np.random.default_rng(439720)
    refprob = abs(psi) ** 2
    ref = rng.multinomial(READS, refprob / refprob.sum()).tolist()
    half = unitary(h0, STEP * (1 - DUTY) / 2)
    arms = {}
    rejected = {}
    for name, h in (("blank", h0), ("ordinary", ho), ("written", hw)):
        cycle = half @ unitary(h, STEP * DUTY) @ half
        prob = time_bank(cycle, psi)
        arms[name] = rng.multinomial(READS, prob).tolist()
        rejected[name] = int(rng.negative_binomial(READS, success))
    out = dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (reference, HERE / "POISSON_RESPONSE.md", HERE / "WEAK_SOURCE.md")
        },
        cycles=CYCLES,
        width=STEP,
        duty=DUTY,
        accepted=READS,
        filter=dict(
            length=FILTER_LENGTH,
            stages=2,
            target=target,
            pass_probability=success,
            **certified,
            reference_rejects=int(rng.negative_binomial(READS, success)),
            rejected=rejected,
        ),
        clamp=y,
        reference=ref,
        arms=arms,
    )
    out["decoded"] = decode(out, tape)
    return out


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    modes = cli.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    path = HERE / "weak_source.json"
    if args.counts_only:
        out = json.loads(path.read_bytes())
        lab.lab.cold.apparatus.check(
            out["decoded"],
            decode(out, json.loads((HERE / "mediated_release.json").read_bytes())),
        )
    else:
        out = run()
        if args.write:
            path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
        else:
            lab.lab.cold.apparatus.check(json.loads(path.read_bytes()), out)
    print(out["filter"])
    print(
        {
            k: (v["relative_l2_difference"], v["regions"])
            for k, v in out["decoded"].items()
        }
    )
