"""Count two separate ordinary-energy losses and a surviving charged core."""

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh

from ...row_09_neutral.extended_decay.preparation import mount as gates, synthesize
from ...row_09_neutral.joined_decay.measure import (
    core_pulse,
    delay_kernel,
    draw,
    N_CLOCK,
    PERIOD,
)
from ...row_09_neutral.passive_decay.material import columns, edge
from ...row_09_neutral.passive_decay.measure import (
    AcquisitionLaw,
    SHOTS,
    compare_records,
)
from .material import formation, generator, mount, SOURCE, X, Z, IDENTITY

HERE = Path(__file__).resolve().parent
ACQUISITION = HERE / "acquisition.npz"
WAITS = (0, 1, 2, 4, 8, 16, 32, 64, 128)
XFRAME = np.array([[1.0, 1.0], [1.0, -1.0]]) / np.sqrt(2)
SLEVELS, SFRAME = eigh(SOURCE)


def local(operator, site):
    terms = [IDENTITY, IDENTITY, IDENTITY]
    terms[site] = operator
    return np.kron(terms[0], np.kron(*terms[1:]))


def real_counts(rng, p, law=None):
    p = np.real_if_close(p).real.ravel()
    assert p.min() > -2e-9 and abs(sum(p) - 1) < 1e-8, p
    p = np.maximum(0, p) / sum(np.maximum(0, p))
    return draw(rng, p if law is None else law(p))


class Reader:
    def __init__(self, energy, frame, site, cutoff):
        delta = energy[None, :] - energy[:, None]
        assert np.ptp(energy) < PERIOD / 2
        late = np.exp(-2j * np.pi * (N_CLOCK - 1) / PERIOD * energy)
        pulse = (frame.T @ local(Z, site) @ frame) * late[:, None]
        bins = np.arange(N_CLOCK // (4 * PERIOD), 2 * N_CLOCK // PERIOD + 1)
        self.accept = delay_kernel(delta, bins) * pulse
        self.light = self.accept * ((delta > 0) & (delta < cutoff))
        delta_diff = delta[:, :, None, None] - delta[None, None, :, :]
        mean = delay_kernel(delta_diff, np.array([0])).reshape(delta_diff.shape)
        self.total = pulse[:, :, None, None] * pulse.conj()[None, None, :, :] * mean

    def records(self, rho):
        good = np.sum(
            (self.accept @ rho) @ self.accept.conj().transpose(0, 2, 1), axis=0
        )
        total = np.einsum("ijkl,jl->ik", self.total, rho)
        np.testing.assert_allclose(np.trace(total), np.trace(rho), atol=1e-12)
        return [total - good, good]


def pull(left, weight, right):
    return (left.conj().transpose(0, 2, 1) @ (weight @ right)).sum(axis=0)


def null_bound(first, second, effect):
    """Squared norm of the surviving paths if the two light lowerings are removed."""
    a, b = first.accept, second.accept
    u, v = first.light, second.light
    both = pull(a, pull(b, effect, b), a)
    good = pull(u, pull(v, effect, v), u)
    cross = pull(a, pull(b, effect, v), u)
    false = both + good - cross - cross.conj().T
    return max(0.0, float(eigh(false, eigvals_only=True)[-1])) + 1e-12


def final_read(readers, rho, effect):
    weights = np.empty((2, 2, 2))
    for i, first in enumerate(readers[0].records(rho)):
        for j, last in enumerate(readers[1].records(first)):
            one = np.trace(last @ effect).real
            weights[i, j] = [one, np.trace(last).real - one]
    return weights


def bank(name, h, initial, cutoff, rng, address_rng, law, include_parking=True):
    e, frame = eigh(h)
    core = local(np.outer(SFRAME[:, 0], SFRAME[:, 0]), 0)
    address = local(np.diag([1.0, 0.0]), 0)  # vacancy on primary when a=0
    state = frame.T @ initial
    energy_basis = np.kron(SFRAME, np.kron(XFRAME, XFRAME))
    address_basis = np.kron(IDENTITY, np.kron(XFRAME, XFRAME))
    flight = []
    for t in WAITS:
        psi = frame @ (np.exp(-1j * e * t) * state)
        flight.append(
            dict(
                wait=t,
                core_x_x=real_counts(rng, abs(energy_basis.T @ psi) ** 2, law),
                vacancy_x_x=real_counts(rng, abs(address_basis.T @ psi) ** 2, law),
            )
        )
    readers = [Reader(e, frame, site, cutoff) for site in (1, 2)]
    moving = state * np.exp(-16j * e)
    rho = np.outer(moving, moving.conj())
    signed = []
    phases = np.exp(-16j * e)
    for pair in [readers, readers[::-1]]:
        core_effect = frame.T @ core @ frame
        raw = final_read(pair, rho, core_effect)
        result = dict(
            arm_loss_core=real_counts(rng, raw, law),
            false_double_loss_lower_core=null_bound(*pair, core_effect),
        )
        if include_parking:
            conditioned = []
            for projector in (address, np.eye(8) - address):
                post = phases * (frame.T @ projector @ initial)
                conditioned.append(
                    final_read(
                        pair, np.outer(post, post.conj()), frame.T @ address @ frame
                    )
                )
            result["initial_vacancy_arm_loss_final_vacancy"] = real_counts(
                address_rng, np.array(conditioned), law
            )
            result["false_double_loss_primary_vacancy"] = null_bound(
                *pair, frame.T @ address @ frame
            )
        signed.append(result)
    return dict(name=name, flight=flight, signed=signed)


def gap_clock(h, rng, law):
    e, v = eigh(h)
    gap = e[1] - e[0]
    rows = []
    prior = 0.0
    for j, t in enumerate([1 / 16, 1 / 4, 1, 4, 16, 64]):
        # Prepare and recombine equal amplitudes of the two local energy modes.
        psi = v @ (np.exp(-1j * e * t) / np.sqrt(2))
        z = (v[:, 0] @ psi) * (v[:, 1] @ psi).conj() * 2
        counts = [int(rng.binomial(SHOTS, law((1 + s) / 2))) for s in (z.real, z.imag)]
        c, s = 2 * np.array(counts) / SHOTS - 1
        angle = np.arctan2(s, c)
        prior = angle + 2 * np.pi * round(
            ((4 * prior if j else 0) - angle) / (2 * np.pi)
        )
        rows.append([t, *counts])
    observed = float(prior / 64)
    assert abs(gap - observed) < 0.00007
    return dict(time_cos_sin=rows, per_native_time=observed)


def run(law=None):
    """Acquire every bank and clock. With an archived law, replay its draws."""
    # Flights of at most 128 units and two 16,384-label clocks on eight
    # levels of magnitude below 3.2 units accumulate eigenvalue roundoff below
    # 1e-11 in any detector probability, far below one count in 2**22.
    law = AcquisitionLaw() if law is None else law
    board, dimers, heads, stock = mount()
    form = formation(board, dimers, heads)
    levels = eigh(form, eigvals_only=True)
    cutoff = levels[-1] / 10
    pulse = local(core_pulse(SFRAME), 0)
    rng = np.random.default_rng(241109)
    arng = np.random.default_rng(241110)
    energies = []
    banks = []
    stocks = {}
    for label, cuts in [
        ("native", ()),
        ("cut_left", (1,)),
        ("cut_right", (2,)),
        ("cut_both", (1, 2)),
    ]:
        h = generator(board, dimers, heads, cuts)
        ground = eigh(h)[1][:, 0]
        hh = dict(heads)
        for site in cuts:
            hh[edge(dimers[0][site - 1], dimers[site][0])] = 32
        links, phase = gates(board, dimers, hh)
        err, counts = synthesize(ground, links)
        energies.append(
            dict(
                name=label,
                native_ground_error=float(err),
                gate_counts=counts,
                phases=phase,
            )
        )
        for arm, initial in [("pulse", pulse @ ground), ("idle", ground)]:
            if not cuts:
                stocks[arm] = initial
            banks.append(bank(label + "_" + arm, h, initial, cutoff, rng, arng, law))
    h = generator(board, dimers, heads)
    for arm, initial in stocks.items():
        banks.append(
            bank(
                "changed_" + arm,
                generator(board, dimers, heads, changed=True),
                initial,
                cutoff,
                rng,
                arng,
                law,
            )
        )
        banks.append(
            bank(
                "parked_" + arm, h + local(X, 0) / 2, initial, cutoff, rng, arng, law
            )
        )
    source_clock = gap_clock(form, rng, law)
    receiver_clocks = [gap_clock(-X / 2, rng, law) for _ in (1, 2)]
    law.complete()
    return dict(
        shots=SHOTS,
        dimers=dimers,
        heads=[[p, q, k] for (p, q), k in sorted(heads.items())],
        stock=stock,
        formation=levels.tolist(),
        cutoff=cutoff,
        collar_residual=columns(),
        changed_residual=columns(True),
        isolated_rail_probabilities=(SFRAME**2).T.tolist(),
        density_transition=float(abs(SFRAME[1, 0] * SFRAME[1, 1])),
        source_clock=source_clock,
        receiver_clocks=receiver_clocks,
        ground_preparations=energies,
        banks=banks,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()
    target = HERE / "records.json"
    if args.save:
        law = AcquisitionLaw()
        result = run(law)
        target.write_text(json.dumps(result, indent=1) + "\n")
        law.save(ACQUISITION)
    else:
        result = run(AcquisitionLaw(ACQUISITION))
        compare_records(json.loads(target.read_text()), json.loads(json.dumps(result)))
    print("Verified native charged fork and all signed-energy and vacancy banks.")
