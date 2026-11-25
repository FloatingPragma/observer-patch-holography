"""Does one negative-charge source select three independently distinguishable species?"""

import json
from pathlib import Path

import numpy as np
from scipy.stats import beta as binomial_beta

from ..measure import BASE, CY, ONE, build, scan
from .certificate import certify

HERE = Path(__file__).resolve().parent
LAB = HERE.parents[5] / "code/particles/data/particle_reference_values.json"
DELAYS = 2**31
TICK = 1 / 32
POSITION_SHOTS = 2**32
RETURN_SHOTS = 2**45


def lab_comparison(floor, ceiling):
    fixture = json.loads(LAB.read_text())["entries"]
    bounds = {}
    ids = {}
    for name in ("electron", "muon", "tau"):
        entry = fixture[name]
        value = entry["value_gev"]
        bounds[name] = [
            value - entry["error_minus_gev"],
            value + entry["error_plus_gev"],
        ]
        ids[name] = entry["source"]["summary_id"]
    targets = {
        name: [v[0] / bounds["electron"][1], v[1] / bounds["electron"][0]]
        for name, v in bounds.items()
        if name != "electron"
    }
    assert all(v[0] > ceiling / floor for v in targets.values())
    lows, highs = np.array(list(bounds.values())).T
    koide = [sum(lows) / sum(np.sqrt(highs)) ** 2, sum(highs) / sum(np.sqrt(lows)) ** 2]
    return dict(
        source="code/particles/data/particle_reference_values.json",
        lab_koide_interval=koide,
        any_electron_anchor_heaviest_MeV_upper=float(
            bounds["electron"][1] * 1000 * ceiling / floor
        ),
        pdg_summary_ids=ids,
        outside_after_one_electron_anchor={
            k: [round(z, 8) for z in v] for k, v in targets.items()
        },
    )


def panel():
    energy, amplitudes, beta, q, v, h, perms = build()
    reached = scan(energy, amplitudes, beta)
    # Fixed selection before external mass comparison: strongest negative
    # original-source PT at an isolated positive energy in each comparison.
    chosen = [i for minimum, _, i, j in reached if minimum < 0 and j == i + 1][:3]
    assert len(chosen) == 3
    cyc = perms[CY]

    def lift(weights):
        x = q @ (v @ weights)
        return np.column_stack((x, (x[cyc] - x[cyc[cyc]]) / np.sqrt(3)))

    src = lift(beta)
    # Orientation gauge of the second row must agree with the actual emitter.
    sign = np.sign(BASE[:, 1] @ src[ONE, 1])
    src[:, 1] *= sign
    expected = np.zeros_like(src)
    expected[ONE] = BASE
    np.testing.assert_allclose(src, expected, atol=1e-12)
    forward = []
    returns = []
    models = []
    rows = []
    priors = []
    for i in chosen:
        center = round(float(energy[i]), 11)
        phase = (energy - center) * TICK / (2 * np.pi)
        window = (
            np.exp(-1j * np.pi * (DELAYS - 1) * phase)
            * np.sinc(DELAYS * phase)
            / np.sinc(phase)
        )
        ket = lift(beta * window * np.exp(-32j * (energy - center)))
        bra = lift(beta * window.conj())
        ket[:, 1] *= sign
        bra[:, 1] *= sign
        admitted = float(np.linalg.norm(ket) ** 2 / 2)
        forward.append(ket.reshape(-1, 64, 2))
        returns.append(bra.reshape(-1, 64, 2))
        priors.append(admitted)
        pure = lift(np.eye(1, len(energy), i).ravel()).reshape(-1, 64, 2)
        models.append(np.sum(abs(pure) ** 2, axis=(1, 2)) / 2)
        # Forgetting a sharp vacancy address destroys the intersite hopping
        # coherences. Every same-position head block of this Hamiltonian is diagonal.
        loc = np.repeat(np.arange(len(ket) // 64), 64)
        sparse = h.tocoo()
        assert all(a == b or loc[a] != loc[b] for a, b in zip(sparse.row, sparse.col))
        before = np.trace(ket.conj().T @ (h @ ket)).real / (2 * admitted)
        after = (h.diagonal() @ np.sum(abs(ket) ** 2, axis=1)) / (2 * admitted)
        rows.append(
            dict(
                center=center,
                clock_admission=admitted,
                component_weight=float(abs(beta[i] * window[i]) ** 2 / admitted),
                before_location_read=float(before),
                after_location_read=float(after),
            )
        )
    # Single-shot position classification on balanced successfully emitted
    # states, with a decoder fixed from ideal poles. Preparation refusals are
    # separately recorded for every original start.
    rng = np.random.default_rng(610734)
    decoder = np.argmax(models, axis=0)
    confusions = []
    for ket, admitted in zip(forward, priors):
        p = np.sum(abs(ket) ** 2, axis=(1, 2)) / (2 * 3)
        assert abs(p.sum() - admitted / 3) < 1e-12
        ends = rng.multinomial(POSITION_SHOTS, np.r_[p, 1 - p.sum()])
        counts = np.bincount(decoder, weights=ends[:-1], minlength=3).astype(int)
        confusions.append(
            dict(assigned_counts=counts.tolist(), preparation_refusals=int(ends[-1]))
        )
    total_conditional = (
        sum(
            row["assigned_counts"][i] / sum(row["assigned_counts"])
            for i, row in enumerate(confusions)
        )
        / 3
    )
    # The source and its inverse each cost 1/3. The second coherent clock is
    # a forward wait (the adjoint filter enters its bra by algebra).
    bands = []
    for i, d in enumerate(confusions):
        n = sum(d["assigned_counts"])
        yes = d["assigned_counts"][i]
        tail = 0.01 / (2 * len(confusions))
        bands.append(
            [
                binomial_beta.ppf(tail, yes, n - yes + 1),
                binomial_beta.isf(tail, yes + 1, n - yes),
            ]
        )
    interventions = []
    for source in forward:
        arms = []
        for detector in returns:
            pieces = np.einsum("xdi,xdj->xij", detector.conj(), source)
            measured = float(np.sum(abs(pieces) ** 2) / 18)
            idle = float(np.sum(abs(pieces.sum(axis=0)) ** 2) / 18)
            arms.append(
                dict(
                    location_return_probability=measured,
                    idle_return_probability=idle,
                    location_hits=int(rng.binomial(RETURN_SHOTS, measured)),
                    idle_hits=int(rng.binomial(RETURN_SHOTS, idle)),
                )
            )
        interventions.append(arms)
    return dict(
        reached_positive_groups=len(reached),
        selection="three most negative standard-source PT eigenvalues on simple standard rows",
        positive_row_span=[float(energy.min()), float(energy.max())],
        centers=rows,
        delays=DELAYS,
        native_delay_spacing=TICK,
        free_hold=32,
        position_starts=POSITION_SHOTS,
        decoder_by_port_index=decoder.tolist(),
        classification=confusions,
        balanced_accuracy_given_preparation=float(total_conditional),
        balanced_accuracy_simultaneous99=np.mean(bands, axis=0).tolist(),
        ideal_position_laws=[np.array(m).round(12).tolist() for m in models],
        ideal_balanced_location_ceiling=float(np.sum(np.max(models, axis=0)) / 3),
        return_starts_per_setting=RETURN_SHOTS,
        second_filter_return=interventions,
    )


def positive_controls():
    # Attach rather than infer three spectator labels and a mass term.
    # Tensor identities verify that H and all charge reads preserve them.
    label = np.diag([-1, 0, 1])
    h = np.array([[0.5, -0.5], [-0.5, 0.5]])
    masses = np.array([1, 4, (6 + np.sqrt(39)) ** 2])
    ext = np.kron(h, np.eye(3)) + np.kron(np.eye(2), np.diag(masses))
    f = np.kron(np.eye(2), label)
    for obs in (ext, np.kron(np.diag([1, 0]), np.eye(3))):
        np.testing.assert_array_equal(f @ obs, obs @ f)
    ground = np.array([1, 1]) / np.sqrt(2)
    copies = np.kron(ground[:, None], np.eye(3))
    np.testing.assert_allclose(ext @ copies, copies @ np.diag(masses), atol=1e-12)
    koide = float(masses.sum() / np.sqrt(masses).sum() ** 2)
    np.testing.assert_allclose(koide, 2 / 3)
    return dict(
        supplied_tag_eigenvalues=label.diagonal().tolist(),
        supplied_formation_masses=masses.tolist(),
        koide=koide,
    )


def run():
    bound = certify()
    readings = panel()
    return dict(
        registered_certificate=bound,
        observer_panel=readings,
        external_targets=lab_comparison(8, 27),
        supplied_control=positive_controls(),
    )


if __name__ == "__main__":
    (HERE / "results.json").write_text(json.dumps(run(), indent=2) + "\n")
