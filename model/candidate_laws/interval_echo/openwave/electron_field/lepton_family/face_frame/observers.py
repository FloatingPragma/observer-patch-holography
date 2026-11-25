"""What two local observers agree on when weak, colour and phase frames differ."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag, expm
from scipy.stats import beta

from .field import port_chart
from .measure import character_frames, faces, native, theorem_projectors

HERE = Path(__file__).resolve().parent
SHOTS = 2**20
SEED = 221568
SIGMA = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])


def su2(axis, angle):
    return expm(
        -0.5j
        * angle
        * np.einsum("a,aij->ij", np.asarray(axis) / np.linalg.norm(axis), SIGMA)
    )


def wedge(u):
    """Exterior action of the actual five-mode block matrix on all 32 occupation masks."""
    index = [[k for k in range(5) if m >> k & 1] for m in range(32)]
    out = np.zeros((32, 32), complex)
    for a, left in enumerate(index):
        for b, right in enumerate(index):
            if len(left) == len(right):
                out[a, b] = np.linalg.det(u[np.ix_(left, right)])
    np.testing.assert_allclose(out.conj().T @ out, np.eye(32), atol=5e-14)
    return out


def frames():
    ca = expm(1j * np.array([[0.2, 0.6j, 0.4], [-0.6j, -0.1, 0.3], [0.4, 0.3, -0.1]]))
    cb = expm(1j * np.array([[0.5, 0.7, 0.1j], [0.7, -0.2, -0.2], [-0.1j, -0.2, -0.3]]))
    out = []
    for color, weak, angle in [
        (ca, su2([1, 2, 3], 0.8), 0.9),
        (cb, su2([2, -1, 0], 1.9), -0.7),
    ]:
        fundamental = block_diag(
            np.exp(-1j * angle / 3) * color, np.exp(1j * angle / 2) * weak
        )
        out.append((fundamental, wedge(fundamental)))
    return out


def rotation(u):
    return np.array(
        [[np.trace(a @ u @ b @ u.conj().T).real / 2 for b in SIGMA] for a in SIGMA]
    )


def observed_frame(ideal, rng):
    # Each of the six local input states is read along all three remote axes.
    n = SHOTS
    intervals = []
    counts = []
    midpoint = np.zeros((3, 3))
    radius = np.zeros_like(midpoint)
    tail = 0.01 / ((18 + 24) * 2)
    for i in range(3):
        for j in range(3):
            draws = rng.binomial(
                n, np.clip((1 + ideal[i, j] * np.array([1, -1])) / 2, 0, 1)
            )
            bounds = np.array([binomial(c, n, tail) for c in draws])
            lo, hi = bounds[0, 0] - bounds[1, 1], bounds[0, 1] - bounds[1, 0]
            midpoint[i, j] = (hi + lo) / 2
            radius[i, j] = (hi - lo) / 2
            counts.append(draws.tolist())
            intervals.append([lo, hi])
    tetra = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3)
    axes = np.array([[1, 2, 3], [2, -1, 2], [-2, 2, 1]], dtype=float)
    axes /= np.linalg.norm(axes, axis=1, keepdims=True)
    held = []
    for matter in (1, -1):
        # Complex conjugation reverses the Pauli Y coordinate on both ends.
        conjugate = np.array([1, matter, 1])
        for u in tetra:
            for a in axes:
                source, read = u * conjugate, a * conjugate
                center = (1 + read @ midpoint @ source) / 2
                spread = abs(read) @ radius @ abs(source) / 2
                predicted = [max(0, center - spread), min(1, center + spread)]
                correct = (1 + read @ ideal @ source) / 2
                count = int(rng.binomial(n, correct))
                lo, hi = binomial(count, n, tail)
                naive = (1 + read @ source) / 2
                held.append(
                    dict(
                        matter=matter,
                        source=u.tolist(),
                        analyzer=a.tolist(),
                        count=count,
                        predicted=predicted,
                        measured=[lo, hi],
                        overlap=bool(lo <= predicted[1] and predicted[0] <= hi),
                        untransported_frame_fits=bool(lo <= naive <= hi),
                    )
                )
    return dict(
        training_counts_plus_minus=counts,
        rotation_intervals=intervals,
        held=held,
        common_confidence=0.99,
    )


def binomial(c, n, tail):
    return [
        0 if c == 0 else float(beta.ppf(tail, c, n - c + 1)),
        1 if c == n else float(beta.isf(tail, c + 1, n - c)),
    ]


def joint_read(ga, gb, rng):
    """A's charged rest preparation, B's energy flag and weak conversion on the same start."""
    p = theorem_projectors()[0][1][0]
    basis = port_chart(p)
    ps = np.array([basis.conj().T @ q @ basis for q in character_frames(faces()[0], p)])
    # Independent family coordinate bases. These are book-keeping choices, not gauge fields.
    fa = expm(1j * np.array([[0.3, 0.2j, 0.4], [-0.2j, -0.1, -0.5], [0.4, -0.5, -0.2]]))
    fb = expm(1j * np.array([[-0.2, 0.8, 0], [0.8, 0.2, 0.2j], [0, -0.2j, 0]]))
    rows = []
    for site in (0, 5):
        x = basis.conj().T[:, site]
        x = x / np.linalg.norm(x)
        for charge in (-1, 1):
            tx, rx, projectors = (
                (fa, fb, ps) if charge < 0 else (fa.conj(), fb.conj(), ps.conj())
            )
            a = tx @ x
            b = rx @ tx.conj().T @ a
            family_lab = np.array([np.vdot(x, proj @ x).real for proj in projectors])
            family_b = np.array(
                [np.vdot(b, rx @ proj @ rx.conj().T @ b).real for proj in projectors]
            )
            np.testing.assert_allclose(family_b, family_lab, atol=3e-14)
            left, right = (ga, gb) if charge < 0 else (ga.conj(), gb.conj())
            lab_e = np.eye(32)[:, 23]
            local_e = right @ left.conj().T @ (left @ lab_e)
            ladder = right @ np.outer(np.eye(32)[:, 15], lab_e) @ right.conj().T
            accepted = np.linalg.norm(ladder @ local_e) ** 2 / 2
            # At positive rest energy a Dirac spinor has equal L and R norms.
            # The conjugate lane tests the conjugate doublet as usual.
            joint = np.column_stack(
                [family_b * accepted, family_b * (1 - accepted)]
            ).ravel()
            np.testing.assert_allclose(joint, np.repeat(family_lab / 2, 2), atol=3e-14)
            ignored = np.array([np.vdot(b, proj @ b).real for proj in projectors])
            wrong_weak = abs(local_e[23]) ** 2 / 2
            naive = np.column_stack(
                [ignored * wrong_weak, ignored * (1 - wrong_weak)]
            ).ravel()
            rows.append(
                dict(
                    charge=charge,
                    source=site,
                    counts_order="character 0,1,2; weak click then refusal",
                    agreed_probabilities=joint.tolist(),
                    two_receiver_record_banks=rng.multinomial(
                        SHOTS, joint, size=2
                    ).tolist(),
                    untransported_total_variation=float(abs(naive - joint).sum() / 2),
                )
            )
    return rows


def field_records(rng):
    (a, ga), (b, gb) = frames()
    link = gb @ ga.conj().T
    doublet = [15, 23]
    block = link[np.ix_(doublet, doublet)]
    calibration = observed_frame(rotation(block), rng)
    reference = np.eye(32)
    charged = reference[:, 23]
    produce = ga @ charged
    receive = link @ produce
    physical = gb @ charged
    weak = gb @ np.outer(reference[:, 15], charged) @ gb.conj().T
    np.testing.assert_allclose(receive, physical, atol=3e-14)
    np.testing.assert_allclose(weak @ receive, gb @ reference[:, 15], atol=3e-14)
    assert np.linalg.norm(weak @ gb @ reference[:, 24]) < 1e-13
    wrong = float(abs(np.vdot(charged, receive)) ** 2)
    # Colour acts nontrivially on the one-colour part of the same exterior.
    triplet = [1, 2, 4]
    colour = link[np.ix_(triplet, triplet)]
    assert np.linalg.norm(colour - np.eye(3) * np.trace(colour) / 3) > 0.1
    # A central representative of the Z6 quotient is the identity on the five fundamental
    # modes when the integer U(1) generator is 6Y.
    center = block_diag(
        (np.exp(2j * np.pi / 3) * np.exp(-2j * np.pi / 3)) * np.eye(3),
        (-np.exp(1j * np.pi)) * np.eye(2),
    )
    np.testing.assert_allclose(wedge(center), np.eye(32), atol=2e-14)
    # Two paths contain identity or a declared electric holonomy; both map the
    # same charged state from the same source to the same receiving reference.
    Q = np.diag([-(m & 7).bit_count() / 3 + ((m >> 3) & 1) for m in range(32)])
    loop = []
    for sign, index in [(-1, 23), (1, 24)]:
        ket = ga @ reference[:, index]
        for flux in (0, 0.7):
            first = link @ ket
            second = gb @ expm(1j * flux * Q) @ ga.conj().T @ ket
            contrast = np.vdot(first, second)
            np.testing.assert_allclose(contrast, np.exp(1j * flux * sign), atol=3e-14)
            draw = rng.binomial(
                SHOTS, np.clip((1 + np.array([contrast.real, contrast.imag])) / 2, 0, 1)
            )
            loop.append(
                dict(charge=sign, inserted_loop=flux, real_imag_counts=draw.tolist())
            )
    return dict(
        calibration=calibration,
        joint_rest_energy_and_weak_read=joint_read(ga, gb, rng),
        charged_read_with_matched_Higgs_axis=float(
            abs(np.vdot(physical, receive)) ** 2
        ),
        charged_read_using_fixed_coordinate_mask=wrong,
        declaring_colour_controls=triplet,
        closed_path_counts=loop,
    )


def native_records(rng):
    p = theorem_projectors()[0][1][0]
    face = faces()[0]
    chars = character_frames(face, p)
    spins = np.array([1j ** i.bit_count() for i in range(8)]) / np.sqrt(8)
    rows = []
    phases_a = 0.23 * np.arange(12)
    phases_b = 1.6 * np.cos(np.arange(12))
    for charge in (-1, 1):
        h = native(face, charge)[0] / 8
        free = expm(-1j * h)
        a = np.repeat(np.exp(1j * charge * phases_a), 8)
        b = np.repeat(np.exp(1j * charge * phases_b), 8)
        native_detector = np.kron(p, np.outer(spins, spins.conj()))
        receiver = b[:, None] * native_detector * b.conj()[None, :]
        link = b[:, None] * free * a.conj()[None, :]
        for j, pro in enumerate(chars):
            src = np.kron(pro[:, face[0]], spins)
            src /= np.linalg.norm(src)
            out = link @ (a * src)
            projected = np.vdot(out, receiver @ out).real
            raw = free @ src
            truth = np.vdot(raw, native_detector @ raw).real
            wrong = np.vdot(out, native_detector @ out).real
            np.testing.assert_allclose(projected, truth, atol=3e-14)
            # Independent original-start count banks in the two coordinate descriptions.
            count = rng.binomial(SHOTS, [truth, projected]).tolist()
            rows.append(
                dict(
                    charge=charge,
                    character=j,
                    agreed_probability=float(truth),
                    two_banks=count,
                    untransported_detector_probability=float(wrong),
                )
            )
    return rows


def run():
    rng = np.random.default_rng(SEED)
    return dict(
        seed=SEED,
        shots_per_bank=SHOTS,
        declared_electroweak_and_colour=field_records(rng),
        native_local_port_frames=native_records(rng),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    r = run()
    if args.write:
        (HERE / "observers.json").write_text(json.dumps(r, indent=2) + "\n")
    held = r["declared_electroweak_and_colour"]["calibration"]["held"]
    print(
        "Transported / unchanged coordinate successes:",
        sum(x["overlap"] for x in held),
        sum(x["untransported_frame_fits"] for x in held),
        "of",
        len(held),
    )
