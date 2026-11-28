"""Record material inside a ball while withholding the source and clamp labels."""

import numpy as np
from scipy.special import entr

from . import regional_tomography as meter


def local_matrix(rho, mask):
    inside = np.flatnonzero(mask)
    block = rho[np.ix_(inside, inside)]
    if mask.all():
        return block
    matrix = np.zeros((len(block) + 1, len(block) + 1))
    matrix[:-1, :-1] = block
    matrix[-1, -1] = 1 - block.trace()
    return matrix


def frozen_read(
    rows, laws, masks, members, active, cells, native_logs, delta, mix, rng
):
    # Each member uses the SAME analyser in every source preparation. Randomly
    # select a position or a quadrature setting with probability proportional to
    # its absolute coefficient. Idle and unused weight have mean sign zero.
    offset = min(k.diagonal().real.min() for k in native_logs)
    by_mask = []
    for mat in native_logs:
        diag = mat.diagonal().real - offset
        i, j = np.triu_indices(len(mat), 1)
        z = mat[i, j]
        by_mask.append((diag, i, j, z))
    norm = max(
        diag.max() + np.sum(abs(z.real) + abs(z.imag)) for diag, _, _, z in by_mask
    )
    probabilities = np.empty((len(rows), 3))
    for k, ((c, _, s, _), m) in enumerate(zip(rows, members, strict=True)):
        mask = masks[m][active]
        r = local_matrix(cells[c, s][3], mask)
        diag, i, j, z = by_mask[m]
        pos = r.diagonal()
        # A position reading with nonnegative weight contributes a positive sign.
        # A quadrature reading contributes its signed port or an idle.
        hit = pos @ diag + (pos[i] + pos[j]) @ (abs(z.real) + abs(z.imag))
        signal = pos @ diag + 2 * (r[i, j].real @ z.real + r[i, j].imag @ z.imag)
        np.testing.assert_allclose(
            signal, np.sum(r * native_logs[m].T).real - offset, atol=8e-13
        )
        probabilities[k] = [
            (hit + signal) / (2 * norm),
            (hit - signal) / (2 * norm),
            1 - hit / norm,
        ]
    assert probabilities.min() >= -1e-14
    sharp = []
    for head in range(4):
        selected = rows[:, 1] == head
        weights = laws[1, selected]
        sharp.append((weights / weights.sum()) @ probabilities[selected])
    preparations = np.vstack((sharp, laws[1:] @ probabilities))
    counts = np.array([rng.multinomial(meter.N, p / p.sum()) for p in preparations])
    freq = norm * (counts[:, 0] - counts[:, 1]) / meter.N
    predicted = meter.product_interval(
        delta, mix @ freq[:4] + np.array([-8, 8]) * norm * meter.HOEFFDING
    )
    held = freq[5] - freq[4] + np.array([-4, 4]) * norm * meter.HOEFFDING
    assert max(held[0], predicted[0]) <= min(held[1], predicted[1])
    return dict(
        norm=norm,
        offset=offset,
        sign_and_idle_counts=counts.tolist(),
        predicted=predicted.tolist(),
        held=held.tolist(),
    )


def entropy_without_labels(rows, laws, masks, members, active, cells, rng, delta, mix):
    """Bound the tangent of S(material | ball) with two measured concavity chords."""
    # Only the source preparation uses the head selector. The regional observer
    # knows the marker membership and the material detector banks for that ball.
    assert np.min(laws) >= 0
    records, entropies, native_logs = [], [], []
    for arm, law in enumerate(laws):
        prior = np.bincount(members, law)
        counts = rng.multinomial(meter.N, prior / prior.sum())
        frequency = counts / meter.N
        local, errors, banks = [], [], []
        for m, marker_mask in enumerate(masks):
            mask = marker_mask[active]
            dim = int(mask.sum()) + int(not mask.all())
            matrix = np.zeros((dim, dim))
            take = np.flatnonzero(members == m)
            for i in take:
                c, _, s, _ = rows[i]
                matrix += law[i] / prior[m] * local_matrix(cells[c, s][3], mask)
            rec, _, receipt = meter.acquire(
                matrix, rng, homogeneous_outside=not mask.all()
            )
            vals, vecs = np.linalg.eigh(rec)
            if arm == 1:
                score = (vecs * -np.log(vals)) @ vecs.conj().T
                if not mask.all():
                    # The homogeneous outside outcome is a classical flag.
                    score[-1, :-1] = score[:-1, -1] = 0
                native_logs.append(score)
            local.append(sum(entr(vals)))
            errors.append(receipt["entropy_error"])
            banks.append(
                dict(marker_ids=np.flatnonzero(marker_mask).tolist(), **receipt)
            )
        err = meter.frequency_error(frequency)
        middle = frequency @ local
        width = (frequency + err) @ errors + err @ abs(np.array(local) - min(local))
        entropies.append([middle - width, middle + width])
        records.append(dict(membership_counts=counts.tolist(), material_banks=banks))
    neg, ref, pos = entropies
    bounds = [pos[0] - ref[1], ref[1] - neg[0]]
    assert bounds[0] <= bounds[1]
    return dict(
        arm_order=["opposite", "native", "repaired"],
        arms=records,
        entropy=entropies,
        concavity_slope_bounds=bounds,
        material_read=frozen_read(
            rows, laws, masks, members, active, cells, native_logs, delta, mix, rng
        ),
    )


def decompose(rows, w, error, material, entropy_error, cells, beta, members):
    """Which part comes from the written labels, ball identity and ordinary material?"""
    group = np.bincount(members, w[0])[members]
    group_error = np.bincount(members, error[0])[members]
    energy, uncertainty = [], []
    for c, _, s, _ in rows:
        h, _, _, _, rec, e, _ = cells[c, s]
        h = h.copy()
        offset = np.mean(h.diagonal())
        h.flat[:: len(h) + 1] -= offset
        energy.append(offset + np.sum(h * rec.real))
        uncertainty.append(np.sum(abs(h) * e))
    # Include +log p_B in the conditional entropy score.
    scores = np.column_stack(
        (-np.log(w[0]), np.log(group), material[:, 0], beta * np.array(energy))
    )
    radii = np.column_stack(
        (
            -np.log1p(-error[0] / w[0]),
            -np.log1p(-group_error / group),
            entropy_error,
            beta * np.array(uncertainty),
        )
    )
    bounds = meter.shift_interval(
        w[2] - w[0], error[2] + error[0], scores - scores.min(0), radii
    )
    return dict(
        zip(
            (
                "sharp_records",
                "conditioning_on_ball",
                "conditional_material",
                "beta_material_energy",
            ),
            bounds.tolist(),
            strict=True,
        )
    )
