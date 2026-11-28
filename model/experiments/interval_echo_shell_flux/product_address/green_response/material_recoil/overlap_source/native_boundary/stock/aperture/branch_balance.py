"""Identify the entropy beyond the ordinary area and volume of native message channels."""

from fractions import Fraction as Q

import numpy as np

from . import moving_service as lab

N = 2**34


def classify():
    f, h = lab.gates()
    flag = np.zeros(99, int)
    for j in range(99):
        full = f[:, j]
        if np.count_nonzero(full) == 4:
            flag[j] = 1
            assert full[j] == 1 and h[j, j] == 5
            assert sorted(full[full != 0]) == [1] * 4
            assert sorted(h[h[:, j] != 0, j]) == [1, 1, 1, 5]
        else:
            assert full.max() == 4 and (
                full[j] == 4 or sorted(h[h[:, j] != 0, j]) == [4, 4]
            )
        # On two fresh copies with the same sharp start, require two different
        # registered outcomes which both differ from that start.
        distinct = sum(
            int(full[x]) * int(full[y])
            for x in np.flatnonzero(full)
            for y in np.flatnonzero(full)
            if x != y and x != j and y != j
        )
        assert distinct == 6 * flag[j]
    return flag


def pair_reads(pre, rng):
    # Directly sample both physical full pulses; mark a trial only by
    # comparing its two observed outcomes and common preparatory record.
    f, _ = lab.gates()
    count = 0
    for j, n in enumerate(rng.multinomial(N, np.asarray(pre, float))):
        if not n:
            continue
        p = f[:, j] / 4
        first = rng.multinomial(int(n), p)
        for x in np.flatnonzero(first):
            if x != j:
                second = rng.multinomial(int(first[x]), p)
                count += int(first[x] - second[x] - second[j])
    return count


def matched_meter(reads, weights, kind, rng, trials):
    # Choose an address proportional to its multiplicity inside the volume
    # or cut. All observations use fresh copies of the prescribed mixture.
    pre = weights @ reads / int(weights.sum())
    draws = rng.multinomial(trials, np.asarray(pre, float))
    full, half = lab.gates()
    counts = np.zeros(8, dtype=np.int64)
    for j, n in enumerate(draws):
        if not n:
            continue
        if kind == "branch":
            first = rng.multinomial(int(n), full[:, j] / 4)
            for x in np.flatnonzero(first):
                if x != j:
                    second = rng.multinomial(int(first[x]), full[:, j] / 4)
                    counts[0] += int(first[x] - second[x] - second[j])
        else:
            kernel = half if kind == "history" else full
            out = rng.multinomial(int(n), kernel[:, j] / kernel[:, j].sum())
            if kind == "history":
                keep = kernel[:, j] > 0
                np.add.at(counts, kernel[keep, j] - 1, out[keep])
            else:
                counts[0] += int(n - out[j])  # Wait is independently deterministic.
    mean = lab.bits(counts / trials) if kind == "history" else counts[0] / trials
    return counts.tolist(), float(weights.sum() * mean)


def matched_experiment(arms, interior, cut, mix, gamma, n, seed):
    rng = np.random.default_rng(seed)
    old = lab.before(arms["original"])
    donor = lab.before(arms["ordinary"])
    source = lab.before(arms["source"])
    prep = [(1 - mix) * old + mix * donor, source]
    rows = {}
    for tag, w, range_ in (
        ("area", cut, 1),
        ("volume", 2 * interior + cut, 1),
        ("history", interior, 3),
        ("branch", interior, 1),
    ):
        values = [matched_meter(arm, w, tag, rng, n) for arm in prep]
        diff = values[1][1] - values[0][1]
        err = float(w.sum() * range_ * np.sqrt(np.log(2 * 8 / 0.01) / n))
        rows[tag] = dict(counts=[v[0] for v in values], change=diff, halfwidth=err)
    # Forecast the boundary with history, volume and the duplicate-pulse
    # branch read taken on three fresh banks. Hold out the boundary bank.
    coeff = dict(volume=1, history=-2, branch=2 * gamma * 8 / 3)
    forecast = sum(v * rows[k]["change"] for k, v in coeff.items())
    ferr = sum(abs(v) * rows[k]["halfwidth"] for k, v in coeff.items())
    area = rows["area"]
    assert abs(rows["volume"]["change"]) <= rows["volume"]["halfwidth"]
    assert area["change"] + area["halfwidth"] < 0
    assert abs(forecast - area["change"]) <= ferr + area["halfwidth"]
    if n == 2**47:
        assert forecast + ferr < 0
    return dict(shots=n, meters=rows, forecast=[forecast - ferr, forecast + ferr])


def isolate_branch(arms, i, cut, gamma, exact, certs):
    # Match both measured incident volume and interior history with independent
    # ordinary single and double additions. Both physical arms are mixtures.
    v = 2 * i + cut
    response = {}
    for tag in ("source", "ordinary", "single"):
        change = certs["original"]["signals"][:, 2] - certs[tag]["signals"][:, 2]
        response[tag] = (
            v @ change,
            i @ change,
            i @ (exact[tag] - exact["original"]),
            cut @ change,
        )
    s, d, u = (response[k] for k in ("source", "ordinary", "single"))
    # First fix the single weight to a rational approximation of the history
    # match. Then fix the double weight by the exact incident-volume equation.
    matrix = np.array(
        [
            [float(d[0]), float(u[0])],
            [float(d[1]) + gamma * float(d[2]), float(u[1]) + gamma * float(u[2])],
        ]
    )
    b = Q(
        round(
            np.linalg.solve(matrix, [float(s[0]), float(s[1]) + gamma * float(s[2])])[1]
            * 1e9
        ),
        10**9,
    )
    a = (s[0] - b * u[0]) / d[0]
    assert a < 0 < b and a + b < 1 and d[2] == 0
    norm = 1 - a
    # Orient the intervention by increasing the measured four-way frequency.
    weights = (
        {"original": (1 - a - b) / norm, "single": b / norm},
        {"source": 1 / norm, "ordinary": -a / norm},
    )
    assert all(sum(p.values()) == 1 and min(p.values()) >= 0 for p in weights)
    delta = [(a * x + b * y - z) / norm for z, x, y in zip(s, d, u, strict=True)]
    vol, internal, branch, area = delta
    assert vol == 0 and branch > 0 and area > 0
    from mpmath import iv

    iv.dps = 40

    def exact_iv(x):
        return iv.mpf(int(x.numerator)) / int(x.denominator)

    bits = exact_iv(internal) + (
        iv.mpf(9) / 4 - iv.mpf(5) / 8 * iv.log(5) / iv.log(2)
    ) * exact_iv(branch)
    history = [
        float(np.nextafter(float(bits.a), -np.inf)),
        float(np.nextafter(float(bits.b), np.inf)),
    ]
    assert max(map(abs, history)) < 1e-10
    read = [sum((lab.before(arms[k]) * v for k, v in p.items())) for p in weights]
    rng = np.random.default_rng(550340)
    n = 2**40
    # The only enrolled ambiguous seam is the third one. Sample it directly
    # and verify the boundary on separate copies selected by their cut address.
    assert i[1:].tolist() == [0, 0, 1]
    local = np.array([0, 0, 0, 1])
    out = {}
    for tag, where in (("branch", local), ("area", cut)):
        values = [matched_meter(r, where, tag, rng, n) for r in read]
        meas = values[0][1] - values[1][1]
        err = float(where.sum() * np.sqrt(np.log(4 / 0.01) / n))
        out[tag] = dict(counts=[v[0] for v in values], change=meas, halfwidth=err)
    prediction = 2 * gamma * 8 / 3 * out["branch"]["change"]
    error = 2 * gamma * 8 / 3 * out["branch"]["halfwidth"] + 2 * max(map(abs, history))
    observed = out["area"]
    assert prediction - error > 0 and observed["change"] - observed["halfwidth"] > 0
    assert abs(prediction - observed["change"]) <= error + observed["halfwidth"]
    return dict(
        arm_weights=[{k: str(v) for k, v in p.items()} for p in weights],
        shots=n,
        volume_change=str(vol),
        history_change_bits=history,
        branch_change=str(branch),
        area_change=str(area),
        meters=out,
        predicted_area=[prediction - error, prediction + error],
    )


def run(arms, numbers, edges, types, hops, certs, meters, radii):
    flag = classify()
    gamma = 2.25 - 5 / 8 * np.log2(5)
    # H(5/8,1/8,1/8,1/8)-3/4. All other columns have h=p(change).
    bulk, _, tv = lab.readers()
    np.testing.assert_allclose(
        lab.bits(bulk.T / 8) - (1 - tv[2] / 16), gamma * flag, atol=2e-15
    )
    exact = {k: lab.before(v) @ flag for k, v in arms.items()}
    rng = np.random.default_rng(790133)
    banks = {
        k: [pair_reads(row, rng) for row in lab.before(v)] for k, v in arms.items()
    }
    reports = {}
    for k in ("source", "ordinary", "single"):
        change = exact[k] - exact["original"]
        gate = (np.array(banks[k]) - banks["original"]) / N
        # Three arms times four types, with a simultaneous 99% difference band.
        err = np.sqrt(np.log(2 * 12 / 0.01) / N)
        measured = gate * 8 / 3 * gamma
        error = err * 8 / 3 * gamma
        expect = gamma * np.array(change, float)
        assert np.max(abs(expect - measured)) <= error
        responses = meters[k]
        # Independently obtained half-pulse history and full-pulse return.
        sums = [row[0] + row[3] for row in responses["held"]]
        residual_error = sum(responses["measured_halfwidth"][i] for i in (0, 3))
        assert np.all(abs(np.asarray(sums) - expect) < residual_error)
        reports[k] = dict(
            four_way_frequency_change=list(map(str, change)),
            pulse_pair_counts=banks[k],
            residual_bits_per_edge=expect.tolist(),
            predicted_by_pair=measured.tolist(),
            prediction_halfwidth=error,
            independent_history_plus_return=sums,
            independent_halfwidth=residual_error,
            whole_sweep_residual_bits=float(numbers @ expect),
        )
    # A voltage-free, native binary message is gate / wait. The jointly
    # readable endpoint and head determine whether the original word changed.
    # A station counts its incident message advantages, a boundary counts
    # each crossing advantage once. V=2*interior_advantage+A on every region.
    survey = {}
    ordinary_rows = []
    p0 = 1 - certs["original"]["signals"][:, 2]
    ps = 1 - certs["source"]["signals"][:, 2] - p0
    po = 1 - certs["ordinary"]["signals"][:, 2] - p0
    chosen = None
    for r in radii:
        positive = 0
        signs = np.zeros(3, int)
        enriched = np.zeros(3, int)
        max_needed = 0.0
        for station, h in enumerate(hops):
            ends = (h < r)[edges]
            i = np.bincount(types[ends.all(1)], minlength=4)
            c = np.bincount(types[ends[:, 0] ^ ends[:, 1]], minlength=4)
            # Calibrate using a blank ordinary double. Every contact has the
            # same measured nonzero dP(change); cancel this common factor.
            ordinary_rows.append(
                (int(c.sum()), -int(2 * i.sum() + c.sum()), int(-i.sum()))
            )
            delta = exact["source"] - exact["original"]
            positive += int(i @ delta > 0)
            need = 2 * gamma * float(i @ delta)
            max_needed = max(max_needed, need)
            v = 2 * i + c
            alpha = (v @ ps) / (v @ po)
            assert 0 <= alpha <= 1 and v @ (ps - alpha * po) == 0
            area_change = c @ (ps - alpha * po)
            idx = 1 + int(area_change > 0) - int(area_change < 0)
            signs[idx] += 1
            if i @ delta > 0:
                enriched[idx] += 1
            if chosen is None and i[3]:
                chosen = dict(
                    center=station,
                    radius=r,
                    ordinary_mixture=str(alpha),
                    interior=i.tolist(),
                    boundary=c.tolist(),
                    area_change=str(area_change),
                    branched_residual=gamma * float(i @ delta),
                    history_and_volume_matched=isolate_branch(
                        arms, i, c, gamma, exact, certs
                    ),
                    pilot=matched_experiment(arms, i, c, alpha, gamma, 2**41, 52083390),
                    independent=matched_experiment(
                        arms, i, c, alpha, gamma, 2**47, 62843222
                    ),
                )
        survey[str(r)] = dict(
            positive_regions=positive,
            largest_further_area_loss=max_needed,
            volume_matched_area_signs=signs.tolist(),
            positive_source_area_signs=enriched.tolist(),
        )
    first = ordinary_rows[0]
    second = next(
        row for row in ordinary_rows if first[0] * row[1] != row[0] * first[1]
    )
    assert all(2 * s == a + b for a, b, s in ordinary_rows)
    assert not np.any(exact["ordinary"] - exact["original"])
    assert np.all((exact["source"] - exact["original"])[1:] > 0)
    return dict(
        four_way_columns=np.flatnonzero(flag).tolist(),
        bits_per_four_way=gamma,
        fixed_four_way_pair_miss_after_ten=(5 / 8) ** 10,
        ideal_direct_pairs_per_arm=int(
            np.ceil(
                4
                * np.log(100)
                / (3 / 8 * float((exact["source"] - exact["original"])[3])) ** 2
            )
        ),
        paired_shots_per_arm_type=N,
        reference_pulse_pair_counts=banks["original"],
        ordinary_area_volume_equations=[first, second],
        surveys=survey,
        first_ball_enclosing_third=chosen,
        changes=reports,
    )
