"""Momentum and recoil readings on one conjugate charged band at two lengths."""

import argparse
import json
from math import log, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.stats import beta

from .apparatus import Guide, MOTIFS, WIDTHS, geometry, certify_chart
from ..conjugate_flight.measure import meter

HERE = Path(__file__).resolve().parent
SHOTS = 2**26
ENERGY_SHOTS = 2**38
CURRENT_SHOTS = 2**34
WAITS = tuple(2**k for k in range(0, 22, 3))
# Simultaneous bounds allow 10,000 appointments, categories and quadratures.
TAIL = 0.01 / 10000
PROB_ERROR = sqrt(log(2 / TAIL) / (2 * SHOTS))
PHASE_ERROR = np.arcsin(2 * np.sqrt(2) * PROB_ERROR)


def interval(count, n):
    return np.array(
        [
            beta.ppf(TAIL / 2, count, n - count + 1) if count else 0,
            beta.ppf(1 - TAIL / 2, count + 1, n - count) if count < n else 1,
        ]
    )


def roll(prob, n, rng):
    return int(rng.binomial(n, np.clip(prob, 0, 1)))


def energy(g, y, rng, shots=ENERGY_SHOTS):
    # H=sum t|x-y|^2 + V|x|^2. All read weights are positive, retain every launch.
    banks = []
    for first in (0, 1):
        probs = {1: 0.0, 8: 0.0}
        for j in range(first, g.n - 1, 2):
            probs[int(2 * g.w[j])] += g.read_pair(y, j, j + 1)[1]
        banks.append(probs)
    banks.append({k: float(sum(abs(y[g.potential == k]) ** 2)) for k in (8, 16)})
    bounds = np.zeros(2)
    estimate = 0.0
    records = []
    for probs in banks:
        values = list(probs)
        p = list(probs.values())
        counts = rng.multinomial(shots, [*p, max(0.0, 1 - sum(p))])
        records.append(dict(values=values, counts=counts.tolist()))
        for w, v in zip(values, counts):
            bounds += w * interval(v, shots)
            estimate += w * v / shots
    return dict(
        launches=shots,
        banks=records,
        estimate=estimate,
        interval99=bounds.tolist(),
    )


def momentum(g, y, rng):
    p = g.period
    mid = roll(sum(abs(y[p:-p]) ** 2), SHOTS, rng)
    raw = 0.0
    records = []
    for parity in (0, 1):
        prob = np.zeros(5)
        for c in range(parity, g.n // p - 1, 2):
            mult = 1 if c in (0, g.n // p - 2) else 2
            for j in range(c * p, (c + 1) * p):
                plus, minus = g.read_pair(y, j, j + p)
                prob[2 + mult] += plus
                prob[2 - mult] += minus
        prob[2] = max(0.0, 1 - sum(prob))
        counts = rng.multinomial(SHOTS, prob / prob.sum())
        records.append(counts.tolist())
        raw += counts @ np.arange(-2, 3) / SHOTS
    pint = mid / SHOTS
    if pint < 0.25:
        return dict(accepted=False, interior_count=mid)
    c = raw / (4 * pint)
    # Two bank means on [-2,2], one interior probability on [0,1].
    candidates = [
        (raw + dr) / (4 * (pint + dp))
        for dr in (-8 * PROB_ERROR, 8 * PROB_ERROR)
        for dp in (-PROB_ERROR, PROB_ERROR)
    ]
    lo, hi = min(candidates), max(candidates)
    assert -1 < lo < hi < 1
    return dict(
        accepted=True,
        interior_count=mid,
        weighted_counts=records,
        phase=float(np.arccos(c)),
        phase_interval99=[float(np.arccos(hi)), float(np.arccos(lo))],
    )


def signature(g, r):
    # Same Bloch recurrence in every physical sublattice, independent of the energy fit.
    y = g.columns[:, r]
    p = g.period
    mid = np.vdot(y[p:-p], y[p:-p]).real
    c = (y[p:-p] @ (y[: -2 * p] + y[2 * p :])) / (2 * mid)
    return float(mid), float(c)


def selected(g):
    # One mode per period is the first lobe. Boundary modes fail a material width gate.
    return [r for r in range(g.n // g.period) if signature(g, r)[0] >= 0.25]


def ramsey(g, r0, r, rng, coarse):
    x, y = g.columns[:, r0], g.columns[:, r]
    state = g.prepare((x + y) / sqrt(2))
    det = [g.prepare((x + s * y) / sqrt(2)) for s in (1, -1j)]
    initial = g.columns.T @ state
    assert coarse[1] - coarse[0] < np.pi
    prev = sum(coarse) / 2
    previous_radius = (coarse[1] - coarse[0]) / 2
    rows = []
    for t in WAITS:
        # Shift by a common energy only to keep phases well conditioned.
        v = g.columns @ (initial * np.exp(-1j * (g.levels - g.levels[r0]) * t))
        probs = [abs(np.vdot(z, v)) ** 2 for z in det]
        # <(x-i y)/sqrt(2)|v> yields (1+i exp(-i delta t))/2.
        cnt = [roll(q, SHOTS, rng) for q in probs]
        arg = np.arctan2(2 * cnt[1] / SHOTS - 1, 2 * cnt[0] / SHOTS - 1)
        assert previous_radius * t + PHASE_ERROR < np.pi
        correction = round((prev * t - arg) / (2 * pi))
        prev = (arg + correction * 2 * pi) / t
        previous_radius = PHASE_ERROR / t
        rows.append(cnt)
    eps = PHASE_ERROR / WAITS[-1]
    return dict(counts=rows, delta=prev, interval99=[prev - eps, prev + eps])


def reciprocity(g, r0, r, sign, ports, rng, beat, beat_interval):
    p = g.period
    # Launch the opposite real phases at the reflected calibration address.
    y0, y1 = g.columns[:, r0], -g.columns[:, r]
    packet = g.prepare((y0 + y1) / sqrt(2))
    # Two material appointments in separate interior wells, both adjacent to a
    # written mark. Source reversal uses an independently sharp material state.
    banks = []
    local_indices = [p + p // 2, (g.n // p - 2) * p + p // 2]
    assert len(set(ports[j] // 12 for j in local_indices)) == 2
    for detector in local_indices:
        assert g.word[detector - 1] == "1"
        sharp = np.eye(1, g.n, detector).ravel()
        counts = []
        incoherent = roll(
            (abs(y0[detector]) ** 2 + abs(y1[detector]) ** 2) / 2, SHOTS, rng
        )
        for t in (0, 32768, 131072, 524288):
            evolve = np.exp(-1j * (g.levels - g.levels[r0]) * t)
            forward = g.columns @ ((g.columns.T @ packet) * evolve)
            back = g.columns @ ((g.columns.T @ sharp) * evolve)
            # Instrument A is the bare terminal load. B reverses the native
            # narrow source and reads a sharp port on a distinct broad launch.
            pa = abs(forward[detector]) ** 2
            pb = abs(np.vdot(packet, back)) ** 2
            np.testing.assert_allclose(pa, pb, atol=2e-12)
            a, b = roll(pa, SHOTS, rng), roll(pb, SHOTS, rng)
            assert abs(a - b) / SHOTS < 2 * PROB_ERROR
            triggered = roll(pa, SHOTS, rng)
            assert abs(triggered - a) / SHOTS < 2 * PROB_ERROR
            counts.append(
                dict(
                    t=t,
                    local_load=a,
                    reciprocal_analyzer=b,
                    arrival_and_native_head=triggered,
                )
            )
        held_check = counts[2]["local_load"]  # fixed second positive wait
        assert max(
            interval(held_check, SHOTS)[0], interval(incoherent, SHOTS)[0]
        ) > min(interval(held_check, SHOTS)[1], interval(incoherent, SHOTS)[1])
        initial_count = counts[0]["local_load"]
        for row in counts[1:]:
            # Calibrate only the DC level (incoherent launch) and fringe at t=0.
            # Both wait frequencies and subsequent local readings are held out here.
            contrast = (initial_count - incoherent) / SHOTS
            row["predicted_local"] = incoherent / SHOTS + contrast * np.cos(
                beat * row["t"]
            )
            eps = (
                4 * PROB_ERROR
                + (abs(contrast) + 2 * PROB_ERROR)
                * row["t"]
                * (beat_interval[1] - beat_interval[0])
                / 2
            )
            row["prediction_and_count_bound99"] = eps
            assert abs(row["local_load"] / SHOTS - row["predicted_local"]) < eps
            assert (
                abs(row["reciprocal_analyzer"] / SHOTS - row["predicted_local"]) < eps
            )
        banks.append(
            dict(
                incoherent_source=incoherent,
                port=int(ports[detector]),
                owner=int(ports[detector] // 12),
                index=detector,
                signed_load=int(meter()[1 + sign] - 1),
                incoming_head=int(sign < 0),
                readings=counts,
            )
        )
    # A releasable pair of standing modes delivers opposite signed material
    # currents when its relative quarter phase is reversed. The controls use
    # the literal, native half-contact splitter at a physical cut.
    currents = []
    j = 2 * p - 1  # Counter at a physical period boundary.
    for direction in (-1, 1):
        src = g.prepare((y0 + direction * 1j * y1) / sqrt(2))
        left, right = src[j : j + 2]
        after = (left + 1j * right) / sqrt(2)
        p0 = abs(after) ** 2
        p1 = abs((right + 1j * left) / sqrt(2)) ** 2
        tallies = rng.multinomial(CURRENT_SHOTS, [p0, p1, 1 - p0 - p1])
        # Position-to-the-right current for real hopping -w.
        val = g.w[j] * (int(tallies[1]) - int(tallies[0])) / CURRENT_SHOTS
        lci, rci = (
            interval(tallies[0], CURRENT_SHOTS),
            interval(tallies[1], CURRENT_SHOTS),
        )
        bound = g.w[j] * np.array([rci[0] - lci[1], rci[1] - lci[0]])
        evolved = g.columns @ (
            (g.columns.T @ src) * np.exp(-1j * (g.levels - g.levels[r0]) * 8192)
        )
        region = roll(sum(abs(evolved[j + 1 :]) ** 2), SHOTS, rng)
        currents.append(
            dict(
                counts=tallies.tolist(),
                coordinate_current=val,
                signed_current=sign * val,
                coordinate_interval99=bound.tolist(),
                right_region_after_8192=region,
            )
        )
    # Time reversal relates the two quadrature launches. Their right-region
    # difference equals twice the integrated initial current times sinc(delta*T).
    t = 8192
    dt = (
        currents[1]["right_region_after_8192"] - currents[0]["right_region_after_8192"]
    ) / SHOTS
    sinc = np.sinc(beat * t / np.pi)
    flow = dt / (2 * t * sinc)
    flow_err = (
        2 * PROB_ERROR + abs(dt) * t * (beat_interval[1] - beat_interval[0])
    ) / (2 * t * abs(sinc))
    left, right = currents[1]["coordinate_interval99"]
    assert max(left, flow - flow_err) <= min(right, flow + flow_err)
    assert np.sign(currents[1]["coordinate_current"]) != np.sign(
        currents[0]["coordinate_current"]
    )
    assert left * right > 0 and abs(flow) > flow_err
    return dict(
        mode=r,
        receivers=banks,
        currents=currents,
        cut_port_indices=[j, j + 1],
        independent_timed_current=flow,
        timed_radius99=flow_err,
        flow_shots=CURRENT_SHOTS,
    )


def acquire(n, p, sign, ports, board, rng):
    g = Guide(n, p)
    candidates = []
    for i in range(n // p):
        state = g.prepare(g.columns[:, i])
        read = momentum(g, state, rng)
        candidates.append(dict(index=i, **read))
    indices = [r["index"] for r in candidates if r["accepted"]]
    assert indices == selected(g)
    reference = indices[0]
    anchor = energy(g, g.prepare(g.columns[:, reference]), rng)
    assert anchor["interval99"][0] <= g.levels[reference] <= anchor["interval99"][1]
    modes = []
    for i in indices:
        m = candidates[i]
        coarse_energy = energy(
            g,
            g.prepare(g.columns[:, i]),
            np.random.default_rng(63210 + n * 100 + sign * 17 + i),
            2**20,
        )
        coarse = (
            coarse_energy["interval99"][0] - anchor["interval99"][1],
            coarse_energy["interval99"][1] - anchor["interval99"][0],
        )
        delta = (
            dict(delta=0.0, interval99=[0.0, 0.0], counts=[])
            if i == reference
            else ramsey(g, reference, i, rng, coarse)
        )
        assert (
            delta["interval99"][0] - 1e-13
            <= g.levels[i] - g.levels[reference]
            <= delta["interval99"][1] + 1e-13
        )
        modes.append(
            dict(
                index=i,
                phase=m["phase"],
                phase_interval99=m["phase_interval99"],
                recoil=delta,
                coarse_energy=coarse_energy,
            )
        )
    # Faster signed-current control uses the two most separated tested momenta.
    cross = reciprocity(
        g,
        reference,
        indices[-1],
        sign,
        ports,
        rng,
        modes[-1]["recoil"]["delta"],
        modes[-1]["recoil"]["interval99"],
    )
    chart = certify_chart(ports, board, g, sign)
    return dict(
        n=n,
        period=p,
        charge=sign,
        anchor=anchor,
        momentum=candidates,
        modes=modes,
        source_receiver=cross,
        reference_index=reference,
        charged_chart=chart,
    )


def fit(train, held):
    # Train only on the first two interior modes at the shorter physical width.
    a, b = train["modes"][:2]
    ea = train["anchor"]["estimate"]
    db = b["recoil"]["delta"]

    def predictions(theta, ref, pa, pb, d, erest, kind):
        if kind == "quadratic":
            return d * (theta**2 - ref**2) / (pb**2 - pa**2)
        if kind == "lattice":
            return d * (np.cos(ref) - np.cos(theta)) / (np.cos(pa) - np.cos(pb))
        if kind == "lattice_shell":
            theta, ref, pa, pb = 2 * np.sin(np.array([theta, ref, pa, pb]) / 2)
        # The first two measured energies set c and the zero-phase mass.
        c2 = (2 * erest * d + d * d) / (pb**2 - pa**2)
        mass2 = erest**2 - c2 * pa**2
        assert mass2 > 0
        return (
            c2
            * (theta**2 - ref**2)
            / (sqrt(mass2 + c2 * theta**2) + sqrt(mass2 + c2 * ref**2))
        )

    out = {}
    for model in ("quadratic", "relativistic", "lattice", "lattice_shell"):
        samples = []
        for source in (train, held):
            ref = source["modes"][0]
            for row in source["modes"][1:]:
                if source is train and row is b:
                    continue
                expected = predictions(
                    row["phase"], ref["phase"], a["phase"], b["phase"], db, ea, model
                )
                # Outward natural interval extension on positive, disjoint phases.
                tl, th = row["phase_interval99"]
                rl, rh = ref["phase_interval99"]
                al, ah = a["phase_interval99"]
                bl, bh = b["phase_interval99"]
                dl, dh = b["recoil"]["interval99"]
                el, eh = train["anchor"]["interval99"]
                if model == "lattice_shell":
                    tl, th, rl, rh, al, ah, bl, bh = 2 * np.sin(
                        np.array([tl, th, rl, rh, al, ah, bl, bh]) / 2
                    )
                numerator = (tl**2 - rh**2, th**2 - rl**2)
                den = (bl**2 - ah**2, bh**2 - al**2)
                assert min(*numerator, *den, dl) > 0
                if model == "quadratic":
                    lower, upper = (
                        dl * numerator[0] / den[1],
                        dh * numerator[1] / den[0],
                    )
                elif model == "lattice":
                    lower = dl * (np.cos(rh) - np.cos(tl)) / (np.cos(al) - np.cos(bh))
                    upper = dh * (np.cos(rl) - np.cos(th)) / (np.cos(ah) - np.cos(bl))
                else:
                    sl, sh = (
                        (2 * el * dl + dl**2) / den[1],
                        (2 * eh * dh + dh**2) / den[0],
                    )
                    ml, mh = el**2 - sh * ah**2, eh**2 - sl * al**2
                    lower = (
                        sl
                        * numerator[0]
                        / (sqrt(mh + sh * th**2) + sqrt(mh + sh * rh**2))
                    )
                    upper = (
                        sh
                        * numerator[1]
                        / (sqrt(ml + sl * tl**2) + sqrt(ml + sl * rl**2))
                    )
                measured = row["recoil"]["interval99"]
                samples.append(
                    dict(
                        length=source["n"],
                        charge=source["charge"],
                        independent_dataset=source is held,
                        mode=row["index"],
                        predicted=expected,
                        predicted_interval=[lower, upper],
                        observed=row["recoil"]["delta"],
                        residual=row["recoil"]["delta"] - expected,
                        compatible=bool(
                            max(lower, measured[0]) <= min(upper, measured[1])
                        ),
                    )
                )
        out[model] = dict(
            compatible=sum(s["compatible"] for s in samples),
            appointments=len(samples),
            held=samples,
        )
    c2 = (2 * ea * db + db * db) / (b["phase"] ** 2 - a["phase"] ** 2)
    out["relativistic_calibration"] = dict(
        rest_native=sqrt(ea**2 - c2 * a["phase"] ** 2),
        fitted_speed_in_port_spacings_per_native_time=train["period"] * sqrt(c2),
    )
    return out


def run():
    ports, board = geometry()
    assert meter() == [0, 1, 2]
    # Exposure and geometric screen depend on carrier lengths, before spectroscopy.
    # Of the eligible motif family, select the lowest interior line on 96 ports.
    screen = []
    for period in MOTIFS:
        guide = Guide(WIDTHS[0], period)
        i = selected(guide)[0]
        z = guide.prepare(guide.columns[:, i])
        lottery = np.random.default_rng(952822 + period)
        counter = energy(guide, z, lottery)
        interior = roll(sum(abs(z[period:-period]) ** 2), SHOTS, lottery)
        assert (
            interior > SHOTS / 4
            and counter["interval99"][0] < guide.levels[i] < counter["interval99"][1]
        )
        screen.append(
            dict(
                period=period,
                phase_interior=signature(guide, i)[0],
                material_interior_count=interior,
                energy_counter=counter,
                line=guide.levels[i],
            )
        )
    chosen = min(screen, key=lambda z: z["energy_counter"]["estimate"])
    assert all(
        chosen["energy_counter"]["interval99"][1]
        < row["energy_counter"]["interval99"][0]
        for row in screen
        if row is not chosen
    )
    p = chosen["period"]
    rows = [
        acquire(n, p, sign, ports[:n], board, np.random.default_rng(830021 + n + sign))
        for sign in (-1, 1)
        for n in WIDTHS
    ]
    fits = {str(sign): fit(rows[j], rows[j + 1]) for sign, j in [(-1, 0), (1, 2)]}
    # Train once on the negative 96-port source, including a blind change of sign.
    joint = {}
    for other in rows[1:]:
        answer = fit(rows[0], other)
        for method in ("quadratic", "relativistic", "lattice", "lattice_shell"):
            if method not in joint:
                joint[method] = [
                    g for g in answer[method]["held"] if not g["independent_dataset"]
                ]
            joint[method].extend(
                g for g in answer[method]["held"] if g["independent_dataset"]
            )
    return dict(
        reference_clock="one native time",
        position_unit="one consecutive physical route port",
        phase_unit="radians per repeated motif",
        sites=[int(p) for p in ports],
        shots=SHOTS,
        waits=WAITS,
        screen=screen,
        selected_period=p,
        observations=rows,
        fits=fits,
        negative_training_joint_holdouts=joint,
    )


if __name__ == "__main__":
    arg = argparse.ArgumentParser()
    arg.add_argument("--write", action="store_true")
    out = json.dumps(run(), indent=2, default=lambda v: v.item()) + "\n"
    if arg.parse_args().write:
        (HERE / "records.json").write_text(out)
    else:
        print(out, end="")
