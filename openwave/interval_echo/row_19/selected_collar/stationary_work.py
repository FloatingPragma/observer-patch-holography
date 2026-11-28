"""Read a contracted equal-volume ball after ordinary work and a finite free wait."""

from fractions import Fraction as F
from math import log, sqrt

from mpmath import iv
import numpy as np

from . import coherent_balls as coherent
from .phase_and_flux import edge_read

local = coherent.local
steady = coherent.steady
bounds = coherent.balls.bounds
N = coherent.N
SEED = 28122926
DOSE = F(1, 4)


def preparations(board, setup):
    _, cells, ports, gates, partner, deg, edges, marks = setup
    ops = steady.native_energies(deg, edges, marks)
    free = ops[0, 0]
    n = len(ports)
    init = (np.ones((n, n), dtype=object), n)
    sign = np.ones(n, int)
    sign[gates[0]] = -1
    deposit = np.zeros((n, n), dtype=object)
    deposit[partner[0], partner[0]] = 1

    def floats(z):
        m, den = z
        return np.vectorize(lambda v: float(F(v, den)))(m)

    ordinary = floats(free.pinch((deposit, 1)))
    ref = {}
    work = {}
    for f, op in ops.items():
        num, den = free.pinch(op.pinch(init))
        assert np.array_equal(free.h @ num, num @ free.h)
        twisted = sign[:, None] * num * sign[None, :]
        out = free.pinch((twisted, den))
        # Exact preservation of supplied kinetic work under the last free wait.
        assert np.trace(free.h @ out[0]) * den == np.trace(free.h @ twisted) * out[1]
        ref[f], work[f] = floats((num, den)), floats(out)
    fwd, rev, _ = steady.parity(board, cells)
    arms = local.inputs("reference_arrival")["arms"]
    for instrument, base, shared in (
        ("serial", 25 / 36, 101 / 144),
        (
            "reference",
            arms["fresh"]["collision_acceptance"],
            arms["shared"]["collision_acceptance"],
        ),
    ):
        assert 0 < shared - base <= 1 / 144 + 1e-16
        for direction, law in (("forward", fwd), ("reversed", rev)):

            def mix(table, p):
                return (1 - p) * table[0, 0] + p * sum(
                    float(w) * table[f] for f, w in law.items()
                )

            yield (
                instrument,
                direction,
                {
                    "fresh": mix(ref, base),
                    "work": mix(work, base),
                    "ordinary": ordinary,
                    "overlap": mix(ref, shared),
                },
            )


def position_bounds(matrix, rng, half):
    k = coherent.draws(matrix.diagonal(), rng)
    return k, [iv.mpf([n / N - half, n / N + half]) for n in k]


def survey_probabilities(pos, sites, neighbors, beta):
    # Public radius, center-clock class and ordinary load inside the ball.
    probabilities = []
    for j, ball in enumerate(sites):
        for slow in (0, 1):
            for hit in (0, 1):
                probabilities.append(
                    (beta if j else 1 - beta)
                    * sum(
                        p
                        for i, p in enumerate(pos)
                        if (i in neighbors) == slow and (i in ball) == hit
                    )
                )
    return probabilities


def survey(pos, sites, neighbors, beta, rng):
    return coherent.draws(survey_probabilities(pos, sites, neighbors, beta), rng)


def held_differences(counts, baseline, volumes, areas, rate, half):
    # One clock calibration is shared by the two independent acquisitions.
    # Subtract counts before propagating its uncertainty. Hoeffding bounds
    # for two held means each contribute the full observable range times half.
    sharp = [
        (volumes[j] * rate ** (3 * slow), areas[j] * rate ** (2 * slow), hit)
        for j in (0, 1)
        for slow in (0, 1)
        for hit in (0, 1)
    ]
    result = []
    for values in zip(*sharp, strict=True):
        mean_difference = sum(
            iv.mpf(n - m) * v / N
            for n, m, v in zip(counts, baseline, values, strict=True)
        )
        span = max(float(iv.mpf(v).b) for v in values) - min(
            float(iv.mpf(v).a) for v in values
        )
        error = 2 * half * (span + 1e-12)
        result.append(mean_difference + iv.mpf([-error, error]))
    return result


def run():
    iv.dps = 45
    rng = np.random.default_rng(SEED)
    board = local.survey.bulk.federation(0)
    setup, _ = steady.symmetric_instrument(board)
    edges = setup[6]
    ports = setup[2]
    gates = setup[3]
    neighbors = {g: {a + b - g for a, b in edges if g in (a, b)} for g in gates}
    score = [
        int(i in neighbors[gates[0]] and i not in gates)
        - int(i in neighbors[gates[1]] and i not in gates)
        for i in range(len(ports))
    ]
    saved = local.inputs("restored_balls")["surveys"]
    previous = local.inputs("stationary_polarization")["mirror_matched"]
    wait = F(previous["combined_trace_error_upper"])
    # One more independent restored wait contracts trace norm. The previous
    # two-stage bound also bounds the error of that restored stage alone.
    wait_error = float(2 * wait) + 4e-13
    half = sqrt(log(2 * 100000 / 0.01) / (2 * N)) + wait_error + 1e-13
    clocks, rate = coherent.balls.proper_clock(rng)
    source = {}
    surveys = []
    regions = {}
    for instrument, label, states in preparations(board, setup):
        key_prefix = instrument + "_" + label
        pos = {}
        raw = {}
        for name, m in states.items():
            raw[name], pos[name] = position_bounds(m, rng, half)
        phase = np.ones(len(ports))
        phase[gates[0]] = -1
        raw_work = {}
        energies = {}
        for name, state in (
            ("fresh", states["fresh"]),
            ("phase", phase[:, None] * states["fresh"] * phase[None, :]),
        ):
            energies[name] = iv.mpf(0)
            raw_work[name] = {}
            for i, j in edges:
                if gates[0] in (i, j):
                    k, probs = edge_read(state, i, j, rng, half)
                    raw_work[name][f"{ports[i]}_{ports[j]}"] = k
                    energies[name] += probs[1]
        source[key_prefix] = {
            "work_seam_counts": raw_work,
            "kinetic_work": bounds(energies["phase"] - energies["fresh"]),
            "positions": raw,
            "available_changes": {
                k: bounds(
                    sum(
                        s * (v - u)
                        for s, v, u in zip(score, p, pos["fresh"], strict=True)
                    )
                )
                for k, p in pos.items()
                if k != "fresh"
            },
        }
        cache = {}
        for row in saved:
            if row["direction"] != label:
                continue
            g = ports.index(row["center"])
            nb = neighbors[g]
            bs = row["ball_sites"]
            V = row["volumes"]
            A = row["cuts"]
            for B in bs:
                key = tuple(B)
                if key not in cache:
                    r, ent, slopes = coherent.read_region(
                        B,
                        states,
                        edges,
                        rng,
                        wait_error,
                        {"work": 2, "ordinary": 2, "overlap": F(1, 72)},
                    )
                    cache[key] = ent, slopes
                    regions[key_prefix + ":" + str(B)] = r
            e1, e2 = (cache[tuple(B)] for B in bs)
            base = sum(pos["fresh"][i] for i in nb)
            c30 = 1 + (rate**3 - 1) * base
            q0 = states["fresh"].diagonal()
            fresh_counts = survey(q0, bs, nb, 0.5, rng)
            changes = {}
            for name in ("work", "ordinary", "overlap"):
                p = pos[name]
                a, _, coin = coherent.balls.analyze(
                    pos["fresh"], p, set(bs[0]), set(bs[1]), nb, V, A, rate
                )
                ds = (e1[1][name] + e2[1][name]) / 2 + coin * (e2[0] - e1[0])
                enclosed = [sum(p[i] - pos["fresh"][i] for i in B) for B in bs]
                dp = sum(enclosed) / 2 + coin * sum(
                    pos["fresh"][i] for i in set(bs[1]) - set(bs[0])
                )
                c3 = 1 + (rate**3 - 1) * (
                    (1 - float(DOSE)) * base + float(DOSE) * sum(p[i] for i in nb)
                )
                c2 = 1 + (rate**2 - 1) * (
                    (1 - float(DOSE)) * base + float(DOSE) * sum(p[i] for i in nb)
                )
                intended = ((V[0] + V[1]) * c30 / (2 * c3) - V[0]) / (V[1] - V[0])
                assert float(intended.a) > 0 and float(intended.b) < 1
                beta = float(intended.mid)
                q = (1 - float(DOSE)) * q0 + float(DOSE) * states[name].diagonal()
                counts = survey(q, bs, nb, beta, rng)
                dv, measured_area, measured_mass = held_differences(
                    counts, fresh_counts, V, A, rate, half
                )
                # Correct the held result to exactly the measured reference
                # volume, including the residual from calibrating the public coin.
                slope_a = (A[1] - A[0]) * c2 / ((V[1] - V[0]) * c3)
                da = measured_area - dv * slope_a
                extra = set(bs[1]) - set(bs[0])
                slope_m = sum(
                    (1 - float(DOSE)) * pos["fresh"][i] + float(DOSE) * p[i]
                    for i in extra
                ) / ((V[1] - V[0]) * c3)
                dm = measured_mass - dv * slope_m
                changes[name] = {
                    "area_tangent": bounds(a),
                    "entropy_tangent": bounds(ds),
                    "ordinary_mass_tangent": bounds(dp),
                    "required_gain": bounds(-ds / a) if 0 not in a else None,
                    "public_probability": beta,
                    "held_counts": counts,
                    "held_volume_difference": bounds(dv),
                    "held_area_difference": bounds(measured_area),
                    "held_mass_difference": bounds(measured_mass),
                    "held_equal_volume_area_change": bounds(da),
                    "held_equal_volume_mass_change": bounds(dm),
                }

            def det(x, y):
                a, b = [changes[k] for k in (x, y)]
                return bounds(
                    iv.mpf(a["area_tangent"]) * iv.mpf(b["entropy_tangent"])
                    - iv.mpf(b["area_tangent"]) * iv.mpf(a["entropy_tangent"])
                )

            ordinary_gain = iv.mpf(changes["ordinary"]["required_gain"])
            # Leave the independently measured ordinary localization calibration
            # fixed and ask how much additional entropy a new measured sector
            # would need to supply. Only positive ordinary gains are admissible.
            missing = (
                {
                    k: bounds(
                        -iv.mpf(v["entropy_tangent"])
                        - ordinary_gain * iv.mpf(v["area_tangent"])
                    )
                    for k, v in changes.items()
                    if k != "ordinary"
                }
                if ordinary_gain.a > 0
                else None
            )
            surveys.append(
                {
                    "center": row["center"],
                    "direction": label,
                    "instrument": instrument,
                    "radii": row["hop_radii"],
                    "fresh_counts": fresh_counts,
                    "changes": changes,
                    "additional_entropy_at_localized_ordinary_calibration": missing,
                    "work_ordinary_determinant": det("work", "ordinary"),
                    "overlap_ordinary_determinant": det("overlap", "ordinary"),
                }
            )
    return {
        "seed": SEED,
        "starts": N,
        "dose": str(DOSE),
        "finite_wait_trace_error_upper": wait_error,
        "clock_counts": clocks,
        "clock_ratio": bounds(rate),
        "regions": regions,
        "source": source,
        "surveys": surveys,
        "inputs": steady.exp.save_inputs(
            "STATIONARY_WORK.md",
            "phase_and_flux.json",
            "reference_arrival.json",
            "restored_balls.json",
            "stationary_polarization.json",
        ),
    }


if __name__ == "__main__":
    z = run()
    steady.exp.emit("stationary_work", z)
    for r in z["surveys"]:
        print(
            r["instrument"],
            r["direction"],
            r["center"],
            r["radii"],
            [
                (v["held_equal_volume_area_change"], v["held_equal_volume_mass_change"])
                for k, v in r["changes"].items()
                if k != "overlap"
            ],
        )
