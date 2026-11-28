"""Follow an unreinitialized source under dense and sparse repeated service."""

import argparse
import json
from hashlib import sha256

import numpy as np
from scipy.sparse.linalg import expm_multiply

from . import observation_density as obs

pos, c, lab, reg, iv = obs.pos, obs.c, obs.lab, obs.reg, obs.iv
HERE = pos.HERE
SEED = 28112106
EXPOSURE = 2**56


def superoperator(kraus):
    k = kraus.reshape(-1, 16)
    gram = (k.T @ k.conj()).reshape(4, 4, 4, 4)
    out = np.einsum("aibj,kij->kab", gram, pos.AXES)
    return np.einsum("lab,kba->lk", pos.AXES, out).real


def flights(ports):
    native = lab.six.survey.native
    board = native.iso.pairing.formation.federation(1)
    region = native.iso.pairing.formation.region(board)
    bank = np.zeros((2, 2, 12, 16, 16))
    for b in (0, 1):
        for sign in (-1, 1):
            h = native.generator(board, region, sign, bool(b)) * (-0.25j)
            for home in ports:
                start = np.zeros((h.shape[0], 4))
                start[4 * home : 4 * home + 4] = np.eye(4)
                tube = expm_multiply(
                    h, start, start=0, stop=24, num=25, traceA=h.diagonal().sum()
                ).reshape(25, -1, 4, 4)
                for t in range(1, 13):
                    free = tube[2 * t]
                    click = tube[t] @ tube[t, home]
                    bank[0, b, t - 1] += superoperator(free) / (2 * len(ports))
                    bank[1, b, t - 1] += (
                        superoperator(click) + superoperator(free - click)
                    ) / (2 * len(ports))
    np.testing.assert_allclose(
        bank[:, :, :, 0, :],
        np.broadcast_to(np.eye(16)[0], bank[:, :, :, 0, :].shape),
        atol=2e-12,
    )
    return bank


def schedules(record):
    ticks = np.zeros(12, dtype=object)
    ticks[1:] = np.array(record["times"]["reference"], dtype=object).sum(0)
    out = {"uniform": [1] * 12}
    for cap in (4, 8, 12):
        hist = ticks.copy()
        hist[cap - 1] = sum(hist[cap - 1 :])
        hist[cap:] = 0
        out[str(cap)] = list(hist)
    return out


def calibrate(M, rng):
    # A known parity can be selected on kept rails or heralded on dephased
    # rails. Each Pauli eigenspace mixture uses local product preparations.
    inputs = [np.r_[0.5, np.zeros(15)]]
    for j in range(15):
        for s in (-1, 1):
            z = inputs[0].copy()
            z[j + 1] = s / 2
            inputs.append(z)
    values = (np.asarray(inputs) @ M.T)[:, 1:]
    hits = rng.binomial(
        2**52, np.clip(0.5 + values, 0, 1), size=(EXPOSURE // 2**52, *values.shape)
    ).sum(0)
    estimates = hits / EXPOSURE - 0.5
    A = (estimates[2::2] - estimates[1::2]).T
    offset = estimates[0]
    delta = 0.5 * np.sqrt(2 * np.log(2 * 100000 / 0.01) / EXPOSURE) + 2e-12
    eps = 30 * delta + 1e-12
    inverse = np.linalg.inv(np.eye(15) - A)
    defect = (
        np.linalg.norm(np.eye(15) - inverse @ (np.eye(15) - A))
        + np.linalg.norm(inverse) * eps
        + 2e-11
    )
    out = {
        "plus_outcomes": hits.tolist(),
        "traceless_operator_error": eps,
        "fixed_inverse_defect": float(defect),
    }
    if defect < 1:
        z = inverse @ offset
        residual = np.linalg.norm(offset - (np.eye(15) - A) @ z) + 2e-12
        out["stationary_coordinate_radius"] = float(
            np.linalg.norm(inverse)
            * (np.sqrt(15) * delta + eps * np.linalg.norm(z) + residual)
            / (1 - defect)
        )
    # On Hermitian trace-zero 4x4 matrices every power has HS operator
    # norm <=2 by trace-distance contractivity and rank four.
    power = A.copy()
    error = eps
    for j in range(21):
        length = 2**j
        upper = np.linalg.norm(power, ord=2) + error + 1e-10
        if upper < 0.9:
            out["power"] = length
            out["power_norm_bound"] = float(upper)
            blocks = max(1, int(np.ceil(np.log(1e-11 / 2) / np.log(upper))))
            out["rounds_for_stationary_read"] = blocks * length
            return out
        error = (2 + np.linalg.norm(power, ord=2) + 1e-10) * error + 1e-11
        power = power @ power
    return out


def read_area(means, hi, lo, lam, eta):
    a, b = (np.asarray(means[k], dtype=object) for k in (hi, lo))
    dr = a - b
    box = [iv.mpf([min(u.a, v.a), max(u.b, v.b)]) for u, v in zip(a, b, strict=True)]
    A, V = pos.geometry_jet(box)
    # Straight integration can lose precision for long paths. Intersect
    # with independently bounded endpoints before reporting the sign.
    Aa, Va = c.dose.geometry(a)
    Ab, Vb = c.dose.geometry(b)

    def intersect(p, q):
        assert p.a <= q.b and q.a <= p.b
        return iv.mpf([max(p.a, q.a), min(p.b, q.b)])

    area = intersect(A.g @ dr, Aa - Ab)
    vol = intersect(V.g @ dr, Va - Vb)
    adjusted = intersect((V.g * lam - A.g * eta) @ dr, lam * vol - eta * area)
    return {
        k: lab.six.bounds(v)
        for k, v in (("area", area), ("volume", vol), ("adjusted_nats", adjusted))
    }


def run():
    iv.dps = 35
    rng = np.random.default_rng(SEED)
    _, sources, *_ = c.trials()
    x = sources["reference"]
    rec = json.loads((HERE / "finite_dose_curvature.json").read_bytes())
    prior = json.loads((HERE / "link_source_contraction.json").read_bytes())
    eta = iv.mpf([23.1042361, 23.1042363])
    lam = iv.mpf(prior["ordinary_volume_multiplier"])
    ports, bins = lab.six.survey.return_bins()
    ids = [i for i, v in enumerate(lab.six.survey.ROUTES) if max(v) < 4]
    bins = bins[ids]
    bank = flights(ports)
    ref_tape = lab.pairs(x, rng, False)
    ref = reg.src.observed(ref_tape, lab.SHOTS, False)
    out = {}
    for name, counts in schedules(rec).items():
        weights = np.asarray(counts, dtype=float) / sum(counts)
        matrices = np.einsum("abtkj,t->abkj", bank, weights)
        calibration = [[calibrate(M, rng) for M in arm] for arm in matrices]
        setting = {"clock_counts": counts, "process_calibration": calibration}
        if all("rounds_for_stationary_read" in j for arm in calibration for j in arm):
            # Use a common multiple of all certified powers and exceed the
            # number of rounds required by each sector and arm.
            period = max(j["power"] for arm in calibration for j in arm)
            steps = max(
                j["rounds_for_stationary_read"] for arm in calibration for j in arm
            )
            rounds = ((steps + period - 1) // period) * period
            result = {}
            for label, arm in zip(("sparse", "dense"), matrices, strict=True):
                terminal = []
                for R, M in zip(x, arm, strict=True):
                    coeff = np.einsum("kij,ji->k", pos.AXES, R).real
                    M = M.copy()
                    M[0] = np.eye(16)[0]
                    z = np.linalg.matrix_power(M, rounds) @ coeff
                    terminal.append(np.tensordot(z, pos.AXES, axes=1))
                result[label] = np.array(terminal)
            heads = {k: lab.pairs(z, rng, False) for k, z in result.items()}
            src = {k: reg.src.observed(z, lab.SHOTS, False) for k, z in heads.items()}
            # Transfer terminal head intervals to the stationary limit using
            # the certified 1e-11 HS error; any Pauli coordinate changes by <=2e-11.
            stable = {k: (z, e + 2e-11) for k, (z, e) in src.items()}
            scores = {
                k: lab.six.bounds(
                    obs.price_contrast(
                        refval, stable["sparse"], stable["dense"], "connected"
                    )
                )
                for k, refval in (
                    ("ordinary_reference", ref),
                    ("dense_reference", stable["dense"]),
                )
            }
            banks = {
                **result,
                **{
                    k + "_pulse": (1 - obs.DOSE) * x + obs.DOSE * z
                    for k, z in result.items()
                },
            }
            time = {
                k: [
                    lab.counts(np.r_[p, 1 - p.sum()], rng)
                    for p in c.station.pair(bins, z)
                ]
                for k, z in banks.items()
            }
            means, _ = reg.timing_slopes(time, lab.SHOTS, (("sparse", "dense", 1),))
            # Each capped clock has operator norm at most three. Even the
            # larger two-sector bound on its limiting error fits this pad.
            pad = iv.mpf([-2e-10, 2e-10])
            means = {k: [v + pad for v in row] for k, row in means.items()}
            setting.update(
                stationary_clock_allowance=2e-10,
                rounds=rounds,
                terminal_heads=heads,
                terminal_radar=time,
                stationary_connected=scores,
                terminal_geometry=read_area(means, "sparse", "dense", lam, eta),
                pulse_geometry=read_area(
                    means, "sparse_pulse", "dense_pulse", lam, eta
                ),
            )
        out[name] = setting
    return {
        "seed": SEED,
        "head_and_clock_starts": lab.SHOTS,
        "process_setting_starts": EXPOSURE,
        "reference_head_counts": ref_tape,
        "schedules": out,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in (
                "PERSISTENT_DENSITY_PROTOCOL.md",
                "finite_dose_curvature.json",
                "link_source_contraction.json",
                "finite_pairing_dose.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    path = HERE / "persistent_density.json"
    if cli.parse_args().write:
        path.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    print(
        {
            k: [v.get("rounds"), v.get("stationary_connected"), v.get("pulse_geometry")]
            for k, v in result["schedules"].items()
        }
    )
