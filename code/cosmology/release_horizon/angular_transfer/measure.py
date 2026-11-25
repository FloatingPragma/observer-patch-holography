"""Read fixed-frame angular contrast, polarization, loss and finite screen transfer."""

import argparse
import json
from math import log, sqrt

import numpy as np
from scipy.special import eval_legendre

if __package__:
    from .rays import Rays, HERE
    from .polarimeter import Polarimeter
else:
    from rays import Rays, HERE
    from polarimeter import Polarimeter

N = 65536
SCALE = 10**9


def angular(directions, record, arrivals, guard, proposals):
    bins = len(directions)
    dots = directions @ directions.T
    contrast = np.r_[np.full(12, 1 / 12), np.full(20, -1 / 20)]
    h6 = eval_legendre(6, dots[:12]).mean(axis=0)
    q = eval_legendre(2, dots[0])
    # Rows are incident directions. The Z analyzer labels are two accepted
    # transverse rails and one rejected rail. Keep original-denominator rates.
    z = record[:, :, :2, 0, :]
    p = z[:, :, :, :2].sum(axis=-1).mean(axis=2)
    np.testing.assert_allclose(p, (1 + dots**2) / 2, rtol=0, atol=2e-13)
    kernel = p * proposals[None, :]
    # In this bank the actual arriving owner is the directional detector.
    # Address aliases would invalidate the identity response.
    stations = np.array([b for _, b in arrivals])
    bypass = (stations[:, None] == stations[None, :]).astype(float)
    np.testing.assert_array_equal(bypass, np.eye(bins))
    stats = {}
    for name, transport in (
        ("ballistic", bypass),
        ("collision", kernel),
        ("half_bypass", guard[0] * bypass + guard[1] * kernel),
    ):
        survival = transport.sum(axis=1)
        assert np.ptp(survival) < 2e-13
        answer = dict(survival=round(float(survival.mean()), 12))
        for label, launch, read in (("six", contrast, h6), ("two", q / bins, q)):
            raw = launch @ transport @ read / (launch @ read)
            answer[label] = [
                round(float(raw), 12) + 0.0,
                round(float(raw / survival.mean()), 12) + 0.0,
            ]
        stats[name] = answer

    # Signed quadrupole response per outgoing proposal in its local frame.
    # The three analysis settings read Stokes Q,U,V by complementary rails.
    quad = np.tensordot(q / bins, record[:, :, :2].mean(axis=2), axes=1)
    stokes = np.column_stack(
        (quad[:, 0, :2].sum(axis=1), quad[:, :, 0] - quad[:, :, 1])
    )
    np.testing.assert_allclose(stokes[:, 0], q / 15, rtol=0, atol=2e-13)
    np.testing.assert_allclose(
        np.hypot(stokes[:, 1], stokes[:, 2]),
        (1 - dots[0] ** 2) / 10,
        rtol=0,
        atol=2e-13,
    )
    np.testing.assert_allclose(stokes[:, 3], 0, rtol=0, atol=2e-13)

    # Separate integer cohorts for each of the two incoherent input states.
    rng = np.random.default_rng(145116)
    counts = np.empty(z.shape, dtype=np.int64)
    for idx in np.ndindex(z.shape[:-1]):
        probabilities = z[idx]
        # At sharp endpoints the exactly impossible branch has zero probability.
        probabilities = np.clip(np.round(probabilities, 12), 0, 1)
        counts[idx] = rng.multinomial(N, probabilities / probabilities.sum())
    hits = counts[:, :, :, :2].sum(axis=-1)
    bound = sqrt(log(2 * hits.size / 0.005) / (2 * N))
    deviation = np.max(abs(hits / N - z[:, :, :, :2].sum(axis=-1)))
    assert deviation < bound and np.all(counts.sum(axis=-1) == N)
    estimate = contrast @ (hits.mean(axis=2) / (N * bins)) @ h6 / (contrast @ h6)
    squared_ranges = (
        (contrast @ contrast) * (h6 @ h6) / (2 * N * bins**2 * (contrast @ h6) ** 2)
    )
    width = sqrt(squared_ranges * log(2 / 0.005) / 2)
    assert abs(estimate - stats["collision"]["six"][0]) < width
    # Compare the *unnormalized* outgoing density across the two source orbits
    # in each analyzer, preserving loss and checking linear/circular components.
    full_difference = np.tensordot(contrast, record[:, :, :2].mean(axis=2), axes=1)
    assert np.max(abs(full_difference)) < 5e-13
    return dict(
        h6_by_orbit=[round(float(h6[0]), 12), round(float(h6[12]), 12)],
        source_power=[
            round(float(contrast @ eval_legendre(l, dots) @ contrast), 12) + 0.0
            for l in range(7)
        ],
        transfer=stats,
        bypass_weights=(np.rint(guard * 1e12).astype(np.int64) / 1e12).tolist(),
        quadrupole_stokes=(np.rint(stokes * 1e12).astype(np.int64) / 1e12).tolist(),
        originals_per_stratum=N,
        max_hit_error=round(float(deviation), 12),
        simultaneous_max_error=round(bound, 12),
        counted_raw_six_gain=round(float(estimate), 12),
        simultaneous_gain_halfwidth=round(width, 12),
    ), counts


def repeat_encounter(rays, frames, pol, proposals):
    """A secondary two-encounter bank with a tagged intermediate station."""
    n = len(frames)
    native_x = pol.givens(0, 1, np.pi / 2) @ pol.phase(1, np.pi)
    native_z = pol.phase(1, np.pi)
    flips = [np.eye(3), native_x, native_z, native_x @ native_z]
    # Each pair of local bases maps amplitudes across a completed direction
    # step. Its angular index is read at the corresponding receiving station.
    basis = np.stack(frames)
    jones = (basis.conj().transpose(0, 2, 1)[None, :, :, :] @ basis[:, None, :, :])[
        :, :, :2, :2
    ]
    emitted = (jones @ jones.conj().transpose(0, 1, 3, 2)) / 2
    depolarized = sum(f[:2, :2] @ emitted @ f[:2, :2].conj().T for f in flips) / 4
    np.testing.assert_allclose(
        depolarized,
        emitted.trace(axis1=2, axis2=3)[:, :, None, None] * np.eye(2) / 2,
        atol=2e-13,
    )
    q = eval_legendre(2, rays.directions @ rays.directions[0])
    contrast = np.r_[np.full(12, 1 / 12), np.full(20, -1 / 20)]
    h6 = eval_legendre(6, rays.directions[:12] @ rays.directions.T).mean(axis=0)
    stations = []
    for word in rays.words:
        middle = rays.destination(0, word)
        destinations = [rays.destination(middle, w) for w in rays.words]
        assert len(set(destinations)) == n
        stations.append(destinations)
    result = {}
    kernels = []
    for name, density in (
        ("polarization_retained", emitted),
        ("polarization_randomized", depolarized),
    ):
        raw = np.zeros((n, n))
        for middle in range(n):
            j = jones[middle]
            propagated = (
                j[None, ...]
                @ density[:, middle, None, :, :]
                @ j.conj().transpose(0, 2, 1)[None, ...]
            )
            raw += (
                proposals[middle]
                * proposals[None, :]
                * propagated.trace(axis1=2, axis2=3).real
            )
        survival = raw.sum(axis=1)
        np.testing.assert_allclose(survival, 4 / 9, atol=3e-13, rtol=0)
        gain = q @ raw @ q / (q @ q)
        six = contrast @ raw @ h6 / (contrast @ h6)
        result[name] = dict(
            survival=round(float(survival.mean()), 12),
            raw_two=round(float(gain), 12),
            accepted_two=round(float(gain / survival.mean()), 12),
            raw_six=round(float(six), 12) + 0.0,
        )
        kernels.append(raw)
    return (
        result,
        np.rint(np.stack(kernels) * SCALE).astype(np.uint32),
        np.array(stations, dtype=np.uint8),
    )


def run():
    rays = Rays()
    endpoints = rays.audit()
    pol = Polarimeter()
    frames = [pol.frame(n) for n in rays.directions]
    readings = pol.reads(frames, rays.net, endpoints["arrivals"])
    # One native half-splitting and a sharp rail read choose bypass or collision.
    guard = abs(pol.givens(0, 1, np.pi / 4)[:2, 0]) ** 2
    proposals = np.array(
        [np.prod(guard[[(i >> b) & 1 for b in range(5)]]) for i in range(32)]
    )
    np.testing.assert_allclose(proposals, 1 / 32, rtol=0, atol=2e-13)
    angular_data, counts = angular(
        rays.directions, readings, endpoints["arrivals"], guard, proposals
    )
    repeat, kernels, stations = repeat_encounter(rays, frames, pol, proposals)
    result = dict(
        second_encounter=repeat,
        endpoints=endpoints,
        directions=rays.directions.round(12).tolist(),
        polarization_probability_denominator=SCALE,
        angular=angular_data,
    )
    arrays = dict(
        tomography=np.rint(readings * SCALE).astype(np.uint32),
        counts=counts.astype(np.uint32),
        second_kernels=kernels,
        second_stations=stations,
    )
    return result, arrays


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result, arrays = run()
    body = json.dumps(result, separators=(",", ":")) + "\n"
    path = HERE / "records.json"
    data = HERE / "detectors.npz"
    if args.write:
        path.write_text(body)
        np.savez_compressed(data, **arrays)
    else:
        if path.read_text() != body:
            raise SystemExit("Angular readout summary differs")
        with np.load(data, allow_pickle=False) as saved:
            for k, v in arrays.items():
                np.testing.assert_array_equal(saved[k], v)
    print(result["angular"]["transfer"])
