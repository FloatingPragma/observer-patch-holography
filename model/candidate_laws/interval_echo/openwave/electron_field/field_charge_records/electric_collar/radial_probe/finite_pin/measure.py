"""Certify connected static energies with finite pins in the written radial collar."""

import argparse
from functools import lru_cache
from hashlib import sha256
import json
from math import cos, inf, nextafter, sin
from pathlib import Path
from runpy import run_path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import connected_components

from .....fan_path.fixture import FanPath, archived_fixture
from ..apparatus import RadialCollar
from ..measure import LABEL, READS, late_arrival


HERE = Path(__file__).resolve().parent
PINS = ((64, 64), (128, 128), (256, 256), (512, 512), (128, 256), (256, 128))
PRESENCE = ((0, 0), (1, 0), (0, 1), (1, 1))
CONNECTORS = (1, 0.5, 0)
WEAK_RATES = (0.5, 0.25, 1 / 32, 1 / 256, 1 / 4096, 1 / 65536)


def norm_inf(a):
    return float(np.max(np.sum(abs(a), axis=1)))


def spectrum(h):
    """Bound the smallest eigenvalue using a checked real spectral factorization."""
    n = len(h)
    np.testing.assert_array_equal(h, h.T)
    scaled = 2**16 * h
    np.testing.assert_array_equal(scaled, np.rint(scaled))
    assert n <= 512 and norm_inf(h) < 2048
    values, columns = eigh(h)
    absolute = abs(columns)
    guard = 16 * n * 2.0**-53
    gram = (1 + guard) * (
        norm_inf(columns.T @ columns - np.eye(n))
        + guard * (1 + max(absolute.T @ absolute.sum(axis=1)))
    )
    reconstructed = norm_inf((columns * values) @ columns.T - h)
    scale = norm_inf(h) + max(absolute @ (abs(values) * absolute.sum(axis=0)))
    assert gram < 1e-8
    # Polar correction turns the represented columns into an orthogonal
    # basis with norm distance at most gram. The conjugated diagonal's
    # displacement is then at most three times the largest |energy| times gram.
    radius = (1 + guard) * (reconstructed + guard * scale) + 3 * max(abs(values)) * gram
    return values, columns[:, 0], nextafter(float(radius + 2e-14), inf)


class Enrolled:
    def __init__(self, collar, connector, present, density):
        self.collar = collar
        self.present = present
        self.density = density
        terms = collar.terms
        self.ring = (
            terms["fan_rings"] + terms["outer_ring"] + connector * terms["connector"]
        )
        self.h = (
            self.ring
            + present[0] * (terms["source_hop"] + terms["source_degree"])
            + present[1] * (terms["probe_hop"] + terms["probe_degree"])
        )
        self.h.eliminate_zeros()
        self.labels = connected_components(self.h, directed=False)[1]
        self.source_pos = np.array([s for s, _, _ in collar.basis])
        self.probe_pos = np.array([p for _, p, _ in collar.basis])

    def index(self, radius, b):
        i = 2 * (radius * self.collar.backgrounds + b)
        assert self.labels[i] == self.labels[i + 1]
        return i

    @lru_cache(maxsize=None)
    def sector(self, label, radius):
        scope = np.flatnonzero(self.labels == label)
        source_away = self.source_pos[scope] != self.collar.hub
        probe_away = self.probe_pos[scope] != self.collar.radial[radius]
        if not self.present[0]:
            assert not source_away.any()
        if not self.present[1]:
            assert not probe_away.any()
        core = ~(source_away | probe_away)
        h = self.h[scope][:, scope].toarray()
        a = h[np.ix_(core, core)]
        k_source = h[np.ix_(core, source_away & ~probe_away)]
        k_probe = h[np.ix_(core, probe_away & ~source_away)]
        np.testing.assert_array_equal(h[np.ix_(core, source_away & probe_away)], 0)
        source_diagonal = k_source @ k_source.T
        probe_diagonal = k_probe @ k_probe.T
        np.testing.assert_array_equal(
            source_diagonal, self.present[0] * 48 * np.eye(len(a))
        )
        np.testing.assert_array_equal(
            probe_diagonal, np.diag(probe_diagonal.diagonal())
        )
        assert np.all(np.isin(probe_diagonal.diagonal(), [0, 16]))
        e_trial = float(min(a.diagonal()))
        high = norm_inf(h[np.ix_(~core, ~core)]) if (~core).any() else 0.0
        k_squared = norm_inf(source_diagonal + probe_diagonal)
        return (
            scope,
            source_away,
            probe_away,
            core,
            h,
            a,
            probe_diagonal,
            e_trial,
            high,
            k_squared,
        )

    @lru_cache(maxsize=None)
    def simulation(self, label, radius, pins):
        scope, source, probe, _, h, _, _, e_trial, d_norm, k2 = self.sector(
            label, radius
        )
        pin_floor = min(pins)
        delta = pin_floor - e_trial
        assert delta > 0
        matrix = h + np.diag(pins[0] * source + pins[1] * probe)
        values, ground, error = spectrum(matrix)
        energy = float(values[0])
        return (
            energy,
            error,
            ground,
            nextafter(k2 * (d_norm + e_trial) / (pin_floor * delta), inf),
            nextafter(k2 / delta**2, inf),
        )

    def readout(self, radius, b, pins):
        start = self.index(radius, b)
        label = int(self.labels[start])
        raw, error, ground, schur, outside_bound = self.simulation(label, radius, pins)
        scope = self.sector(label, radius)[0]
        indices = np.searchsorted(scope, [start, start + 1])
        wave = ground[indices]
        weight = np.vdot(wave, self.density @ wave).real
        assert weight > 0
        return [raw, error, float(weight), schur, outside_bound]


def physical_spectrometer(p, label, inherited):
    spectra, admissions, _, _ = late_arrival(p, label, inherited)
    density = spectra[0].sum(axis=0)
    assert np.min(np.linalg.eigvalsh(density)) > 1e-5
    return density, float(sum(admissions))


def check_removal(p, collar):
    known = set(port for edge in collar.edges for port in edge)
    reservoir = (
        set.intersection(*(set(p.board.neighbours[x]) for x in collar.radial)) - known
    )
    assert reservoir == {1055}
    assert set(int(p.board.degree[x]) for x in (*collar.radial, 1055)) == {5}
    incoming, outer = p.bridges[-1]
    assert incoming == 1010 and outer == collar.hub
    assert incoming not in known and int(p.board.degree[incoming]) == int(
        p.board.degree[outer]
    )
    # Same-rate selected feeder: the short two-port motion keeps all other
    # written heads, and its own reserved head is thereafter excluded.
    assert all(
        tuple(sorted((x, 1055))) in set(tuple(sorted(e)) for e in p.board.seams)
        for x in collar.radial
    )


def match_pinned_words(quartet, radius, b):
    pinned = []
    for arm in quartet:
        i = arm.index(radius, b)
        scope, _, _, core, _, a, bp, *_ = arm.sector(int(arm.labels[i]), radius)
        pinned.append((scope[core], a, bp))
    (w00, a00, b00), (w10, a10, b10), (w01, a01, b01), (w11, a11, b11) = pinned
    for w in (w10, w01, w11):
        np.testing.assert_array_equal(w, w00)
    np.testing.assert_array_equal(a10 - a00, 12 * np.eye(len(a00)))
    np.testing.assert_array_equal(a11 - a01, 12 * np.eye(len(a01)))
    np.testing.assert_array_equal(b11, b01)
    np.testing.assert_array_equal(b00, b10)


def calculation(collar, density, connector):
    quartet = [Enrolled(collar, connector, setting, density) for setting in PRESENCE]
    output = np.zeros((2, collar.backgrounds, len(PINS), 4, 5))
    for radius, b in np.ndindex(output.shape[:2]):
        match_pinned_words(quartet, radius, b)
        for i, penalty in enumerate(PINS):
            output[radius, b, i] = [arm.readout(radius, b, penalty) for arm in quartet]
    weights = np.array([1, -1, -1, 1])
    connected = output[..., 0] @ weights
    error = np.nextafter(output[..., 1].sum(axis=-1) + 2e-12, np.inf)
    schur = np.nextafter(output[..., 3].sum(axis=-1) * (1 + 2**-48), np.inf)
    if connector == 0:
        assert np.all(abs(connected) <= error)
    else:
        assert np.all(abs(connected) <= schur + error)
    for arm in quartet:
        arm.sector.cache_clear()
        arm.simulation.cache_clear()
    return output, connected, error, schur


def extrema(values, widths, addresses):
    entries = []
    for k in range(len(PINS)):
        vv, ee = values[:, k], widths[:, k]
        result = {}
        for key, choice in (("minimum", np.argmin(vv)), ("maximum", np.argmax(vv))):
            x = int(choice)
            result[key] = {
                "near_and_far_prepared_head_words": list(addresses[x]),
                "energy": float(vv[x]),
                "interval": [
                    nextafter(float(vv[x] - ee[x]), -inf),
                    nextafter(float(vv[x] + ee[x]), inf),
                ],
            }
        entries.append(result)
    return entries


def summarize(collar, computed):
    output, connected, err, schur = computed
    probe_bit = collar.probe_mask >> collar.inside_bits
    accepted = [
        bg
        for bg in range(collar.backgrounds)
        if bool(bg & probe_bit) == (collar.probe_sign == -1)
    ]
    opposite = [bg ^ probe_bit for bg in accepted]
    for a, b in zip(accepted, opposite, strict=True):
        i, j = 2 * a, 2 * (collar.backgrounds + b)
        for word in range(2):
            assert collar.terms["probe_hop"][i + word, j + word] == -4
    identical = connected[1] - connected[0]
    same_sector = connected[1, opposite] - connected[0, accepted]
    return {
        "presence_raw_records": output.tolist(),
        "same_exterior_radial_extrema": extrema(
            identical, err[1] + err[0], [(b, b) for b in range(collar.backgrounds)]
        ),
        "same_probe_sector_radial_extrema": extrema(
            same_sector,
            err[1, opposite] + err[0, accepted],
            list(zip(accepted, opposite, strict=True)),
        ),
        "uniform_second_order_radial_upper": (
            schur[1].max(axis=0) + schur[0].max(axis=0)
        ).tolist(),
        "bare_connected_energy_max_abs": abs(connected).max(axis=(0, 1)).tolist(),
        "computed_max_outside_pin_probability_bound": output[..., 4]
        .max(axis=(0, 1, 3))
        .tolist(),
        "minimum_ground_outcome_probability_per_launch": output[..., 2]
        .min(axis=(0, 1, 3))
        .tolist(),
    }


def conjugate_probe(a, b):
    """Reversing only probe sign and label exchanges the two pin addresses."""
    assert a.radial == b.radial and len(a.basis) == len(b.basis)
    other = dict(zip(a.radial, a.radial[::-1], strict=True))
    permutation = [b.ids[s, other[p], w] for s, p, w in a.basis]
    assert len(set(permutation)) == len(a.basis)
    for label, left in a.terms.items():
        assert (b.terms[label][permutation][:, permutation] != left).nnz == 0


def charge_conjugacy(a, b):
    mask = (1 << len(a.edges)) - 1
    permutation = [b.ids[s, p, w ^ mask] for s, p, w in a.basis]
    assert len(set(permutation)) == len(a.basis) == len(b.basis)
    for label, left in a.terms.items():
        assert (b.terms[label][permutation][:, permutation] != left).nnz == 0


def weak_connector(collar, density, acceptance):
    """Compare the limiting open-chart floor and the actual cut-accessible line."""
    word_pair = (45, 47)
    pins = (64, 64)
    open_chart = [Enrolled(collar, 1, occupation, density) for occupation in PRESENCE]
    cut_chart = [Enrolled(collar, 0, occupation, density) for occupation in PRESENCE]
    connector = collar.terms["connector"]
    # A disjoint family of two-word ring swaps is a projector at native rate.
    assert (connector @ connector != connector).nnz == 0
    visibility = []
    for radius, b in enumerate(word_pair):
        entries = []
        for a, c in zip(open_chart, cut_chart, strict=True):
            i = a.index(radius, b)
            indices, sa, pa, *_ = a.sector(int(a.labels[i]), radius)
            h = c.h[indices][:, indices].toarray() + np.diag(
                pins[0] * sa + pins[1] * pa
            )
            floor, _, interval_error = spectrum(h)
            reached, reached_error, *_ = c.readout(radius, b, pins)
            gap = reached - reached_error - float(floor[0]) - interval_error
            entries.append(
                {
                    "open_chart_limit_energy_and_error": [
                        float(floor[0]),
                        interval_error,
                    ],
                    "cut_chart_reachable_energy_and_error": [reached, reached_error],
                    "gap_above_open_limit_lower": max(
                        0.0, nextafter(gap - 2e-12, -inf)
                    ),
                }
            )
        visibility.append(entries)
    assert visibility[0][2]["gap_above_open_limit_lower"] > 0.26
    assert visibility[0][3]["gap_above_open_limit_lower"] > 0.25
    rows = []
    for g in WEAK_RATES:
        arms = [Enrolled(collar, 2 * g, occupancy, density) for occupancy in PRESENCE]
        values = []
        for radius, b in enumerate(word_pair):
            profile = []
            for j, arm in enumerate(arms):
                i = arm.index(radius, b)
                open_indices = open_chart[j].sector(
                    int(open_chart[j].labels[i]), radius
                )[0]
                actual_indices = arm.sector(int(arm.labels[i]), radius)[0]
                np.testing.assert_array_equal(actual_indices, open_indices)
                energy, error, probability, *_ = arm.readout(radius, b, pins)
                gap = visibility[radius][j]["gap_above_open_limit_lower"] - 2 * g
                # Project the eigenvalue equation into the cut seed component.
                # It is separated by gap from the weakly continued eigenvalue.
                upper = acceptance * ((2 * g) / gap) ** 2 if gap > 0 else acceptance
                profile.append(
                    {
                        "energy_and_error": [energy, error],
                        "nominal_ground_probability_per_launch": probability,
                        "ground_probability_per_launch_upper": nextafter(
                            min(acceptance, upper), inf
                        ),
                    }
                )
            values.append(profile)
        rows.append({"connector_rate": g, "two_radial_presence_banks": values})
    final = rows[-1]["two_radial_presence_banks"][0]
    assert final[2]["ground_probability_per_launch_upper"] < 3e-10
    assert final[3]["ground_probability_per_launch_upper"] < 3e-10
    return {
        "near_and_far_exterior_words": list(word_pair),
        "pins": list(pins),
        "energy_floor_limit_under_cut_in_open_chart_versus_reachable_cut_chart": visibility,
        "rates": rows,
    }


def run():
    parent = json.loads(READS.read_text())
    flight = next(
        f
        for c in parent["cases"]
        if c["fan_count"] == 4
        for s in c["stock_arms"]
        if s["preparation"] == "clock_bypassed"
        for f in s["flights"]
        if f["material_bridge_rate"] == 4 and f["appointment"] == 30.8125
    )
    witness = next(
        w for w in flight["pair_readings"] if w["port_face_and_words"] == LABEL
    )
    inherited_error = flight["observable_absolute_error_upper_per_pair"]
    acceptance_cap = nextafter(
        witness["receiver_reference_pauli_moments"][0][0] + inherited_error + 2e-12,
        inf,
    )
    with archived_fixture() as archive:
        positive, negative = [FanPath(archive, 4, s, asymmetric=True) for s in (1, -1)]
        code_mask = (1 << len(positive.field.edges)) - 1
        gauges = [
            negative.ids[x, tuple(w ^ code_mask for w in words)]
            for x, words in positive.basis
        ]
        assert (negative.h[gauges][:, gauges] != positive.h).nnz == 0
        np.testing.assert_array_equal(negative.source()[gauges], -positive.source())
        cases = []
        collars = {}
        for path in (positive, negative):
            label = (
                LABEL
                if path.sign == 1
                else [*LABEL[:2], *(w ^ code_mask for w in LABEL[-2:])]
            )
            density, acceptance = physical_spectrometer(path, label, witness)
            assert (
                abs(acceptance - witness["receiver_reference_pauli_moments"][0][0])
                < 2e-12
            )
            for probe_sign in (1, -1):
                collar = RadialCollar(path, probe_sign, label[-2:])
                collars[path.sign, probe_sign] = collar
                check_removal(path, collar)
                comparison = []
                for g in CONNECTORS:
                    result = summarize(collar, calculation(collar, density, g))
                    comparison.append({"connector_rate": g / 2, **result})
                cases.append(
                    {
                        "source_sign": path.sign,
                        "probe_sign": probe_sign,
                        "original_face_admission": acceptance,
                        "original_face_admission_upper": acceptance_cap,
                        "connectors": comparison,
                    }
                )
                print(
                    f"source {path.sign:+}, probe {probe_sign:+}, all pin banks complete",
                    flush=True,
                )
                if path.sign == 1 and probe_sign == 1:
                    weak = weak_connector(collar, density, acceptance_cap)
        for sign in (1, -1):
            conjugate_probe(collars[sign, 1], collars[sign, -1])
            charge_conjugacy(collars[1, sign], collars[-1, -sign])
        plus, minus = [
            np.array(case["connectors"][0]["presence_raw_records"])[0, 45, 0]
            for case in cases[:2]
        ]
        weights = np.array([1, -1, -1, 1])
        sign_sum = float((plus[:, 0] + minus[:, 0]) @ weights)
        sign_error = float(np.sum(plus[:, 1] + minus[:, 1]) + 2e-12)
        assert sign_sum - sign_error > 0
        return {
            "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
            "word_reads_sha256": sha256(READS.read_bytes()).hexdigest(),
            "probe_ports_near_far": list(collar.radial),
            "intact_board_hops": collar.ranges,
            "probe_removed_to": 1055,
            "source_removed_to": 1010,
            "zero_connected_energy_in_unfiltered_prepared_mean": True,
            "exact_probe_sign_position_interchange": True,
            "exact_joint_sign_full_head_conjugacy": True,
            "fixed_near_background_probe_sign_sum_example": {
                "source_sign": 1,
                "near_port": 1048,
                "word": 45,
                "pin_strengths": [64, 64],
                "connector_rate": 0.5,
                "connected_energy_probe_sign_plus_and_minus": [
                    float(plus[:, 0] @ weights),
                    float(minus[:, 0] @ weights),
                ],
                "sum_interval": [
                    nextafter(sign_sum - sign_error, -inf),
                    nextafter(sign_sum + sign_error, inf),
                ],
            },
            "source_and_probe_pin_settings": [list(pair) for pair in PINS],
            "bare_pin_phase_calibration_exposure": 1 / 1024,
            "ramsey_calibration_cos_sin_by_penalty": [
                [x, cos(x / 1024), -sin(x / 1024)] for x in (64, 128, 256, 512)
            ],
            "presence_order": [list(x) for x in PRESENCE],
            "raw_record_axes": [
                "probe_radius",
                "exterior_head_word",
                "pin_setting",
                "presence",
            ],
            "raw_record_fields": [
                "ground_energy",
                "absolute_energy_error",
                "nominal_ground_probability_per_launch",
                "second_order_error_upper",
                "outside_pin_probability_upper",
            ],
            "cases": cases,
            "selected_weak_connector_study": weak,
        }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path)
    args = cli.parse_args()
    data = run()
    encoded = json.dumps(data, separators=(",", ":"), allow_nan=False)
    if args.output:
        with args.output.open("x") as out:
            out.write(encoded + "\n")
    else:
        compare = run_path(str(HERE.parents[3] / "replay.py"))["compare"]
        compare(
            json.loads(encoded),
            json.loads((HERE / "results.json").read_text()),
            "finite_pin",
        )
        print(
            "Verified finite pin, matched presence and both radial energy comparisons."
        )


if __name__ == "__main__":
    main()
