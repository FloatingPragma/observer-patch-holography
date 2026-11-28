"""Read the native parked source through resolved seam energies and test its stress."""

from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from os.path import relpath

from mpmath.ctx_iv import MPIntervalContext
import numpy as np
from scipy.linalg import expm
import sympy as sp

from . import reference_arrival as reference


def station_records():
    routes = reference.local.survey.routes
    volume = routes.volume
    law = volume.candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + volume.candidate.KAPPA * contrast
    vacuum = np.eye(law.N_SRC)[law.SRC.index((0, 0))]
    memories = [volume.label(0, 1), volume.label(2, 1), volume.label(2, -1)]
    states = np.column_stack(
        [np.kron(memory, vacuum) for memory in [*memories, np.eye(11)[0]]]
    )
    energies = np.array([0, 32, 33, 0])
    np.testing.assert_allclose(h @ states, states * energies, rtol=0, atol=3e-13)
    np.testing.assert_allclose(
        contrast @ states, states * [0, 4, 4, 0], rtol=0, atol=3e-13
    )
    ordinary = np.kron(np.eye(11)[0], np.eye(law.N_SRC)[law.SRC.index((1, 0))])
    np.testing.assert_allclose(contrast @ ordinary, ordinary, rtol=0, atol=3e-13)

    # Ramsey outcomes of an external controlled-evolution meter. On an
    # eigenpreparation, the minus-port amplitude is (I - exp(-ith))/2.
    def minus(t):
        return np.sum(abs((states - expm(-1j * t * h) @ states) / 2) ** 2, axis=0)

    parity = minus(np.pi)
    np.testing.assert_allclose(parity, [0, 0, 1, 0], rtol=0, atol=1e-12)
    even = np.flatnonzero(energies % 2 == 0)
    half_period = minus(np.pi / 32)[even]
    np.testing.assert_allclose(half_period, [0, 1, 0], rtol=0, atol=1e-12)
    sharp = np.where(np.rint(parity).astype(bool), 33, 0)
    sharp[even] = 32 * np.rint(half_period).astype(int)
    # Separate preparations insert ordinary loads at neither, just one,
    # or both endpoints. Every other term sees at most one inserted load.
    expected = []
    for loads in ((0, 0), (1, 0), (0, 1), (1, 1)):
        endpoints = np.eye(law.N_SRC)[law.SRC.index(loads)]
        probes = np.column_stack(
            [np.kron(m, endpoints) for m in [*memories, np.eye(11)[0]]]
        )
        expected.append(np.sum(probes.conj() * (h @ probes), axis=0).real)
    mixed = expected[3] - expected[1] - expected[2] + expected[0]
    np.testing.assert_allclose(mixed, [0, -32, -33, -17], rtol=0, atol=1e-12)
    mixed = np.rint(mixed).astype(int)
    edges = sorted(routes.cell().seams)
    reports, insertions = [], []
    for word, matching in enumerate(routes.codebook()[0][1:]):
        digits = divmod(word, 3)
        marked = {edge: digits[k] if k < 2 else 0 for k, edge in enumerate(matching)}
        reports.append([sharp[marked.get(edge, 3)] for edge in edges])
        insertions.append([mixed[marked.get(edge, 3)] for edge in edges])
    return edges, np.array(reports), np.array(insertions)


def decks(reports, classes):
    # A constant column cannot separate two words. Test every subset of
    # the remaining seams, including subsets that separate classes only.
    cols = np.flatnonzero(np.ptp(reports, axis=0))
    for size in range(len(cols) + 1):
        valid = []
        for subset in combinations(cols, size):
            letters = [tuple(row) for row in reports[:, subset]]
            if len(set(zip(letters, classes, strict=True))) == len(set(letters)):
                valid.append(subset)
        if valid:
            return valid
    raise ValueError("The local energy deck failed to distinguish source classes")


def controls(fresh):
    # The signed contrasts refer to the existing nine calibration words.
    u = np.array([1, 0, -1, 1, -66, 65, 0, 0, 0])
    v = np.eye(9, dtype=int)[1] - np.eye(9, dtype=int)[3]
    delta = np.outer(u, v) + np.outer(v, u)
    base = 8 * fresh
    return u, v, {"fresh": base, "plus": base + delta, "minus": base - delta}


def categorical(x):
    return (x[:, None] == np.unique(x)[None, :]).astype(int)


def quadratic_class_read(symbols, classes):
    pairs = list(combinations_with_replacement(range(symbols.shape[1]), 2))
    features = sp.Matrix(
        np.column_stack(
            [np.ones(len(symbols), dtype=int), symbols]
            + [symbols[:, i] * symbols[:, j] for i, j in pairs]
        )
    )
    pivots = list(features.rref()[1])
    independent = features[:, pivots]
    indicators = sp.Matrix(categorical(classes)[:, 1:])
    weights = independent.inv() * indicators
    assert independent * weights == indicators
    names = [()] + [(i,) for i in range(symbols.shape[1])] + pairs
    return [
        {"factors": list(names[i]), "weights": [str(w) for w in row]}
        for i, row in zip(pivots, weights.tolist(), strict=True)
        if any(row)
    ]


def run():
    iv = MPIntervalContext()
    iv.dps = 50
    edges, readings, insertions = station_records()
    survey = reference.local.survey
    ref, fresh, shared = survey.standards(survey.bulk.federation(0))[0][:3]
    classes = np.array([(32, 48, 16).index(x) for x in ref[0]])
    stations = reference.local.inputs("pair_clock_information")["local_reader"][
        "calibrations"
    ]["fixed_contact"]["sharp_records_in_calibration_order"]
    energy, bits = np.array(stations).T
    np.testing.assert_array_equal(energy, readings.sum(axis=1))
    energy_bins, bit_bins = categorical(energy), categorical(bits)
    choices = decks(readings, classes)
    assert len(choices) == 1
    deck = choices[0]
    symbols = readings[:, deck]
    assert len(np.unique(symbols, axis=0)) == 9
    quadratics = quadratic_class_read(symbols, classes)
    # Rational coefficients are calibrated on the nine sharp words alone.
    # Read the two reference-class indicators using ordinary material work.
    contacts = [(2, 3), (3, 6), (3, 8), (3, 9)]
    b23, b36, b38, b39 = insertions[:, [edges.index(e) for e in contacts]].T
    high = -19 * 17 - 17 * (b23 + b36 + b38) + 15 * b39
    low = 17 + b36
    work_flags = np.column_stack([high // 544, low // 17])
    np.testing.assert_array_equal(high, 544 * work_flags[:, 0])
    np.testing.assert_array_equal(low, 17 * work_flags[:, 1])
    np.testing.assert_array_equal(work_flags, categorical(classes)[:, 1:])
    block = [[iv.mpf(int(ref[i, j])) for j in (0, 1, 2)] for i in (0, 1, 2)]
    interaction = [
        [
            iv.log(block[i][0] * block[0][j] / (block[i][j] * block[0][0]))
            for j in (1, 2)
        ]
        for i in (1, 2)
    ]
    determinant = (
        interaction[0][0] * interaction[1][1] - interaction[0][1] * interaction[1][0]
    )
    assert interaction[0][0] > 0 and determinant > 0
    u, v, stocks = controls(fresh)
    assert u.sum() == v.sum() == 0
    assert not np.any(u @ readings) and not np.any(u @ bit_bins)
    assert not np.any(v @ energy_bins)
    base = stocks["fresh"]
    for matrix in stocks.values():
        assert matrix.min() > 0
        delta = matrix - base
        for axis in (0, 1):
            assert not np.any(delta.sum(axis=axis))
        for features in (readings, energy_bins, bit_bins):
            assert not np.any(features.T @ delta @ features)

    # Each calibrated meter has a deterministic output on this input deck.
    # Enumerate every trial and publish only the pair of meter outcomes.
    _, order = np.unique(symbols, axis=0, return_index=True)
    counts = {}
    # A second contrast kills even the cross-station cubics. Any cubic
    # involving BOTH stations has degree at most one on at least one side.
    w = np.array([1, 0, -1, -1, 0, 1, 0, 0, 0])
    assert w.sum() == 0 and not np.any(w @ readings)
    stocks["fourth_order"] = base + 32 * np.outer(w, w)
    assert stocks["fourth_order"].min() > 0
    stocks["shared"] = 8 * shared
    quench = 4 * np.count_nonzero(readings, axis=1)
    work = {}
    for name, table in {"reference": 8 * ref, **stocks}.items():
        left, right = np.unravel_index(
            np.repeat(np.arange(81), table.ravel()), table.shape
        )
        report_pairs = np.column_stack([symbols[left], symbols[right]])
        alphabet, histogram = np.unique(report_pairs, axis=0, return_counts=True)
        assert len(alphabet) == 81
        counts[name] = histogram.tolist()
        work[name] = str(Fraction(int(np.sum(quench[left] + quench[right])), len(left)))
    assert len(set(work.values())) == 1
    cost_span = iv.log(iv.mpf(int(ref.max())) / int(ref.min()))
    changes = {}
    for name in ("plus", "minus", "fourth_order", "shared"):
        delta = np.subtract(counts[name], counts["fresh"])
        direct = -sum(
            int(d) * iv.log(iv.mpf(q) / int(base.sum()))
            for d, q in zip(delta, counts["reference"], strict=True)
        ) / int(base.sum())
        exact = direct
        cov = work_flags.T @ (stocks[name] - base) @ work_flags
        from_work = sum(
            int(cov[a, b]) * interaction[a][b] for a, b in np.ndindex(2, 2)
        ) / int(base.sum())
        assert 0 in direct - from_work
        if name in ("plus", "minus"):
            sign = 1 if name == "plus" else -1
            exact = sign * 130 * iv.log(iv.mpf(39) / 35) / int(base.sum())
            assert 0 in exact - direct
        if name == "fourth_order":
            exact = 32 * iv.log(iv.mpf(114244) / 25515) / int(base.sum())
            assert 0 in exact - direct
        changes[name] = {
            "exact": reference.bounds(exact),
            "mixed_work_flag_covariance_delta": [
                [str(Fraction(int(v), int(base.sum()))) for v in row] for row in cov
            ],
            "max_total_variation_per_arm_for_sign_fixed_reference": reference.bounds(
                abs(exact) / (2 * cost_span)
            ),
        }

    galaxies = reference.local.inputs("modular_stock_sky")["galaxies"]
    sky = {
        allowance: sum(
            Fraction(g["bounds"][allowance]["required_fraction_lower"]) > 0
            for g in galaxies.values()
        )
        for allowance in ("catalogue", "expanded_HI")
    }
    return {
        "energy_seams": [list(edges[i]) for i in deck],
        "sharp_spectra_in_calibration_order": symbols.tolist(),
        "classes_in_calibration_order": classes.tolist(),
        "station_energy": energy.tolist(),
        "quadratic_class_read": quadratics,
        "mixed_work_contacts": [list(e) for e in contacts],
        "mixed_work_on_nine_calibrations": insertions[
            :, [edges.index(e) for e in contacts]
        ].tolist(),
        "mixed_work_response_determinant": reference.bounds(determinant),
        "arrival_bit": bits.tolist(),
        "u": u.tolist(),
        "v": v.tolist(),
        "w": w.tolist(),
        "fresh_counts": base.tolist(),
        "num_seam_product_moments_held": readings.shape[1] ** 2,
        "denominator": int(base.sum()),
        "meter_order": order.tolist(),
        "counts": counts,
        "native_work_per_unit_coupling_change": work,
        "charge_change_nats": changes,
        "sparc_galaxies": len(galaxies),
        "sparc_positive_extra_requirement": sky,
        "inputs": reference.coherent.steady.exp.save_inputs(
            "pair_clock_information.json",
            "modular_stock_sky.json",
            relpath(survey.__file__, reference.local.HERE),
            relpath(survey.routes.__file__, reference.local.HERE),
            relpath(survey.routes.volume.__file__, reference.local.HERE),
            relpath(
                survey.routes.volume.candidate.native_law.__file__, reference.local.HERE
            ),
        ),
    }


if __name__ == "__main__":
    result = run()
    reference.coherent.steady.exp.emit("parked_energy_charge", result)
    print(result["charge_change_nats"], result["sparc_positive_extra_requirement"])
