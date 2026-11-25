import json

import numpy as np
from scipy.sparse import csr_matrix, diags

from . import compact, measure, readout

RECORDS = json.loads((measure.HERE / "records.json").read_text())


def test_independent_full_operator_and_charge_conjugation():
    mask, comparison = measure.native_geometry()
    assert comparison == 3 and np.flatnonzero(~mask).tolist() == [1, 3, 10]
    words, swaps = measure.sector(2, 2)
    kinetic, contrast, references = measure.operators(words, swaps, mask)
    k = measure.candidate.KAPPA
    shut = csr_matrix(kinetic) + diags(k * contrast.sum(axis=1))
    active = shut + diags(k * references.sum(axis=1))
    h, c, w = compact.construct()
    ids = {v.tobytes(): i for i, v in enumerate(words)}
    ix = np.array([ids[v.astype(words.dtype).tobytes()] for v in w])
    residual = (h - active[ix][:, ix]).tocsr()
    assert residual.nnz == 0
    med, shells, groups = measure.graph_read(words)
    np.testing.assert_equal(c, np.array(groups)[ix])
    reverse = np.array([ids[(-v).tobytes()] for v in words])
    assert (shut - shut[reverse][:, reverse]).nnz == 0
    odd_term = diags(4 * k * words[:, mask].sum(axis=1))
    assert (active - active[reverse][:, reverse] - odd_term).nnz == 0
    np.testing.assert_equal(med, med[reverse])
    np.testing.assert_equal(shells, -shells[reverse])
    assert (np.ptp(words.sum(axis=1)) == 0) and words[0].sum() == 0


def test_exact_connected_obstruction_and_sparse_floor():
    result = compact.measure()
    assert result["connected_words"] == 1260
    assert result["connected_exact_lower"] == 155
    assert result["full_variational_upper"] == 144
    spec = RECORDS["panels"]["neutral"]["native"]["spectrum"]
    assert abs(result["all"][0] - spec["floor"]) < 1e-8
    assert all(result[key][1] < 1e-8 for key in ("all", "connected_only", "other_only"))
    assert 158.26 < result["connected_only"][0] < 158.27


def test_trial_accounting_and_independent_energy_meter():
    starts = RECORDS["starts"]
    for name, panels in RECORDS["panels"].items():
        charge = sum(x * y for x, y in zip(measure.NAMES[name], (1, -1)))
        for tag, row in panels.items():
            profiles = row if tag == "unfiltered_profile" else row["profile"]
            for item in profiles.values():
                labels = np.array(item["keys"])
                assert np.all(labels[:, :4].sum(axis=1) == charge)
                assert item["accepted"] == sum(item["counts"])
                assert item["accepted"] + item["rejected"] == starts
            if tag == "unfiltered_profile":
                continue
            instrument = row["local_energy"]
            seam, reference = instrument["seam"], instrument["reference"]
            for item in (seam, reference):
                assert item["accepted"] == sum(item["counts"])
                assert item["accepted"] + item["rejected"] == starts
            k = 0 if tag == "exchange" else measure.candidate.KAPPA
            seam_means = (
                np.array(seam["counts"]) @ np.array(seam["keys"]) / seam["accepted"]
            )
            ref_mean = (
                np.dot(reference["counts"], reference["keys"]) / reference["accepted"]
            )
            energy = 30 * (seam_means[0] + k * seam_means[1])
            if tag != "no_reference":
                energy += 9 * k * ref_mean
            assert abs(energy - instrument["energy"]) < 1e-10
            assert abs(energy - row["spectrum"]["mean_energy"]) <= instrument["radius"]
    endpoints = RECORDS["accounting"]["native"]
    assert (
        endpoints["neutral_minus_positive"] - endpoints["energy_difference_radius"] > 0
    )
    assert endpoints["work_to_transfer_vacancy"] + endpoints["work_radius"] < 0


def test_record_only_centre_and_controls():
    native = RECORDS["panels"]["neutral"]["native"]["profile"]
    intrinsic = readout.read(native["intrinsic"])
    assert intrinsic["shells"][0][1] < 0 and intrinsic["distant"][0] > 0
    assert intrinsic["connected"][1] < 0.00241
    for name, row in native.items():
        if name.startswith("gateway"):
            result = readout.read(row)
            assert result["shells"][0][0] > 0.33
            assert result["distant"][1] < -0.09 and result["shells"][3][1] < -0.21
    for tag in ("unfiltered_profile", "no_reference", "exchange"):
        row = RECORDS["panels"]["neutral"][tag]
        for panel in row.get("profile", row).values():
            assert all(lo <= 0 <= hi for lo, hi in readout.read(panel)["shells"])
