from fractions import Fraction
import json

import numpy as np
import pytest
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.executable_model.carrier import HERE as CARRIER_HOME
from ..submission import row_04_check
from . import open_pair as opened


def test_exact_open_pair_and_native_binary_columns():
    setup = opened.seal_from_sharp()
    dark = opened.dark_certificate()
    cut, home = opened.geometry()
    assert dark["dimension"] == 28680 and tuple(setup["final_pair"]) == (156, 157)
    assert len(cut) == len(setup["neutral_edges"]) == len(set(setup["neutral_edges"])) == 11
    assert set(cut) == set(setup["neutral_edges"]) and len(setup["live_swaps"]) == 4
    assert sum(dark["dark_bins"]) == 1 and dark["dark_bins"][0] == dark["dark_bins"][3] == .5
    board = federation(0)
    edges, index, _, _ = row_04_check.screen(CARRIER_HOME/"carriers"/"FEDERATION_LEVEL0.npz")
    _, exchange, contrast, detector, source = opened.full_parts()
    assert [tuple(e) for e in edges if e[0]//12 == e[1]//12 == 13] == list(home)
    assert list(detector) == [index[edge] for edge in home]
    # Directly assemble the three full flight columns from physical seam
    # operations; look for all destinations, including those outside the source.
    for k in (-2.5, 0., 8.):
        direct = np.zeros((exchange.shape[0], 3))
        term = candidate.native_law.seam_term(k)
        alphabet = candidate.native_law.SRC
        for i, pair in enumerate(opened.BASIS):
            for edge in board.seams:
                if edge in cut:
                    continue
                word = tuple(int(v in pair) for v in edge)
                column = term[:, alphabet.index(word)]
                for j in np.flatnonzero(column):
                    result = [v for v in pair if v not in edge]
                    for v, count in zip(edge, alphabet[j], strict=True):
                        result.extend([v]*count)
                    direct[index[tuple(sorted(result))], i] += column[j]
            direct[index[pair], i] += k*sum(6-board.degree[v] for v in pair)
        native = exchange + sparse.diags(k*contrast)
        np.testing.assert_array_equal(native[:, source].toarray(), direct)


def test_symmetric_packet_can_reach_outside_the_home():
    length, numerator, denominator = opened.bright_escape()
    assert length == 4 and Fraction(numerator, denominator) == Fraction(1, 18432)
    sector, exchange, _, _, source = opened.full_parts()
    bright = np.zeros(sector.P)
    bright[source[[0, 2]]] = 1/np.sqrt(2)
    t = .05
    evolved = expm_multiply((-1j*t)*exchange, bright)
    outside = (sector.aa//12 != 13) | (sector.bb//12 != 13)
    ratio = np.linalg.norm(evolved[outside])**2 / t**(2*length)
    assert abs(ratio/(numerator/denominator)-1) < .005
    dephased = np.diag([.5, 0., .5])
    bins = opened.read_release(8, 1, dephased, flight="stochastic_swap")
    assert abs(bins[:-1].sum()-.7560097781598543) < 2e-12
    assert bins[-1] > .24 and abs(bins.sum()-1) < 2e-13


def test_readouts_keep_coherence_and_source_failure():
    dark = np.array([1, 0, -1])/np.sqrt(2)
    rho = .7*np.outer(dark, dark)
    expected = .7*np.asarray(opened.dark_certificate()["dark_bins"])
    for k in (0., .5, 8.):
        observed = opened.read_release(k, 1.5, rho)
        np.testing.assert_allclose(observed, expected, rtol=0, atol=2e-12)
    wrong = opened.read_release(0, 1.5, rho, flight="dephased")
    assert abs(wrong.sum()-.7) < 2e-13 and wrong[-1] > .1


def test_different_contact_rejects_the_exact_dark_claim(monkeypatch):
    sector, exchange, contrast, detector, source = opened.full_parts()
    shifted = contrast.copy()
    shifted[source[0]] += 1
    monkeypatch.setattr(opened, "full_parts", lambda: (sector, exchange, shifted, detector, source))
    with pytest.raises(ArithmeticError, match="dark state"):
        opened.dark_certificate()


@pytest.mark.parametrize("rho", [np.diag([.8, .2, -.1]), np.diag([1,1,0]), np.diag([0,0,np.nan]),
                                 np.ones(3), np.triu(np.ones((3, 3)))])
def test_unphysical_release_rejected(rho):
    with pytest.raises(ValueError, match="accepted density"):
        opened.read_release(8, 1, rho)


def test_raw_board_replay_and_quantized_receipt():
    from . import open_pair_check as check, open_pair_release as release
    rho, report = release.source()
    report["home_seams_in_order"] = opened.geometry()[1]
    report["appointments"] = [{"kappa": .5, "exposure": 1,
                               "bins": release.terminal(rho, report["exhaustion"], .5, 1)}]
    archive = json.loads(check.ARCHIVE.read_text())
    result = check.verify(report, (.5,), (1,), archive)
    assert result["checked_bins"] == 124 and result["max_archive_error"] < check.TOLERANCE
    raw = report["appointments"][0]["bins"]["coherent"]
    position = int(np.argmin(raw))
    raw[position] = False  # numerically indistinguishable from this near-zero cell
    with pytest.raises(ValueError, match="numeric terminal"):
        check.verify(report, (.5,), (1,), archive)
    raw[position] = 0.
    edges = list(report["setup"]["neutral_edges"])
    edges[0] = (156, 157)  # a live path seam substituted for an actual cut
    report["setup"]["neutral_edges"] = edges
    with pytest.raises(ValueError, match="contact installation"):
        check.verify(report, (.5,), (1,), archive)
    report["setup"] = opened.seal_from_sharp()
    report["setup"]["final_pair"] = [156., 157]
    with pytest.raises(ValueError, match="contact installation"):
        check.verify(report, (.5,), (1,), archive)


def test_whole_archive_agrees_with_the_raw_board_replay():
    from . import open_pair_check as check
    archive = json.loads(check.ARCHIVE.read_text())
    _, home, expected = check.independently_replay(check.STRENGTHS, check.TIMES)
    assert len(expected) == 9
    assert check.max_bin_error(check.archive_bins(archive, home), expected) < check.TOLERANCE


def test_external_raw_replay_rejects_duplicate_json_keys(tmp_path, monkeypatch):
    import sys
    from . import open_pair_check as check
    receipt = tmp_path/"receipt.json"
    receipt.write_text('{"setup": {}, "setup": {}}', encoding="utf-8")
    monkeypatch.setattr(sys, "argv", ["check", "--receipt", str(receipt)])
    with pytest.raises(ValueError, match="duplicate JSON key"):
        check.main()


@pytest.mark.parametrize("defect", ["missing arm", "noninteger", "duplicate", "reordered receiver", "wrong bins"])
def test_raw_checker_rejects_corrupted_receipt(defect):
    from . import open_pair_check as check
    archive = json.loads(check.ARCHIVE.read_text())
    row = archive["appointments"][0]
    if defect == "missing arm":
        row["bins"].pop("stochastic_swap")
    elif defect == "noninteger":
        row["bins"]["dephased"][3] = True
    elif defect == "duplicate":
        archive["appointments"][-1] = row
    elif defect == "reordered receiver":
        archive["outcome_order"].reverse()
    else:
        row["bins"]["coherent"].reverse()  # normalized and physical, wrong allocation
    _, home, expected = check.independently_replay((0,), (1,))
    with pytest.raises((ValueError, ArithmeticError)):
        check.max_bin_error(check.archive_bins(archive, home), expected)
