"""Exercise the publication identity, complete count words and conditional hand limits."""
import io
import json
import zipfile

import numpy as np
import pytest

from . import loaded_clock as clock


def _archive(**columns):
    out = io.BytesIO()
    np.savez_compressed(out, **columns)
    return out.getvalue()


def test_published_source_and_signed_counts_replay():
    result = clock.experiment(compare_seed=True)
    assert set(result["same_seed_replay"]) == {"differing_rows", "max_bin_difference"}
    assert result["archive_sha256"] == clock.PUBLISHED_SHA256
    assert result["count_source"] == "published simulated acquisition"
    assert result["launches"] == 1_188 * (1 << 24)
    assert result["informative_count_hands"] == 330 and result["full_circle_count_hands"] == 0
    assert result["dephased_disks"] == result["dephased_disks_including_zero"] == 66
    assert not result["failures"]


def test_one_transferred_count_and_alternate_reencoding(tmp_path):
    protocol = json.loads(clock.BOOK.read_text(encoding="utf-8"))
    keys, words, _ = clock.design(protocol)
    weights = np.stack([clock._grid(word) for word in words])
    cols = dict(keys=keys, weights=weights, counts=clock.counts_for(keys, weights, protocol))
    local = tmp_path / "counts.npz"
    # A reencoded acquisition on this backend must identify itself as local.
    local.write_bytes(_archive(**dict(reversed(cols.items()))))
    with pytest.raises(ValueError, match="SHA-256"):
        clock.experiment(local)
    with pytest.raises(ValueError, match="SHA-256"):
        clock.experiment(local, compare_seed=True)
    local_result = clock.experiment(local, seed_replay=True)
    assert local_result["count_source"] == "local simulated acquisition"
    assert local_result["same_seed_replay"] == {"differing_rows": 0, "max_bin_difference": 0}
    cols["counts"][0, 0] -= 1
    cols["counts"][0, 1] += 1
    local.write_bytes(_archive(**cols))
    with pytest.raises(ValueError, match="seeded acquisition"):
        clock.experiment(local, seed_replay=True)
    # A nearby integer source may pass native floating-point replay. Authenticating
    # an unpublished acquisition requires the exact grid of the asserted backend.
    cols["weights"][0, 0] -= 1
    cols["weights"][0, 1] += 1
    local.write_bytes(_archive(**cols))
    with pytest.raises(ValueError, match="local acquisition grid"):
        clock.experiment(local, seed_replay=True)


def test_ambiguous_or_impossible_source_words_are_rejected():
    keys = np.zeros((1, 5), dtype="<i4")
    words = np.array([[0., .5, 0., 0., .5, 0., 0.]])
    weights = np.array(words * clock.GRID_SIZE, dtype="<i8")
    counts = np.array(words * 100, dtype="<i8")
    def check(k=keys, w=weights, c=counts):
        return clock._records(_archive(keys=k, weights=w, counts=c), keys, words, 100)
    check()
    for field, value in (("k", keys.astype(float)), ("w", weights.astype(float)), ("c", counts.astype(float))):
        with pytest.raises(ValueError):
            check(**{field: value})
    impossible = counts.copy()
    impossible[0, :2] += [1, -1]
    with pytest.raises(ValueError, match="counts"):
        check(c=impossible)
    # Nine-quintillion positive counts can wrap a signed row sum to 100.
    overflow = np.array([[2**62, 2**62, 2**62, 2**62, 100, 0, 0]], dtype="<i8")
    assert overflow.sum() == 100
    with pytest.raises(ValueError, match="counts"):
        check(c=overflow)
    out = io.BytesIO(_archive(keys=keys, weights=weights, counts=counts))
    with zipfile.ZipFile(out, "a") as zipped, pytest.warns(UserWarning, match="Duplicate"):
        zipped.writestr("counts.npy", zipped.read("counts.npy"))
    with pytest.raises(ValueError, match="columns"):
        clock._records(out.getvalue(), keys, words, 100)


def test_count_only_interval_handles_empty_unresolved_and_small_apertures():
    for counts in ((0, 0, 0, 0), (100, 100, 100, 100)):
        assert clock.hand_interval(counts, 10.)[1] == np.pi
    counts = (10_000, 0, 0, 0)
    angle, ideal = clock.hand_interval(counts, 10.)
    assert angle == 0
    _, native = clock.hand_interval(counts, 10., trials=100_000, cell_radius=.001)
    assert 0 < ideal < native < np.pi / 2
    assert clock.hand_interval(counts, 10., trials=10**8, cell_radius=.001)[1] == np.pi
