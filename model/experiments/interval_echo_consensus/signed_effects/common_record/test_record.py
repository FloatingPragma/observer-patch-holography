import gzip
import hashlib
import json

import numpy as np
import pytest

from . import verify


def test_independent_reader_rejects_spliced_tag_with_corrected_digests():
    base = verify.HERE
    receipt = json.loads((base / "results.json").read_bytes())
    stream = gzip.decompress((base / "records.bin.gz").read_bytes())
    assert verify.verify(stream, receipt) == 65536
    rows = np.frombuffer(stream, dtype=np.uint8).reshape(-1, 6).copy()
    admitted = rows[:, 1] != 255
    tags = rows[admitted, 1]
    rows[admitted, 1] = np.roll(tags, 1)
    assert sorted(tags) == sorted(rows[admitted, 1])
    assert np.count_nonzero(tags != rows[admitted, 1]) == receipt["shifted_tag_disagreements_preserving_marginal"]
    edited = rows.tobytes()
    receipt["uncompressed_trace_sha256"] = hashlib.sha256(edited).hexdigest()
    del receipt["run_sha256"]
    receipt["run_sha256"] = hashlib.sha256(json.dumps(receipt, separators=(",", ":"), sort_keys=True).encode()).hexdigest()
    with pytest.raises(ValueError, match="independent source/operation replay"):
        verify.verify(edited, receipt)


@pytest.mark.parametrize("change", ["count", "operation"])
def test_recomputed_manifest_cannot_relabel_the_archived_run(change):
    base = verify.HERE
    receipt = json.loads((base / "results.json").read_bytes())
    stream = gzip.decompress((base / "records.bin.gz").read_bytes())
    if change == "count":
        receipt["refused_launches"] = 1
    else:
        receipt["word_sha256_by_context_including_all_raw_branches"]["phase"] = "0" * 64
    del receipt["run_sha256"]
    receipt["run_sha256"] = hashlib.sha256(json.dumps(receipt, separators=(",", ":"), sort_keys=True).encode()).hexdigest()
    with pytest.raises(ValueError, match="refusal count differs|archived native run"):
        verify.verify(stream, receipt)
