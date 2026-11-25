import gzip
import json

import numpy as np

from . import measure as m


def test_full_native_detector_on_coherent_and_classical_memory_inputs():
    for t, p in ((0, 0), (np.pi / 2, 0.5), (np.pi, 1)):
        head = m.detector(t)
        np.testing.assert_allclose(head[1], [[p, 0], [0, 0]], atol=2e-13)
        for phase in (0, 0.3, 1.7):
            z = np.array([1, np.exp(1j * phase)]) / np.sqrt(2)
            np.testing.assert_allclose(z.conj() @ head[1] @ z, p / 2, atol=2e-13)


def test_all_feedforward_paths_and_public_counters():
    receipt = json.loads(gzip.decompress(m.RECEIPT.read_bytes()))
    for level, board in enumerate(receipt["boards"]):
        real = m.release.federation(level)
        a, b = board["seam"]
        tree = board["next_hop_to_meter_after_public_address"]
        for start in range(real.n):
            seen = set()
            x = start
            while x != a:
                y = tree[x]
                assert x not in seen and y in real.neighbours[x] and {x, y} != {a, b}
                seen.add(x)
                x = y
        for p in board["panels"]:
            for label, r in p["measurements"].items():
                hist = r["measured_words_per_position"]
                assert (
                    sum(v for _, _, v in hist) + r["spectral_refusal"]
                    == receipt["shots"]
                )
                assert sum(v for _, w, v in hist if w == 1) == r["flip_count"]
                if label == "native_pi" and p["operation"] == "one_native_echo":
                    assert 0.018 < r["flip_count"] / r["admitted"] < 0.0191
