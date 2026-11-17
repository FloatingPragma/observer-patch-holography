"""Native blank-bit routing realization and charged endpoint-work audit."""
import json
from fractions import Fraction as F
from pathlib import Path
from native import CALLS, EDGES, endpoint, logprob

base = Path(__file__).resolve().parent
frozen = json.loads((base / "CONFIRM_RECEIPT.json").read_text())
selectors, blank_calls, target_calls = [], {"native": 0, "swap": 0}, {"native": 0, "swap": 0}
for m in (4, 8):
    bits = m.bit_length()
    for law in ("native", "swap"):
        catalog = frozen["instrument_catalog"][f"{m}:{law}"]["target"]
        for slot in range(3 * m):
            v, k = divmod(slot, m)
            words, fine = [], []
            for word in range(2 * m):
                tape = [(word >> shift) & 1 for shift in reversed(range(bits))]
                currents = []
                for coin in tape:
                    out = endpoint(0, 0, coin, law)
                    blank_calls[law] += 1
                    assert out == (0, 0)
                    currents.append(-out[0])
                ticket, edge_bit = word // 2, word % 2
                star = v if ticket < m - k else (v + 1) % 3
                a, b = [edge for edge in EDGES if star in edge][edge_bit]
                words.append({"native_coin_word": tape, "vacuum_current_word": currents,
                              "probability": str(F(1, 2 * m)), "ticket": ticket,
                              "star": star, "edge_bit": edge_bit, "decoded_address": [a, b]})
                for coin in (0, 1):
                    flags = [[[0] * 3 for _ in range(3)] for _ in range(3)]
                    for source in range(3):
                        state = [int(i == source) for i in range(3)]
                        before = state[a]
                        state[a], state[b] = endpoint(state[a], state[b], coin, law)
                        target_calls[law] += 1
                        flags[before - state[a] + 1][state.index(1)][source] = 1
                    fine.extend(flags)
            assert sum(F(w["probability"]) for w in words) == 1
            prefix, ops = [[0] * 3 for _ in range(3)], []
            for i, matrix in enumerate(fine):
                if i:
                    ops.append(([[prefix[r][c] - i * matrix[r][c] for c in range(3)] for r in range(3)],
                                4 * m * i * (i + 1)))
                prefix = [[prefix[r][c] + matrix[r][c] for c in range(3)] for r in range(3)]
            ops.insert(0, (prefix, 48 * m * m))
            assert len(ops) == len(catalog[slot])
            for (matrix, denominator), prior in zip(ops, catalog[slot]):
                assert matrix == [list(map(int, row)) for row in prior["numerator"]]
                assert denominator == prior["squared_denominator"]
            selectors.append({"M": m, "law": law, "slot": slot, "blank_events_per_source_event": bits,
                              "all_words": words, "full_catalog_equality": True,
                              "selector_rejection_probability": "0"})
costs = []
for case in frozen["cases"]:
    m, law, turns = case["M"], case["law"], case["turns"]
    bits = m.bit_length()
    expected = F(2)
    for row in case["records"]:
        reach = F(row["reach"])
        if row["kind"] == "marker":
            expected += reach * (F(255, 128) if law == "native" else 1)
        elif row["kind"] == "target":
            expected += reach * (bits + 1)
        else:
            expected += reach
    success = F(case["success"])
    frames = 3 * m * turns
    nonmarker = 2 + 4 * frames * (bits + 1) + frames
    costs.append({"M": m, "law": law, "turns": turns, "direction": case["direction"],
                  "expected_native_endpoint_opportunities_per_attempt": str(expected),
                  "expected_per_attempt_float": float(expected),
                  "log10_expected_native_endpoint_opportunities_per_herald": logprob(expected / success),
                  "successful_path_endpoint_count_bounds": [nonmarker + 2 * frames, nonmarker + 16 * frames],
                  "declared_non_native_operations": "uniform3-way selector states, coherent erasers/path gates/read/reset; no energy conversion"})
assert all(CALLS[law] == blank_calls[law] + target_calls[law] for law in CALLS)
print(json.dumps({"candidate": "A-native-subdivided-observer-rotation-v3", "passed": True,
                  "selector_tables": selectors, "endpoint_costs": costs,
                  "literal_blank_endpoint_census_calls": blank_calls,
                  "literal_target_endpoint_census_calls": target_calls,
                  "literal_total_endpoint_calls": CALLS,
                  "original_run_counts_unchanged": True}, indent=2))
