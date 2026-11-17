#!/usr/bin/env python3
"""Bind a native finite recorder to the existing Pro5 capacity evaluator."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PUBLIC_SOURCES = (
    "code/capacity_readback/correctable_public_record_capacity.py",
    "code/capacity_readback/checkpoint_channels.py",
    "code/capacity_readback/public_record_csp.py",
)


def dependencies():
    from . import native as source
    capacity_path = ROOT / "code/capacity_readback"
    if str(capacity_path) not in sys.path:
        sys.path.insert(0, str(capacity_path))
    import correctable_public_record_capacity as capacity
    return source, capacity



def word_name(word):
    return json.dumps(word, separators=(",", ":"))


def enumerate_source(source, binary):
    states = source.states(binary=binary)
    words = Counter()
    witnesses = {}
    prefixes = {()}
    transitions = set()
    histories = []
    for initial in states:
        for _, first in source.tick(source.State(initial)):
            prefixes.add(first.prefix)
            transitions.add(((), first.prefix[0], first.prefix))
            for _, last in source.tick(first):
                if source.tick(last) != ((source.Q(1), last),):
                    raise ValueError("The declared terminal continuation does not halt")
                words[last.prefix] += 1
                prefixes.add(last.prefix)
                transitions.add((first.prefix, last.prefix[1], last.prefix))
                witnesses.setdefault(last.prefix, [
                    f"initial_loads:{word_name(initial)}",
                    f"event_1:{word_name(first.prefix[0])}",
                    f"event_2_seam:{word_name(source.selected_seam(first))}",
                    f"event_2:{word_name(last.prefix[1])}",
                ])
                histories.append((initial, first, last))
    return words, witnesses, prefixes, transitions, histories


def full_tape_packet(source, capacity, words, witnesses, name):
    # All 36 terminal tape atoms exist. Reachability selects 11 or 17 of them.
    atoms = [word_name(w) for w in product(source.LETTERS, repeat=2)]
    ids = {w: capacity.section_id({"observer": word_name(w)}) for w in words}
    letters = {r: i for i, r in enumerate(source.LETTERS)}
    projections = {}
    for word, record_id in ids.items():
        sector = letters[word[0]] * 6 + letters[word[1]]
        start = 78 * (1 + 6 + sector)
        projections[record_id] = list(range(start, start + 78))
    packet = {
        "terminal_id": name,
        "capacity_dimension": 3354,
        "observers": {"observer": atoms},
        "interfaces": [],
        "reachability_witnesses": {ids[w]: witnesses[w] for w in words},
        "publicness_policy": [["observer"]],
        "continuation_manifest_complete": True,
        "global_checkpoint_kernels": [{
            "authorized_observers": ["observer"],
            "continuation_id": "terminal_halt_then_exact_full_tape_read",
            "rows": {ids[w]: {word_name(w): 1} for w in words},
        }],
        # There is one observer, so the joint channel is its sole marginal.
        "local_marginal_consistency_passed": True,
        "projection_supports": projections,
    }
    return packet


def readout_capacity(capacity, words, output):
    alphabet = [word_name(w) for w in sorted(words)]
    identity = {"rows": {x: {x: 1} for x in alphabet}}
    readout = {"rows": {word_name(w): {output(w): 1} for w in words}}
    graph = capacity.compound_confusability_graph(alphabet, [identity, readout])
    code = capacity.maximum_independent_set(graph)
    return {"M0": len(code), "code": code,
            "confusability_edges": sum(map(len, graph.values())) // 2}


def encode_prefixes(prefixes, transitions):
    ordered = sorted(prefixes, key=lambda p: (len(p), p))
    codes = {prefix: i for i, prefix in enumerate(ordered)}
    # This is an encoding of reached classical histories, not a quantum map.
    table = [{"from": codes[before], "letter": list(letter), "to": codes[after]}
             for before, letter, after in sorted(transitions)]
    return {"states": len(ordered), "minimum_fixed_width_bits": (len(ordered) - 1).bit_length(),
            "codes": [{"code": codes[p], "prefix": [list(r) for r in p]} for p in ordered],
            "append_table": table,
            "halt_codes": [codes[p] for p in ordered if len(p) == 2]}


def compile_controller(source, histories, encoding):
    """Compile a reached classical support; no unseen state is admitted."""
    decode = {row["code"]: tuple(tuple(r) for r in row["prefix"])
              for row in encoding["codes"]}
    encode = {prefix: code for code, prefix in decode.items()}
    admitted = set()
    for initial, first, last in histories:
        admitted.update(((encode[()], initial), (encode[first.prefix], first.loads),
                         (encode[last.prefix], last.loads)))
    return {"source": source, "decode": decode, "encode": encode, "admitted": admitted}


def compact_view(compiled, code):
    if type(code) is not int or code not in compiled["decode"]:
        raise ValueError("Unknown compact controller code")
    prefix = compiled["decode"][code]
    # Deliberately no loads are available to either control or forecast.
    class PrefixOnly:
        pass
    state = PrefixOnly()
    state.prefix, state.phase = prefix, len(prefix)
    source = compiled["source"]
    return {"phase": state.phase, "prefix": prefix,
            "seam": source.selected_seam(state),
            "forecast": source.issue_forecast(state) if state.phase == 1 else None}


def compact_tick(compiled, code, unchanged_source_loads):
    view = compact_view(compiled, code)
    loads = tuple(unchanged_source_loads)
    if (len(loads) != 12 or any(type(x) is not int for x in loads)
            or (code, loads) not in compiled["admitted"]):
        raise ValueError("Outside the admitted reachable classical support")
    source = compiled["source"]
    state = source.State(loads, view["prefix"])
    rows = []
    for weight, after in source.tick(state):
        after_code = compiled["encode"].get(after.prefix)
        if after_code is None or (after_code, after.loads) not in compiled["admitted"]:
            raise ValueError("Continuation left the compiled classical support")
        rows.append((weight, after_code, after.loads))
    return tuple(rows)


def run():
    source, capacity = dependencies()
    families = {}
    for binary, name, expected_words, expected_prefixes in (
            (True, "binary_66", 11, 16), (False, "full_78", 17, 24)):
        words, witnesses, prefixes, transitions, histories = enumerate_source(source, binary)
        if len(words) != expected_words or len(prefixes) != expected_prefixes:
            raise ValueError("Reachable record inventory changed")
        packet = full_tape_packet(source, capacity, words, witnesses, name)
        result = capacity.evaluate_terminal(packet)
        if (result["status"] != "PASS" or result["exact_zero_error_capacity"] != len(words)
                or result["saturation_passed"]):
            raise ValueError("Unexpected existing-evaluator capacity result")
        encoding = encode_prefixes(prefixes, transitions)
        compiled = compile_controller(source, histories, encoding)
        codes = {tuple(tuple(r) for r in row["prefix"]): row["code"] for row in encoding["codes"]}
        append = {(row["from"], tuple(row["letter"])): row["to"] for row in encoding["append_table"]}
        for initial, first, last in histories:
            if (append[codes[()], first.prefix[0]] != codes[first.prefix]
                    or append[codes[first.prefix], last.prefix[1]] != codes[last.prefix]):
                raise ValueError("Compressed controller history fails append replay")
            restored = source.State.restore_controller(
                {"phase": 1, "prefix": [list(r) for r in first.prefix]}, first.loads)
            if source.tick(restored) != source.tick(first):
                raise ValueError("Decoded controller fails source continuation")
            for state in (source.State(initial), first, last):
                view = compact_view(compiled, codes[state.prefix])
                expected_forecast = source.issue_forecast(state) if state.phase == 1 else None
                expected = tuple((weight, codes[after.prefix], after.loads)
                                 for weight, after in source.tick(state))
                if (view["phase"] != state.phase or view["seam"] != source.selected_seam(state)
                        or view["forecast"] != expected_forecast
                        or compact_tick(compiled, codes[state.prefix], state.loads) != expected):
                    raise ValueError("Compact execution does not preserve the native history")
        views = {
            "full_tape": lambda w: word_name(w),
            "first_letter_only": lambda w: word_name(w[0]),
            "second_letter_only": lambda w: word_name(w[1]),
            "chosen_seam_only": lambda w: "F" if w[0][1] > 0 else "E",
            "complete_erasure": lambda w: "blank",
        }
        if binary:
            views["forecast_only"] = lambda w: str(source.FORECAST[w[0]])
            views["forecast_and_chosen_seam"] = lambda w: str((source.FORECAST[w[0]], w[0][1] > 0))
        controls = {view: readout_capacity(capacity, words, function) for view, function in views.items()}
        first_words = sorted({(w[0],) for w in words})
        first_capacity = readout_capacity(capacity, first_words, word_name)
        families[name] = {
            "initial_states": 66 if binary else 78,
            "positive_complete_histories": len(histories),
            "reachable_first_records": len(first_words),
            "reachable_terminal_records": len(words),
            "terminal_word_counts": [{"word": [list(r) for r in w], "histories": words[w]}
                                     for w in sorted(words)],
            "first_record_capacity": first_capacity["M0"],
            "terminal_capacity_nats": f"ln({len(words)})",
            "existing_evaluator_result": result,
            "projection_rank_each": 78,
            "orthogonal_projection_rank_sum": 78 * len(words),
            "readout_controls": controls,
            "classical_prefix_encoding": encoding,
            "compact_execution_checked_histories": len(histories),
            "admitted_controller_source_pairs": len(compiled["admitted"]),
        }
    return {
        "schema": "oph.native_observer_capacity_join.v1",
        "source_pins": bound_pins({name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                        for name in PUBLIC_SOURCES}),
        "allocation": {"working_dimension": 78, "terminal_tape_words": 36,
                       "terminal_sector_dimension": 2808, "all_phase_record_states": 43,
                       "total_hilbert_dimension": 3354,
                       "naive_separate_phase_and_slot_bits": 8,
                       "allocated_classical_record_index_bits": 6},
        "families": families,
        "capacity_contract": {
            "publicness": "One observer reads the complete classical tape.",
            "continuations": "Any number of terminal halt ticks and exact classical rereads.",
            "coarse_controls": "Diagnostic channels that discard distinctions, not new native operations.",
            "closure_scope": "One bounded finite packet, not a complete OPH terminal fiber or cosmic fixed point.",
            "units": "N_record = ln M0 in nats; slot count, raw dimension and actual Shannon entropy differ.",
            "compression": "Reached classical controller histories only; source state, quantum law and arbitrary-input isometry unchanged.",
        },
    }


def canonical(result):
    return json.dumps(result, indent=2, sort_keys=True) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = run()
    encoded = canonical(result)
    if args.write:
        with args.write.open("x") as output:
            output.write(encoded)
    if args.check and args.check.read_text() != encoded:
        raise SystemExit("Capacity join receipt mismatch")
    print(json.dumps({"M0_binary": result["families"]["binary_66"]["reachable_terminal_records"],
                      "M0_full": result["families"]["full_78"]["reachable_terminal_records"],
                      "receipt_sha256": hashlib.sha256(encoded.encode()).hexdigest()}))
