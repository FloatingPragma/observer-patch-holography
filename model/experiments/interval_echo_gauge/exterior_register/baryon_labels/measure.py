"""Read two charge phases and compare all four-field operators at both stations."""

import argparse
import json
from pathlib import Path

import numpy as np

from .. import material
from . import apparatus, grammar

HERE = Path(__file__).resolve().parent


def run():
    geometry, signals = material.flights()
    # Separate sharp frames carry the internal occupation and the independently
    # initialized formal two-component index. The factorized, deterministic
    # channel supplies the complete unequal as well as equal frame pairs.
    _, headers = material.flights()
    packets = [
        {s: [row[s], headers[spin][s]] for s in ("sender", "receiver", "blocked")}
        for row in signals
        for spin in (0, 1)
    ]
    assert all(
        material.decode(frame) is None for p in packets for frame in p["blocked"]
    )
    connectors = apparatus.connectors(geometry)
    gates = apparatus.pair_gates()
    delivery = np.eye(32)[:, [material.decode(r["receiver"]) for r in signals]]
    np.testing.assert_array_equal(delivery, np.eye(32))
    observations, programs = apparatus.phase_banks(signals, gates, delivery)
    phases, radius = apparatus.phase_read(observations)
    signs, duplicates = apparatus.sorting_trials(signals, gates, delivery)
    exterior = grammar.exterior_signs(signs, True)
    unsigned = grammar.exterior_signs(signs, False)
    stations = []
    for side in ("sender", "receiver"):
        fields, components = grammar.field_census(packets, phases, side)
        votes, admitted = grammar.ballots(fields)
        assert len(votes) == 715 and len(admitted) == 8
        ranks = [
            grammar.serialize_contractions(
                row, grammar.contractions(row, fields, components, exterior, unsigned)
            )
            for row in admitted
        ]
        assert all(result["signed_rank"] == 1 for result in ranks)
        changed = {}
        for i, row in fields.items():
            words = {token[0] // 2 for token in components[i]}
            read = {phases[word, "opposite_L", side] for word in words}
            assert read == {-row[2]}
            changed[i] = (*row[:2], read.pop(), *row[3:])
        control_votes, same = grammar.ballots(changed)
        assert same == admitted
        laws = []
        for actual, rival in zip(votes, control_votes, strict=True):
            if all(actual[2]):
                b, ell = actual[1][1:3]
                assert abs(ell) == 1 and b == 3 * ell and rival[1][1] != 3 * rival[1][2]
                laws.append([actual[0], b, ell, rival[1][2]])
        stations.append(
            dict(
                fields=[list(fields[i]) for i in range(10)],
                ballots=votes,
                contractions=ranks,
                charge_law_control=laws,
            )
        )
    assert stations[0] == stations[1]
    return dict(
        schema="Two-station Weyl field labels 1",
        owner_cells=geometry["owners"],
        shots=apparatus.SHOTS,
        seed=apparatus.SEED,
        simultaneous_probability_radius=radius,
        controlled_phase_routes=connectors,
        parity_programs=programs,
        sharp_packets=packets,
        phase_banks=observations,
        phase_columns=[
            "word",
            "phase program",
            "station",
            "axis",
            "nonzero raw counts",
        ],
        labels=[[w, tag, side, value] for (w, tag, side), value in phases.items()],
        wedge_reference_readings=signs,
        duplicate_source_reads=duplicates,
        names=grammar.NAMES,
        field_columns=["6Y", "3B", "L", "triality", "weak", "left"],
        both_stations=stations[0],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "results.json"
    if args.write:
        with path.open("x") as stream:
            # Put each acquisition or ballot on one line.
            for i, (key, value) in enumerate(result.items()):
                if i == 0:
                    stream.write("{\n")
                if isinstance(value, list) and len(value) > 9:
                    body = (
                        "[\n"
                        + ",\n".join("    " + json.dumps(v) for v in value)
                        + "\n  ]"
                    )
                elif key == "both_stations":
                    body = (
                        "{\n"
                        + ",\n".join(
                            "    "
                            + json.dumps(k)
                            + ": [\n"
                            + ",\n".join("      " + json.dumps(v) for v in group)
                            + "\n    ]"
                            for k, group in value.items()
                        )
                        + "\n  }"
                    )
                else:
                    body = json.dumps(value)
                stream.write(
                    "  "
                    + json.dumps(key)
                    + ": "
                    + body
                    + (",\n" if i < len(result) - 1 else "\n}\n")
                )
    elif json.loads(path.read_text()) != json.loads(json.dumps(result)):
        raise ValueError("The signed charge and operator bank differs")
    print(
        "Two stations, ten signed fields, 715 ballots and eight rank-one contraction spaces"
    )


if __name__ == "__main__":
    main()
