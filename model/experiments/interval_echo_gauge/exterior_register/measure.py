"""Replay the two-station five-mode exterior grammar and every held control."""

import argparse
import json
from pathlib import Path

from . import coherent, material, selection


HERE = Path(__file__).resolve().parent


def run():
    geometry, records = material.flights()
    rows = {}
    ballots = {}
    for owner in ("sender", "receiver"):
        measured = [material.decode(r[owner]) for r in records]
        assert measured == list(range(32))
        rows[owner] = selection.census(measured)
        ballots[owner] = selection.scan(rows[owner])
    assert all(material.decode(r["blocked"]) is None for r in records)
    assert (
        rows["sender"] == rows["receiver"] and ballots["sender"] == ballots["receiver"]
    )
    winners = selection.survivors(ballots["receiver"])
    assert winners == [227, 796]
    challenges = [[m, i, m ^ (1 << i)] for m in winners for i in range(10)]
    assert all(changed not in winners for _, _, changed in challenges)
    signs, refusals = coherent.insertion_trials(records)
    for i, j, base, clear, with_sign, ordinary in signs:
        assert clear == "0100"
        assert material.decode(ordinary) == base | (1 << j)
        assert material.decode(with_sign) == base | (1 << j) | (1 << i)
    # Charge uniquely keys the ten rows in the theorem's independent enumeration.
    lean_charges = [-2, 3, -4, 1, 6, 4, -1, -6, 2, -3]
    component_map = [lean_charges.index(row["charge"]) for row in rows["sender"]]
    assert sorted(component_map) == list(range(10))
    return dict(
        schema="Two-station exterior register 1",
        geometry=geometry,
        sharp_flights=records,
        meter="two agreement bits for each rail against loads 0,2, in park/marked dimer order",
        local_census_at_both_stations=rows["sender"],
        station_to_lean_component_indices=component_map,
        ballots_at_both_stations=ballots["sender"],
        ballot_columns=["mask", "nonempty", "chiral", "gravity", "su3", "su2", "cube"],
        survivors=winners,
        held_one_row_changes=challenges,
        insertion_columns=[
            "i",
            "j",
            "spectator_word",
            "insertion_source",
            "receiver_exterior",
            "receiver_ordinary",
        ],
        unit_probability_insertion_records=signs,
        refused_additions=refusals,
        coherent_block_connectors=coherent.connector_routes(geometry),
        rotations_on_each_connector=coherent.rotations(),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "results.json"
    if args.write:
        with path.open("x") as out:
            out.write("{\n")
            for i, (key, value) in enumerate(result.items()):
                suffix = "," if i < len(result) - 1 else ""
                if isinstance(value, list) and len(value) > 9:
                    serialized = (
                        "[\n"
                        + ",\n".join("    " + json.dumps(v) for v in value)
                        + "\n  ]"
                    )
                else:
                    serialized = json.dumps(value)
                out.write("  " + json.dumps(key) + ": " + serialized + suffix + "\n")
            out.write("}\n")
    elif json.loads(path.read_text()) != json.loads(json.dumps(result)):
        raise ValueError("Measured records differ from the archived bank")
    print(
        "32 delivered, 32 blocked, 1024 ballots, 80 paired sign reads, 80 refused insertions"
    )


if __name__ == "__main__":
    main()
