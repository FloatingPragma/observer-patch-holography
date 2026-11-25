import json
from itertools import product

from .measure import HERE, charge_power, describe, run, signature


def test_complete_charge_profile_and_orientation():
    for h, t, s in product(range(3), repeat=3):
        net, direction, stock = signature((h, t, s))
        assert signature((t, h, s)) == (net, -direction, stock)
        assert charge_power(net, direction) == charge_power(net, -direction)
        if not net:
            assert charge_power(net, direction)[-1] == 4 * (h != t)
            assert (describe((net, direction, stock)) == "resolved_pair") == (h != t)


def test_archived_banks_and_missed_launches():
    current = json.loads(json.dumps(run()))
    assert current == json.loads((HERE / "results.json").read_text())
    assert (current["full_signed_classes"], current["unoriented_classes"]) == (27, 18)
    for row in current["detected_banks"]:
        c = row["classifications"]
        assert (
            sum(c.values()) + row["source_refused"] + row["no_arrival"] == row["starts"]
        )
        assert c["locally_filled"] == 0
        if row["setting"] == "full":
            assert c["resolved_pair"] > 0 and c["charged"] == 0
        if row["setting"] == "off":
            assert c["charged"] > 0 and c["resolved_pair"] == 0
    for sign in (-1, 1):
        blocked = current["independently_prepared_controls"][f"{sign}/depleted_donor"]
        assert blocked["classifications"]["resolved_pair"] == 0
