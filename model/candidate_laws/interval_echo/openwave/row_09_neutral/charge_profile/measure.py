"""Resolve neutral two-port records by their measured charge density at two spatial scales."""

import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path

from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge import (
    neutral_ladder,
)

HERE = Path(__file__).resolve().parent
SOURCE = Path(neutral_ladder.__file__).parent / "results.json"
FILLED = 1
TAIL = ("no arrival", "source refused")


def signature(word):
    hub, target, donor = (n - FILLED for n in word)
    return (hub + target, target - hub, donor)


def describe(signature_value):
    total, difference, _ = signature_value
    if total:
        return "charged"
    return "resolved_pair" if difference else "locally_filled"


def charge_power(total, difference):
    # |h + t exp(i k a)|^2, at k a = 0, pi/2 and pi, for real signed loads.
    return [total**2, (total**2 + difference**2) / 2, difference**2]


def observed(counts, words, starts):
    if len(counts) != len(words) + len(TAIL) or sum(counts) != starts:
        raise ValueError("A complete original-start bank is required")
    classes = Counter()
    orientations = Counter()
    for n, word in zip(counts[: len(words)], words, strict=True):
        total, diff, stock = signature(word)
        if n:
            classes[describe((total, diff, stock))] += n
            orientations[(total, diff, stock)] += n
    return dict(
        classifications={
            k: classes[k] for k in ("charged", "locally_filled", "resolved_pair")
        },
        signed_records=[[*s, n] for s, n in sorted(orientations.items())],
        no_arrival=counts[-2],
        source_refused=counts[-1],
        starts=starts,
    )


def run():
    source = json.loads(SOURCE.read_text())
    words = [tuple(w) for w in source["ternary_word_order"]]
    assert set(words) == set(product(range(3), repeat=3))
    assert list(TAIL) == source["load_count_tail"]
    complete = []
    for word in words:
        total, diff, donor = signature(word)
        assert (total - diff, total + diff) == tuple(2 * (v - 1) for v in word[:2])
        complete.append(
            dict(
                word=word,
                net=total,
                signed_difference=diff,
                donor=donor,
                charge_power=charge_power(total, diff),
                outcome=describe((total, diff, donor)),
            )
        )
    banks = []
    totals = {}
    for experiment in source["experiments"]:
        combined = Counter()
        for basis, column, counts in experiment["load_banks"]:
            data = observed(counts, words, source["starts_per_bank"])
            banks.append(
                dict(
                    sign=experiment["sign"],
                    setting=experiment["setting"],
                    basis=basis,
                    column=column,
                    **data,
                )
            )
            combined.update(data["classifications"])
            combined.update(
                {key: data[key] for key in ("no_arrival", "source_refused", "starts")}
            )
        totals[f"{experiment['sign']}/{experiment['setting']}"] = dict(combined)
    local_controls = {}
    for row in source["local_controls"]:
        # These are separately prepared local controls, which start at the apparatus.
        data = observed(
            [*row["ternary_counts"], 0, 0], words, sum(row["ternary_counts"])
        )
        local_controls[f"{row['sign']}/{row['setting']}"] = data
    return dict(
        basis_order=["net excess", "target minus hub excess", "donor excess"],
        phase_step_radians=["0", "pi/2", "pi"],
        complete_endpoint_catalogue=complete,
        full_signed_classes=len(set(signature(w) for w in words)),
        unoriented_classes=len({(q, abs(d), s) for q, d, s in map(signature, words)}),
        detected_banks=banks,
        by_sign_and_setting=totals,
        independently_prepared_controls=local_controls,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    result = json.dumps(run(), indent=2) + "\n"
    if parser.parse_args().write:
        (HERE / "results.json").write_text(result)
    else:
        print(result, end="")
