from fractions import Fraction
import json

import numpy as np

from .measure import HERE, acquire, twice_generator


def test_all_histories_and_independent_factorization():
    result = acquire()
    assert json.loads(json.dumps(result)) == json.loads(
        (HERE / "frequencies.json").read_text()
    )
    connected = 0
    for bank in result["histories"]:
        if bank["receiver_factorizes"]:
            x, y = bank["diagonal"][2], bank["diagonal"][1]
            two = np.kron(twice_generator([0, x]), np.eye(2)) + np.kron(
                np.eye(2), twice_generator([0, y])
            )
            np.testing.assert_array_equal(two, twice_generator(bank["diagonal"]))
            assert bank["excited_rest_candidates"]
        else:
            connected += bank["both_links_conduct"] + bank["other"]
            assert Fraction(bank["detuning_lower_rational"]) > Fraction("0.19393656")
    assert connected == 256
    assert sum(b["both_links_conduct"] for b in result["histories"]) == 128
    assert sum(b["other"] for b in result["histories"]) == 384
    assert Fraction(result["crossing"]["detuning_lower_rational"]) > Fraction(
        "0.00644841"
    )
