"""Compare joint addressed messages across the same 25 native written seams."""

import argparse
from collections import Counter
from itertools import combinations
import json

from .cut_code import verify as verify_native_columns
from .experiment import HERE

PIVOTS = (0, 1, 2, 4, 8, 16)
PAYLOAD = tuple(i for i in range(25) if i not in PIVOTS)


def affine_mask(seed):
    """Evaluate one affine function on the 25 distinct public five-bit labels."""
    return sum((((seed >> 5) ^ (i & seed).bit_count()) & 1) << i for i in range(25))


def decode(word):
    """Identify and remove the affine pad using the six reference sign readings."""
    offset = word & 1
    slope = sum((((word >> (1 << j)) ^ offset) & 1) << j for j in range(5))
    clean = word ^ affine_mask((offset << 5) | slope)
    return sum(((clean >> i) & 1) << j for j, i in enumerate(PAYLOAD))


def run():
    count, _, full = verify_native_columns()
    assert count == 25 and full == 2**25
    survey = json.loads((HERE / "capacity.json").read_text())
    panel = next(p for p in survey["panels"]
                 if (p["level"], p["seed"], p["dilute"]) == (1, 280222, False))
    station = next(s for s in panel["stations"] if s["origin"] == 0)
    seams = sorted(panel["components"][str(station["component"])]["seams"])

    masks = [affine_mask(seed) for seed in range(64)]
    assert len(set(masks)) == 64 and all(decode(w) == 0 for w in masks)
    # Both encoder and decoder use XORs of linear parity checks. Testing
    # every payload generator and every pad verifies decoding on all messages.
    for j, i in enumerate(PAYLOAD):
        assert all(decode((1 << i) ^ w) == 1 << j for w in masks)

    triples = list(combinations(range(count), 3))
    for addresses in triples:
        hist = Counter(tuple((w >> i) & 1 for i in addresses) for w in masks)
        assert len(hist) == 8 and set(hist.values()) == {8}

    # Four collected records can already carry one source-dependent bit.
    first_parity = sum(1 << i for i in (0, 1, 2, 3))
    assert all((w & first_parity).bit_count() % 2 == 0 for w in masks)
    assert all(((w ^ (1 << 3)) & first_parity).bit_count() % 2 == 1 for w in masks)

    # In the second protocol the sender selects the parity of all 25 signs.
    # Removing any one site gives a bijection from each parity class onto
    # all 24-bit strings, since precisely one missing sign completes the word.
    for missing in range(count):
        visible = ((1 << count)-1) ^ (1 << missing)
        for i in range(count):
            if i != missing:
                even = (1 << i) | (1 << missing)
                odd = even ^ (1 << missing)
                assert even.bit_count() % 2 == 0 and odd.bit_count() % 2 == 1
                assert even & visible == odd & visible == 1 << i

    return dict(cut_seams=seams, affine_reference_seams=[seams[i] for i in PIVOTS],
                fully_checked_triples=len(triples), pad_words_per_message=len(masks),
                pad_words_per_three_sign_outcome=8, first_revealing_parity_seams=seams[:4],
                addressable_messages=1 << len(PAYLOAD), sender_decoupled_messages=1,
                parity_channel=dict(messages=2, proper_subcut_messages=1,
                                    sign_words_per_message=1 << (count-1)),
                scope="externally prepared messages at the write, common joint sign read and independent coins")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    opt = parser.parse_args()
    result = run()
    path = HERE / "cut_messages.json"
    if opt.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    elif result != json.loads(path.read_text()):
        raise ValueError("Joint message check differs from the archived receipt")
    print(f"{result['fully_checked_triples']} triple surveys agree; joint message counts "
          f"{result['addressable_messages']} and {result['sender_decoupled_messages']}.")


if __name__ == "__main__":
    main()
