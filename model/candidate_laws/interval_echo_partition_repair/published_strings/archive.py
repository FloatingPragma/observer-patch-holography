"""Read-only source preparation from the original, already posted clock record."""
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PUBLISHED = HERE.parents[2]/"experiments"/"interval_echo_clock"/"publication.json"
PAIRS = ((0, 0), (0, 1), (1, 0), (1, 1))


def posted_words():
    """Read one pair of actual posted bits at each service of the archived journey."""
    bank = json.loads(PUBLISHED.read_text())
    row, = [s for s in bank["specimens"]
            if (s["cells"], s["width"], s["seed"]) == (80, 3, 31)]
    words = tuple(zip(row["moving_public"], row["home_public"], strict=True))
    if len(words) != sum(bank["schedule"]) or set(words) != set(PAIRS):
        raise ValueError("expected all four posted classes on the registered journey")
    return words


def address_lotteries(words):
    """Return exact new-sector probabilities for each *original service address*."""
    n = len(words)
    old = Counter(words)
    if set(old) != set(PAIRS):
        raise ValueError("this matched-singles control requires all four observed classes")
    traveller = old[1, 0]+old[1, 1]
    home = old[0, 1]+old[1, 1]
    if traveller+home > n:
        raise ValueError("zero coincidence at these marginal rates is impossible")
    # The zero-coincidence coupling at these singles is unique. Each required
    # class is implemented using same-service words, without permuting partners.
    classes = {"blank": (n, 0, 0, 0),
               "ordinary": (n-traveller-home, home, traveller, 0),
               "held": tuple(old[k] for k in PAIRS)}
    out = {}
    for name, masses in classes.items():
        weights = {k: Fraction(c, n*old[k]) for k, c in zip(PAIRS, masses, strict=True)}
        out[name] = tuple(weights[pair] for pair in words)
    return out


def source_densities():
    """Mean traveller post and same-service double post in each read-only lottery."""
    words = posted_words()
    return {k: (float(sum(w*t for w, (t, _) in zip(law, words, strict=True))),
                float(sum(w*t*h for w, (t, h) in zip(law, words, strict=True))))
            for k, law in address_lotteries(words).items()}


def refine_published_blocks(blocks, kernel):
    """Split already-classical transcript blocks, including any quantum reference.

    kernel[transcript] is a normalized distribution of the added record.
    An old classical/quantum observer that discards the new label recovers the
    original unnormalized operator for each full transcript. No native call
    occurs here. This helper is also valid for arbitrary reference dimensions.
    """
    refined = {}
    for transcript, block in blocks.items():
        weights = kernel[transcript]
        if min(weights.values()) < 0 or sum(weights.values()) != 1:
            raise ValueError("expected a normalized nonnegative archive kernel")
        for label, weight in weights.items():
            refined[transcript, label] = float(weight)*block
    return refined
