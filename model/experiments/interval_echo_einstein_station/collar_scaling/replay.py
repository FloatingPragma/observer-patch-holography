"""Reproduce collar samples and bound the effect of finite sampler arithmetic."""
from hashlib import sha256
import json

import numpy as np

from .. import assert_same_replay


# Hash the ordered list of complete count histograms as compact sorted-key JSON.
# Entries bind, respectively, the deposited experiment and the independent replay.
STREAMS = {
    "reads": ("08c34c43e7df6694b07158fc30ce1562343d4033cca1cdbe323b5b3091f73733",
              "8207805a6c02544150247b72afe98701a31ac8eff230a8278953830d07490451"),
    "readings": ("5abc89c820d5396abeb5674cadba336ffe9b0385762f6966e5d73716deeca285",
                 "8094225676fbeb468cfd4fce5ad6249546d30fcfa04f1111cd98c329f7543377"),
}


def _extract(record, name, histograms):
    if isinstance(record, dict):
        other = {}
        for key, value in sorted(record.items()):
            if key == name:
                histograms.extend(value)
                other[key] = len(value)
            else:
                other[key] = _extract(value, name, histograms)
        return other
    if isinstance(record, list):
        return [_extract(value, name, histograms) for value in record]
    return record


def _check_counts(record, name, laws):
    from . import ancestor, measure

    def check(histogram, weight, words, a, b):
        if a == b == 2:
            conditional = np.ones(1)
        elif a == 2 or b == 2:
            conditional = words.sum(axis=0 if a == 2 else 1)
        else:
            conditional = words.ravel()
        counts = np.asarray(histogram)
        chances = weight * conditional
        if counts.shape != chances.shape or np.any(counts[chances == 0] != 0):
            raise ArithmeticError("Outcome outside the collar measurement law")
        # Retrospective Bernstein diagnostic for each bin in both full records.
        # Even 10^6 bins have combined tail probability at most 2*10^6*exp(-32).
        expected = measure.SHOTS * chances
        if np.any(np.abs(counts-expected) > np.sqrt(64*expected) + 64/3 + 1):
            raise ArithmeticError("Collar counts disagree with native placement and word probabilities")

    for panel in record["panels"]:
        if name == "reads":
            for setting, sizes in zip(panel["settings"], panel["geometry"]["annulus_counts"], strict=True):
                for mode, arm in setting["arms"].items():
                    placement = measure.slot_membership(sizes, mode == "colocated")
                    for hist, law in zip(arm[name], laws, strict=True):
                        for a, b in np.ndindex(3, 3):
                            check(hist[f"{a}{b}"], placement[a, b], law, a, b)
        else:
            for arm, setting in zip(panel["arms"], panel["geometry"]["settings"], strict=True):
                chances = ancestor.tagged_slots(setting["admissions"])
                for hist, law in zip(arm[name], laws, strict=True):
                    for tag, a, b in np.ndindex(2, 3, 3):
                        check(hist[f"{tag}{a}{b}"], chances[tag][a, b], law if tag == 0 else laws[1], a, b)


def verify_replay(actual, stored, *, name):
    """Certify raw outcomes exactly and compare the geometry, scores and controls."""
    archived, rerun = STREAMS[name]
    left, right = [], []
    assert_same_replay(_extract(actual, name, left), _extract(stored, name, right))
    def digest(histograms):
        return sha256(json.dumps(histograms, sort_keys=True, separators=(",", ":"),
                                 allow_nan=False).encode()).hexdigest()
    if digest(right) != archived or digest(left) not in (archived, rerun):
        raise ValueError("Unrecognized complete collar count stream")
    from .measure import source_reader
    laws, _, _ = source_reader()
    for record in (actual, stored):
        _check_counts(record, name, laws)
    for observed, retained in zip(left, right, strict=True):
        # At most three trial outcomes change bins among 2**48 trials per arm.
        # This bounds the change in the empirical probability of *any* event.
        differences = [n-m for key in retained for n, m in zip(observed[key], retained[key], strict=True)]
        if (observed.keys() != retained.keys() or sum(map(sum, observed.values())) != 2**48
                or sum(differences) or sum(map(abs, differences)) > 6):
            raise ArithmeticError("Collar trials lost or event-probability budget exceeded")
