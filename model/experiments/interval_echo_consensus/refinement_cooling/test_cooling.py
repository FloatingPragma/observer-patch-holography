from fractions import Fraction
import json
from math import log
from pathlib import Path

import pytest

from . import measure


def test_complete_native_counts_and_flight_spectrum():
    report = measure.run()
    held = json.loads(Path(measure.__file__).with_name("results.json").read_text())
    assert report["record_panels"] == held["record_panels"]
    for fresh, saved in zip(report["doubled_flight"], held["doubled_flight"], strict=True):
        assert fresh["excited_mean_energy"] == saved["excited_mean_energy"]
        assert fresh["ports"] == saved["ports"]
        for a, b in zip(fresh.get("spectrum_evaluation", []), saved.get("spectrum_evaluation", []), strict=True):
            assert a == pytest.approx(b, abs=1e-11)


def test_refinement_cost_and_uniform_mass_and_entropy():
    for depth in (2, 3, 7, 16, 64):
        words = measure.words(depth)
        assert len(words) == len(set(words)) == depth+1
        assert sum(Fraction(1, 2**len(h)) for h in words) == 1
        for h in words:
            for coarse in range(2, depth+1):
                assert h[:coarse] in measure.words(coarse)
                assert h[:coarse][:2] == h[:2]
                assert len(h[:coarse]) <= len(h)
        for beta in (1, 1.5, 2, 4, 10):
            r = 2**(-beta)
            weights = [r**(len(h)-1) for h in words]
            z = sum(weights)
            # Direct entropy expression, stable under underflow for long words.
            s = log(z) + beta*log(2)*sum((len(h)-1)*w/z for h, w in zip(words, weights))
            assert (z-1)/z <= 2*r/(1+2*r) + 1e-14
            assert s <= log(1+2*r) + 2*beta*log(2)*r/(1-r)**2 + 1e-14
        for beta in (1, 2, 3):
            acceptance = sum(Fraction(1, 2**(beta*len(h))) for h in words)
            assert acceptance >= Fraction(1, 2**beta)
            mean_cost = sum(Fraction(len(h), 2**len(h)) for h in words)
            assert mean_cost == 2*(1-Fraction(1, 2**depth))
            assert beta*mean_cost/acceptance <= 2*beta*2**beta
