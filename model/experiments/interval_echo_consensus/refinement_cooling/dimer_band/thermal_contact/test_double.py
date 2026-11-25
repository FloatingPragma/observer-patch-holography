import json
from math import cos, sqrt, log
import numpy as np
from .measure import HERE
from .test_contact import read_summary


def test_independent_double_energy_frequency_and_energy_balanced_contact():
    j = json.loads((HERE / "double.json").read_text())
    signed = np.array([1, -1, -1, 1])

    def eps(n, width=1):
        return width * sqrt(log(2 * 2048 / 0.01) / (2 * n))

    for row in j["clocks"]:
        c = np.array(row["counts"])
        assert (c.sum(axis=1) == j["clock_shots"]).all()
        observed = (c[0] - c[1]) @ signed / j["clock_shots"]
        assert abs(observed - cos(2 * row["time"])) < 2 * eps(j["clock_shots"], 2)
    _, temps, _ = read_summary()
    p, n = temps["B", 1]

    def paired(x):
        return x * x / (x * x + (1 - x) ** 2)

    lo, hi = paired(p - eps(n)), paired(p + eps(n))
    for row in j["contacts"]:
        events = row["events"]
        good = sum(e[-1] for e in events)
        assert good + sum(row["refusal_counts"]) == j["shots"]
        high = sum(c * k for a, b, c, o, d, e, f, k in events) / good
        flow = sum((d + e - a - b) * k for a, b, c, o, d, e, f, k in events) / good
        assert all(a + b + 2 * c == d + e + 2 * f for a, b, c, o, d, e, f, k in events)
        if not row["copy"]:
            assert lo - eps(good) < high < hi + eps(good)
            assert abs(flow) < eps(good, 4)
        else:
            assert high > hi + eps(good)
            if row["parked"]:
                assert flow == 0
            else:
                assert abs(flow - 2 / 27) < eps(good, 4) and flow > eps(good, 4)
