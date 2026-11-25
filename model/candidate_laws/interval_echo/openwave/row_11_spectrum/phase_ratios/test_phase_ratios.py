import json
from cmath import exp, pi
from math import cos, sqrt
from pathlib import Path

from model.experiments.interval_echo_lepton_register import measure as register

from . import certificate as run


def test_exact_minors_and_independent_reader():
    assert run.certificate() == json.loads(
        Path(run.__file__).with_name("results.json").read_text()
    )
    root = exp(2j * pi / 24)
    ports = register.addresses((0, 1))
    prepared = register.prepare(1)
    for theta in (2 / 9, 0, pi / 3, -0.37):
        z = exp(1j * theta)
        weights = []
        for k in range(3):
            p = run.channel(k)
            observed = sum(
                complex(v.subs(run.T, root)) * z**n for n, v in enumerate(p)
            ) / (36 * z**2)
            weight = (sqrt(2) + 2 * cos(theta + 2 * pi * k / 3)) ** 2 / 36
            assert abs(weight - observed) < 3e-14
            weights.append(weight)
        material = abs(register.reader(ports, theta, 1) @ prepared) ** 2
        terminal = sorted(material[[ports[0], ports[2], ports[3]]])
        assert (
            max(abs(x - y) for x, y in zip(sorted(weights), terminal, strict=True))
            < 5e-13
        )
        if theta == 2 / 9:
            light, middle, heavy = sorted(weights)
            assert abs(middle / light - 206.770316) < 5e-7
            assert abs(heavy / light - 3477.472837) < 5e-7
        if theta == 0:
            assert abs(weights[1] - weights[2]) < 3e-15
