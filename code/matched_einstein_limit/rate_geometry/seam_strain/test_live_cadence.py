from fractions import Fraction
import json

import numpy as np
from scipy.integrate import dblquad
from scipy.linalg import expm
from scipy.sparse import csr_matrix

from . import live_cadence as cadence


def direct_branches(hop, launch, start, relay, p, steps):
    live = launch.copy()
    returning = []
    mass = np.zeros(steps+1)
    for t in range(1, steps+1):
        live = hop @ live
        relay_branch = np.zeros_like(live)
        relay_branch[relay] = np.sqrt(p)*live[relay]
        live[relay] *= np.sqrt(1-p)
        for i, branch in enumerate(returning):
            after = hop @ branch
            mass[t] += p*np.linalg.norm(after[start])**2
            after[start] *= np.sqrt(1-p)
            returning[i] = after
        returning.append(relay_branch)
        survivor = np.linalg.norm(live)**2 + sum(np.linalg.norm(v)**2 for v in returning)
        assert abs(1-survivor-mass.sum()) < 2e-14
    return mass


def test_live_head_renewal_against_all_explicit_detector_branches():
    rng = np.random.default_rng(26507)
    values = rng.integers(-15, 16, size=(6, 6))
    h40 = values+values.T
    rows = np.array([[0, 1], [2, 3]])
    step, steps = .25, 13
    hop = expm(-1j * step * h40/40)
    kernels = cadence.kernel_grid(csr_matrix(h40), rows, step*np.arange(steps+1))
    initial = np.array([1, 1j])/np.sqrt(2)
    for p in (1., -np.expm1(-step)):
        explicit = np.zeros(steps+1)
        for start, relay in ((rows[0], rows[1]), (rows[1], rows[0])):
            launch = np.zeros(6, complex)
            launch[start] = initial
            explicit += direct_branches(hop, launch, start, relay, p, steps)/2
        means = cadence.effects(kernels, 0, 1, p)
        np.testing.assert_allclose(np.einsum("i,tij,j->t", initial.conj(), means, initial).real,
                                   explicit, atol=2e-14, rtol=0)


def test_rate_one_rule_approaches_exact_damped_dimer():
    # H exchanges two sites at unit amplitude. The continuum miss
    # generator has eigenvalues -1/4 ± i sqrt(15)/4.
    omega = np.sqrt(15)/4
    def density(t):
        return np.exp(-t/2)*(np.sin(omega*t)/omega)**2
    exact = dblquad(lambda t,s: density(t)*density(s), 0, 3., 0, lambda s:3-s)[0]
    errors = []
    for width in (8, 16, 32, 64):
        t = np.arange(3*width+1)/width
        kernel = np.zeros((len(t), 2, 2, 1, 1), dtype=complex)
        kernel[:, 0, 0, 0, 0] = kernel[:, 1, 1, 0, 0] = np.cos(t)
        kernel[:, 0, 1, 0, 0] = kernel[:, 1, 0, 0, 0] = -1j*np.sin(t)
        discrete = np.sum(cadence.effects(kernel, 0, 1, -np.expm1(-1/width))).real
        errors.append(abs(discrete-exact))
    assert all(b < .55*a for a,b in zip(errors[:-1], errors[1:], strict=True))
    assert errors[-1] < 0.006


def test_native_reference_ranges_and_checked_cadence_receipt():
    fresh = cadence.run()
    logged = json.loads((cadence.HERE/"live_cadence.json").read_text())
    survey = json.loads((cadence.HERE/"live_survey.json").read_text())
    assert fresh["continuous_near_return_lower"] == logged["continuous_near_return_lower"]
    assert all(Fraction(x) > Fraction(18, 10**18) for x in fresh["continuous_near_return_lower"].values())
    assert fresh["leakage_squared_upper"] == logged["leakage_squared_upper"]
    for a,b in zip(fresh["panels"], logged["panels"], strict=True):
        for k in ("detector", "step", "ports", "sharp_uniform_complete_upper"):
            assert a[k] == b[k]
        for charge in ("-1","1"):
            for arm, result in a["by_charge"][charge].items():
                for key, value in result.items():
                    np.testing.assert_allclose(value, b["by_charge"][charge][arm][key], rtol=1e-9, atol=2e-13)
                if arm != "hidden":
                    assert abs(1-result["complete"]+result["at_deadline"]-result["saturated"]) < 1e-13
                    if a["detector"] == "sharp":
                        assert result["complete"] <= a["sharp_uniform_complete_upper"] + 1e-13
        if a["step"] == .5 and a["detector"] == "sharp":
            source, = (s for s in survey["panels"] if s["deadline"] == 12 and s["coupled"] and s["read"] == "pooled")
            i, j = (survey["stations"].index(port) for port in a["ports"])
            assert abs(a["pooled"]["reference"]["capped_range"]-source["lengths"][i][j]) < 2e-13
    near = [row for row in fresh["panels"] if row["ports"] == [0,1] and row["detector"] == "rate_one"]
    assert abs(near[-1]["pooled"]["reference"]["capped_range"]-2.96828756) < 1e-8
    remote, = (r for r in fresh["panels"] if r["ports"] == [5,40] and r["step"] == 1/64 and r["detector"] == "sharp")
    assert 2.7e-10 < remote["pooled"]["reference"]["complete"] < 2.9e-10
