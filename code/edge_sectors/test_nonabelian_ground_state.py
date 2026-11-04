"""Independent original-matrix controls for the finite S3 sector readout.

The oracle diagonalizes the specified character Hamiltonian at 350 digits.
It does not call the producer's secular root, normalization or log formula.
"""

from decimal import Decimal
from fractions import Fraction
from numbers import Integral
from pathlib import Path
import re
import subprocess
import sys

import mpmath
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from edge_sectors.nonabelian_ground_state import (
    ground_root_bracket, s3_diagnostics, s3_edge_distribution,
)


def _exact_scalar(value):
    if isinstance(value, Integral):
        return Fraction(int(value))
    if isinstance(value, Fraction):
        return value
    numerator, denominator = value.as_integer_ratio()
    return Fraction(int(numerator), int(denominator))


def _matrix_oracle(h):
    ctx = mpmath.mp.clone()
    ctx.dps = 350
    original = _exact_scalar(h)
    h = ctx.mpf(original.numerator)/original.denominator
    matrix = ctx.matrix([[0, 0, -1], [0, 24*h, -1], [-1, -1, 12*h-1]])
    eigenvalues, eigenvectors = ctx.eigsy(matrix)
    p = {name: eigenvectors[row, 0]**2
         for row, name in enumerate(("triv", "sign", "std"))}
    time = ctx.log(2*p["triv"]/p["std"])/3
    predicted = p["triv"]*ctx.exp(-6*time)
    discrepancy = ctx.log(p["sign"]/predicted)
    return ctx, -eigenvalues[0], {
        "probabilities": p,
        "fit_time": time,
        "diffusion_fit": time > 0,
        "predicted_sign": predicted,
        "log_discrepancy": discrepancy,
        "normalized_log_residual": discrepancy/abs(ctx.log(predicted/p["triv"])),
        "log_ratio_excess": ctx.log(p["sign"]/p["triv"])
                            /ctx.log(p["std"]/(2*p["triv"]))-2,
    }


def _assert_readout_matches_oracle(h):
    _, _, expected = _matrix_oracle(h)
    measured = s3_diagnostics(h)
    probabilities = s3_edge_distribution(h)
    assert measured["diffusion_fit"] == expected["diffusion_fit"]
    for name, value in expected["probabilities"].items():
        assert probabilities[name] == pytest.approx(float(value), rel=5e-13, abs=0)
        assert measured["probabilities"][name] == pytest.approx(float(value), rel=5e-13, abs=0)
    for name in ("fit_time", "predicted_sign", "log_discrepancy",
                 "normalized_log_residual", "log_ratio_excess"):
        assert measured[name] == pytest.approx(float(expected[name]), rel=5e-13, abs=0)
    # Every finite coupling misses the d_R ansatz on the held-out sector.
    # Rounded equality of probabilities must not turn this into exact fit.
    assert measured["log_discrepancy"] < 0
    assert measured["normalized_log_residual"] < 0
    assert measured["log_ratio_excess"] != 0


@pytest.mark.parametrize("h", [0., .05, .5, 100., 1e15, 1e40, 1e70, 1e75, 1e76])
def test_original_character_hamiltonian_agrees_in_each_component(h):
    _assert_readout_matches_oracle(h)


_ZERO_TIME = float("0.054909885216490926863573234054007262015547224101393")


@pytest.mark.parametrize("h", [np.nextafter(_ZERO_TIME, -np.inf),
                               _ZERO_TIME, np.nextafter(_ZERO_TIME, np.inf)])
def test_neighbors_of_zero_diffusion_keep_the_original_input_sign(h):
    _assert_readout_matches_oracle(h)
    # The exactly supplied middle float lies above the irrational boundary.
    assert s3_diagnostics(h)["diffusion_fit"] == (h >= _ZERO_TIME)


def test_unresolved_zero_time_normalization_preserves_probabilities():
    ctx = mpmath.mp.clone()
    ctx.dps = 150
    boundary = (1+ctx.sqrt(3)-ctx.sqrt(2))/24
    original = Fraction(int(ctx.floor(boundary*10**100)), 10**100)
    _, _, expected = _matrix_oracle(original)
    measured = s3_edge_distribution(original)
    for name in measured:
        assert measured[name] == pytest.approx(float(expected["probabilities"][name]), rel=5e-13)
    with pytest.raises(ValueError, match="unresolved.*root precision"):
        s3_diagnostics(original)


@pytest.mark.parametrize("h", [0., .05, .5, 100., 1e15, 1e40, 1e76])
@pytest.mark.parametrize("bits", [64, 256])
def test_exact_polynomial_bracket_contains_independent_matrix_eigenvalue(h, bits):
    lower, upper = ground_root_bracket(h, bits=bits)
    original = _exact_scalar(h)
    # Independently expanded characteristic polynomial of H at E=-x.
    def polynomial(x):
        return (x**3+(36*original-1)*x*x
                +(288*original*original-24*original-2)*x-24*original)
    assert isinstance(lower, Fraction) and isinstance(upper, Fraction)
    assert 0 < lower <= upper
    assert polynomial(lower) <= 0 <= polynomial(upper)
    assert upper-lower <= 6*lower/Fraction(2**bits)
    ctx, root, _ = _matrix_oracle(h)
    lo = ctx.mpf(lower.numerator)/lower.denominator
    hi = ctx.mpf(upper.numerator)/upper.denominator
    if lower == upper:
        assert abs(root-lo) <= ctx.mpf("1e-330")
    else:
        assert lo < root < hi


@pytest.mark.parametrize("h", [2**53+1, np.int64(2**53+1),
                               Fraction(2**53+1), Decimal(2**53+1)])
def test_original_scalar_is_not_rounded_before_the_root_certificate(h):
    lower, upper = ground_root_bracket(h)
    rounded_lower, _ = ground_root_bracket(float(h))
    assert upper < rounded_lower
    _assert_readout_matches_oracle(h)


@pytest.mark.parametrize("h", [1e80, Fraction(10**100), Decimal("1e100")])
@pytest.mark.parametrize("evaluate", [s3_edge_distribution, s3_diagnostics])
def test_unrepresentable_positive_sector_is_refused(h, evaluate):
    _, _, expected = _matrix_oracle(h)
    assert 0 < expected["probabilities"]["sign"] < float.fromhex("0x0.0000000000001p-1022")
    with pytest.raises(ValueError, match="probability.*reporting range"):
        evaluate(h)


@pytest.mark.parametrize("bad", [True, np.bool_(False), 1+0j, np.complex128(.5),
                                 np.nan, np.inf, -np.inf, -1., Fraction(-1, 3),
                                 Decimal("NaN"), np.ma.masked,
                                 np.ma.array(.5, mask=True), np.ma.array(.5, mask=False)])
@pytest.mark.parametrize("evaluate", [ground_root_bracket, s3_edge_distribution, s3_diagnostics])
def test_invalid_original_coupling_is_refused(bad, evaluate):
    with pytest.raises(ValueError, match="coupling"):
        evaluate(bad)


def test_paper_s3_identity_replays_from_original_character_hamiltonian():
    paper = (Path(__file__).resolve().parents[2]/"paper"/
             "deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex")
    text = " ".join(paper.read_text(encoding="utf-8").split())
    assert r"x+12h-1-1/x-1/(24h+x)=0" in text
    assert r"=\left(\frac{z-1}{z}\right)^2<1" in text
    for h in (.5, 1., 2., 5., 12., 100.):
        ctx, x, reference = _matrix_oracle(h)
        p = reference["probabilities"]
        assert abs(x+12*h-1-1/x-1/(24*h+x)) < ctx.mpf(10)**-300
        z = x+12*h
        ratio = p["sign"]/(p["std"]**2/(4*p["triv"]))
        assert abs(ratio-((z-1)/z)**2) < ctx.mpf(10)**-300
        assert ratio < 1


def test_public_s3_cli_keeps_signed_fit_and_tiny_nonzero_discrepancies():
    script = Path(__file__).with_name("heat_kernel_holdout_validation.py")
    couplings = [.05, .5, 1e15, 1e70]
    run = subprocess.run([sys.executable, "-W", "error", str(script), "--groups", "S3",
                          "--h", *map(str, couplings)], capture_output=True, text=True, check=True)
    rows = [line.split() for line in run.stdout.splitlines() if re.match(r"^\d", line)]
    details = re.findall(r"log\(measured/predicted\)=([+-][\d.e+-]+); "
                         r"positive diffusion fit: (True|False)", run.stdout)
    assert len(rows) == len(details) == len(couplings)
    for h, row, (discrepancy, diffusion_fit) in zip(couplings, rows, details):
        _, _, expected = _matrix_oracle(h)
        for cell, name in zip(row[1:4], ("triv", "sign", "std")):
            assert float(cell) == pytest.approx(float(expected["probabilities"][name]), rel=5e-5, abs=0)
        assert float(row[4]) == pytest.approx(float(expected["fit_time"]), abs=5e-5, rel=0)
        for cell, name, tolerance in ((row[5], "predicted_sign", 5e-4),
                                      (row[6], "normalized_log_residual", 5e-4),
                                      (row[7], "log_ratio_excess", 5e-7),
                                      (discrepancy, "log_discrepancy", 5e-7)):
            assert float(cell) == pytest.approx(float(expected[name]), rel=tolerance, abs=0)
        assert float(discrepancy) < 0
        assert (diffusion_fit == "True") == expected["diffusion_fit"]
