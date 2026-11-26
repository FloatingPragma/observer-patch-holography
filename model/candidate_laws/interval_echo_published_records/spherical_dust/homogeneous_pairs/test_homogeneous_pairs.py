from fractions import Fraction as F
import math

from .experiment import (CENTERS, DIRECTIONS, FRAMES, REST, UNENROLLED, calibration,
                         complete_pairs, gate, lattice_points, read_bank, relative_factor, verify)


def test_ten_named_interpolations_and_held_orientation():
    determinant, weights = calibration()
    assert F(6994924, 10**8) < determinant[0] < determinant[1] < F(6994925, 10**8)
    for frame in FRAMES:
        lo, hi = relative_factor(frame, weights)
        assert lo <= frame[0]**2 <= hi and hi-lo < F(1, 10**18)
    lo, hi = relative_factor(UNENROLLED, weights)
    assert F(-462, 1000) < lo < hi < F(-461, 1000)


def test_native_spectra_allow_both_signed_exposures():
    for e in (-0.001, 1.e-12, 0.003):
        pair = gate(e)
        assert math.isclose(pair.cut_entropy_change, e, rel_tol=1e-12)
        assert 0 < pair.cut_relative_entropy <= 30*e*e
        assert 0 < pair.whole_modular_change <= 40*e*e


def test_overlapping_displaced_regions_independent_of_preparation():
    den = 256
    radius = F(1, 64)
    reference = {tuple(p) for p in lattice_points(den, CENTERS[0], radius, REST)}
    moved = lattice_points(den, CENTERS[1], radius, REST)
    points = {tuple(p) for p in moved}
    assert reference & points and reference - points and points - reference
    for group in DIRECTIONS:
        for b in group:
            steps = math.floor(den/(8*math.sqrt(sum(x*x for x in b))))
            brute = [sum(tuple(x+j*d for x, d in zip(p, b)) in points for p in points)
                     for j in range(1, steps+1)]
            assert complete_pairs(moved, b, steps) == (sum(brute), brute[0])
    blank = read_bank(den, CENTERS[1], radius, REST, 0, [0, 0, 0])
    assert abs(blank["entropy_over_radius4"]) < 1e-14
    assert blank["relative_entropy_over_radius4"] < 1e-30


def test_archive_and_control():
    results = verify()
    assert math.isclose(results["held_rest_density"], 19/(44*math.pi))
    separate = results["separately_enrolled_blank_control"]
    assert separate["flat_projected_radius"] < 0.0156252 < separate["held_projected_radius"]
    assert separate["blank_entropy_increase_lower"] > 100
    assert len(results["surveys"]) == 18
    last = next(s for s in results["surveys"] if s["mesh_denominator"] == 1024
                and s["radius"] == "1/64" and s["center"][0] == "81/1021")
    assert last["frames"][-1]["entropy_over_radius4"] < 0
    assert max(row["relative_entropy_over_radius4"] for row in last["frames"]) < 2.64e-7
