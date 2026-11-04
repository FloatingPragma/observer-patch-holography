"""Independent spectral and charge controls for the matching lemma."""
import pytest
import sympy as s

from .matching import check_intertwiner, example_certificate


def test_rotated_rectangular_exact_embedding_preserves_joint_spectrum():
    j = s.Matrix([[s.Rational(3, 5), 0], [0, 1], [s.Rational(4, 5), 0]])
    complement = s.Matrix([-s.Rational(4, 5), 0, s.Rational(3, 5)])
    he, qe = s.diag(2, 3), s.diag(-1, 1)
    hs = j*he*j.H + 7*complement*complement.H
    qs = j*qe*j.H
    report = check_intertwiner(hs, he, j, charges={"electric": (qs, qe)})
    assert report["generator_intertwining_exact"]
    assert report["spectral_distance_upper"] == "0"
    assert set(hs.eigenvals()) == {2, 3, 7}
    assert not report["source_spectrum_exhausted"]


def test_compression_alone_misses_source_leakage():
    report = example_certificate()
    assert report["leakage_control_compression_equals_field_generator"]
    leaking = report["leakage_control"]
    assert not leaking["generator_intertwining_exact"]
    assert leaking["generator_residual_frobenius_squared"] == "1/100"
    assert leaking["spectral_distance_upper"] == "1/10"
    # Directly diagonalize the nontrivial 2x2 block to check the bound.
    closest = 3-s.sqrt(s.Rational(401, 100))
    assert abs(closest-1) < s.Rational(1, 10)


def test_matching_rejects_charge_relabelling_and_nonconservation():
    h = s.diag(1, 2)
    with pytest.raises(ValueError, match="does not intertwine"):
        check_intertwiner(h, h, s.eye(2), charges={"q": (s.diag(1, -1), s.diag(-1, 1))})
    with pytest.raises(ValueError, match="not a conserved"):
        check_intertwiner(h, h, s.eye(2), charges={"q": (s.Matrix([[0, 1], [1, 0]]), s.eye(2))})


@pytest.mark.parametrize("embedding", [s.diag(1, 2), s.ones(3, 1)])
def test_nonisometric_or_wrong_size_maps_are_rejected(embedding):
    with pytest.raises(ValueError, match="isometry|dimensions"):
        check_intertwiner(s.eye(2), s.eye(2), embedding)


def test_generator_self_adjointness_and_exactness_are_required():
    with pytest.raises(ValueError, match="Hermitian"):
        check_intertwiner(s.Matrix([[1, 1], [0, 2]]), s.eye(2), s.eye(2))
    with pytest.raises(ValueError, match="finite exact"):
        check_intertwiner(s.diag(1.0, 2), s.eye(2), s.eye(2))
    with pytest.raises(ValueError, match="finite exact"):
        check_intertwiner(s.diag(s.oo, 2), s.eye(2), s.eye(2))
