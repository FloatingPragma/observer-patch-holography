from itertools import product
from math import atanh, pi

import numpy as np
import pytest

from model.executable_model import federation
from . import admission, boundary_sources, transport


@pytest.fixture(scope="module")
def panel():
    return admission.run()


def test_fresh_reference_comparator_and_visible_tangent(panel):
    source = panel["material"][3]
    q, d = (np.array(source["records"][-1][key]) for key in ("q", "d"))
    assert d[-1] == pytest.approx(-0.014370355, abs=2e-8)
    for f in (.001, .05, .5, 1.):
        shift = f*d[-1]
        probs = np.array([1-admission.compare_returns(q[-1]+shift, q[-1]),
                          admission.compare_returns(q[-1]+shift, q[-1])])
        tangent = np.array([shift, -shift])/2
        assert -tangent @ np.log(probs) == pytest.approx(admission.balanced_charge(q, d, f), abs=2e-16)
        assert 0 <= shift*atanh(shift) - shift**2 <= shift**4/(3*(1-shift**2)) + 1e-18
        expected = -f*np.dot(d, np.log(q+f*d))
        assert expected == pytest.approx(f*source["records"][-1]["linear_modular"]-admission.weaken(q, d, f))
    changed = panel["material"][-1]["records"][-1]
    assert abs(changed["d"][-1]) < 1e-13


def test_information_ceiling_and_admission_scaling(panel):
    assert all(arm["one_copy_classical_ceiling"] < arm["quantum_log"] for arm in panel["material"])
    arm = panel["material"][0]
    assert arm["quantum_log"]-arm["one_copy_classical_ceiling"] > 4.8e-5
    assert admission.sphere_count(2) == sum(i*i+j*j+k*k <= 4 for i in range(-2,3)
                                          for j in range(-2,3) for k in range(-2,3))
    flat = panel["flat_admission"]
    final = flat["rows"][-1]
    assert abs(final["balanced_growing_bank"]/flat["limiting_balanced_growing_bank"] - 1) < .001
    assert final["balanced_one_mark"]/final["balanced_growing_bank"] < .0001
    assert flat["limiting_balanced_growing_bank"] == pytest.approx(
        (4*pi/3 * panel["material"][3]["records"][-1]["d"][-1])**2
    )


def test_point_injection_and_absorption():
    result = transport.run()
    conservative, lossy = result["complete_box"]
    for box in [conservative, lossy]:
        assert box["maximum_cell_residual"] < 1e-13
        assert box["total_sink"] == pytest.approx(1)
    assert lossy["rows"][1]["local_density_ratio"] > 1000
    near_asymptote = result["larger_boxes"]
    errors = [row["density_ratio"]-4 for row in near_asymptote]
    assert all(a>b>0 for a, b in zip(errors, errors[1:]))
    assert errors[-1] < .13


def test_ancestral_cut_and_oriented_response():
    rows = boundary_sources.run()["rows"]
    assert [row["boundary"] for row in rows] == [3, 3, 6, 12, 12, 48]
    for row in rows:
        assert row["observed_per_u_squared"] == pytest.approx(row["coefficient"], rel=.004)
        assert row["one_bit_observed_per_u_squared"] == pytest.approx(row["one_bit_coefficient"], rel=.004)
        assert 0 < row["one_bit_coefficient"] < row["coefficient"]
    assert boundary_sources.absolute_flux(12) == 693/256
    assert boundary_sources.absolute_flux(0) == 0
    assert boundary_sources.absolute_flux(1024)**2/1024 == pytest.approx(2/pi, rel=.001)
    signs = np.array(list(product((-1, 1), repeat=6)))
    full_flux = signs.sum(axis=1)
    for available in range(7):
        observed_flux = signs[:, :available].sum(axis=1)
        assert np.mean(np.sign(observed_flux)*full_flux) == pytest.approx(
            boundary_sources.absolute_flux(available)
        )


def test_explicit_sign_records_for_unequal_neighboring_sources():
    board = federation(0)
    sites = [port//12 for port in board.external[0]]
    cut = boundary_sources.boundary_coefficients(board, sites, [1., 1.])
    field = boundary_sources.boundary_coefficients(board, sites, [.3, .8])
    active = (field != 0) | (cut != 0)
    words = np.array(list(product((-1, 1), repeat=sum(active))))
    flux = words @ cut[active]
    signed = words @ field[active]
    exposure = .015
    y = exposure*signed
    reference = .5
    p = np.column_stack((1-reference*y, 1+reference*y))/(2*len(words))
    weakening = np.column_stack((y, -y))/(2*len(words))
    assert -np.sum(weakening * np.log(p)) == pytest.approx(np.mean(y*np.arctanh(reference*y)))
    assert np.mean(signed**2) == pytest.approx(field @ field)
    np.testing.assert_allclose(p.sum(axis=0), [1/2, 1/2])
    # A fair sign on zero flux has conditional mean zero.
    aligned = np.mean(np.sign(flux)*signed)
    b = int(cut @ cut)
    assert aligned == pytest.approx(boundary_sources.absolute_flux(b)/b*(cut @ field))
    folded_yes = (1+np.sign(flux)*reference*y)/2
    assert folded_yes.mean() == pytest.approx((1+reference*exposure*aligned)/2)
    pair_agreement = np.mean(((1+reference*y)/2)**2 + ((1-reference*y)/2)**2)
    assert pair_agreement == pytest.approx((1+reference**2*np.mean(y**2))/2)


def test_joint_limit_admission_saturation():
    contrast = -.014370365536930434
    reference = .5
    n = 256**2
    boundary = 3*256
    for c in (.1, 1., 10.):
        u = c/n
        gain = contrast*admission.any_hit(u, n)/(3*n)
        aligned = gain*boundary_sources.absolute_flux(boundary)
        sparse = reference*(contrast*u/3)**2
        flagged_ratio = boundary_sources.conditional_charge(boundary, gain, reference)/(sparse*boundary)
        oriented_ratio = aligned*atanh(reference*aligned)/(sparse*boundary_sources.absolute_flux(boundary)**2)
        target = (np.expm1(-c)/c)**2
        assert flagged_ratio == pytest.approx(target, rel=1e-5)
        assert oriented_ratio == pytest.approx(target, rel=1e-5)
