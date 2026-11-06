"""Check preparation/readout composition beyond the stationary native input."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import certificate as c


def half_trace(matrix):
    return np.linalg.svd(matrix, compute_uv=False).sum()/2


def reference_state(rho, dim, ref):
    return np.trace(rho.reshape(dim, ref, dim, ref), axis1=0, axis2=2)


def cq_distance(left, right):
    assert tuple(left) == tuple(right)
    return sum(half_trace(left[record]-right[record]) for record in left)


def test_exact_prefix_uses_full_carrier_and_all_law_envelope():
    report = c.preparation_certificate()
    assert report["full_carrier"] and report["law"] == "all"
    assert report["sweeps"] == 19 and report["native_attempts"] == 570
    assert report["collision_absorption_sweeps"] == 1
    assert report["trace_distance_upper"] < F(1, 2**24)
    assert c.sweeps.certificate(18, "all", full=True)["trace_distance_upper"] > F(1, 2**24)


@pytest.mark.parametrize("law", ("C", "E"))
def test_arbitrary_input_executor_matches_original_stationary_all_records(law):
    actual = c.native_outputs(np.eye(66)/66, law, 4)
    parent = dict(c.readout.execute_native(law, 4).outputs)
    assert tuple(actual) == c.readout.records(4)
    assert cq_distance(actual, parent) < 2e-15
    assert abs(sum(np.trace(x) for x in actual.values())-1) < 2e-15


@pytest.mark.parametrize("law", ("C", "E"))
def test_entangled_collision_input_full_suffix_and_reference(law):
    # Every occupied amplitude is in a collision state. Distinct reference
    # amplitudes and their phases test the original marginal and coherences.
    collisions = [i for i, x in enumerate(c.sweeps.FULL_STATES) if max(x) == 2]
    psi = np.zeros(78*2, complex)
    psi[2*collisions[0]] = 1/2
    psi[2*collisions[3]+1] = 1j/2
    psi[2*collisions[8]] = np.sqrt(.5)
    rho = np.outer(psi, psi.conj())
    rho_r = reference_state(rho, 78, 2)
    actual = c.execute(rho, law, 4, reference_dimension=2)
    target = c.native_outputs(np.kron(np.eye(66)/66, rho_r), law, 4, reference_dimension=2)
    assert cq_distance(actual, target) <= float(c.preparation_certificate()["trace_distance_upper"])
    aggregate = sum(actual.values())
    np.testing.assert_allclose(reference_state(aggregate, 66, 2), rho_r, atol=5e-15)
    assert abs(np.trace(aggregate)-1) < 3e-15
    assert len(actual) == 8*4+3
    assert (0, 0, 0, -1) in actual and (1, 4, 0, -1) in actual and (2, 4, 0, -1) in actual


def test_without_preparation_the_stationary_herald_rate_is_not_general():
    # A definite configuration with equal loads at the herald seam never
    # heralds; substituting the stationary5/33 rate would be an error.
    seam = c.readout.interface.quantum.mixing.seams()[0]
    index = next(i for i, x in enumerate(c.readout.clock.quantum.STATES) if x[seam[0]] == x[seam[1]])
    rho = np.zeros((66, 66), complex)
    rho[index, index] = 1
    output = c.native_outputs(rho, "C", 4)
    assert abs(np.trace(output[0, 0, 0, -1])-1) < 5e-16
    assert abs(float(1-c.readout.YIELD)-np.trace(output[0, 0, 0, -1])) > .1


def test_vacuous_tolerance_preserves_required_collision_absorption():
    report = c.preparation_certificate(F(1))
    assert report["sweeps"] == report["collision_absorption_sweeps"] == 1
    rho = np.zeros((78, 78), complex)
    collision = next(i for i, x in enumerate(c.sweeps.FULL_STATES) if max(x) == 2)
    rho[collision, collision] = 1
    outputs = c.execute(rho, "C", 1, tolerance=F(1))
    assert abs(sum(np.trace(x) for x in outputs.values())-1) < 2e-15


@pytest.mark.parametrize("mass", (0., 2**-30, .375))
def test_no_conditional_renormalization(mass):
    rho = np.eye(78, dtype=complex)/78
    output = c.execute(mass*rho, "E", 1, tolerance=F(1, 1000))
    assert abs(sum(np.trace(x) for x in output.values())-mass) < 2e-15*max(mass, 1e-30)


@pytest.mark.parametrize("bad", (True, 0, -1, .001, F(2)))
def test_bad_precision_rejected(bad):
    with pytest.raises(ValueError):
        c.preparation_certificate(bad)


def test_parent_binding_and_exact_completed_event_intervals():
    receipt = c.run()
    assert receipt["parent"]["sha256"] == c.PARENT_SHA256
    assert len(receipt["histories"]) == 6
    epsilon = receipt["preparation"]["trace_distance_upper"]
    delta = receipt["source_mixer_two_call_half_norm_error_upper"]
    for row in receipt["histories"]:
        lo, hi = map(F, row["ideal_phase_match_probability_interval"])
        exact = c.experiment.exact_native_statistics(row["law"], row["cap"])["phase_match_and_success"]
        assert (lo+hi)/2 == exact and hi-lo == 2*epsilon
        slo, shi = map(F, row["source_mixer_phase_match_weight_interval"])
        assert shi-slo == 2*epsilon+4*delta
        assert row["maximum_native_attempts_including_prefix_herald_and_pointer"] == 572+row["cap"]
        assert row["complete_joint_suffix_records"] == 128*row["cap"]+3
    assert F(receipt["cap64_source_mixer_C_minus_E_phase_match_weight_interval"][0]) > 0


def test_changed_parent_receipt_is_rejected(monkeypatch):
    monkeypatch.setattr(c, "PARENT_SHA256", "0"*64)
    with pytest.raises(ValueError, match="parent receipt"):
        c.parent_receipt()
