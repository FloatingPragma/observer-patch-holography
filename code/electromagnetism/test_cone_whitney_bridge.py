"""Semantic mutation controls for the finite volume interpolation."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_cone_whitney_bridge as verifier


def test_independent_volume_replay():
    result = verifier.verify(verifier.load())
    assert result["field_variations_per_history"] == 68
    assert result["ampere_hat_max_abs"] > 1


def test_fresh_producer_passes_independent_quadrature():
    import cone_whitney_bridge as producer
    fresh = json.loads(producer.canonical(producer.build()))
    assert verifier.verify(fresh)["field_variations_per_history"] == 68
    stored = verifier.load()
    assert fresh["cochain_maps"] == stored["cochain_maps"]
    for current, previous in zip(fresh["executions"], stored["executions"], strict=True):
        assert current["cone_cochains"] == previous["cone_cochains"]


def replace(path, value):
    def mutation(packet):
        target = packet
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
    return mutation


@pytest.mark.parametrize("mutation", [
    replace(["scope"], "PHYSICAL_CONTINUUM_CONFIRMED"),
    lambda p: p["pins"].pop("Lean/Screen/ConeCochainBridge.lean"),
    replace(["pins", "Lean/Screen/WhitneyTimeBridge.lean"], "0"*64),
    replace(["mesh", "apex_vertex"], False),
    replace(["mesh", "port_vertex_map", 0], True),
    replace(["mesh", "vertices", 0], [0, 0, 0.1]),
    lambda p: p["mesh"]["tetrahedra"].pop(),
    lambda p: p["mesh"]["faces"][0].reverse(),
    replace(["cochain_maps", "P0", 0, 0], "0"),
    replace(["cochain_maps", "P1", 0, 0], "0"),
    replace(["cochain_maps", "P2", 20, 0], "0"),
    replace(["cochain_contract", "magnetic_domain"], "all boundary fluxes"),
    replace(["numeric_policy", "comparison_atol"], 10.0),
    replace(["whitney", "M1"], np.eye(42).tolist()),
    replace(["whitney", "tetrahedra", 0, "one_coefficients", 0, 0, 0], 100.0),
    replace(["time_source_contract", "current"], "ordinary constant current on each slab"),
    lambda p: p["executions"].pop(),
    replace(["executions", 0, "h"], "180"),
    replace(["executions", 0, "decode_event_ids", 0], 0),
    replace(["executions", 0, "source_history", "A", 0, 0], "999"),
    replace(["executions", 0, "cone_cochains", "E", 0, 0], "0"),
    replace(["executions", 0, "action", "volume_prism"], 0.0),
    replace(["executions", 0, "action", "temporal_defect"], 0.0),
    replace(["executions", 0, "residuals", "ampere_hat"], [0.0]*42),
    lambda p: p["executions"][0]["residuals"].__setitem__(
        "ampere_hat", p["executions"][0]["residuals"]["ampere_hat"][12:]),
    replace(["executions", 0, "residuals", "gauss"], [[0.0]*13]*2),
    replace(["executions", 0, "whitney_fields", "B_barycentric_coefficients", 0, 0, 0, 0], 100.0),
    replace(["executions", 0, "counting_temporal_control", "temporal_defect"], "0"),
    replace(["controls", "full_volume_stationarity"], True),
    lambda p: p["nonclaims"].pop(),
])
def test_false_green_mutations_fail(mutation):
    original = verifier.load()
    packet = deepcopy(original)
    mutation(packet)
    # JSON types are part of the mutation, including bool-versus-int changes.
    assert json.dumps(packet, sort_keys=True) != json.dumps(original, sort_keys=True)
    with pytest.raises(ValueError):
        verifier.verify(packet)


@pytest.mark.parametrize("raw", ['{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}'])
def test_ambiguous_json_rejected(tmp_path, raw):
    path = tmp_path / "bad.json"
    path.write_text(raw, encoding="utf-8")
    with pytest.raises(ValueError):
        verifier.load(path)


def test_exact_time_action_control_and_gauge_failure_of_naive_current():
    a, b = [row["counting_temporal_control"] for row in verifier.load()["executions"]]
    assert Q(a["temporal_defect"]) == Q(663629599, 63369648)
    assert Q(a["hat_residual"][0]) == Q(44749, 55152)
    assert Q(a["hat_residual_squared_norm"]) == Q(10737664307, 3041743104)
    assert a["temporal_defect"] == b["temporal_defect"]
    assert Q(a["ordinary_constant_current_action_defect"]) == Q(4001, 1532)
    assert Q(b["ordinary_constant_current_action_defect"]) == Q(37965, 10724)


def test_temporal_refinement_bound_is_not_total_action_bound():
    # B(t)=t on fixed [0,1]. Direct exact integration and endpoint quadrature:
    # the total action defect contains h/4, while its interior correction is h²/12.
    for count in (4, 8, 16, 32):
        h = Q(1, count)
        right_energy = h*sum(Q(i, count)**2 for i in range(1, count+1))/2
        integrated_energy = Q(1, 6)
        defect = right_energy-integrated_energy
        assert defect == h/4+h*h/12
        assert defect-h/4 == h*h/12
        assert defect > h*h


def test_nonuniform_time_steps_need_additional_endpoint_terms():
    # Two slabs with B=0,1,2 and lengths 1,2: no common-h telescoping formula.
    h = [Q(1), Q(2)]
    b = [Q(0), Q(1), Q(2)]
    direct = sum(h[n]*(b[n+1]**2/2-(b[n]**2+b[n]*b[n+1]+b[n+1]**2)/6)
                 for n in range(2))
    wrong = h[0]*(b[-1]**2-b[0]**2)/4+sum(h[n]*(b[n+1]-b[n])**2/12 for n in range(2))
    assert direct != wrong


def test_cross_platform_explicit_utf8(monkeypatch):
    original = Path.read_text
    def windows_default(self, *args, **kwargs):
        if not args and "encoding" not in kwargs:
            kwargs["encoding"] = "cp1252"
        return original(self, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", windows_default)
    assert len(verifier.source_mesh()[0]) == 12
    assert verifier.verify(verifier.load())["gauge_histories"] == 2


def test_physical_charge_sign_from_weak_gradient_adjoint():
    # One oriented edge: a negative coefficient source charge moves with J=-q/h.
    # The action's +rho*phi has physical charge -rho under E=-d_t A-grad phi.
    d = np.array([[-1, 1]])
    rho0, rho1, current = np.array([1, 0]), np.array([0, 1]), np.array([-2])
    assert np.array_equal(rho1-rho0+d.T@current/2, np.zeros(2))
    physical_rho0, physical_rho1 = -rho0, -rho1
    weak_divergence = -d.T
    assert np.array_equal(physical_rho1-physical_rho0+weak_divergence@current/2, np.zeros(2))


def test_transitive_source_pins_rechecked_after_prior_success(monkeypatch):
    packet = verifier.load()
    verifier.verify(packet)
    original = Path.read_bytes
    target = verifier.ROOT / "Lean/Screen/NeutralPairJointStationaryWitness.lean"
    def changed_provider(path):
        data = original(path)
        return data+b"\n" if path == target else data
    monkeypatch.setattr(Path, "read_bytes", changed_provider)
    with pytest.raises(ValueError, match="source digest"):
        verifier.verify(packet)
