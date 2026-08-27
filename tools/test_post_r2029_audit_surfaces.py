"""Regression gates for the post-r2029 adversarial audit.

These checks keep theorem hypotheses, no-go scope, numerical-certificate
status, and simulation custody aligned across the scientific surfaces.
"""

from __future__ import annotations

import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def _claims() -> dict[str, dict]:
    rows = yaml.safe_load(
        (ROOT / "claims/claim_registry.yaml").read_text(encoding="utf-8")
    )["claims"]
    return {row["claim_id"]: row for row in rows}


def _text(relative_path: str) -> str:
    return " ".join(
        (ROOT / relative_path).read_text(encoding="utf-8").split()
    )


def _json(relative_path: str) -> dict:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8"))


def test_clock_claims_preserve_nonzero_and_live_route_scope() -> None:
    claims = _claims()
    recurrence = claims["OPH-QFT-SOURCE-RECURRENCE-CLOCK-PRECURSOR"]["statement"]
    assert "standard Markov-chain period is one" in recurrence
    assert "no deterministic return interval exists" in recurrence
    assert "no single period exists" not in recurrence

    action = claims["OPH-GEOMETRY-PROPER-TIME-INTERNAL-ACTION"]["statement"]
    assert action.count("nonzero ledger") >= 2
    accrual = claims["OPH-GEOMETRY-SOURCE-CLOCK-RATE-ALONG-WORLDLINES"]["statement"]
    assert "For E nonzero" in accrual
    chain = claims["OPH-GEOMETRY-PROPER-LENGTH-CLOCKED-CHAIN"]
    combined = f'{chain["statement"]} {chain["oph_specific_delta"]}'
    assert "at least one rest between crossings" in combined
    assert "all-crossing case remains matched by a constant duration" in combined


def test_force_transfer_keeps_all_converse_hypotheses() -> None:
    statement = _claims()["OPH-GEOMETRY-TIMELIKE-CLASS-FORCE-LAW"]["statement"]
    for token in (
        "q and h are nonzero",
        "respective scaled Ampere evolutions",
        "original and delayed hopping currents",
        "original field static across the crossing step",
    ):
        assert token in statement


def test_divided_step_charge_and_transport_claims_require_nonzero_h() -> None:
    claims = _claims()
    charge_fixed = claims["OPH-GEOMETRY-CHARGE-FIXED-INTERACTION"]["statement"]
    coupling = claims["OPH-GEOMETRY-PORT-CHARGE-MINIMAL-COUPLING"]["statement"]
    transport = claims["OPH-GEOMETRY-WORLDLINE-HOP-TRANSPORT"]["statement"]
    assert "Field sector at a declared nonzero step h" in charge_fixed
    assert "For a declared nonzero step h" in coupling
    assert "For a declared nonzero step h" in transport
    assert "any divided-step conclusion is asserted at h = 0" in claims[
        "OPH-GEOMETRY-PORT-CHARGE-MINIMAL-COUPLING"
    ]["falsifier"]
    assert "a divided-step conclusion is asserted at h = 0" in claims[
        "OPH-GEOMETRY-WORLDLINE-HOP-TRANSPORT"
    ]["falsifier"]

    charge_lean = _text("Lean/Geometry/ChargeFixedInteraction.lean")
    for token in (
        "theorem induced_continuity (κ h : ℝ) (hh : h ≠ 0)",
        "theorem coupledAction_reduction (κ h : ℝ) (hh : h ≠ 0)",
        "theorem coupled_field_equations (κ h : ℝ) (hh : h ≠ 0)",
        "theorem identification_iff_kappa (κ h q : ℝ) (hh : h ≠ 0)",
    ):
        assert token in charge_lean

    coupling_lean = _text("Lean/Geometry/PortChargeMinimalCoupling.lean")
    for token in (
        "theorem hopping_continuity (q h : ℝ) (hh : h ≠ 0)",
        "theorem sourcePairing_gauge_invariant_iff (h : ℝ) (hh : h ≠ 0)",
        "theorem hopping_work_energy (q h : ℝ) (hh : h ≠ 0)",
    ):
        assert token in coupling_lean
    transport_lean = _text("Lean/Geometry/WorldlineHopTransport.lean")
    for token in (
        "theorem transport_load_charge_fixed (q h τ : ℝ) (hh : h ≠ 0)",
        "theorem transported_field_equations (q h τ : ℝ) (hh : h ≠ 0)",
        "theorem unit_unconstrained_by_field_sector (q h : ℝ) (hh : h ≠ 0)",
    ):
        assert token in transport_lean

    flagship = _text("flagship/from_observer_consensus_to_standard_physics.tex")
    observers = _text("paper/observers_are_all_you_need.tex")
    assert "at a declared nonzero step $h$" in flagship
    assert "At a declared nonzero step \\(h\\)" in observers
    assert "declared nonzero field step \\(h\\)" in observers


def test_generic_mode_oscillator_keeps_lower_courant_and_nonzero_h() -> None:
    claim = _claims()["OPH-EM-CARRIER-MODE-OSCILLATORS"]
    assert "declared nonzero step h with 0 <= h^2 lam <= 4" in claim["statement"]
    assert "Under h != 0 and 0 <= h^2 lam <= 4" in claim["falsifier"]
    lean = _text("Lean/Screen/CarrierModeOscillators.lean")
    assert (
        "theorem modeOscillator (h : ℝ) (hh : h ≠ 0) (lam : ℝ) "
        "(h0 : 0 ≤ h ^ 2 * lam) (h4 : h ^ 2 * lam ≤ 4)"
    ) in lean
    observers = _text("paper/observers_are_all_you_need.tex")
    assert "declared nonzero step \\(h\\) with \\(0\\le h^2\\lambda\\le 4\\)" in observers


def test_native_decision_and_coefficient_flow_trust_boundaries_are_public() -> None:
    repair = (ROOT / "Lean/Screen/PortGramRepairCovariance.lean").read_text(
        encoding="utf-8"
    )
    isometry = (ROOT / "Lean/Screen/PortGramA5Isometry.lean").read_text(
        encoding="utf-8"
    )
    assert sum(line.strip() == "native_decide" for line in repair.splitlines()) == 7
    assert sum(line.strip() == "native_decide" for line in isometry.splitlines()) == 6

    completion = _claims()["OPH-GR-PORT-RESPONSE-COMPLETION"]
    assert "Thirteen finite covariance and isometry table identities" in completion[
        "statement"
    ]
    assert "not kernel-only proofs" in completion["statement"]
    flow = _claims()["OPH-EM-CARRIER-EVOLUTION-FLOW"]
    assert "coefficient-state group" in flow["statement"]
    assert "coefficient-to-field map need not be injective" in flow["statement"]
    assert "presented as faithful without nonzero independence or a quotient" in flow[
        "falsifier"
    ]

    for relative_path in (
        "flagship/from_observer_consensus_to_standard_physics.tex",
        "paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex",
        "paper/observers_are_all_you_need.tex",
        "paper/screen_microphysics_and_observer_synchronization.tex",
    ):
        surface = _text(relative_path)
        assert "native" in surface and "kernel-only" in surface


def test_flagship_uses_mixed_verification_and_declared_closure_map_status() -> None:
    flagship = _text("flagship/from_observer_consensus_to_standard_physics.tex")
    assert "its proof machine-checked" not in flagship
    assert "Lean theorem, exact executable certificate" in flagship
    assert "conditional paper argument, numerical receipt, or explicit" in flagship
    assert "unproved premise" in flagship
    assert "come from the closure hypothesis directly" not in flagship
    assert "Neither clears the source-only selection" in flagship
    assert "the choice of map is itself declared" in flagship


def test_exact_courant_boundary_is_field_energy_not_potential_boundedness() -> None:
    lean = _text("Lean/Screen/ScaledMaxwellStability.lean")
    for token in (
        "theorem critical_mode",
        "theorem instability_at_golden_boundary",
        "theorem zeroCurrentElectricEnergyBounded_iff_window",
        "No bound on the gauge potential",
    ):
        assert token in lean

    claim = _claims()["OPH-EM-SCALED-MAXWELL-STABILITY-ACTION"]["statement"]
    assert "at equality" in claim
    assert "electric energy proportional to (2n+1)^2" in claim
    assert "no gauge-potential bound" in claim


def test_d10_germ_cannot_promote_and_exact_slice_stays_diagnostic() -> None:
    payload = _json(
        "code/particles/runs/calibration/"
        "d10_ew_tau2_current_carrier_obstruction.json"
    )
    assert payload["diagnostic_only"] is True
    direction = payload["direction_obstruction"]
    assert direction["global_no_go_inferred_from_germ"] is False
    nonlinear = payload["closed_form_nonlinear_point_test"]
    assert nonlinear["obstruction_established"] is True
    assert nonlinear["exact_interval_certificate"]["zero_excluded"] is True
    assert nonlinear["nonlinear_routes_left_open"]


def test_w5_classifies_quartic_minima_without_closing_boundary_route() -> None:
    cubic = _json(
        "code/particles/runs/leptons/charged_entropic_branch_no_go.json"
    )
    assert cubic["promotion_allowed"] is False
    gates = cubic["epistemic_gates"]
    assert gates["cubic_global_extremum_proved"] is True
    assert gates["finite_seed_search_can_close_route"] is False
    assert gates["quartic_global_optimality_proved"] is True
    assert gates["strict_full_support_quartic_no_go_conditional"] is True
    assert gates["quartic_packet_globally_excluded"] is False
    assert gates["full_entropic_mechanism_no_go_certified"] is False

    quartic = _json(
        "code/particles/runs/flavor/entropy_w5_shape_certificate.json"
    )
    assert quartic["schema"] == "oph.entropy_w5_shape_certificate.v3"
    boundary = quartic["exhaustiveness_boundary"]
    assert boundary["global_minimizer_classification_proved"] is True
    assert boundary["full_critical_orbit_classification_proved"] is False
    assert boundary["quartic_packet_globally_excluded"] is False
    assert "zero-weight closed-simplex boundary states" in boundary[
        "viable_routes_not_excluded"
    ]
    assert quartic["comparison"]["boundary_reentry_compare"]["verdict"] == (
        "BOUNDARY_BRANCH_CAN_MATCH_CENTRAL_SHAPE_COMPARE_ONLY"
    )


def test_flagship_simulation_language_matches_actual_custody() -> None:
    flagship = _text(
        "flagship/from_observer_consensus_to_standard_physics.tex"
    )
    assert "$65{,}536$ finite patch rows" in flagship
    assert "$2{,}048$ patch-observer neighborhoods" in flagship
    assert "finite-consensus theorem receipt nevertheless fails" in flagship
    assert "source-to-observer contract is also false" in flagship
    assert "b52196b296435d704b14d005d1f69caaaa662f97" in flagship
    assert "evidence/e6_64k_dense_20260820" in flagship
    assert "An unpublished local simulator bundle" not in flagship


def test_mandatory_suite_executes_the_affected_receipt_gates() -> None:
    suite = _text("tools/run_mandatory_suite.py")
    for relative_path in (
        "code/particles/calibration/test_d10_ew_tau2_current_carrier_obstruction.py",
        "code/particles/leptons/test_charged_entropic_branch_no_go.py",
        "code/particles/flavor/test_entropy_w5_shape_certificate.py",
        "tools/test_post_r2029_audit_surfaces.py",
    ):
        assert relative_path in suite


def test_measurement_claims_require_a_declared_instrument_and_state_support() -> None:
    flagship = _text("flagship/from_observer_consensus_to_standard_physics.tex")
    reality = _text("paper/reality_as_consensus_protocol.tex")
    screen = _text("paper/screen_microphysics_and_observer_synchronization.tex")
    assert "swap-twisted instrument has the same effects and different updates" in flagship
    assert "declare the selective outcome instrument to be the L\\\"uders map" in reality
    assert "an immediate reread with the same projector" in reality
    assert "explicit persistence hypothesis" in reality
    assert "absence of a touched-support marker" in reality
    assert "effect projectors alone do not select this update" in screen
    assert "same-effect non-L\\\"uders" in screen


def test_thermodynamic_scope_separates_stationary_kernels_from_normalizers() -> None:
    claim = _claims()["OPH-THERMO-FOUR-LAW-PACKAGE"]["statement"]
    assert "five typed source and physical receipts" in claim
    assert "They are not a source-realized theorem of the three axioms" in claim
    assert "stationary law" in claim
    assert "deterministic strict-descent normalizer is a distinct map with no entropy inequality" in claim
    observers = _text("paper/observers_are_all_you_need.tex")
    assert "-\\Delta D" in observers
    assert "data processing gives \\(\\Delta D\\le0\\)" in observers


def test_wz_and_quark_surfaces_do_not_promote_declared_maps_or_close_future_routes() -> None:
    particle = _text("paper/deriving_the_particle_zoo_from_observer_consistency.tex")
    assert "target-free declared-map running/chart coordinates" in particle
    assert "Conditional on the current registered source signature" in particle
    assert "does not define a unique quark mass spectrum" in particle
    quark = _claims()["OPH-QUARK-YUKAWA-TWO-MODULUS-NONIDENTIFIABILITY"]
    assert "current registered A1--A3 source signature" in quark["statement"]
    assert "not a machine proof of registry completeness" in quark["statement"]
    assert "generic_six_scalar_selector_open" in quark["status"]
    artifact = _json(
        "code/particles/runs/flavor/"
        "quark_axiom_level_yukawa_moduli_nonidentifiability.json"
    )
    assert artifact["machine_proof"] is False
    assert "current registered source signature" in " ".join(artifact["notes"])
    assert artifact["registered_signature_completeness_verified_by_emitter"] is False


def test_heat_kernel_formula_and_finite_diagnostics_keep_their_conditions() -> None:
    paper = _text("paper/tex_fragments/PAPER.tex")
    assert "log(p_0/d_0)-\\log(p_j/d_j)" in paper
    assert "H = -\\frac{g^2}{2} \\Delta_G" in paper
    assert "H_{\\mathrm{edge}} = \\Delta_G" in paper
    assert "separate Gibbs parameter is \\(t\\)" in paper
    assert "conjugacy-invariant and inverse-symmetric" in paper
    assert "character average is only its normalized trace" in paper
    claim = _claims()["OPH-EDGE-HEAT-KERNEL-FINITE-DIAGNOSTICS"]
    assert "Z2 has only one nontrivial eigenvalue" in claim["statement"]
    assert "No limit theorem" in claim["statement"]
    assert "no frozen scientific receipt" in claim["statement"]
    assert "declared_finite_diagnostic_without_frozen_receipt" in claim["status"]
    assert "earlier truncated SU3 numerical table" in claim["statement"]
    assert "withheld" in claim["statement"]


def test_maxent_locality_and_markov_alignment_hypotheses_are_not_inferred() -> None:
    paper = _text("paper/tex_fragments/PAPER.tex")
    claim = _claims()["OPH-GR-D4-NULLNET-STANDARDNESS"]["statement"]
    assert "relative interior of the attainable moment set" in paper
    assert "Boundary optima require restriction to their support" in paper
    assert "Modular support-envelope inclusion" in paper
    assert "need not be inner in the algebra itself" in paper
    assert "permute isomorphic central blocks" in paper
    assert "faithful/interior stagewise MaxEnt states" in claim
    assert "specified endpoint-to-cell map" in claim
    assert "finite-range locality alone is insufficient" in claim
    assert "boundary maxent optima do not supply" in claim.lower()
    assert "required uniformly for every half-line" in claim
    assert "cyclic for their relative commutant" in claim
    assert "weak additivity is also a separate cofinal receipt" in claim
    assert "not derived from one half-line receipt" in claim
    assert "finite one-particle modular compression" in claim
    assert "do not instantiate those scaling-limit receipts" in claim


def test_field_subset_search_records_nonlumpability_and_keeps_other_maps_open() -> None:
    receipt = _json("code/thermodynamics/runtime/collar_matrix_realization_probe.json")
    audit = receipt["raw_coarsening_audit"]
    assert audit["distinct_induced_partition_count"] == 4
    scope = audit["enumeration_scope"]
    assert scope["arbitrary_set_partitions_enumerated"] is False
    assert scope["weak_lumpability_tested"] is False
    repair_rows = [
        row for row in audit["rows"]
        if row["packet_fields"] == ["repair_load_bucket"]
    ]
    assert len(repair_rows) == 1
    repair = repair_rows[0]
    assert repair["state_count"] == 8
    assert repair["fine_chain_strongly_lumpable_at_tolerance"] is False
    assert repair["fine_chain_strong_lumpability_max_err"] > 0.9
    claim = _claims()["OPH-THERMO-FOUR-LAW-PACKAGE"]["statement"]
    assert "only four distinct partitions" in claim
    assert "not a certified Markov quotient" in claim
    assert "weakly lumpable" in claim


def test_particle_exact_formulae_preserve_degeneracy_and_domain_boundaries() -> None:
    particle = _text("paper/deriving_the_particle_zoo_from_observer_consistency.tex")
    assert "distinguished calibrated positive gap vector" in particle
    assert "simple and strictly positive" in particle
    assert "contains one zero singular value, its zero column retains" in particle
    assert "If \\(P=0\\), all three off-diagonal entries vanish" in particle
    assert "strictly positive Thomson limit" in particle
    assert "local existence and uniqueness" in particle
    assert "full-interval non-blowup and domain hypotheses" in particle
    assert "matching output lies in the domain of the next beta system" in particle


def test_yang_mills_legacy_gap_does_not_fake_the_missing_bridge() -> None:
    gap = _text("Lean/ObserverPatchHolography/YangMillsGap.lean")
    witness = _text("Lean/ObserverPatchHolography/YangMillsGapWitness.lean")
    assert "positive rates with the uniform floor and `prop_8_1`" in gap
    assert "does **not** consume `lemma_7_2`" in gap
    assert "source theorem connecting the uniform-fiber relaxation data" in gap
    assert "fixed space is `{0}`" in witness
    assert "space of constant functions is all of `W`, not `{0}`" in witness


def test_screen_variational_and_finite_descent_theorems_have_complete_hypotheses() -> None:
    screen = _text("paper/screen_microphysics_and_observer_synchronization.tex")
    assert "rate support is bidirectional" in screen
    assert "modulo that represented commutant" in screen
    assert "K,K'\\in\\mathcal A" in screen
    assert "minimum-cost locus contains a six-axis configuration" in screen
    assert "source selector first minimizes the additive cost" in screen
    assert "Treat a repair of a locally repaired datum" in screen
    assert "as the identity, not as a step" in screen
    assert "every accepted nonidentity repair strictly lowers" in screen
    assert "every state outside the declared visible local normal forms admits such a repair" in screen
    assert "one common finite commutative event algebra" in screen
    assert "Pointwise arrows on regulator-dependent finite sets do not by themselves" in screen
    assert "continuum convergence receipt" in screen
    assert "d=q_i-q_j\\ge2" in screen
    assert "uphill and equalizing-reversal moves are not repair steps" in screen
    assert "every accepted settling path" in screen


def test_stitch_symmetry_no_go_requires_no_fixed_admissible_stitch() -> None:
    particle = _text("paper/deriving_the_particle_zoo_from_observer_consistency.tex")
    assert "its induced action fixes no admissible stitch" in particle
    assert "Without the fixed-point-free clause" in particle
    assert "a unique natural stitch need only be invariant" in particle


def test_consensus_confluence_is_quotient_not_representative_uniqueness() -> None:
    reality = _text("paper/reality_as_consensus_protocol.tex")
    assert "unique \\emph{quotient} normal form" in reality
    assert "Microscopic terminal representatives need not be equal" in reality
    assert "only its quotient is canonical" in reality
    assert "no representative equality follows from quotient confluence" in reality
    assert "Let \\(\\overline M:Q\\to Y\\) be a quotient observable" in reality
    assert "quotient-level confluence theorem" in reality
    assert "supremum product-pseudometric distance" in reality
    assert "does not follow from vanishing pseudometric defects alone" in reality


def test_qumond_crosscheck_remains_uncertified_and_nonpromoting() -> None:
    receipt = _json(
        "code/cosmology/rar_deep_regime/receipts/"
        "qumond_quadrupole_crosscheck.json"
    )
    assert receipt["physical_claim"] is False
    assert receipt["promotion_allowed"] is False
    rows = [
        row for key, row in receipt["field_law"].items()
        if key.endswith("_inputs")
    ]
    assert len(rows) == 3
    assert all(row["numerically_certified"] is False for row in rows)
    assert sum(row["integration_warning_count"] for row in rows) == 4903


def test_registry_epistemic_classes_and_geometry_assumptions_match_theorems() -> None:
    claims = _claims()
    assert claims["OPH-GR-E2E-BRANCH-ENTRY"]["claim_class"] == "conditional_implication"
    assert claims["OPH-EDGE-HEAT-KERNEL-FINITE-DIAGNOSTICS"]["claim_class"] == "declared_structure"
    for claim_id in (
        "OPH-DM-MATCHED-OBSERVABLE-DIAGNOSTIC",
        "OPH-DM-ROTATION-CURVE-PROFILE-LIKELIHOOD",
    ):
        assert claims[claim_id]["claim_class"] == "empirical_implementation"

    geometry = claims["OPH-GR-D3H-GEOMETRY-PRODUCER"]
    assert "Refinement generally changes the literal complex" in geometry["statement"]
    assert "certified_simplicial_refinement_and_topology_receipt" in geometry["assumptions"]

    nullnet = claims["OPH-GR-D4-NULLNET-STANDARDNESS"]
    for token in (
        "EC_aligned_Markov_replacement_control",
        "endpoint_cell_uniform_locality_control",
        "null_net_weak_additivity_receipt",
        "standard_HSMI_relative_commutant_cyclicity_receipt",
    ):
        assert token in nullnet["assumptions"]
    for claim_id in ("OPH-GR-D5", "OPH-GR-D5A-ABSOLUTE-EINSTEIN"):
        assert "EC_aligned_Markov_replacement_control" in claims[claim_id]["assumptions"]

    events = claims["OPH-GR-D4B-EVENT-MANIFOLD"]["statement"]
    assert "second-countable nonseparated uniform space" in events
    consensus = claims["OPH-CONS-D1"]
    assert "unique schedule-independent quotient normal form" in consensus["statement"]
    assert "distinct declared quotient observables" in consensus["falsifier"]


def test_cosmology_and_source_no_gos_leave_alternative_typed_routes_open() -> None:
    radial = _text("paper/tex_fragments/RADIAL_LIFT_THEOREMS_330.tex")
    assert "not asserted to be exhaustive or mutually exclusive" in radial
    assert "another fully typed orbit-separating construction remains admissible" in radial

    scale = _text("paper/tex_fragments/PHYSICAL_SCALE_BRIDGE_THEOREMS.tex")
    assert "angular label \\(\\ell_{\\rm src}\\) alone does not determine" in scale
    assert "bounded Lipschitz" in scale
    assert "uniform tightness" in scale

    supplement = _text("paper/tex_fragments/DERIVATION_TECHNICAL_SUPPLEMENT_PORT.tex")
    assert "The unnormalized holonomy itself tends to the identity for every finite" in supplement
    assert "No gauge- or matter-sector result supplied in this packet selects" in supplement

    paper = _text("paper/tex_fragments/PAPER.tex")
    assert "this null-data route cannot fix" in paper
    assert "Bound neutral composites, dipole or higher-multipole sources" in paper
