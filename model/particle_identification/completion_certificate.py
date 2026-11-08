"""Replay the finite interfaces supporting the particle-identification route."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import capacity, dense_effects, native_capacity, passive_clock, uncertainty
from . import native_field_response
from ..native_causal_dataflow import dataflow
from ..native_lorentz_transition import certificate as lorentz
from ..field_cone_compatibility import certificate as field_cone
from ..source_state_bridge import bridge as state_bridge
from ..protected_reference import reference as protected_reference
from ..gauge_defect_transport import transport as gauge_transport
from ..native_capacity_restriction import packet as capacity_restriction
from ..native_record_response import response as record_response
from ..native_record_causality import packet as record_causality
from ..native_gauge_strings import certificate as gauge_strings
from ..native_clock_observables import packet as clock_observables
from ..weighted_source_family import packet as weighted_family
from ..native_effect_coverage import packet as effect_coverage
from ..native_clock_comparison import packet as clock_comparison
from ..weighted_local_limit import packet as weighted_limit
from ..native_local_response import packet as local_response
from ..native_effect_synthesis import packet as effect_synthesis
from ..native_sector_readout import packet as sector_readout
from ..native_correlation_sectors import packet as correlation_sectors
from ..native_spatial_scaling import packet as spatial_scaling
from ..native_correlation_preparation import packet as correlation_preparation
from ..native_pair_response import packet as pair_response
from ..native_exchange_clock import packet as exchange_clock
from ..native_local_history_limit import packet as local_history_limit
from ..native_coexistent_effects import packet as coexistent_effects
from ..native_contact_response import packet as contact_response
from ..native_vacuum_response import packet as vacuum_response
from ..native_transfer_control import packet as transfer_control
from ..native_connected_contact import packet as connected_contact
from ..native_rooted_refinement import packet as rooted_refinement

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "completion_results.json"


def _same_json(left, right):
    """Preserve JSON scalar types when authenticating a scientific receipt."""
    encode = lambda value: json.dumps(value, sort_keys=True, allow_nan=False,
                                      separators=(",", ":"))
    return encode(left) == encode(right)


def source_bindings():
    paths = integration_source_paths(ROOT, {
        "domain_and_proof": "model/particle_identification/CAPACITY.md",
        "source_receipts": ["model/interacting_extension/higgs_join/results.json"],
        "independent_checks": list(capacity.SOURCE_FILES + passive_clock.SOURCE_FILES
                                   + native_capacity.SOURCE_FILES
                                   + state_bridge.SOURCE_FILES
                                   + protected_reference.SOURCE_FILES
                                   + gauge_transport.SOURCE_FILES
                                   + dataflow.SOURCE_FILES + lorentz.SOURCE_FILES
                                   + field_cone.SOURCE_FILES),
    })
    paths.update(uncertainty.source_bindings())
    for producer in (capacity_restriction, record_response, record_causality,
                     gauge_strings, clock_observables, weighted_family,
                     effect_coverage, clock_comparison, weighted_limit,
                     local_response, effect_synthesis, sector_readout,
                     correlation_sectors, spatial_scaling,
                     correlation_preparation, pair_response, exchange_clock,
                     local_history_limit, coexistent_effects, contact_response,
                     vacuum_response, transfer_control, connected_contact,
                     rooted_refinement):
        paths.update(producer.source_bindings())
        paths.add(str(producer.OUTPUT.relative_to(ROOT)))
    paths.update("model/particle_identification/" + name for name in (
        "DENSE_EFFECTS.md", "UNCERTAINTY.md", "SCATTERING.md", "ROADMAP.md",
        "NATIVE_FIELD_RESPONSE.md", "native_field_response.py",
        "test_native_field_response.py"))
    paths.update(("Lean/EventAlgebra/FiniteBuschGleason.lean",
                  "Lean/EventAlgebra/OperationalBornBridge.lean",
                  "Lean/EventAlgebra/Robertson.lean",
                  "model/native_lorentz_transition/results.json",
                  "model/gauge_defect_transport/results.json",
                  "model/protected_reference/results.json",
                  "model/field_cone_compatibility/results.json"))
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in sorted(paths)}


def report():
    before = source_bindings()
    frames, cones = lorentz.build(), field_cone.build()
    gauge = gauge_transport.certificate()
    protected = protected_reference.certificate()
    restrictions = capacity_restriction.produce()
    responses = record_response.report()
    causal_records = record_causality.produce()
    strings = gauge_strings.produce()
    observables = clock_observables.produce()
    weighted = weighted_family.produce()
    effects = effect_coverage.produce()
    clocks = clock_comparison.produce()
    local_states = weighted_limit.produce()
    local_records = local_response.produce()
    synthesized_effects = effect_synthesis.produce()
    sector_records = sector_readout.produce()
    correlation_records = correlation_sectors.produce()
    spatial_records = spatial_scaling.produce()
    prepared_correlations = correlation_preparation.produce()
    spatial_pairs = pair_response.produce()
    exchange_responses = exchange_clock.produce()
    local_history = local_history_limit.produce()
    coexistent_menu = coexistent_effects.produce()
    contact_records = contact_response.produce()
    vacuum_records = vacuum_response.produce()
    transfer_experiment = transfer_control.produce()
    connected_response = connected_contact.produce()
    rooted_windows = rooted_refinement.produce()
    for producer, value in ((capacity_restriction, restrictions),
                            (record_response, responses),
                            (record_causality, causal_records),
                            (gauge_strings, strings),
                            (clock_observables, observables),
                            (weighted_family, weighted),
                            (effect_coverage, effects),
                            (clock_comparison, clocks),
                            (weighted_limit, local_states),
                            (local_response, local_records),
                            (effect_synthesis, synthesized_effects),
                            (sector_readout, sector_records),
                            (correlation_sectors, correlation_records),
                            (spatial_scaling, spatial_records),
                            (correlation_preparation, prepared_correlations),
                            (pair_response, spatial_pairs),
                            (exchange_clock, exchange_responses),
                            (local_history_limit, local_history),
                            (coexistent_effects, coexistent_menu),
                            (contact_response, contact_records),
                            (vacuum_response, vacuum_records),
                            (transfer_control, transfer_experiment),
                            (connected_contact, connected_response),
                            (rooted_refinement, rooted_windows)):
        if not _same_json(value, strict_json.load(producer.OUTPUT)):
            raise ValueError(f"retained {producer.OUTPUT.parent.name} receipt differs")
    if not _same_json(frames, strict_json.load(lorentz.OUTPUT)):
        raise ValueError("retained Lorentz transition receipt differs")
    if not _same_json(cones, strict_json.load(field_cone.OUTPUT)):
        raise ValueError("retained field-cone receipt differs")
    if not _same_json(gauge, strict_json.load(ROOT / "model/gauge_defect_transport/results.json")):
        raise ValueError("retained gauge-defect transport receipt differs")
    if not _same_json(protected, strict_json.load(ROOT / "model/protected_reference/results.json")):
        raise ValueError("retained protected-reference receipt differs")
    result = {
        "schema": "oph.particle_identification.finite_interfaces.v1",
        "dense_rational_effects": dense_effects.build(),
        "compensated_probe_uncertainty": uncertainty.certificate(),
        "local_capacity_comparison": capacity.certificate(),
        "native_FIFO_capacity_equivalence": native_capacity.certificate(),
        "direct_native_state_and_carrier_bounds": state_bridge.certificate(),
        "protected_native_reference_and_gauge_join": protected,
        "complete_gauge_defect_transport": gauge,
        "fair_source_canonical_restriction": restrictions,
        "weighted_native_complete_record_response": responses,
        "disjoint_native_complete_record_equivalence": causal_records,
        "declared_gauge_tree_propagation": strings,
        "native_coherent_clock_observables": observables,
        "weighted_connected_source_extension": weighted,
        "native_operational_effect_coverage": effects,
        "native_clock_comparison": clocks,
        "weighted_canonical_local_limit": local_states,
        "native_local_record_response": local_records,
        "native_rank_one_effect_synthesis": synthesized_effects,
        "native_occupation_sector_readout": sector_records,
        "native_correlation_degree_responses": correlation_records,
        "native_routed_graph_motion": spatial_records,
        "native_correlation_preparation": prepared_correlations,
        "native_spatial_pair_response": spatial_pairs,
        "native_exchange_clock_response": exchange_responses,
        "native_complete_local_history_limit": local_history,
        "native_coexistent_effect_instruments": coexistent_menu,
        "native_contact_spectral_response": contact_records,
        "native_positive_energy_record_reconstruction": vacuum_records,
        "native_controlled_transfer_experiment": transfer_experiment,
        "native_connected_contact_response": connected_response,
        "native_rooted_refinement_windows": rooted_windows,
        "common_native_field_response": native_field_response.compose(
            prepared_correlations, exchange_responses),
        "passive_source_clock_boundary": passive_clock.certificate(),
        "causal_dataflow": dataflow.certificate(),
        "native_lorentz_transitions": frames,
        "field_cone_compatibility": cones,
        "scope": {
            "mathematical_prose_proofs": ["dense rational-effect extension",
                                           "complete local-history error composition",
                                           "passive source unitary obstruction",
                                           "native state and retained-carrier bounds",
                                           "protected stationary states and modular repair covariance",
                                           "complete finite gauge transport and conditional common reference",
                                           "fair canonical restriction and exact nesting",
                                           "actual weighted-source complete record response",
                                           "complete disjoint-command exchange with references",
                                           "declared tree Gauss-sector exhaustion and propagation",
                                           "reflected native record transfer and clock distinction",
                                           "adaptive independent-product dynamic restriction",
                                           "native coherent Pauli instruments and finite uncertainty",
                                           "weighted connected repair references and boundary comparison",
                                           "complete native configuration readout and coherent effect coverage",
                                           "operational separation of working and record-transfer clocks",
                                           "weighted canonical local limit without global enumeration",
                                           "native dependency cones and complete local response bounds",
                                           "bounded-precision native rank-one effect synthesis",
                                           "complete native tree readout across occupation sectors",
                                           "native correlation-degree decomposition and reflected record probes",
                                           "routed graph halo and discrete native displacement bounds",
                                           "arbitrary-input low-degree native response preparation",
                                           "same-sector spatial pair responses and quadratic diffusive limit",
                                           "native attempt comparison with the local exchange action",
                                           "common-source preparation and exchange-response error composition",
                                           "complete local record capacity limit with declared common clock bins",
                                           "bounded noncommuting coexistent native effect instruments",
                                           "native contact resolvent and observable finite correlated pole",
                                           "controlled positive-energy continuation of prepared native records",
                                           "controlled native transfer-ground preparation, time and phase response",
                                           "connected native contact dressing and comparison threshold",
                                           "exact corner topology and certified finite routed record windows",
                                           "finite null-transition and field-cone criteria"],
            "finite_fixtures_are_not_all_domain_proofs": True,
            "native_dense_effect_menu_derived": False,
            "active_record_assisted_encodings_excluded": False,
            "physical_capacity_to_local_action_map_derived": False,
            "fair_federation_and_weighted_clock_are_one_source_family": False,
            "native_record_relaxation_factors_are_particle_masses": False,
            "record_transfer_and_working_modular_clocks_identified": False,
            "coherent_clock_services_lifted_to_growing_family": False,
            "weighted_connected_extension_admitted_by_native_program": False,
            "physical_occupation_density_selected": False,
            "full_global_record_archive_identified_as_local": False,
            "record_correlation_degree_is_particle_species": False,
            "stationary_single_load_motion_is_physical_light_propagation": False,
            "bounded_repair_only_exact_reset_derived": False,
            "observable_preparation_is_full_state_preparation": False,
            "exchange_endpoint_comparison_is_complete_history_equivalence": False,
            "quadratic_diffusive_limit_is_interacting_relativistic_theory": False,
            "bounded_coexistent_menu_is_universal_effect_coverage": False,
            "local_history_limit_selects_physical_capacity_or_geometry": False,
            "finite_contact_mode_is_a_measured_particle": False,
            "reconstructed_transfer_vacuum_is_native_operational_vacuum": False,
            "controlled_transfer_time_is_autonomous_repair_time": False,
            "finite_rooted_windows_prove_all_radius_numeric_routing": False,
            "negative_parameter_contact_resolvent_identifies_a_positive_energy_particle": False,
            "common_native_to_gauge_Higgs_intertwiner_derived": False,
            "physical_position_momentum_identified": False,
            "relativistic_interacting_particle_spectrum_identified": False,
            "physical_collider_predictions_added": 0,
            "canonical_physical_premises_discharged": [],
        },
        "source_sha256": before,
    }
    if source_bindings() != before:
        raise RuntimeError("particle interface sources changed during replay")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", nargs="?", const=OUTPUT, type=Path)
    mode.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    value = report()
    if args.write is not None:
        args.write.write_text(json.dumps(value, sort_keys=True, indent=2,
                                        allow_nan=False) + "\n")
    elif not _same_json(strict_json.load(args.check), value):
        raise ValueError("particle interface receipt differs from complete replay")
    print("Verified finite effect, uncertainty and local-history interfaces; physical joins remain explicit.")


if __name__ == "__main__":
    main()
