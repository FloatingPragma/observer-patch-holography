"""Verify artifact bindings and replay the runnable model receipts."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

from tools import strict_json
from .replay_empirical import EMPIRICAL_REPLAYS

ROOT = Path(__file__).resolve().parents[1]
PACKETS = ("model/observer", "model/observer/p_response", "model/local_tick",
           "model/observer/scalar_response", "model/observer/source_scalar",
           "model/predictions", "model/observer/source_discrimination",
           "model/common_source", "model/spatial_interaction", "model/relativistic_scalar",
           "model/scalar_observations", "model/observer_reconstruction",
           "model/interacting_reconstruction", "model/coupled_reconstruction",
           "model/anharmonic_observations", "model/anharmonic_transfer",
           "model/identifiability_checks", "model/occupation_control",
           "model/identification_design", "model/primitive_observer",
           "model/native_quantum_history", "model/native_clock",
           "model/native_fluctuations", "model/empirical_diagnostics",
           "model/native_relational_readout", "model/native_circulation",
           "model/native_interface", "model/measurement_bounds", "model/clocked_common_observer",
           "model/event_circulation", "model/compiled_native_interface",
           "model/calibrated_measurement", "model/clocked_quantum_readout",
           "model/native_modular_reference", "model/native_effective_fields",
           "model/relative_observer_geometry", "model/relative_modular_dynamics",
           "model/transfer_modular_probe", "model/federated_observer",
           "model/native_equilibration", "model/relational_event_clock",
           "model/coherent_screen_observer",
           "model/federated_quantum_sector", "model/readout_controls",
           "model/observer_time_fields",
           "model/native_action", "model/joint_clock_process",
           "model/closed_prediction_protocol", "model/characteristic_constraints",
           "model/record_laws", "model/physical_space_time",
           "model/repair_action_audit", "model/laboratory_clock",
           "model/observer_correspondence",
           "model/quantum_observer_correspondence",
           "model/interacting_band_extension", "model/repair_identification", "model/compiled_clock_history",
           "model/repair_scheduler")
OCCUPATION_TRANSFER_REPLAY = ("model.occupation_control.transfer", None,
                              "model/occupation_control/transfer_results.json.gz")
REPLAYS = (
    ("model.common_source.native_family_protocol", None, "model/common_source/native_family_protocol_results.json"),
    ("model.common_source.native_repair_mass", None, "model/common_source/native_repair_mass_results.json"),
    ("model.common_source.native_gauge_matter_transport", None, "model/common_source/native_gauge_matter_transport_results.json"),
    ("model.common_source.native_sm_current", None, "model/common_source/native_sm_current_results.json"),
    ("model.common_source.native_matter_attachment", None, "model/common_source/native_matter_attachment_results.json"),
    ("model.common_source.native_particle_records", None, "model/common_source/native_particle_records_results.json"),
    ("model.common_source.native_gauge_object", None, "model/common_source/native_gauge_object_results.json"),
    ("model.common_source.native_persistent_response", None, "model/common_source/native_persistent_response_results.json"),
    ("model.common_source.native_virtual_exchange", None, "model/common_source/native_virtual_exchange_results.json"),
    ("model.common_source.native_code_transport", None, "model/common_source/native_code_transport_results.json"),
    ("model.common_source.native_charged_spectroscopy", None, "model/common_source/native_charged_spectroscopy_results.json"),
    ("model.observer_time_fields.preparation", None, "model/observer_time_fields/preparation_results.json"),
    ("model.native_clock.passive_age", None, "model/native_clock/passive_age_results.json"),
    ("model.repair_scheduler", None, "model/repair_scheduler/results.json"),
    ("model.predictions.prepared_clock", None, "model/predictions/prepared_clock_results.json"),
    ("model.repair_identification", None, "model/repair_identification/results.json"),
    ("model.observer_correspondence.quantum_experiment", None, "model/observer_correspondence/quantum_results.json"),
    ("model.observer_correspondence", None, "model/observer_correspondence/results.json"),
    ("model.record_laws.instrument_contract", None, "model/record_laws/instrument_results.json"),
    ("model.physical_space_time", None, "model/physical_space_time/results.json"),
    ("model.physical_space_time.occupation_field", None, "model/physical_space_time/occupation_field_results.json"),
    ("model.repair_action_audit", None, "model/repair_action_audit/results.json"),
    ("model.laboratory_clock", None, "model/laboratory_clock/results.json"),
    ("model.compiled_clock_history", None, "model/compiled_clock_history/results.json"),
    ("model.compiled_clock_history.law_calibration", None, "model/compiled_clock_history/law_calibration_results.json"),
    ("model.relative_observer_geometry.temporal_records", None, "model/relative_observer_geometry/temporal_results.json"),
    ("model.native_clock.clock_bridge", None, "model/native_clock/clock_bridge_results.json"),
    ("model.native_clock.bit_readout", None, "model/native_clock/bit_readout_results.json"),
    ("model.native_clock.configuration", None, "model/native_clock/configuration_results.json"),
    ("model.native_clock.trajectory_readout", None, "model/native_clock/trajectory_results.json"),
    ("model.native_clock.regional_preparation", None, "model/native_clock/regional_preparation_results.json"),
    ("model.observer_time_fields.canonical_sector", None, "model/observer_time_fields/canonical_sector_results.json"),
    ("model.predictions.bell_records", None, "model/predictions/bell_results.json"),
    ("model.relative_observer_geometry.binary_dynamics", None, "model/relative_observer_geometry/binary_results.json"),
    ("model.native_clock.host_clock", None, "model/native_clock/host_clock_results.json"),
    ("model.native_action.current_clock_action", None, "model/native_action/current_clock_action_results.json"),
    ("model.native_action.history_action", None, "model/native_action/history_action_results.json"),
    ("model.relative_observer_geometry.intrinsic_geometry", None, "model/relative_observer_geometry/intrinsic_results.json"),
    ("model.native_clock.comparison", None, "model/native_clock/comparison_results.json"),
    ("model.observer_time_fields.refinement", None, "model/observer_time_fields/refinement_results.json"),
    ("model.observer_time_fields", None, "model/observer_time_fields/results.json"),
    ("model.native_action.incident_refinement", None, "model/native_action/incident_results.json"),
    ("model.native_action.screen_compatibility", None, "model/native_action/screen_results.json"),
    ("model.native_action", None, "model/native_action/results.json"),
    ("model.joint_clock_process", None, "model/joint_clock_process/results.json"),
    ("model.closed_prediction_protocol", None, "model/closed_prediction_protocol/results.json"),
    ("model.characteristic_constraints", None, "model/characteristic_constraints/results.json"),
    ("model.native_equilibration", None, "model/native_equilibration/results.json"),
    ("model.native_equilibration.entropy_refinement", None, "model/native_equilibration/entropy_results.json"),
    ("model.native_equilibration.cut_completion", None, "model/native_equilibration/cut_results.json"),
    ("model.native_equilibration.record_thermodynamics", None, "model/native_equilibration/record_results.json"),
    ("model.native_equilibration.common_variation", None, "model/native_equilibration/variation_results.json"),
    ("model.relational_event_clock", None, "model/relational_event_clock/results.json"),
    ("model.federated_quantum_sector", None, "model/federated_quantum_sector/results.json"),
    ("model.coherent_screen_observer", None, "model/coherent_screen_observer/results.json"),
    ("model.readout_controls", None, "model/readout_controls/results.json"),
    ("model.federated_observer", None, "model/federated_observer/results.json"),
    ("model.native_effective_fields", None, "model/native_effective_fields/results.json"),
    ("model.native_effective_fields.neutral_transport", None, "model/native_effective_fields/neutral_transport_results.json"),
    ("model.native_effective_fields.neutral_instrument", None, "model/native_effective_fields/neutral_instrument_results.json"),
    ("model.relative_observer_geometry", None, "model/relative_observer_geometry/results.json"),
    ("model.relative_modular_dynamics", None, "model/relative_modular_dynamics/results.json"),
    ("model.transfer_modular_probe", None, "model/transfer_modular_probe/results.json"),
    ("model.native_modular_reference", None, "model/native_modular_reference/results.json"),
    ("model.native_modular_reference.preparation", None, "model/native_modular_reference/preparation_results.json"),
    ("model.native_modular_reference.full_flow", None, "model/native_modular_reference/full_flow_results.json"),
    ("model.native_modular_reference.regional", None, "model/native_modular_reference/regional_results.json"),
    ("model.native_modular_reference.completed_clock", None, "model/native_modular_reference/completed_clock_results.json"),
    ("model.native_modular_reference.common_caps", None, "model/native_modular_reference/common_caps_results.json"),
    ("model.native_modular_reference.local_preparation", None, "model/native_modular_reference/local_preparation_results.json"),
    ("model.native_modular_reference.local_history", None, "model/native_modular_reference/local_history_results.json"),
    ("model.native_modular_reference.cap_local_clock", None, "model/native_modular_reference/cap_local_clock_results.json"),
    ("model.native_modular_reference.reference_tomography", None, "model/native_modular_reference/reference_tomography_results.json"),
    ("model.native_modular_reference.learned_clock", None, "model/native_modular_reference/learned_clock_results.json"),
    ("model.native_modular_reference.phase_clock", None, "model/native_modular_reference/phase_clock_results.json"),
    ("model.native_modular_reference.native_census", None, "model/native_modular_reference/native_census_results.json"),
    ("model.native_modular_reference.census_clock", None, "model/native_modular_reference/census_clock_results.json"),
    ("model.native_modular_reference.history_clock", None, "model/native_modular_reference/history_clock_results.json.gz"),
    ("model.event_circulation", None, "model/event_circulation/results.json"),
    ("model.compiled_native_interface", None, "model/compiled_native_interface/results.json"),
    ("model.calibrated_measurement", None, "model/calibrated_measurement/results.json"),
    ("model.clocked_quantum_readout", None, "model/clocked_quantum_readout/results.json"),
    ("model.primitive_observer.complete_certificate", None, "model/primitive_observer/complete_certificate_results.json"),
    ("model.native_circulation", None, "model/native_circulation/results.json"),
    ("model.native_interface", None, "model/native_interface/results.json"),
    ("model.native_interface.coherence_budget", None, "model/native_interface/coherence_results.json"),
    ("model.native_interface.pointer", None, "model/native_interface/pointer_results.json"),
    ("model.native_interface.tomography", None, "model/native_interface/tomography_results.json"),
    ("model.measurement_bounds", None, "model/measurement_bounds/results.json"),
    ("model.native_relational_readout", None, "model/native_relational_readout/results.json"),
    ("model.native_quantum_history", None, "model/native_quantum_history/results.json"),
    ("model.native_clock", None, "model/native_clock/results.json"),
    ("model.native_fluctuations", None, "model/native_fluctuations/results.json"),
    ("model.empirical_diagnostics", None, "model/empirical_diagnostics/results.json"),
    ("model.primitive_observer.transfer", None, "model/primitive_observer/transfer_results.json"),
    ("model.common_source.open_observer", None, "model/common_source/open_observer_results.json"),
    ("model.occupation_control", None, "model/occupation_control/results.json"),
    ("model.occupation_control.shared_calibration", None, "model/occupation_control/shared_calibration_results.json"),
    OCCUPATION_TRANSFER_REPLAY,
    ("model.identification_design", None, "model/identification_design/results.json"),
    ("model.primitive_observer", None, "model/primitive_observer/results.json"),
    ("model.anharmonic_transfer", None, "model/anharmonic_transfer/results.json"),
    ("model.anharmonic_observations", None, "model/anharmonic_observations/results.json"),
    ("model.coupled_reconstruction", None, "model/coupled_reconstruction/results.json"),
    ("model.interacting_reconstruction", None, "model/interacting_reconstruction/results.json"),
    ("model.observer_reconstruction", None, "model/observer_reconstruction/results.json"),
    ("model.scalar_observations", None, "model/scalar_observations/results.json"),
    ("model.scalar_observations.occupation", None, "model/scalar_observations/occupation_results.json"),
    ("model.scalar_observations.experimental_readout", None,
     "model/scalar_observations/experimental_readout/results.json"),
    ("model.scalar_observations.experimental_readout.executable", None,
     "model/scalar_observations/experimental_readout/executable_results.json"),
    ("model.relativistic_scalar", None, "model/relativistic_scalar/results.json"),
    ("model.spatial_interaction", None, "model/spatial_interaction/results.json"),
    ("model.observer.source_selection", None, "model/observer/results/source_selection.json"),
    ("model.observer.family_port_identifiability", None,
     "model/observer/results/family_port_identifiability.json"),
    ("model.observer.native_loop_readout", None,
     "model/observer/results/native_loop_readout.json"),
    ("model.observer.coherent_continuation", None, "model/observer/results/coherent_continuation.json"),
    ("model.observer", "native", "model/observer/results/native.json"),
    ("model.observer.quantum_mixing_receipt", None, "model/observer/results/quantum_mixing.json"),
    ("model.observer", "capacity", "model/observer/results/capacity.json"),
    ("model.observer", "run", "model/observer/results/run.json"),
    ("model.observer.p_response", None, "model/observer/p_response/results.json"),
    ("model.local_tick", None, "model/local_tick/results.json"),
    ("model.observer.scalar_response", None, "model/observer/scalar_response/results.json"),
    ("model.observer.scalar_response", "--continue-field", "model/observer/scalar_response/continuation_results.json"),
    ("model.observer.scalar_response.local_evolution", None,
     "model/observer/scalar_response/local_results.json"),
    ("model.observer.source_scalar", None, "model/observer/source_scalar/results.json"),
    ("model.predictions", None, "model/predictions/results.json"),
    ("model.observer.source_discrimination", None,
     "model/observer/source_discrimination/results.json"),
    ("model.common_source", None, "model/common_source/results.json"),
    ("model.common_source.composition", None, "model/common_source/composition_results.json"),
    ("model.common_source.smooth_history", None, "model/common_source/smooth_history_results.json"),
    ("model.common_source.output_certificate", None, "model/common_source/output_certificate_results.json"),
    ("model.common_source.quantum_native", None, "model/common_source/quantum_native_results.json"),
    ("model.common_source.encoded_observer", None, "model/common_source/encoded_observer_results.json"),
    ("model.common_source.record_observer", None, "model/common_source/record_observer_results.json"),
    ("model.common_source.infinite_history", None, "model/common_source/infinite_results.json"),
    ("model.common_source.sharp_output", None, "model/common_source/sharp_output_results.json"),
    ("model.common_source.bridge", None, "model/common_source/bridge_results.json"),
    ("model.common_source.smooth_bridge", None, "model/common_source/smooth_results.json"),
    ("model.common_source.thermal_bridge", None, "model/common_source/thermal_results.json"),
    ("model.common_source.causality", None, "model/common_source/causality_results.json"),
    ("model.common_source.native_fields", None, "model/common_source/native_fields_results.json"),
    ("model.common_source.native_defect_decision", None,
     "model/common_source/native_defect_decision_results.json"),
    ("model.observer.record_recovery_receipt", None, "model/observer/results/record_recovery.json"),
    ("model.clocked_common_observer", None, "model/clocked_common_observer/results.json"),
    ("model.identifiability_checks.rule_observability", None, "model/identifiability_checks/rule_observability_results.json"),
    ("model.identifiability_checks.dephasing_observability", None, "model/identifiability_checks/dephasing_observability_results.json"),
    ("model.identifiability_checks.consensus_observability", None, "model/identifiability_checks/consensus_observability_results.json"),
    ("model.record_laws", None, "model/record_laws/results.json"),
    ("model.quantum_observer_correspondence", None,
     "model/quantum_observer_correspondence/results.json"),
    ("model.interacting_band_extension", None,
     "model/interacting_band_extension/results.json"),
    ("model.interacting_band_extension.thermal_certificate", None,
     "model/interacting_band_extension/thermal_results.json"),
    ("model.predictions.optical_clock", None, "model/predictions/optical_clock_results.json"),
    ("model.relative_observer_geometry.clock_reconstruction", None,
     "model/relative_observer_geometry/clock_reconstruction_results.json"),
    ("model.record_laws.burn_in", None, "model/record_laws/burn_in_results.json"),
    ("model.predictions.regional_optical", None, "model/predictions/regional_optical_results.json"),
    ("model.predictions.modular_optical", None, "model/predictions/modular_optical_results.json"),
    ("model.physical_space_time.continuum_field", None,
     "model/physical_space_time/continuum_field_results.json"),
    ("model.relative_modular_dynamics.passive", None, "model/relative_modular_dynamics/passive_results.json"),
    ("model.observer_time_fields.thermal", None, "model/observer_time_fields/thermal_results.json"),
    ("model.transfer_modular_probe.event_clock", None, "model/transfer_modular_probe/event_clock_results.json"),
    ("model.physical_space_time.observation_process", None,
     "model/physical_space_time/observation_process_results.json"),
    ("model.observer_time_fields.spatial_dynamics", None,
     "model/observer_time_fields/spatial_dynamics_results.json"),
    ("model.observer_time_fields.wave_reference", None,
     "model/observer_time_fields/wave_reference_results.json"),
    ("model.predictions.inverse_clock", None, "model/predictions/inverse_clock_results.json"),
    ("model.physical_space_time.local_continuum", None,
     "model/physical_space_time/local_continuum_results.json"),
    ("model.observer_time_fields.continuum_dynamics", None,
     "model/observer_time_fields/continuum_dynamics_results.json"),
    ("model.predictions.repair_response", None, "model/predictions/repair_response_results.json"),
    ("model.common_source.replay", "mass_response", "model/common_source/mass_response/results.json"),
    ("model.common_source.replay", "particle_identification", "model/common_source/particle_identification/results.json"),
    ("model.common_source.replay", "chiral_candidate", "model/common_source/chiral_candidate/results.json"),
    ("model.common_source.hierarchy_diagnosis.audit", None, "model/common_source/hierarchy_diagnosis/results.json"),
    ("model.common_source.full_repair_response.run", None, "model/common_source/full_repair_response/results.json"),
    ("model.empirical_diagnostics.cross_scale", None, "model/empirical_diagnostics/cross_scale_results.json"),
    ("model.scalar_observations.experimental_readout.shared_calibration", None,
     "model/scalar_observations/experimental_readout/shared_calibration_results.json"),
    ("model.common_source.chiral_candidate.protected_momentum", None,
     "model/common_source/chiral_candidate/protected_momentum_results.json"),
    ("model.common_source.chiral_candidate.record_two_point", None,
     "model/common_source/chiral_candidate/record_two_point_results.json"),
    ("model.native_equilibration.vacuum_geometry", None,
     "model/native_equilibration/vacuum_geometry_results.json"),
    ("model.common_source.full_repair_response.stationary_selection", None,
     "model/common_source/full_repair_response/stationary_selection_results.json"),
    ("model.common_source.chiral_candidate.dispersion_identifiability", None,
     "model/common_source/chiral_candidate/dispersion_identifiability_results.json"),
    ("model.common_source.overlap_regulator", None,
     "model/common_source/overlap_regulator_results.json"),
    ("model.common_source.calibrated_matter", None,
     "model/common_source/calibrated_matter_results.json"),
    ("model.common_source.calibrated_matching", None,
     "model/common_source/calibrated_matching_results.json"),
    ("model.common_source.calibrated_spectrum_observer", None,
     "model/common_source/calibrated_spectrum_observer_results.json"),
    ("model.native_equilibration.ball_modular", None,
     "model/native_equilibration/ball_modular_results.json"),
    ("model.common_source.parameter_selection", None,
     "model/common_source/parameter_selection_results.json"),
    ("model.common_source.fermion_ball", None,
     "model/common_source/fermion_ball_results.json"),
    ("model.common_source.calibrated_transport", None,
     "model/common_source/calibrated_transport_results.json"),
    ("model.native_equilibration.vacuum_repair", None,
     "model/native_equilibration/vacuum_repair_results.json"),
    ("model.common_source.repair_continuation", None,
     "model/common_source/repair_continuation_results.json"),
    ("model.common_source.repair_field_attachment", None,
     "model/common_source/repair_field_attachment_results.json"),
    ("model.native_equilibration.reciprocal_vacuum", None,
     "model/native_equilibration/reciprocal_vacuum_results.json"),
    ("model.common_source.native_staggered_field", None,
     "model/common_source/native_staggered_field_results.json"),
    ("model.common_source.native_field_history", None,
     "model/common_source/native_field_history_results.json"),
    ("model.common_source.native_charged_spectrum", None,
     "model/common_source/native_charged_spectrum_results.json"),
    ("model.native_equilibration.wedge_modular", None,
     "model/native_equilibration/wedge_modular_results.json"),
    ("model.common_source.native_response_compatibility", None,
     "model/common_source/native_response_compatibility_results.json"),
    ("model.particle_screening.screening", None,
     "model/particle_screening/screening_results.json"),
    ("model.active_code_energy.energy_bridge", None,
     "model/active_code_energy/energy_bridge_results.json"),
    ("model.clock_defect_energy.clock_count_inference", None,
     "model/clock_defect_energy/clock_count_inference_results.json"),
    ("model.local_repair_family.observer", None,
     "model/local_repair_family/observer_results.json"),
    ("model.local_repair_family.record_clock", None,
     "model/local_repair_family/record_clock_results.json"),
    ("model.screened_recovery.experiment", None,
     "model/screened_recovery/screened_recovery_results.json"),
    ("model.conserved_sector_action.sector_action", None,
     "model/conserved_sector_action/sector_action_results.json"),
    ("model.clock_defect_energy.shared_event_clock.run", None,
     "model/clock_defect_energy/shared_event_clock/results.json"),
    ("model.screened_transport.transport", None,
     "model/screened_transport/transport_results.json"),
    ("model.exchange_sector.statistics_packet", None,
     "model/exchange_sector/statistics_results.json"),
    ("model.exchange_sector.many_body_packet", None,
     "model/exchange_sector/many_body_results.json"),
    ("model.bound_pair", None, "model/bound_pair/results.json"),
    ("model.charged_dynamics.packet", None, "model/charged_dynamics/results.json"),
    ("model.particle_clock", None, "model/particle_clock/results.json"),
    ("model.field_regions", None, "model/field_regions/results.json"),
    ("model.field_stress", None, "model/field_stress/results.json"),
    ("model.field_tastes", None, "model/field_tastes/results.json"),
    ("model.support_refinement", None, "model/support_refinement/results.json"),
    ("model.particle_comparison.compare", None, "model/particle_comparison/results.json"),
    ("model.biased_reset", None, "model/biased_reset/results.json"),
    ("model.band_clock", None, "model/band_clock/results.json"),
    ("model.band_clock_budget", None, "model/band_clock_budget/results.json"),
    ("model.native_localization.experiment", None, "model/native_localization/results.json"),
    ("model.native_program.receipt", None, "model/native_program/results.json"),
    ("model.particle_screening.mobile_packet", None,
     "model/particle_screening/mobile_results.json"),
)


def replay_groups(replays=None):
    """Partition the declaration without splitting the empirical cache worker."""
    replays = REPLAYS if replays is None else tuple(replays)
    if len(set(replays)) != len(replays):
        raise ValueError("Replay declarations must be unique")
    groups, index, empirical_groups = [], 0, 0
    while index < len(replays):
        if replays[index:index+len(EMPIRICAL_REPLAYS)] == EMPIRICAL_REPLAYS:
            groups.append(EMPIRICAL_REPLAYS)
            empirical_groups += 1
            index += len(EMPIRICAL_REPLAYS)
        else:
            if replays[index] in EMPIRICAL_REPLAYS:
                raise ValueError("Empirical replays must remain contiguous and in their canonical order")
            groups.append((replays[index],))
            index += 1
    if empirical_groups != 1:
        raise ValueError("Exactly one complete empirical replay group is required")
    return tuple(groups)


def selected_replay_groups(shard_index=None, shard_count=None, *, replays=None):
    """Keep the default order; isolate expensive shared-cache work when sharded.

    Shard zero owns the empirical group. With three or more shards, shard one
    owns the all-case occupation transfer; remaining groups retain declaration
    order and are distributed round-robin over shards two onward. Two shards
    keep every non-empirical group together on shard one. Every new declared
    replay is included automatically; no replay is silently skipped.
    """
    groups = replay_groups(replays)
    if shard_index is None and shard_count is None:
        return groups
    if (type(shard_index) is not int or type(shard_count) is not int
            or not 1 <= shard_count <= len(groups)
            or not 0 <= shard_index < shard_count):
        raise ValueError(f"Supply both shard arguments: 1 <= count <= {len(groups)}, 0 <= index < count")
    if shard_count == 1:
        return groups
    transfer = (OCCUPATION_TRANSFER_REPLAY,)
    if shard_count >= 3 and transfer not in groups:
        raise ValueError("Dedicated occupation-transfer replay is missing")
    if shard_index == 0:
        return (EMPIRICAL_REPLAYS,)
    remaining = tuple(group for group in groups if group != EMPIRICAL_REPLAYS)
    if shard_count >= 3:
        if shard_index == 1:
            return (transfer,)
        remaining = tuple(group for group in remaining if group != transfer)
        return tuple(group for index, group in enumerate(remaining)
                     if index % (shard_count-2) == shard_index-2)
    return tuple(group for index, group in enumerate(remaining)
                 if index % (shard_count-1) == shard_index-1)


def verify_receipt_source_bindings(digests=None):
    """Check live replay source maps without recalculating their science.

    Recognized maps use their producers' repository or package-relative
    path convention. Local-only maps reject missing or escaping files.
    The mixed ``sources`` maps also contain external custody metadata and
    aliases; those and frozen captures outside REPLAYS retain their owning
    verifiers. This does not replace schema validation or scientific replay.
    """
    digests = {} if digests is None else digests
    root = ROOT.resolve()
    repository_fields = {"sources", "source_sha256", "source_bindings",
                         "input_sha256", "source_binding_sha256",
                         "source_control_bindings", "source_pins"}
    package_levels = {"implementation_sha256": 0, "source_interface_hashes": 1}
    # The endpoint audit's generic field name is scoped to its actual location.
    endpoint_pins = "root.retrospective_endpoint_audit.source.pins"

    def walk(value, receipt, location="root"):
        if isinstance(value, dict):
            for key, child in value.items():
                here = f"{location}.{key}"
                recognized = key in repository_fields or key in package_levels or here == endpoint_pins
                local_only = recognized and key != "sources"
                if local_only and (not isinstance(child, dict) or not child):
                    raise ValueError(f"Empty or malformed local source map: {receipt} at {here}")
                if recognized and isinstance(child, dict):
                    base = ((root/receipt).parents[package_levels[key]]
                            if key in package_levels else root)
                    for relative, expected in child.items():
                        if not isinstance(relative, str) or "://" in relative:
                            if local_only:
                                raise ValueError(f"Invalid local source path: {receipt} at {here}")
                            continue
                        path = Path(relative)
                        if path.is_absolute():
                            if local_only:
                                raise ValueError(f"Absolute local source path: {receipt} at {here}[{relative}]")
                            continue
                        try:
                            target = (base/path).resolve()
                            local_file = target.is_relative_to(root) and target.is_file()
                        except (OSError, ValueError):
                            local_file = False
                        if not local_file:
                            if local_only:
                                raise ValueError(f"Missing or escaping local source path: {receipt} at {here}[{relative}]")
                            continue
                        if target not in digests:
                            digests[target] = hashlib.sha256(target.read_bytes()).hexdigest()
                        if digests[target] != expected:
                            raise ValueError(
                                f"Receipt source digest mismatch: {receipt} at {here}[{relative}]")
                walk(child, receipt, here)
        elif isinstance(value, list):
            for index, child in enumerate(value):
                walk(child, receipt, f"{location}[{index}]")

    for _, _, receipt in REPLAYS:
        path = ROOT/receipt
        raw = path.read_bytes()
        if path.suffix == ".gz":
            raw = gzip.decompress(raw)
        walk(strict_json.loads(raw.decode("utf-8")), receipt)


def verify_artifacts():
    covered = set()
    digests = {}
    for name in PACKETS:
        folder = ROOT / name
        manifest = json.loads((folder / "manifest.json").read_text())
        files = manifest["files"]
        if not files:
            raise ValueError(f"Empty artifact manifest: {name}")
        for relative, expected in files.items():
            target = (folder / relative).resolve()
            if not target.is_relative_to(folder.resolve()):
                raise ValueError(f"Artifact outside its package: {name}/{relative}")
            actual = hashlib.sha256(target.read_bytes()).hexdigest()
            if actual != expected:
                raise ValueError(f"Artifact digest mismatch: {name}/{relative}")
            covered.add(target)
            digests[target] = actual
    for name in PACKETS:
        for target in (ROOT / name).rglob("*"):
            if (target.is_file() and (target.suffix in {".py", ".md", ".json"}
                                      or target.name.endswith(".json.gz"))
                    and target.name != "manifest.json"
                    and target.resolve() not in covered):
                raise ValueError(f"Unbound artifact: {target.relative_to(ROOT)}")
    verify_receipt_source_bindings(digests)
    return len(covered)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hashes-only", action="store_true",
                        help="check artifacts and live receipt source bindings without recalculation")
    parser.add_argument("--shard-index", type=int,
                        help="zero-based replay shard; requires --shard-count")
    parser.add_argument("--shard-count", type=int,
                        help="number of replay shards; every shard checks all artifact hashes")
    args = parser.parse_args()
    try:
        groups = selected_replay_groups(args.shard_index, args.shard_count)
    except ValueError as error:
        parser.error(str(error))
    count = verify_artifacts()
    print(f"Verified {count} model artifacts", flush=True)
    if args.hashes_only:
        return
    environment = dict(os.environ)
    environment.setdefault("OPENBLAS_NUM_THREADS", "1")
    environment.setdefault("OMP_NUM_THREADS", "1")
    environment.setdefault("MKL_NUM_THREADS", "1")
    environment.setdefault("VECLIB_MAXIMUM_THREADS", "1")
    selected_count = sum(map(len, groups))
    if args.shard_count is not None:
        print(f"Replaying shard {args.shard_index+1}/{args.shard_count}: "
              f"{selected_count} receipts", flush=True)
    for group in groups:
        if group == EMPIRICAL_REPLAYS:
            # Keep every canonical check, sharing only the expensive source
            # compilation caches in this one child. Failures still propagate.
            subprocess.run([sys.executable, "-m", "model.replay_empirical"],
                           cwd=ROOT, env=environment, check=True)
            continue
        (module, command, receipt), = group
        argv = [sys.executable, "-m", module]
        if command is not None:
            argv.append(command)
        argv.extend(["--check", receipt])
        print(f"Replaying {receipt}", flush=True)
        subprocess.run(argv, cwd=ROOT, env=environment, check=True)
    if args.shard_count is None:
        print(f"All {selected_count} model receipts reproduced", flush=True)
    else:
        print(f"All {selected_count} receipts in shard {args.shard_index+1}/{args.shard_count} reproduced",
              flush=True)


if __name__ == "__main__":
    main()
