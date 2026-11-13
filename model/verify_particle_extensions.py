"""Additional live particle producers, independent of retained replay ordering.

Run ``python -m model.verify_particle_extensions --hashes-only`` for complete
live inventories, or omit that flag for actual producer replays. The dedicated
particle-extensions workflow composes this check with the retained observer
workflow without changing the scientific contracts that bind that workflow.
"""
import argparse
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
PACKAGES=(
    'openwave_beta_angular', 'openwave_lepton_widths', 'openwave_neutrino_matter',
    'openwave_spin_gradient_nuclei', 'openwave_rest_spectroscopy', 'openwave_colour_string',
    'openwave_common_electromagnetism', 'openwave_boosted_records',
    'openwave_single_valley_leptons', 'openwave_mass_catalogue',
    'openwave_longitudinal_review', 'openwave_boost_review',
    'openwave_string_positive_cone', 'openwave_common_hyperfine_spectrum',
    'openwave_strong_step_scaling', 'openwave_transverse_completion',
    'openwave_native_running_response', 'openwave_passive_clock_records',
    'openwave_quark_flux_records', 'openwave_universal_energy_geometry',
    'openwave_common_gravity_records', 'openwave_baryon_form_factors',
    'openwave_pure_su2_source', 'openwave_su2_transfer', 'openwave_su2_spectral',
    'openwave_gauge_completion_review',
    'openwave_hadron_refinement_flow', 'openwave_strange_composites',
    'openwave_light_nuclear_curve', 'openwave_native_centre_string',
    'openwave_dressed_core_response', 'openwave_dirac_maxwell_refinement',
    'openwave_autonomous_clock_channel', 'openwave_transverse_gravity',
    'openwave_su2_volume_gap', 'openwave_capacity_stress',
    'openwave_clock_string_review', 'openwave_quantum_observable_review',
    'openwave_baryon_dispersion', 'openwave_cluster_bandwidth',
    'openwave_nuclear_isospin', 'openwave_clock_phase_resource',
    'openwave_quotient_homotopy', 'openwave_resolved_core',
    'openwave_dirac_field_connection', 'openwave_gravity_unit_constraint',
    'openwave_horizon_capacity', 'openwave_su2_propagation',
    'openwave_baryon_poles', 'openwave_pole_weak_current',
    'openwave_strange_shell_current', 'openwave_coupled_haar_source',
    'openwave_coupled_haar_spectrum', 'openwave_euclidean_refinement',
    'openwave_source_ward', 'openwave_einstein_error',
    'openwave_radiation_capacity',
    'openwave_uniform_pole_bands', 'openwave_normalizable_particle_records',
    'openwave_scattering_current_support', 'openwave_haar_coupling_flow',
    'openwave_cubic_haar_source', 'openwave_gauge_refinement_clock',
    'openwave_staggered_ward', 'openwave_replicated_geometry',
    'openwave_capacity_perturbations',
    'openwave_full_strong_current', 'openwave_lambda_wavepacket',
    'openwave_cluster_weyl_sequence', 'openwave_refined_haar_source',
    'openwave_refined_haar_spectrum', 'openwave_haar_memory_closure',
    'openwave_tensor_record_transport', 'openwave_collective_gravity_strength',
    'openwave_capacity_transfer',
    'openwave_lambda_spectral_measure',
    'openwave_finite_phase_clock',
    'openwave_localized_quark_record',
    'openwave_nuclear_em_transfer',
    'openwave_fock_field_transport',
    'openwave_energy_schur',
    'openwave_collective_finite_strength',
    'openwave_gravity_noise_normalization',
    'openwave_capacity_normalization_transfer',
)

def source_bindings():
    """Reconstruct each producer's live import/receipt inventory, fail closed."""
    for name in PACKAGES:
        module=importlib.import_module(f'model.{name}.__main__')
        expected=module.pins()
        receipt=json.loads((ROOT/'model'/name/'results.json').read_text())
        if not expected or receipt.get('source_sha256')!=expected:
            raise ValueError(f'particle extension source inventory changed: {name}')
        if not str(receipt.get('verdict','')).startswith('PASS'):
            raise ValueError(f'particle extension producer failed: {name}')
    return len(PACKAGES)

def selected(index=None,count=None):
    if index is None and count is None:return PACKAGES
    if index is None or count is None or count<1 or not 0<=index<count:
        raise ValueError('both valid shard index and count required')
    return tuple(name for i,name in enumerate(PACKAGES)if i%count==index)

def replay(index=None,count=None):
    names=selected(index,count)
    environment=dict(os.environ)
    for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'):environment.setdefault(k,'1')
    for name in names:
        argv=[sys.executable,'-m',f'model.{name}','--check']
        if name=='openwave_boosted_records':
            upstream=environment.get('OPH_OPENWAVE_ROOT')
            python=environment.get('OPH_TAICHI_PYTHON')
            if not upstream or not python:
                raise ValueError('boosted replay requires pinned OPH_OPENWAVE_ROOT and OPH_TAICHI_PYTHON; no numerical fallback')
            argv+=['--openwave-root',upstream,'--python',python]
        subprocess.run(argv,cwd=ROOT,env=environment,check=True)
    return len(names)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hashes-only',action='store_true')
    parser.add_argument('--shard-index',type=int)
    parser.add_argument('--shard-count',type=int)
    args=parser.parse_args()
    selected(args.shard_index,args.shard_count)
    count=source_bindings()
    print(f'All {count} additional particle source inventories verified',flush=True)
    if not args.hashes_only:
        count=replay(args.shard_index,args.shard_count)
        print(f'All {count} selected additional particle receipts reproduced',flush=True)

if __name__=='__main__':
    main()
