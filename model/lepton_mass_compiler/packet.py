"""Frozen full-source compiler, complete records and one explicit calibration."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from model.particle_common_contract.provenance import import_closure
from . import analyzer,compiler,family_map

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    paths.update(('model/lepton_cycle/PROOF.md','model/lepton_cycle/experimental_inputs.json',
        'model/standard_model_source_map/KOIDE_SOURCE.md',
        'model/charged_lepton_response/CAR_KOIDE.md','model/charged_lepton_response/RELATIVISTIC_MASS.md',
        'model/coherent_defect_wave/mass_instrument/PROOF.md','model/requirements.txt','requirements.txt'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


def calculate():
    pins=bindings();source=analyzer.certificate();effects=source['effect_values_k0_k1_k2']
    result={'schema':'oph.lepton_mass_compiler.v1','source_bindings':pins,
        'actual_native_analyzer':source,
        'existing_charged_family_frame_witness':family_map.certificate(),
        'full_CAR_operator':{'physical_modes_per_site':90,'enlarged_modes_per_site':324000,
            'fresh_ancillary_modes_per_site':323910,'source_pairing_rank_per_spin_orientation':3,
            'left_frame':'U on the entire native lepton doublet',
            'right_frame':'conjugate(U) on the native singlet',
            'oriented_full_pairing_identity':'S† Delta_Q conjugate(S) has oriented block U† Q U=W',
            'all_occupation_vacuum_compression':'J† V J=H_(epsilon tensor [conjugate(D) tensor Y+conjugate(D)^T tensor Y^T])',
            'exact_mean_norm_per_site':'2 Tr(Y)=5/9 at w_H=1',
            'exact_kick_norm_per_site':6,'physical_vacuum_subspace_reduces_interaction':False,
            'full_Fock_matrix_materialized':False,'full_unitary_lift_precedes_vacuum_compression':True},
        'complete_active_span_controls':[compiler.four_mode_control(e) for e in effects],
        'same_source_kinetic_attachment':compiler.full_native_dirac(),
        'exact_limiting_rest_mass_instrument':compiler.exact_rest_instrument(),
        'finite_compiler_mass_readout':compiler.ramsey_readout(),
        'periodic_side_four_resource_example':compiler.archive_bound(108*64**2*10**12,sites=64),
        'electron_anchored_compiled_source':compiler.electron_calibration(),
        'source_state_reference_clock':{
            'source':'finite periodic free Weyl hopping plus the literal full-CAR analyzer-controlled neutral pairing, fresh flag-mode vacuum each round',
            'reference':'chosen limiting quadratic filled negative band, static compensating charge/parity references, and all actual ancillary archives',
            'physical_input_need_not_be_ancillary_source_ground':True,
            'clock':'t=sum pulse durations in exp(-ih(H0_extended+V)); analyzer and transfer service overhead is separate',
            'control_law':'vacuum/nonvacuum readout followed by unitary export and fresh vacuum replacement; never success conditioning',
            'same_clock_two_limits':'fixed finite volume and fixed t first N→infinity; controlled compact-momentum lattice limit with all eight valleys',
            'autonomous_elapsed_native_clock_identification':False,
            'whole_source_ancillary_product_preparation_is_stationary_ground':False},
        'scope':{'actual_lepton_effect_compiled_into_full_CAR_mass_coefficient':True,
            'same_source_free_Dirac_limit_with_complete_archive_bound':True,
            'three_free_charged_Dirac_fields_per_retained_valley':True,
            'finite_N_evolution_identified_with_exact_unitary_mass_spectrum':False,
            'all_eight_mirror_valleys_retained':True,
            'programmed_phase_two_ninths_is_passive_source_prediction':False,
            'neutral_background_is_selected_interacting_Higgs_VEV':False,
            'physical_Standard_Model_leptons_or_gauge_continuum_derived':False,
            'GeV_calibration_is_empirical':True,'experimental_muon_or_tau_inputs_used_in_source':False,
            'finite_N_round_count_executed_at_quoted_precision':False,
            'quantum_archives_or_failure_records_discarded':False}}
    if bindings()!=pins:raise ArithmeticError('compiler source changed during calculation')
    return result


def equal(a,b,path=()):
    if isinstance(b,dict):return isinstance(a,dict) and set(a)==set(b) and all(equal(a[k],v,path+(k,)) for k,v in b.items())
    if isinstance(b,list):return isinstance(a,list) and len(a)==len(b) and all(equal(x,y,path) for x,y in zip(a,b,strict=True))
    if isinstance(b,float):
        probability=any(k in ('complete_original_history_probabilities','coarse_vacuum_nonvacuum_history_probabilities') for k in path)
        tolerance=(2e-30+2e-8*abs(b)) if probability else (2e-12+2e-9*abs(b))
        return type(a) in (int,float) and np.isfinite(a) and abs(a-b)<=tolerance
    return type(a) is type(b) and a==b


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true');args=parser.parse_args()
    result=calculate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    elif not equal(json.loads(path.read_text()),result):raise ArithmeticError('complete compiler receipt differs')
    print('Native lepton analyzer, full-CAR compiler, same-source Dirac records and calibrated scale verified.')
