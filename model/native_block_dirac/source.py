"""Admit the actual native constituent and its complete local reference."""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from model.native_abelian_field_join import certificate as native

PARENT = 'model/native_abelian_field_join/results.json'
PARENT_SHA = 'f6aadcd79eec889de0c3243925c950e482bd456ea5c22fcab6a8406978472327'
PREPARATION = 'model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@lru_cache(None)
def admit():
    if digest(ROOT/PARENT) != PARENT_SHA:
        raise ValueError('the native constituent/reference receipt changed')
    saved = json.loads((ROOT/PARENT).read_text())
    native.validate_pins(saved['input_sha256'])
    if saved != native.certificate():
        raise ValueError('fresh native source/reference controls differ')
    ref = saved['local_reference']
    if not (ref['complete_reference_dimension'] == 2**31+1
            and ref['vacuum_shared_and_unique']
            and ref['proper_cube_covariance_controls'] == 24):
        raise ValueError('the complete old/dual/port reference is required')
    low, high = map(F, saved['intervals']['native_selected_E_interval'])
    if not 0 < low < high:
        raise ValueError('positive actual native spectral channel required')
    prep = json.loads((ROOT/PREPARATION).read_text())['preparation']
    if prep['actual_local_preparation']['preparation_success_probability'] != '1':
        raise ValueError('native algebraic trial must be deterministic')
    return saved, prep


def preparation_parameters():
    _, prep = admit()
    p = F(prep['atomic_cube_trial']['exact_interacting_ground_trial_overlap_lower'])
    epsilon = F(prep['finite_energy_agreement_filter']['mixed_ground_excited_amplitude_upper'])
    length = int(prep['finite_energy_agreement_filter']['explicit_integer_M'])
    odds = (1-p)/p
    contraction = odds+2*epsilon*epsilon
    assert 0 < odds < contraction < 1
    return odds, contraction, length


def constituent_bounds():
    saved, prep = admit()
    from model.coherent_defect_wave.physical_attachment.chiral import representation
    actual = representation()
    basis, currents = actual[0], actual[2]
    assert basis[9] == (3, 4)
    assert currents[-1][9, 9] == 1
    assert all(g[9, 9] == 0 for g in currents[:-1])
    gaps = {
        'complete_neutral_source': F(prep['finite_energy_agreement_filter']['spectral_gap_lower']),
        'both_complete_port_sectors': F(saved['intervals']['complete_port_source_gap_lower']),
        'all_old_and_dual_occupied_reference_sectors': F(saved['intervals']['complete_old_and_dual_reference_gap_lower'])}
    assert min(gaps.values()) > 0
    native_norm = F(prep['finite_energy_agreement_filter']['full_absolute_H_norm_upper'])
    assert native_norm == 8*57+24+F(360, 2**34)
    return {'full_physical_constituent_gap_lower': str(min(gaps.values())),
            'complete_sector_gap_lowers': {k: str(v) for k, v in gaps.items()},
            'gap_domain': 'complete native physical constituent with shared vacuum, old/dual reference and both ports',
            'raw_full_kinematic_native_H_norm_upper': str(native_norm),
            'raw_local_reference_H_norm': '120',
            'selected_source_charge_convention': 'Y=+1 for addition, Y=-1 for removal; port Q is the same source charge',
            'underlying_ec_current_Y': '1', 'underlying_ec_Peter_Weyl_integer_charge_6Y': '6',
            'raw_lepton_current_is_not_used_as_hypercharge_normalization': True,
            'no_affine_outer_charged_ground_transfer': True}
