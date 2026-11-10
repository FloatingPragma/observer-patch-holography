"""All-column compensated swaps and their complete reduced records."""
import sympy as s
from .instrument import (SOURCE_MODES, branch_metadata, ground_preparation,
                         reduced_checks, swap_checks)


def test_all_four_actual_fermionic_swaps_and_complete_reductions():
    for mode in SOURCE_MODES:
        swap = swap_checks(mode)
        reduced = reduced_checks(mode)
        assert swap['full_fermion_columns'] == 16384
        assert swap['actual_annihilator_column_intertwiners'] == 14 * 16384
        assert swap['actual_local_Gauss_current_checks'] == 24
        assert swap['source_plus_apparatus_dimension'] == 2457600
        assert reduced['all_source_columns_in_completeness'] == 4096
        assert reduced['even_reference_branch_columns_checked'] == 8192
        assert reduced['all_source_parity_instrument_complete']
        site = mode // 6
        for setting, outcome, charge in ((1, 0, -1), (1, 1, 0), (0, 1, 1), (0, 0, 0)):
            branch = branch_metadata(mode, setting, outcome)
            assert branch['boundary_hypercharges'][site] == charge
            assert branch['boundary_hypercharges'][1 - site] == 0
            assert branch['source_delta_L'] == charge


def test_physical_ground_preparation_retains_original_failures_and_charges():
    preparation = ground_preparation()
    assert preparation['physical_dimension'] == 1053
    assert preparation['ground_projector_rank'] == 1
    assert s.Rational(preparation['original_ground_success_probability']) + s.Rational(preparation['original_preparation_failure_probability']) == 1
    assert preparation['ground_L_twoSz_parity'] == [2, 0, 'even']
    assert preparation['post_ground_ancilla_L_twoSz_parity'] == [-2, 0, 'even']
    assert not preparation['apparatus_work_calculated']
