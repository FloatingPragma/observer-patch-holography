from fractions import Fraction as F
from . import source as s
from .certificate import calculate, full_fock_local_bound
from .gauss import fusion_certificate


def test_native_pair_is_positive_complete_and_pure_quark_zero():
    b, m = s.components('baryon'), s.components('meson')
    assert b['parts']['Pair4'].nnz == 0
    P = m['parts']['Pair4']
    assert (P@P-8*P).nnz == 0
    assert P.shape == (576,576)


def test_full_spectra_calibrations_failures_and_current_transfer():
    d = calculate()
    assert sum(x['full_magnetic_dimension'] for x in d['complete_nuclear_channels']) == 2520
    assert F(1,100) < F(d['native_axial_ground']['five_input_ground_projector_norm_distance_upper']) < F(1,20)
    assert not d['held_out_neutral_pion_control']['matches']
    assert all(F(x['singlet_preparation_success'])+F(x['coloured_preparation_failure']) == 1 for x in d['complete_preparation_and_spectral_records'])


def test_all_fock_counts_and_unbounded_electric_irrep_cutoff():
    d = calculate()['common_length_attachment']
    bound = full_fock_local_bound(F(d['kappa_required_interval'][1]))
    assert sum(x['dimension'] for x in bound['complete_count_sectors']) == 2**24
    assert F(bound['uniform_number_floor']) > 0
    fusion = fusion_certificate()
    assert fusion['uniform_open_chain_factor'] == '22/7'
    assert fusion['edge_sum_cutoff'] == '56/3'
