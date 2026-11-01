from copy import deepcopy

import pytest

from . import matrix, packet, identification


@pytest.mark.parametrize('field',[
    'all_560_eigenvalues_or_degeneracies_resolved',
    'gap_bound_applied_to_unobserved_sector',
    'infinite_volume_locality_or_gap_inferred_from_cube',
    'common_Lorentz_metric_or_physical_masses_identified',
    'native_measurement_executed'])
def test_unproved_scope_promotions_are_rejected(field):
    result=deepcopy(matrix.certificate())
    result[field]=True
    with pytest.raises(ValueError,match='scope'):
        packet.validate({'schema':'oph.matrix_spectral_response.v1','all_pass':True,
                         'matrix_response':result})


def test_prior_local_spectrum_authentication_retains_entire_source_domain():
    result=packet.prior_certificate()
    assert result['all_prior_response_source_hashes_verified'] is True
    assert result['prior_full_source_dimension']==560
    assert result['prior_local_spectrum_lines']==78


def test_replay_binds_its_proof_and_original_source():
    bindings=packet.source_bindings()
    for path in ('model/matrix_spectral_response/matrix.py',
                 'model/matrix_spectral_response/README.md',
                 'model/retarded_matter_response/results.json',
                 'model/particle_species_detector/physical_targets.json',
                 'paper/tex_fragments/MATRIX_SPECTRAL_RESPONSE.tex'):
        assert path in bindings and len(bindings[path])==64


@pytest.mark.parametrize('field',[None,'unit_factor_or_source_parameters_fitted',
    'finite_gaps_proved_to_be_relativistic_rest_masses',
    'unobserved_sectors_composites_or_continuum_excluded'])
def test_finite_dictionary_is_required_and_cannot_promote_scope(field):
    value={'schema':'oph.matrix_spectral_response.v1','all_pass':True,
           'matrix_response':matrix.certificate(),
           'finite_gap_dictionary':deepcopy(identification.certificate())}
    if field is None:del value['finite_gap_dictionary']
    else:value['finite_gap_dictionary'][field]=True
    with pytest.raises(ValueError):packet.validate(value)
