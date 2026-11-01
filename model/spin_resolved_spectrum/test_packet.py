from functools import lru_cache
import copy

import pytest

from . import packet as p


@lru_cache(maxsize=1)
def sample():
    return {'schema':'oph.spin_resolved_spectrum.v1','all_pass':True,
            'source_spin':p.spin.certificate(),'resolved_spectrum':p.spectrum.certificate()}


def test_complete_packet_and_new_sources_are_bound():
    p.validate(sample())
    bindings=p.source_bindings()
    for name in ('spin.py','spectrum.py','test_spin.py','test_spectrum.py','test_packet.py'):
        assert 'model/spin_resolved_spectrum/'+name in bindings
    assert 'paper/tex_fragments/SPIN_RESOLVED_SPECTRUM.tex' in bindings


@pytest.mark.parametrize('key',[
    'dark_energy_eigenvalues_or_full_energy_degeneracies_computed',
    'arbitrary_ground_selection_adds_frequencies',
    'native_full_quintet_preparation_or_79_line_resolution_executed',
    'total_hidden_spin_identified_with_particle_Lorentz_spin',
    'physical_particle_species_or_relativistic_rest_masses_identified',
])
def test_unsupported_promotions_rejected(key):
    value=copy.deepcopy(sample());value['resolved_spectrum'][key]=True
    with pytest.raises(ValueError):p.validate(value)


def test_ground_completeness_and_dictionary_boundaries_are_mandatory():
    for path,new in ((('source_spin','ground_family','creator_columns_checked'),16),
                     (('resolved_spectrum','finite_gap_dictionary','eligible_physical_comparisons_added'),1),
                     (('resolved_spectrum','finite_gap_dictionary','unseen_composites_other_particle_numbers_or_continuum_excluded'),True)):
        value=copy.deepcopy(sample());target=value
        for key in path[:-1]:target=target[key]
        target[path[-1]]=new
        with pytest.raises(ValueError):p.validate(value)
