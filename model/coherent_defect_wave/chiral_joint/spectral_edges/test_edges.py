"""The finite edge premise comes from spectrum and original response weight."""
from copy import deepcopy
from fractions import Fraction as F
import pytest
from .model import source_receipt,attach_branch,certificate_data


def test_both_source_edges_instantiate_the_existing_original_weight_theorem():
    result=certificate_data()
    for row in result['actual_finite_source_edges']:
        lo,hi=map(F,row['first_atom_source_energy_interval'])
        bound=row['existing_edge_theorem_conclusion']
        assert F(bound['threshold_lower_strict'])<lo<hi==F(bound['threshold_upper'])
        a=F(row['source_tail_premise']['amplitude'])
        assert 0<a<F(1,1053)
        assert not bound['conditional_success_renormalization_used']
    assert not result['one_isolated_pole_required_for_every_physical_species']


def test_unresolved_cluster_weight_cannot_certify_the_first_atom():
    source=deepcopy(source_receipt())
    source['conditional_readout']['records'][0]['complete_energy_windows'][0]['spectral_rank']=2
    with pytest.raises(ValueError,match='first-atom weight'):attach_branch(source,'add')


def test_an_empty_observed_bin_is_not_a_lower_threshold_proof():
    source=deepcopy(source_receipt())
    source['conditional_readout']['records'][0]['unreachable_conserved_blocks_have_exact_zero_effect_on_prepared_ground']=False
    with pytest.raises(ValueError,match='uncontrolled lower'):attach_branch(source,'add')
