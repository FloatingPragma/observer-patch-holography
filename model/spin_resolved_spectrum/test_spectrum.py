import copy
from fractions import Fraction as F
import json

import pytest

from . import spectrum as s


def test_complete_spin_resolved_lines_and_probability_normalization():
    rows=s.resolved_rows()
    high=[row for row in rows if row['total_hidden_spin']=='5/2']
    assert len(rows)==78 and len(high)==1
    assert high[0]['mixed_reference_local_addition_weight']==['3/10','3/10']
    assert high[0]['conditional_on_addition_local_weight']==['2/5','2/5']
    assert sum(row['reachable_N5_dimension'] for row in rows)==1256
    for key,target in (('mixed_reference_local_addition_weight',F(3,4)),
                       ('conditional_on_addition_local_weight',F(1))):
        assert sum(F(row[key][0]) for row in rows)<=target<=sum(F(row[key][1]) for row in rows)
        assert all(F(row[key][0])>0 for row in rows)


def test_lambda_identification_uses_exact_interval_and_band_not_decimal_labels():
    rows=copy.deepcopy(s.prior_spectrum()['matrix_response']['matrix_residue_rows'])
    assert s.lambda_row(rows)==2
    rows[1]['band']=-1
    with pytest.raises(ValueError):s.lambda_row(rows)
    rows[1]['band']=1;rows.append(rows[1])
    with pytest.raises(ValueError):s.lambda_row(rows)


def test_actual_physical_frame_retains_intersite_spin_mixing():
    value=s.frame_certificate()
    assert value['exact_Gaussian_integer_unitarity_and_kinetic_conjugacy']
    assert value['physical_intersite_blocks_with_spin_mixing']>0
    assert value['physical_local_response']=='g_xx(z) I_spin'


def test_prior_source_digest_cannot_be_changed(tmp_path,monkeypatch):
    prior=s.prior_spectrum();prior['source_bindings'][next(iter(prior['source_bindings']))]='0'*64
    path=tmp_path/'bad.json';path.write_text(json.dumps(prior))
    monkeypatch.setattr(s,'PRIOR',str(path))
    with pytest.raises(ValueError,match='source changed'):s.prior_spectrum()


def test_dark_dimensions_and_mixed_GNS_dimension_are_not_conflated():
    value=s.certificate()
    assert sum(row['reachable_dimension'] for row in value['N5_decomposition'])==1256
    assert sum(row['dark_dimension'] for row in value['N5_decomposition'])==3112
    assert value['mixed_retarded_minimal_cyclic_dimension']==632
    assert value['removal']['reachable_N3_dimension']==16
    assert value['dark_energy_eigenvalues_or_full_energy_degeneracies_computed'] is False
    assert value['physical_particle_species_or_relativistic_rest_masses_identified'] is False
