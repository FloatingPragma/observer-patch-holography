"""Independent source columns, complete cut laws and broken attachments."""
import json
import random
import pytest
import sympy as s
from . import algebra,records
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
from model.coherent_defect_wave.chiral_joint.peter_weyl_links import model as pw
from model.coherent_defect_wave.chiral_joint.periodic_family.model import local_yukawa


def test_full_actual_gauss_compiler_and_all_cut_identities():
    census=algebra.gauss_cut_census()
    assert census['independent_actual_Gauss_column_factor_basis_inputs']==3105
    rng=random.Random(918)
    y=representation()[2][11];h=higgs_memory()[1][11]
    for _ in range(12):
        key=(rng.getrandbits(240),tuple(rng.randrange(119) for _ in range(24)),
             tuple(rng.randrange(5) for _ in range(8)))
        density=[];gauss=[]
        for cell in range(8):
            density.append(sum(y[mode%15,mode%15] for mode in range(30)
                if key[0]>>(30*cell+mode)&1)+h[key[2][cell],key[2][cell]])
            actual=pw.gauss_column(key,cell,11)
            assert set(actual)<={key}
            gauss.append(actual.get(key,0))
        for mask in range(256):
            selected=[cell for cell in range(8) if mask>>cell&1]
            assert sum(gauss[x] for x in selected)==sum(density[x] for x in selected)-s.Rational(algebra.cut_flux_six(mask,key[1]),6)


def test_six_actual_routes_fix_unique_flux_with_reference_sign():
    result=algebra.route_certificate()
    assert result['all_route_cut_identities']==1536
    patterns={tuple(row['three_half_cube_fluxes']) for row in result['actual_routes']}
    assert patterns=={('-1','0','0'),('0','-1','0'),('0','0','-1')}
    for row in result['actual_routes']:
        endpoint=row['origin_endpoint'];i=row['forced_full_link_basis_index']
        vector=s.eye(119)[:,i]
        assert (algebra.current(11,endpoint)+s.eye(119))*vector==s.zeros(119,1)
        # A wrong reference sign or wrong outward-flux sign cannot match
        # this same actual route. Defining another source is not a repair.
        assert (algebra.current(11,endpoint)-s.eye(119))*vector!=s.zeros(119,1)
        assert algebra.current(11,1-endpoint)[i,i]==1
        assert row['arrival_enclosing_outward_Gauss_flux']=='-1'


def test_complete_boundary_spectrum_has_no_dropped_link_states():
    values=[int(x) for x in algebra.charge_weights(0)]
    independent={}
    for first in values:
        for second in values:independent[-first-second]=independent.get(-first-second,0)+1
    assert dict(algebra.boundary_spectrum(2))==independent
    for count in (0,1,6,12,24):
        result=algebra.boundary_spectrum(count)
        assert sum(result.values())==119**count
        assert max(result)==6*count and min(result)==-6*count
    for invalid in (-1,25,1.0):
        with pytest.raises(ValueError):algebra.boundary_spectrum(invalid)


def test_scalar_and_weak_gauge_current_omission_controls():
    y=representation()[2][11];ys=higgs_memory()[1][11]
    image=local_yukawa(0,0,(s.Integer(1),)*3)
    assert image
    missing_scalar=[]
    for (mask,scalar),amplitude in image.items():
        charge=sum(y[mode%15,mode%15] for mode in range(30) if mask>>mode&1)
        assert charge+ys[scalar,scalar]==0 and amplitude!=0
        missing_scalar.append(charge)
    assert any(value!=0 for value in missing_scalar)
    full_q=sum((algebra.current(a,end) for a in (10,11) for end in (0,1)),s.zeros(119))
    assert full_q!=s.zeros(119)
    assert algebra.current(11,0)+algebra.current(11,1)==s.zeros(119)


def test_joint_records_keep_actual_probabilities_and_every_failure():
    value=records.certificate()
    rows=value['complete_original_preparation_records']['records']
    assert len(rows)==10
    arrivals=[row for row in rows if 'joint_charge_and_flux_atom' in row]
    assert len(arrivals)==6
    for row in arrivals:
        lo,hi=map(s.Rational,row['probability_interval'])
        assert s.Rational('6.430004784e-38')<lo<hi<s.Rational('1.274554272e-37')
        assert row['joint_atom_probability_equals_actual_arrival_probability']
        assert row['joint_charge_and_flux_atom']['physical_charge_flux_mismatch_probability']=='0'
    assert all('joint_flux_refinement' in row for row in rows if row not in arrivals)
    assert not value['postprocessing_scope']['energy_and_arrival_projectors_assumed_commuting']
    assert not value['physical_scope']['all_failure_flux_probabilities_numerically_evaluated']


def test_custody_rejects_missing_source_and_changed_data(monkeypatch):
    from . import certificate as custody
    value={'input_sha256':{'model/openwave_native_charge_flux/algebra.py':'0'*64}}
    with pytest.raises(ValueError,match='changed'):custody.validate_pins(value)
    monkeypatch.setattr(custody,'DEPENDENCIES',(*custody.DEPENDENCIES,'missing-native-charge-flux.json'))
    with pytest.raises(FileNotFoundError,match='missing native charge/flux'):custody.input_hashes()
