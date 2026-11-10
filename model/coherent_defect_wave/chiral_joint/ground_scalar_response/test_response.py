"""Full native scalar identity, complete poles and hostile source controls."""
import numpy as np
import pytest
import scipy.linalg
import sympy as s
from . import model,spectrum,apparatus
from ..family_cell.quotient import sqrt_interval


def test_full_native_scalar_star_current_and_complete_zero_records():
    source=model.source_identities()
    assert source['full_native_anticommutator_N_Y_equals_Y']
    assert source['full_native_scalar_current_squared_equals_Y_squared']
    assert not source['individual_Yukawa_channels_substituted']
    result=spectrum.certificate()
    assert result['complete_joint_records']==22 and result['positive_probability_records']==6
    assert result['exact_zero_probability_records']==16
    assert not result['retarded_zero_frequency_pole']
    assert not result['transition_frequency_identified_with_relativistic_rest_mass']
    for setting in result['settings']:
        records=setting['complete_energy_records']
        assert len(records)==9 and sum(row['complete_spectral_rank'] for row in records)==32
        assert sum(row['role']=='zero' for row in records)==7
        assert [row['complete_spectral_rank'] for row in records if row['role']!='zero']==[2,2]
        assert records[0]['excitation_energy_interval']==['0','0']
        dark=next(row for row in records if row['target_energy_interval']==['1','1'])
        assert dark['complete_spectral_rank']==16 and dark['joint_probability_interval']==['0','0']


def test_independent_generalized_eigenvector_scalar_weights():
    data=model.matrices();g=np.array(data['G'],dtype=float)
    energies,vectors=scipy.linalg.eigh(np.array(data['G']*data['H'],dtype=float),g)
    assert abs(energies[0]-energies[1])<2e-10
    for setting,name,normalization in zip(spectrum.certificate()['settings'],('N','Y'),(16,1024)):
        amplitude=vectors.T@(g@(np.array(data[name],dtype=float)@vectors[:,:2]))
        for row in setting['complete_energy_records']:
            center=float(sum(map(s.Rational,row['target_energy_interval']))/2)
            ids=np.flatnonzero(np.abs(energies-center)<2e-10)
            assert len(ids)==row['complete_spectral_rank']
            probability=np.sum(amplitude[ids,:]**2)/normalization
            low,high=map(lambda value:float(s.Rational(value)),row['joint_probability_interval'])
            assert low-2e-10<=probability<=high+2e-10


def test_false_scalar_number_and_penalty_rejected():
    data=model.matrices()
    with pytest.raises(ValueError,match='must be a projector'):
        model.validate_scalar_star(data['H'],2*data['N'])
    with pytest.raises(ValueError,match='Yukawa anticommutator'):
        model.validate_scalar_star(data['H']+data['N'],data['N'])


def test_full_square_root_series_bound_and_even_apparatus():
    ratio=s.Rational(165,256);exact=sqrt_interval((1-ratio,1-ratio),bits=256)
    for order in (0,16,64,128):
        coefficients,tail=apparatus.square_root_series(order)
        value=sum(coefficient*ratio**k for k,coefficient in enumerate(coefficients))
        assert value>=exact[1] and value-exact[0]<=tail
    with pytest.raises(ValueError,match='nonnegative exact'):
        apparatus.square_root_series(-1)
    with pytest.raises(ValueError,match='nonnegative exact'):
        apparatus.square_root_series(1.0)
    result=apparatus.certificate()
    assert result['Y_failure_is_defined_on_full_carrier']
    assert result['both_probe_operators_are_even']
    assert all(not row['finite_polynomial_claimed_exact_unitary'] for row in result['uniform_full_carrier_approximation_bounds'])


def test_missing_declared_and_transitive_inputs_rejected(monkeypatch):
    from . import certificate as custody
    with monkeypatch.context() as patch:
        patch.setattr(custody.action,'SOURCE_FILES',(*custody.action.SOURCE_FILES,'missing-ground-scalar.json'))
        with pytest.raises(FileNotFoundError,match='missing declared'):
            custody.input_hashes()
    monkeypatch.setattr(custody,'import_closure',lambda seeds:set(seeds)|{'missing-ground-scalar-import.py'})
    with pytest.raises(FileNotFoundError,match='missing ground scalar import'):
        custody.input_hashes()
