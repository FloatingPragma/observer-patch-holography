from dataclasses import replace
from fractions import Fraction as F
import numpy as np
import pytest
import sympy as s
from . import observables as o


@pytest.fixture(scope='module')
def d():
    return o.data()


@pytest.fixture(scope='module')
def compiled():
    arrays=o.source_arrays()
    return arrays,o.compile_certificate(arrays)


def halfdistance(left,right):
    assert left.keys()==right.keys()
    return sum(np.linalg.norm(left[k]-right[k],'nuc') for k in left)/2


def reference_density(ref=2):
    rng=np.random.default_rng(291)
    v=rng.normal(size=66*ref)+1j*rng.normal(size=66*ref)
    v/=np.linalg.norm(v)
    return np.outer(v,v.conj())


def test_complete_actual_flags_and_full_carrier_identities(d):
    result=o.instrument_certificate(d)
    assert result['raw_slots_per_axis']==25
    assert result['zero_native_slots_per_axis']==19
    assert result['exact_full_carrier_fine_branch_identities']==18
    assert sum(r.weight for r in d['branches'] if r.pairs and r.contrast==0)==1
    assert {r.weight for r in d['branches'] if r.pairs}=={F(8,33),F(25,33)}


def test_every_zero_sector_matrix_unit_is_preserved(d):
    # In particular the empty/double endpoint coherence is not secretly read.
    zero=[i for i in range(66) if d['zero'][i,i]]
    assert len(zero)==46
    i=next(k for k in zero if d['states'][k][1]==d['states'][k][7]==1)
    j=next(k for k in zero if d['states'][k][1]==d['states'][k][7]==0)
    assert i!=j
    for row in d['branches']:
        if row.pairs and row.contrast==0:
            V=o.branch_operator(d,row,'Z')
            assert V==d['zero']
            assert V[i,i]*V[j,j]==1


@pytest.mark.parametrize('selected',o.AXES)
def test_entangled_reference_complete_lueders_maps(d,selected):
    rho=reference_density()
    result=o.ideal_complete_output(d,rho,selected,corrected=True)
    for k,row in enumerate(d['branches']):
        if row.pairs:
            P=np.kron(np.asarray(d['projectors'][selected][row.contrast],complex),np.eye(2))
            expected=float(row.weight)*P@rho@P
        else:
            expected=np.zeros_like(rho)
        assert np.linalg.norm(result[k]-expected,'nuc')<2e-14
    assert abs(sum(np.trace(x) for x in result.values())-1)<2e-14


@pytest.mark.parametrize('selected',o.AXES)
def test_corrected_repeatability_and_raw_destructivity(d,selected):
    i,j,*_=d['pairs'][0]
    P=np.asarray(d['projectors'][selected][1],complex)
    v=P[:,i];v/=np.linalg.norm(v);rho=np.outer(v,v.conj())
    restored=sum(o.ideal_complete_output(d,rho,selected,corrected=True).values())
    assert np.linalg.norm(restored-rho,'nuc')<2e-14
    raw=sum(o.ideal_complete_output(d,rho,selected).values())
    assert np.linalg.norm(raw-rho,'nuc')>0.1


def test_coin_mutation_rejected(d):
    bad=dict(d)
    bad['branches']=tuple(replace(row,weight=row.weight+F(1,100)) if row.pairs else row for row in d['branches'])
    with pytest.raises(ArithmeticError,match='branch labels'):
        o.instrument_certificate(bad)


def test_spectator_dependent_correction_rejected(d):
    bad=dict(d);S=s.MutableSparseMatrix(d['swap'])
    i,j,*_=d['pairs'][0];S[i,j]*=-1
    bad['swap']=s.ImmutableSparseMatrix(S)
    with pytest.raises(ArithmeticError,match='correction'):
        o.instrument_certificate(bad)


def test_zero_sector_refinement_rejected(d):
    bad=dict(d);rows=list(d['branches'])
    k=next(k for k,row in enumerate(rows) if row.pairs and row.contrast==0)
    rows[k]=replace(rows[k],pairs=rows[k].pairs[:-1]);bad['branches']=tuple(rows)
    with pytest.raises(ArithmeticError,match='branch labels'):
        o.instrument_certificate(bad)


def test_y_orientation_mutation_rejected(d):
    bad=dict(d);bad['Y']=-d['Y']
    with pytest.raises(ArithmeticError,match='commutator'):
        o.instrument_certificate(bad)


def test_stationary_robertson_is_not_zero_or_conditioned(d):
    result=o.uncertainty_certificate(d)
    assert result['Robertson_lower']=='2719201/56325025'
    assert result['Robertson_margin']=='301088/2253001'
    assert result['active_probability']=='3201/7505'


def test_same_modular_clock_sign_on_full_carrier(d):
    weights=np.asarray(d['pi'].diagonal(),float).ravel()
    X,Y=(np.asarray(d[name],complex) for name in ('X','Y'))
    for parameter in (0.5,1):
        U=np.diag(np.exp(-1j*parameter*np.log(weights)))
        angle=parameter*np.log(25/8)
        assert np.linalg.norm(U.conj().T@X@U-(np.cos(angle)*X-np.sin(angle)*Y))<1e-14
        assert np.linalg.norm(U.conj().T@Y@U-(np.sin(angle)*X+np.cos(angle)*Y))<1e-14
    assert not o.modular_clock_certificate(d)['physical_time_identified']


def test_preparation_archive_boundary():
    result=o.preparation_certificate()
    assert result['half_diamond_upper']=='542159/1000000000000'
    assert not result['joint_preparation_archive_included']


@pytest.mark.parametrize('name',('X','X_inverse','Y','Y_inverse'))
def test_full512_compiled_word_resources_and_bound(compiled,name):
    arrays,cert=compiled
    row=cert['words'][name]
    assert arrays[name].shape==(512,512)
    assert F(row['operator_error_upper'])<F(1,10**7)
    assert F(row['returned_array_operator_error_upper'])==o.cap.ROTATION_GATE
    assert row['resources']['coherent_work_qubits']==9
    assert row['resources']['native_pulses']>10000
    assert F(row['first_processor_failure_mass_upper'])>0


@pytest.mark.parametrize('selected',o.AXES)
@pytest.mark.parametrize('corrected',(False,True))
def test_literal_complete_outputs_on_entangled_reference(d,compiled,selected,corrected):
    arrays,cert=compiled;rho=reference_density()
    actual=o.compiled_output(d,rho,selected,arrays,corrected=corrected)
    ideal=o.ideal_output_with_same_exports(d,rho,selected,corrected=corrected)
    row=o.protocol_certificate(cert)[selected]['corrected_worst_branch' if corrected else 'raw']
    assert halfdistance(actual,ideal)<=float(F(row['complete_halfdiamond_upper']))+1e-14
    assert len([key for key in actual if key[0]=='native'])==25
    assert abs(sum(np.trace(block).real for block in actual.values())-1)<=2*float(F(row['complete_halfdiamond_upper']))+1e-14


def test_padding_escape_is_not_discarded(d,compiled):
    arrays,_=compiled;mutated={**arrays,'X':arrays['X'].copy()}
    codes=set(o.clock.common.full.support()[0]);outside=next(i for i in range(512) if i not in codes)
    source=next(iter(codes));mutated['X'][outside,source]+=0.1
    rho=np.eye(66,dtype=complex)/66
    actual=o.compiled_output(d,rho,'X',mutated)
    assert np.trace(actual[('analyzer_export_abort',)]).real>1e-5
    with pytest.raises(ArithmeticError,match='full512'):
        o.compile_certificate(mutated)


def test_literal_inverse_word_is_executed(compiled):
    arrays,_=compiled
    assert np.linalg.norm(arrays['Y_inverse']-arrays['Y'].conj().T)<1e-10
    assert len(o.words()['Y_inverse'])==len(o.words()['Y'])


def test_compiled_noncommuting_effect_margin(compiled):
    _,cert=compiled
    assert F(cert['raw_noncommuting_effect_commutator_norm_lower'])>F(49,100)


@pytest.mark.parametrize('invalid',(None,0,True,'x','XX'))
def test_axis_validation(invalid):
    with pytest.raises(ValueError):
        o.axis(invalid)


def test_product_arithmetic_has_real_resource_domain():
    assert F(o.multiplication_error(10416)['operator_error_upper'])<F(1,10**7)
    for invalid in (True,0,12001):
        with pytest.raises(ValueError):
            o.multiplication_error(invalid)
    with pytest.raises(ValueError):
        o.multiplication_error(10416,11)


def test_dropped_zero_address_is_rejected(d):
    bad=dict(d)
    k=next(k for k,row in enumerate(d['branches']) if not row.pairs)
    bad['branches']=d['branches'][:k]+d['branches'][k+1:]
    with pytest.raises(ArithmeticError,match='branch labels'):
        o.instrument_certificate(bad)


def test_empty_reference_carrier_rejected(d):
    with pytest.raises(ValueError,match='nonempty'):
        o.native_output(d['branches'][0],np.zeros((0,0)))
