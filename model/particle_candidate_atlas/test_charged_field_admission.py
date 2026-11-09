"""Exact source/whole-Fock controls and scoped-admission negative checks."""
from copy import deepcopy
from fractions import Fraction as F
import json

import pytest
import sympy as s

if __package__:
    from . import charged_field_admission as a
else:
    import charged_field_admission as a
from model.standard_model_source_map import fermion_attachment as fermion


@pytest.fixture(scope='module')
def report():
    value=a.certificate()
    json.dumps(value,allow_nan=False)
    return value


@pytest.mark.parametrize('value',[
    s.zeros(3),s.diag(1,1,0),s.diag(1,1,-1),s.eye(2),
    s.Matrix([[1,1,0],[0,1,0],[0,0,1]]),s.eye(3)*1.0,
    s.diag(1,1,s.Symbol('x',positive=True)),s.diag(1,1,s.oo)])
def test_positive_matrix_admission_rejects_missing_premises(value):
    with pytest.raises(ValueError):a.positive_family(value)


def test_generic_complex_family_has_two_exact_singular_copies():
    u=s.Matrix([[s.Rational(3,5),4*s.I/5,0],[4*s.I/5,s.Rational(3,5),0],[0,0,1]])
    assert u.H*u==s.eye(3)
    y=u*s.diag(1,2,3)*u.H
    k=a.family_kernel(y)
    active=[i for i in range(45) if any(k[i,j]!=0 for j in range(45))]
    assert len(active)==6
    z=k.extract(active,active)
    t=s.Symbol('lambda')
    assert s.factor((z*z.H).charpoly(t).as_expr())==(t-1)**2*(t-4)**2*(t-9)**2
    assert k==k.T and a.clean(k-k.H)!=s.zeros(45)


def test_full_fock_pairing_norm_from_every_occupation():
    # One family has two independent spin-singlet pair blocks, on all16states.
    eps=s.Matrix([[0,1],[-1,0]])
    delta=s.kronecker_product(eps,s.Matrix([[0,3],[3,0]]))
    pair=s.zeros(16)
    for mask in range(16):
        for i in range(4):
            for j in range(i+1,4):
                first=fermion.creation_action(j,mask,4)
                if first is None:continue
                second=fermion.creation_action(i,first[0],4)
                if second is not None:pair[second[0],mask]+=delta[i,j]*first[1]*second[1]
    full=pair+pair.H
    values=full.eigenvals()
    assert sum(values.values())==16
    assert min(values)==-6 and max(values)==6
    assert values[0]==6


@pytest.mark.parametrize('bad',[True,0.25,float('nan')])
def test_controller_bounds_require_exact_supplied_values(bad):
    with pytest.raises(ValueError):a.physical_controller_budget(10,1,bad)


@pytest.mark.parametrize('args',[dict(steps=0,time=1,p_upper='1/12'),
    dict(steps=2,time=-1,p_upper='1/12'),dict(steps=2,time=1,p_upper=1),
    dict(steps=2,time=1,p_upper='1/12',cells=True),
    dict(steps=2,time=1,p_upper='1/12',gain_upper=-1)])
def test_controller_rejects_invalid_service(args):
    with pytest.raises(ValueError):a.physical_controller_budget(**args)


def test_complete_controller_bound_and_finite_volume_dependence():
    b=a.physical_controller_budget(10000,'3/5','1/12',cells=2,gain_upper='7/3',higgs_norm='2/5')
    p=1-F(11,12)**2;norm=2*2*F(7,3)*F(2,5);t=F(3,5)
    expected=(4*p*norm*t+4*(1+p)*(norm*t)**2)/(10000*(1-p))
    assert F(b['uniform_fixed_memory_bound'])==expected
    assert F(b['archive_distance_squared_upper'])==2*expected
    assert F(b['archive_distance_upper'])**2>=2*expected
    assert b['common_H0_norm_enters_bound'] is False
    assert b['archive_uses_same_record_law_and_initial_correlations'] is True
    assert b['native_attempt_clock_identified'] is False
    assert b['volume_uniform'] is False


def test_actual_complete_native_source_and_spectator_countercontrol(report):
    c=report['source_C_duty_control']
    assert c['fine_records_all_old_inputs']==1536
    assert c['unconditional_trace']=='1/12' and c['physical_pairing_norm']=='1/6'
    assert all(s.sympify(x).is_positive is True for x in c['exact_positive_principal_minors'])
    zero=report['spectators']
    assert (zero['massless_quark_internal_coefficients'],zero['massless_neutral_internal_coefficients'])==(36,3)
    control=zero['free_zero_momentum_control']
    assert control['electric_charge']=='-1' and control['particle_number']==3
    assert control['global_color_generators_annihilate']==8
    assert control['energy_above_a_free_ground_with_empty_spectator_zero_modes']=='0'
    assert control['local_Gauss_constraint_or_interacting_hadron_claim'] is False


def test_current_gated_defaults_remain_distinct_from_coherent_default(report):
    c=report['distinct_default_action_exclusions']
    assert c['coherent_default_is_current_gated_C_action'] is False
    assert c['literal_C_duty_muon_ratio_disjoint_from_reference'] is True
    assert c['distinct_coherent_default_ratios']==['1','1','34+24*sqrt(2)']
    assert [r['duty'] for r in c['current_gated_actions']]==['always','C','E','CE']
    assert all(r['two_thirds_excluded'] for r in c['current_gated_actions'])


def test_forged_parent_koide_scope_is_rejected():
    parent=json.loads((a.ROOT/'model/standard_model_source_map/particle_map.json').read_text())
    data=deepcopy(parent)
    entry=next(r for r in data['finite_source_attachments']['source_action_masses']['policies'][0]['channels'] if r['channel']=='charged_lepton')
    entry['Koide_diagnostic']['positive_root_Koide_ratio_interval']=['0','1']
    with pytest.raises(ValueError):a.default_action_exclusions(data)


def test_all18_targets_and_entire_regulator_spectrum_are_retained(report):
    assert len(report['admission_rows'])==18
    assert len(report['retained_valleys'])==8
    assert sorted(r['winding'] for r in report['retained_valleys'])==[-1]*4+[1]*4
    for row in report['retained_valleys']:
        derivative=s.diag(*[-s.cos(s.pi*b) for b in row['bits']])
        assert row['winding']==derivative.det()
    assert report['configured_family_count']==3
    assert report['families_selected_by_primitive_repair'] is False
    assert report['physical_particle_map']==report['physical_mass_predictions']==[]
    for row in report['admission_rows']:
        assert row['physical_particle_identified'] is False
        assert row['physical_mass'] is None and row['defect_link'] is None
        assert row['amber_unproved_attachments']
        assert all(set(r)=={'claim','reason'} for r in row['red_exactly_excluded_assertions'])


def test_missing_or_duplicated_target_is_not_hidden():
    p=json.loads((a.ROOT/'model/standard_model_source_map/particle_map.json').read_text())
    p['rows'][-1]=p['rows'][0]
    with pytest.raises(ValueError):a.admission_rows(p)
