"""Independent CAR, full-source energy and complete-record checks."""
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import json

import pytest
import sympy as s

from . import model as m


def _subtract(left,right,scale=1):
    out=dict(left);m.add(out,right,-scale);return m.clean_state(out)


def _literal_bond(vector,axis,link):
    terms=[];h=m.transport.spin_blocks()[axis]
    for a in range(4):
        for b in range(4):
            if not h[a,b]:continue
            for i in range(15):
                for j in range(15):
                    c=h[a,b]*link[i,j]
                    if c:terms.extend(((15*a+i,60+15*b+j,c),(60+15*b+j,15*a+i,s.conjugate(c))))
    return m.sectors._act(vector,terms)


def _literal_graded(delta,label,j,axis,link):
    polynomial=next(p for ell,k,p in m.irrep_polynomials(delta) if (ell,k)==(label,j))
    omega=m.spin.OMEGA|(m.spin.OMEGA<<60)
    vector=m.apply_polynomial({omega:s.S.One},polynomial,0)
    def comm(v):
        return _subtract(_literal_bond(m.apply_polynomial(v,polynomial,1),axis,link),
                         m.apply_polynomial(_literal_bond(v,axis,link),polynomial,1))
    return _subtract(m.apply_polynomial(comm(vector),polynomial,0,True),
                     comm(m.apply_polynomial(vector,polynomial,0,True)),(-1)**(sum(map(abs,delta))%2))


def _coefficient_graded(delta,label,j,link):
    polynomial=next(p for ell,k,p in m.irrep_polynomials(delta) if (ell,k)==(label,j))
    omega=m.spin.OMEGA|(m.spin.OMEGA<<60);vector=m.apply_polynomial({omega:s.S.One},polynomial,0)
    blank=m.apply_polynomial(vector,polynomial,0,True);answer={};sign=(-1)**(sum(map(abs,delta))%2)
    for (q,block),increment in zip(m.sectors.native_blocks(),delta):
        if not increment:continue
        for left,i in enumerate(block):
            for right,j in enumerate(block):
                z=s.conjugate(link[i,j]) if increment<0 else link[i,j]
                first=m.apply_polynomial(m.derivative_coefficient(vector,polynomial,q,left,right,increment<0),polynomial,0,True)
                second=m.derivative_coefficient(blank,polynomial,q,left,right,increment<0)
                m.add(answer,first,z);m.add(answer,second,-sign*z)
    return m.clean_state(answer)


def test_all_forty_literal_irrep_columns_and_gauge_currents():
    rows=m.calculated_rows()
    assert len(rows)==40
    assert min(row['contractive_probe_graded_witness_norm_squared'] for row in rows)==F(1,972)
    for _,base in m.identities.FAMILIES:
        for sign in (1,-1):
            delta=tuple(sign*x for x in base);local=m.spin.local_multiplet(delta)
            for gen in local['gauge_matrices']:assert gen*local['kernel_basis']==s.zeros(gen.rows,local['kernel_basis'].cols)
            for label,j,polynomial in m.irrep_polynomials(delta):
                assert sum(s.conjugate(v)*v for _,v in polynomial)==1
                assert all(m.identities.polynomial_delta(word)==delta for word,_ in polynomial)


def test_full_120_mode_Wilson_action_matches_derivative_on_nondiagonal_native_link():
    link=m.transport.nontrivial_link()
    assert link.H!=link and m.spin.clean(link.H*link)==s.eye(15)
    for family,(_,base) in enumerate(m.identities.FAMILIES):
        for sign in (1,-1):
            delta=tuple(sign*x for x in base);label,j,_=m.irrep_polynomials(delta)[0]
            actual=_literal_graded(delta,label,j,family,link)
            coefficient=_coefficient_graded(delta,label,j,link)
            assert actual and actual==coefficient


def test_Haar_norms_and_full_energy_means_include_both_source_signs():
    expected={
        ('weak_pair_charge_six','G1'):{F(-5513,8)},
        ('weak_pair_charge_six','H'):{F(-690)},
        ('color_triple','G1'):{F(-6739,10)},
        ('color_triple','H'):{F(-2021,3)},
        ('even_determinant_hole','E'):{F(-1353,2)},
        ('even_determinant_hole','T1'):{F(-1353,2),F(-14841,22)},
        ('even_determinant_hole','T2'):{F(-1353,2)}}
    seen={key:set() for key in expected}
    for row in m.calculated_rows():
        seen[row['family'],row['probe_irrep']].add(row['normalized_first_action_source_energy'])
        assert row['raw_graded_witness_norm_squared']==sum(q['Haar_norm_squared'] for q in row['source_Haar_coefficients'])
        weighted=sum(q['Haar_norm_squared']*F(105,2)*q['source_link_Casimir'] for q in row['source_Haar_coefficients'])
        assert row['normalized_first_action_source_energy']==-720+3*row['fixed_total_defect_count']+weighted/row['raw_graded_witness_norm_squared']
        assert row['output_defects_at_x_y']==[1,row['fixed_total_defect_count']-1]
    assert seen==expected
    for row in m.calculated_rows():
        other=next(r for r in m.calculated_rows() if r['family']==row['family'] and r['sign']==-row['sign'] and (r['probe_irrep'],r['component'])==(row['probe_irrep'],row['component']))
        for key in ('raw_graded_witness_norm_squared','raw_polynomial_monomials','normalized_first_action_source_energy'):
            assert row[key]==other[key]


def test_spectral_Markov_bound_and_clock_widening_need_sector_floor():
    plan={'energy_resolution':'1/64','clock_plan':{'complete_clock_tail_upper':'1/1024'},'one_complete_clock_error_upper':'1/1024'}
    value=m.spectral_window(-5,-4,-8,-7,plan)
    assert value['ideal_final_source_energy_interval']==[F(-5),F(-2)]
    assert value['first_action_reported_interval_weight_lower']==F(2,3)-F(1,512)
    # A three-atom positive measure checks the bound without assuming an eigenvector.
    energy=[F(-5),F(-4),F(2)];weight=[F(3,4),F(1,5),F(1,20)]
    mean=sum(e*p for e,p in zip(energy,weight));cut=3*mean+10
    assert sum(p for e,p in zip(energy,weight) if e<=cut)>=F(2,3)
    with pytest.raises(ValueError):m.spectral_window(-5,-4,-8,-5,plan)
    with pytest.raises(ValueError):m.spectral_window(-5,-6,-8,-7,plan)


def test_factored_original_contrast_controls_every_component_without_expanding_c():
    for R in (1,F(45,2),10**5):
        for gamma in (F(1),F(1,972),F(1,10000)):
            plan=m.finite_time_plan(R,582407133,gamma)
            time=plan['time']['rational_prefactor']
            assert time==gamma/(64*R**2)
            assert plan['nonneighbor_to_neighbor_bound_ratio_upper']<F(1,2)
            assert plan['neighbor_original_joint_event_probability_lower']['negative_exponent']==1747221399
            assert plan['nonneighbor_original_joint_event_probability_upper']['negative_exponent']==2329628532
            for c in (F(1),F(1,8),F(1,1024)):
                t=c*time;near=c*t*t*gamma/16;far=9*R**4*t**4
                assert far<near/2 and t<=1/R
                # Rational lower 1/(2 sqrt2)>1/3 verifies the amplitude margin.
                assert F(1,3)-3*R**2*t/gamma>=F(1,4)
    with pytest.raises(ValueError):m.finite_time_plan(1,0,0)


def test_complete_contraction_dilation_retains_noncommuting_source_reference():
    # Exact nonnormal K requires the two different square roots in the Julia dilation.
    K=s.Matrix([[0,s.Rational(3,5)],[0,0]])
    left=s.diag(s.Rational(4,5),1);right=s.diag(1,s.Rational(4,5))
    U=K.row_join(left).col_join(right.row_join(-K.H))
    assert U.H*U==s.eye(4) and K.H*K+right*right==s.eye(2)
    psi=s.Matrix([1,0,0,1])/s.sqrt(2);rho=s.eye(2)/2
    p_success=(psi.H*s.kronecker_product(K.H*K,s.eye(2))*psi)[0]
    p_failure=(psi.H*s.kronecker_product(right*right,s.eye(2))*psi)[0]
    assert p_success+p_failure==1
    assert K*rho*K.H+right*rho*right!=rho
    # Marginal positivity is sufficient for the scalar event, never joint-pure domination.
    assert p_success==s.trace(rho*K.H*K)
    assert (psi*psi.H).rank()==1


@lru_cache(maxsize=1)
def _actual():
    member=m.contract.member(0);return member,m.bind_member(member)


def test_actual_complete_energy_bands_source_ids_and_occupation_records():
    member,bound=_actual();result=bound['result']
    assert all(result[key]==value for key,value in member['ids'].items())
    assert result['physical_density_floor_factored']['negative_exponent']==582407133
    assert len(result['component_rows'])==40
    assert result['common_original_spatial_energy_contrast']['common_witness_lower']=='1/972'
    labels=set()
    for row in result['component_rows']:
        event=row['complete_final_energy_event'];labels.add(tuple(row['actual_final_native_occupation']))
        assert F(event['reported_final_excitation_interval'][0])>0
        assert F(event['first_action_reported_interval_weight_lower'])>F(1,2)
        assert not event['energy_label_is_transfer'] and not event['initial_energy_conditioning_used']
    assert len(labels)==6
    menu=result['complete_instrument']
    assert menu['net_native_increment']==[0]*5 and menu['CAR_parity']==0
    assert menu['all_failed_interventions_other_sectors_and_energy_records_retained']
    assert not menu['original_input_density_replaced'] and not menu['nonselective_output_density_claimed_unchanged']
    assert result['mass_predictions_GeV']==[] and result['physical_particle_map']==[]


def test_mutation_of_source_energy_normalization_or_complete_event_is_rejected():
    member,bound=_actual();m.verify_bound(member,bound)
    for path in (('component_rows',0,'raw_polynomial_monomials'),('component_rows',0,'normalized_first_action_source_energy')):
        mutated=deepcopy(bound);mutated['result'][path[0]][path[1]][path[2]]='0'
        with pytest.raises(ValueError):m.verify_bound(member,mutated)
    mutated=deepcopy(bound);mutated['result']['complete_instrument']['original_input_density_replaced']=True
    with pytest.raises(ValueError):m.verify_bound(member,mutated)
    mutated=deepcopy(member);mutated['ids']['clock_id']='another-clock'
    with pytest.raises(ValueError):m.bind_member(mutated)


def test_canonical_receipt_serialization_and_scope_do_not_promote_a_mass():
    _,bound=_actual();result=bound['result']
    assert json.loads(json.dumps(bound,sort_keys=True))==bound
    assert result['nonzero_original_joint_spatial_final_energy_contrast_proved']
    assert not result['original_odd_neighboring_two_point_m1_value_or_sign_computed']
    assert not result['physical_neutrino_or_neutron_identified']
    assert not result['isolated_line_or_relativistic_dispersion_proved']
