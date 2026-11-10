"""Actual compiled flux routes, independent cube moments and full physical law."""
import pytest
import sympy as s
from . import algebra,transport
from ..interacting_charged_response.model import FAMILIES
from ..periodic_family.model import Configuration,kinetic_terms
from ..peter_weyl_links import model as pw


def test_every_actual_incident_kinetic_route_creates_the_retained_flux():
    terms={term.path[0]:term for term in kinetic_terms(Configuration())}
    for route in algebra.routes():
        for family in ('ec','L'):
            index=FAMILIES[family][0]
            # An independent full-source compiler column, with a particle at
            # the arrival cell and a hole at origin, tests the refilling branch.
            key=(1<<(30*route['arrival']+index),(0,)*24,(0,)*8)
            image=pw.hopping_column(key,(terms[route['edge']],))
            assert image and sum(s.conjugate(v)*v for v in image.values())==s.Rational(1,16)
            for (mask,links,scalars),value in image.items():
                assert mask.bit_count()==1 and mask.bit_length()<=30
                assert (mask.bit_length()-1)%15 in FAMILIES[family]
                assert [i for i,x in enumerate(links) if x]==[route['edge']]
                assert scalars==(0,)*8
    result=algebra.certificate()
    assert not result['reference_moved_to_arrival_cell']
    assert result['complete_three_axis_spin_matrix_identities']==3


def test_independent_256_state_cube_projection_of_the_arrival_effect():
    from ..atomic_spin_geometry.cube import full_operators,magnetization_block,pauli
    _,_,unitary=full_operators();basis,h=magnetization_block(4)
    w=s.Matrix.hstack(*(h**3+20*h**2-16*h-576*s.eye(70)).nullspace())
    gram=w.T*w;small=gram.inv()*w.T*h*w
    x=s.Symbol('x');f=s.Poly(x**3+20*x**2-16*x-576,x)
    projector=(small**2+(20+x)*small+(x*x+20*x-16)*s.eye(3))/f.diff().as_expr()
    b,effects=algebra.spin_effects()
    sigma=(s.eye(2),s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))
    for axis,site in enumerate((4,2,1)):
        effect=effects[axis].subs(b,s.Rational(-1,8))
        full=s.SparseMatrix(256,256,{})
        for left in range(4):
            for right in range(4):
                coefficient=s.trace(effect*s.kronecker_product(sigma[left],sigma[right]))/4
                if coefficient:
                    a=s.eye(256) if left==0 else pauli(0,left-1)
                    c=s.eye(256) if right==0 else pauli(site,right-1)
                    full+=coefficient*a*c
        rotated=(unitary.H*full*unitary).extract(basis,basis)
        actual=s.trace(projector*gram.inv()*w.T*rotated*w)
        beta=s.Rational(-1,8)
        expected=(1+3*beta**2)/128+beta**2*(1-beta)*x/768
        numerator=s.together(actual-expected).as_numer_denom()[0]
        assert s.Poly(numerator,x).rem(f).is_zero


def test_preparation_error_preserves_strict_unconditional_small_arrival():
    for family in ('ec','L'):
        result=transport.certificate(family)
        controls=result['full_H_time_and_preparation_control']
        low,high=map(s.Rational,controls['per_route_arrival_interval_from_original_physical_vacuum_input'])
        error=s.Rational(controls['retained_environment_purification_vector_error_upper'])
        assert 0<low<high<error**2
        assert s.Rational(controls['strict_scaled_amplitude_margin'])>0
        assert not controls['preparation_error_subtracted_as_global_probability_error']
        records=result['complete_original_physical_records']
        assert len(records['records'])==records['coarse_record_count']==10
        assert len([row for row in records['records'] if row['record']=='charged_arrival_with_selected_flux'])==6
        assert all(row['static_reference_cell']==0 for row in records['records'] if 'static_reference_cell' in row)
        assert result['same_preparation_energy_setting']['energy_and_arrival_are_distinct_observer_settings']
        assert not result['same_preparation_energy_setting']['joint_commuting_energy_location_measurement_or_energy_resolved_arrival_inferred']
    q,a=s.symbols('q a');p=s.symbols('p:6')
    assert s.expand((1-q)+q*(1-a)+q*a/2+q*sum(p)+q*(a/2-sum(p)))==1


def test_exact_idle_kernel_and_invalid_transfer_controls():
    valid=(s.Rational(1,100),s.Rational(1,90))
    idle=transport.arrival_bounds(valid,s.Rational(1,2**34),0,15,s.Rational(1,1000),360,576)
    assert idle['probability_interval']==(0,0)
    overwhelmed=transport.arrival_bounds(valid,s.Rational(1,2**34),s.Rational(1,2**24),15,0,360,576,1)
    assert overwhelmed['probability_interval'][0]==0 and overwhelmed['strict_scaled_amplitude_margin']<0
    for args in ((valid,0,1,15,0,360,576), (valid,1,-1,15,0,360,576),
                 (valid,1,1,15,0,1,576), ((0,1),1,1,15,0,360,576),
                 (valid,0.1,1,15,0,360,576)):
        with pytest.raises(ValueError):transport.arrival_bounds(*args)


def test_missing_source_cannot_be_pinned(monkeypatch):
    from . import certificate as custody
    monkeypatch.setattr(custody,'DEPENDENCIES',(*custody.DEPENDENCIES,'missing-charge-route.json'))
    with pytest.raises(FileNotFoundError,match='missing declared charged transport'):
        custody.input_hashes()
