"""Full CAR probe checks, exact local support and charged-window controls."""
import pytest
import sympy as s
from . import model,local,transfer
from ..model import sequence
from ..collective_spin_response.bounds import probability_interval


def test_actual_full_CAR_reference_vacuum_contraction_for_all_families():
    result=model.operator_certificate()
    assert result['full_reference_Fock_dimension']==2**30
    assert len(result['all_five_native_elementary_probe_controls'])==5
    # Independent 60-mode signed CAR column, with source and reference both graded.
    mask=sum(1<<i for i in range(30) if i%3!=1)
    for family,indices in model.FAMILIES.items():
        image={}
        occupied=0
        for spin in range(2):
            for a in indices:
                mode=15*spin+a;occupied+=(mask>>mode)&1
                step=sequence(mask,((30+mode,True),(mode,False)))
                if step is not None:image[step[0]]=step[1]/s.sqrt(2*len(indices))
        norm=sum(s.conjugate(value)*value for value in image.values())
        assert norm==s.Rational(occupied,2*len(indices))<=1
        assert all((target>>30).bit_count()==1 for target in image)
        assert all(target.bit_count()%2==mask.bit_count()%2 for target in image)


def test_complete_positive_local_support_and_compensated_ground_comparison():
    for family,count,blocks in (('ec',9,8),('L',17,9)):
        response=local.positive_response(family)
        assert response['q'].degree()==len(response['rows'])==count
        assert all(0<row['weight_interval'][0]<=row['weight_interval'][1]<s.Rational(1,2)
            for row in response['rows'])
        assert all(row['source_energy_difference_interval'][0]>0 for row in response['rows'])
        sector=local.compensated_blocks(family)
        assert len(sector['blocks'])==blocks
        assert sector['complete_compensated_local_gap_lower']>s.Rational(9,100)
    assert not local.compensated_blocks('ec')['global_compensated_ground_identified']
    assert local.compensated_blocks('L')['global_compensated_ground_identified']


def test_positive_resolvent_against_independent_exact_matrix_inverse():
    from ..charged_ground_exchange.resolvent import field_data,field_value
    from ..ground_composite_response.model import source
    e,f,field=field_data();h=source(1)['H'];x=s.Symbol('x')
    p=s.Poly(h.charpoly(x).as_expr(),x).sqf_part()
    for family in ('ec','L'):
        response=local.positive_response(family)
        pullback=sum((b['adjoint']*(32*s.eye(b['H'].rows)-b['H']).inv()*b['map']
            for b in response['blocks']),s.zeros(h.rows))
        actual=field.zero;power=s.eye(h.rows)
        for i in range(p.degree()):
            coefficient=s.Poly(sum(p.nth(j)*x**(j-i-1) for j in range(i+1,p.degree()+1)),x)
            actual+=field_value(coefficient)/field_value(p.diff())*field.convert(s.trace(power*pullback))
            power=h*power
        n=response['numerator'];q=response['denominator']
        numerator=sum((n.rep.nth(k)*field.convert(32**k) for k in range(n.degree()+1)),field.zero)
        denominator=sum((q.rep.nth(k)*field.convert(32**k) for k in range(q.degree()+1)),field.zero)
        assert actual==numerator/denominator


def test_full_window_bound_admits_unseen_intruder_eigenstates():
    # The measured bin has rank two although the atomic probe has one visible atom there.
    c=s.Rational(999999,1000001);t=s.Rational(2000,1000001)
    rotation=s.Matrix([[c,-t,0],[t,c,0],[0,0,1]])
    h0=s.diag(-3,1,s.Rational(-29,10));h=rotation*h0*rotation.T
    projector=rotation*s.diag(1,0,1)*rotation.T
    assert projector**2==projector and projector*h==h*projector and s.trace(projector)==2
    v=4*t
    assert (h-h0)**2==s.diag(v*v,v*v,0)
    trial=s.Matrix([1,1,0])/s.sqrt(2)
    atomic=s.Matrix([1,0,0])/s.sqrt(2)
    window,error,_=transfer.window_error(((-3,-3),(1,1)),0,s.Rational(1,5),v,0)
    difference=projector*trial-atomic
    assert (difference.T*difference)[0]<=error**2
    probability=(trial.T*projector*trial)[0]
    low,high=probability_interval((s.Rational(1,2),)*2,error)
    assert low<=probability<=high


def test_actual_interacting_odd_channel_and_physical_preparation_scope():
    for family in ('ec','L'):
        result=transfer.certificate(family)
        reference=result['reference_sector_comparison']
        assert s.Rational(reference['zero_reference_cost_full_lattice_one_particle_sector_gap_lower'])>s.Rational(4,100)
        assert s.Rational(reference['positive_cost_full_Fock_occupied_reference_sector_gap_lower'])>s.Rational(16,100)
        readout=result['full_interacting_energy_readout']
        assert readout['odd_source_parity_on_every_pump_success']
        for low,high in readout['conditional_probabilities'].values():assert 0<s.Rational(low)<=s.Rational(high)<1
        records=result['original_preparation_records']
        assert records['positive_coarse_record_count']==len(records['records'])==4
        assert not records['initial_failure_source_only_Gauss_constraints_asserted']
        physical=result['constructive_approximate_preparation_records']
        assert physical['all_preparation_success_and_failure_branches_obey_source_Gauss']
        assert physical['all_detector_success_and_failure_branches_obey_joint_source_static_reference_Gauss']
        assert not readout['isolated_charged_eigenvalue_or_unchanged_multiplicity_claimed']


def test_missing_support_wrong_family_and_unseparated_window_are_rejected(monkeypatch):
    for family in ('Q','electron','ec addition'):
        with pytest.raises(ValueError):local.positive_response(family)
    for args in ((((0,0),(1,1)),0,2,0,0), (((0,0),),1,1,0,0),
                 (((1,0),),0,1,0,0), (((0,0),),0,-1,0,0)):
        with pytest.raises(ValueError):transfer.window_error(*args)
    from . import certificate as custody
    monkeypatch.setattr(custody,'DEPENDENCIES',(*custody.DEPENDENCIES,'missing-odd-source.json'))
    with pytest.raises(FileNotFoundError,match='missing declared interacting charged'):
        custody.input_hashes()
