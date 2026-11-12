"""Independent boundary, charge, perturbation and projection controls."""
from itertools import product
import sympy as s
from . import source,threshold


def test_every_one_fermion_translation_square_and_boundary_difference():
    for axis in range(3):
        for mode in range(240):
            key=(1<<mode,(0,)*24,(0,)*8)
            for periodic,expected in ((False,-1),(True,1)):
                one,a=source.translate_key(key,axis,periodic=periodic)
                two,b=source.translate_key(one,axis,periodic=periodic)
                assert two==key and a*b==expected


def test_complete_column_preserves_actual_L_B_and_parity():
    def charges(key):
        mask=key[0]
        l=sum(source.L_WEIGHTS[i%15] for i in range(240) if mask>>i&1)
        b=sum(source.B_WEIGHTS[i%15] for i in range(240) if mask>>i&1)
        return l,b,mask.bit_count()%2
    key=((1<<9)|(1<<11)|(1<<26)|(1<<85),(0,)*24,(0,)*8)
    for periodic in (False,True):
        image=source.full_column(key,periodic=periodic)
        assert image
        assert all(charges(out)==charges(key) for out in image)


def test_K2_sign_cancellation_requires_complete_support_preservation():
    b=s.Matrix([[1,2],[0,1],[2,0],[1,1]])
    d=s.diag(2,3,4,5);sign=s.diag(-1,-1,1,1)
    original=-b.T*d.inv()*b
    assert -(sign*b).T*d.inv()*(sign*b)==original
    changed=d.copy();changed[1,2]=changed[2,1]=s.Rational(1,7)
    assert -(sign*b).T*changed.inv()*(sign*b)!=-b.T*changed.inv()*b


def test_periodic_group_projection_and_antiperiodic_absence():
    # Independent site wave functions, with all eight character effects.
    points=list(product((0,1),repeat=3));projectors=[]
    for momentum in points:
        vector=s.Matrix([(-1)**sum(a*b for a,b in zip(momentum,x)) for x in points])/s.sqrt(8)
        projectors.append(vector*vector.T)
    assert sum(projectors,s.zeros(8))==s.eye(8)
    localized=s.eye(8)[:,0]
    assert all((localized.T*p*localized)[0]==s.Rational(1,8) for p in projectors)
    t=s.Matrix([[0,-1],[1,0]])
    assert t*t==-s.eye(2) and (t-s.eye(2)).det()!=0


def test_positive_threshold_uses_actual_ground_and_calibration_is_input():
    result=threshold.certificate()
    assert result['ground_raw_L']==8 and result['ground_raw_B3']==0
    assert result['calibration']['scale_is_supplied_empirical_calibration']
    for row in result['sectors']:
        low,high=map(s.Rational,row['periodic_variant_zero_momentum_sector_threshold_interval'])
        assert 0<low<high<s.Rational(1,10)
        assert s.Rational(row['remaining_threshold_under_that_perturbation_lower'])==low/2
        assert row['decay_to_L8_ground_forbidden_by_source_L_conservation']
        assert not row['arbitrary_charge_carrying_environment_stability_claimed']


def test_sector_preserving_perturbation_budget_and_symmetry_breaking_control():
    gap=s.Rational(1,10);epsilon=gap/4
    h=s.diag(0,gap);l=s.diag(8,7)
    allowed=s.diag(epsilon,-epsilon)
    assert (h+allowed)[1,1]-(h+allowed)[0,0]==gap/2
    assert allowed*l==l*allowed
    forbidden=s.Matrix([[0,epsilon],[epsilon,0]])
    assert forbidden*l!=l*forbidden


def test_neutrality_is_actual_current_and_dark_control_retains_number():
    result=source.current_certificate()
    assert result['Q_commutator_exactly_zero']
    assert not result['sterile_coordinate_or_zero_charge_added']
    dark=source.dark_number_control()
    assert dark['actual_conserved_number_above_its_vacuum']==4
    assert not dark['number_is_identified_with_R15_electric_current']
