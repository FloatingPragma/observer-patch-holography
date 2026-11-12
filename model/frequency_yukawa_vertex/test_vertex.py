import numpy as np
from scipy import linalg as la
import sympy as s
from . import vertex as v


def test_complete_resolvent_matches_independent_dense_inverse():
    h=np.asarray(v.data()['H'],complex);matrix=1j*np.eye(320)-h-2*np.eye(320)
    assert la.norm(np.asarray(v.resolvent(1),complex)-la.solve(matrix,np.eye(320)),2)<2e-13


def test_all_tensor_components_and_frequency_exchanged_statistics():
    rows=v.complete_tensor()
    assert len(rows)==72 and sum(row['value']['exact']!='0' for row in rows)==8
    value,terms=v.ordered_three_point(0,4,0,1,2)
    reversed_value,_=v.ordered_three_point(4,0,0,2,1)
    assert value+reversed_value==0 and len(terms)==6
    assert len(set(terms))>1 and sum(terms)==value


def test_full_leg_amputation_is_not_the_contact_substitution():
    rows=[v.frequency_row(*pair) for pair in ((1,2),(2,3))]
    a=s.sympify(rows[0]['components'][0]['amputated_coefficient']['exact'])
    b=s.sympify(rows[1]['components'][0]['amputated_coefficient']['exact'])
    assert a!=b and a!=1 and b!=1
    assert rows[0]['components'][0]['connected_three_point']!=rows[0]['components'][0]['amputated_coefficient']


def test_full_pointer_retains_preparation_failures_and_quantum_outputs():
    r=v.complete_pointer();records=r['complete_original_records']
    assert len(records)==17 and sum(s.sympify(x['probability']) for x in records)==1
    assert sum(s.sympify(x['probability']) for x in records if x['stage']=='ground_then_pointer')==s.Rational(1,160)
    assert all(r['total_current_commutators_exact']) and r['four_positive_effects_sum_identity_exact']
