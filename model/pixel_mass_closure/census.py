"""Exact one-loop index arithmetic on the retained native current representation."""
from fractions import Fraction as F

import sympy as s
from model.common_source.chiral_candidate import action


def coefficients():
    operators=action.identification_operators(families=False)
    color,weak,Y=(operators[key] for key in ('color_C2','weak_C2','Y'))
    if any(matrix.shape!=(15,15) or matrix!=s.diag(*matrix.diagonal()) for matrix in (color,weak,Y)):
        raise ArithmeticError('the native fifteen-state joint readout is required')
    # Tr(C2_R)=dim(g)*T(R), with the same normalized source Casimir operators.
    t3=F(str(s.trace(color)/8));t2=F(str(s.trace(weak)/3));yy=F(str(s.trace(Y*Y)))
    groups={}
    for i in range(15):
        key=tuple(str(matrix[i,i]) for matrix in (Y,color,weak))
        groups[key]=groups.get(key,0)+1
    if groups != {('1/6','4/3','3/4'):6,('-2/3','4/3','0'):3,('1/3','4/3','0'):3,
                  ('-1/2','0','3/4'):2,('1','0','0'):1}:
        raise ArithmeticError('native current census differs from the declared minimal matter')
    sm=(F(3,5)*(F(2,3)*3*yy+F(1,3)*F(1,2)),
        -F(22,3)+F(2,3)*3*t2+F(1,3)*F(1,2),-11+F(2,3)*3*t3)
    mssm=(F(3,5)*(3*yy+1),-6+3*t2+1,-9+3*t3)
    if sm != (F(41,10),F(-19,6),F(-7)) or mssm != (F(33,5),F(1),F(-3)):
        raise ArithmeticError('one-loop index arithmetic failed')
    return {'one_copy_indices':{'SU3':str(t3),'SU2':str(t2),'sum_dimension_Y_squared':str(yy)},
            'native_joint_current_blocks':[{'Y':y,'color_C2':c,'weak_C2':w,'dimension':n}
                                           for (y,c,w),n in sorted(groups.items())],
            'indices_computed_from_native_15_state_operators':True,
            'SM_three_families_one_complex_Higgs':[str(x) for x in sm],
            'MSSM_three_chiral_families_two_chiral_Higgs':[str(x) for x in mssm],
            'convention':'d g_i / d log(mu) = b_i g_i^3 / (16 pi^2), g1^2=(5/3)gY^2',
            'four_dimensional_one_loop_functional':'declared conditional input',
            'MSSM_extra_fields_present_in_native_quadratic_CAR_model':False}
