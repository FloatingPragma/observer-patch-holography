"""Exact rational specialization of the transverse spectral-window proof."""
from fractions import Fraction as F


def positive(value, name):
    if type(value) not in (int, F) or value <= 0:
        raise ValueError(name+' must be a positive exact rational')
    return F(value)


def window(momentum, magnetization, exchange=F(1)):
    p=positive(momentum,'momentum'); m=positive(magnetization,'magnetization'); j=positive(exchange,'exchange')
    if p>1 or m>F(1,2):
        raise ValueError('require 0<p<=1 and 0<m<=1/2')
    return {'momentum_radius':str(p),'magnetization_lower_bound':str(m),'exchange':str(j),
            'plateau_energy_lower':str(m*j*p/8),'plateau_energy_upper':str(32*j*p/m**2),
            'filter_energy_lower':str(m*j*p/16),'filter_energy_upper':str(64*j*p/m**2),
            'filtered_weight_lower':str(m*m*p/96),'first_moment_upper':str(j*p*p/2),
            'conditional_mean_energy_upper':str(48*j*p/m**2)}


def selection_defect(field, support, edges, exchange=F(1)):
    b=positive(field,'field'); j=positive(exchange,'exchange')
    if type(support) is not int or support<1 or type(edges) is not int or edges<2:
        raise ValueError('positive support and at least two edges required')
    return {'unit_norm_local_ground_defect_upper':str(b*support+6*j/F(edges)),
            'dressed_undressed_expectation_difference_upper':str(4/F(edges)**3)}


def certificate():
    # E(k)>=p²/20, E(k)<=p²/2, 5<=E(k+Q)<=6 on p/2<=|k|<=p<=1.
    # sqrt(E E_Q)>=p/2 and <=2p.  The exact nonnegative margins below
    # are sufficient for each square-root comparison used in the proof.
    margins={'pi_bound':F(10)-F(22,7)**2,
             'root_two_lower':F(2)-F(7,5)**2,
             'weight_coefficient':F(9,64)*F(7,5)-F(1,16)-F(1,8),
             'root_120_upper':F(12)**2-F(120)}
    if not all(v>0 for v in margins.values()):raise ArithmeticError('spectral constants failed')
    return {'exact_positive_margins':{k:str(v) for k,v in margins.items()},
            'parameter_examples_not_measured_magnetizations':[window(F(1,2**n),F(1,10)) for n in (2,4,8)],
            'selection_examples':[selection_defect(F(1,e),3,e) for e in (12,96,768)],
            'magnetization_numerical_value_computed':False,
            'spectral_weight_is_for_L2_normalized_momentum_packet':True,
            'support_is_strictly_positive_energy':True}


def order_margins():
    # 1/12-2/(5sqrt24)>1/1000 follows by squaring positive numbers.
    margin=(F(1,12)-F(1,1000))**2-F(1,150)
    assert margin>0
    # f(e)=e/9-(2/5)sqrt(e/18) is increasing for e>=3/4.
    derivative_margin=F(1,9)**2-F(2,5)**2/(4*18*F(3,4))
    assert derivative_margin>0
    return {'order_squared_margin':str(margin),'monotonicity_squared_margin':str(derivative_margin),
            'staggered_bragg_weight_lower':'1/1000','selected_magnetization_lower':'1/1000',
            'rigorous_window':window(F(1,1000),F(1,1000))}
