"""Exact elimination and outward bounds excluding one scale/census choice."""
import sympy as s
from .closure import C


def exact_elimination():
    x,y = s.symbols('x y',positive=True)
    A,L_over_two_pi = s.symbols('A L_over_two_pi')
    solution = s.solve((x-A+s.Rational(19,6)*L_over_two_pi,
                        y-A+7*L_over_two_pi),(A,L_over_two_pi))
    if any(s.expand(solution[k]-v) != 0 for k,v in {A:(42*x-19*y)/23,L_over_two_pi:6*(x-y)/23}.items()):
        raise ArithmeticError('running elimination failed')
    inverse_one = s.factor(solution[A]+s.Rational(41,10)*solution[L_over_two_pi])
    scale_linear = s.factor(2*solution[L_over_two_pi]-solution[A]/2)
    if s.expand(inverse_one-(333*x-218*y)/115) != 0 or s.expand(scale_linear+(18*x+5*y)/46) != 0:
        raise ArithmeticError('scale elimination failed')
    return {'inverse_alpha_U':str(solution[A]),'log_MU_over_mZ_over_2pi':str(solution[L_over_two_pi]),
            'inverse_alpha_1':str(inverse_one),'scale_linear_over_pi':str(scale_linear)}


def certificate():
    b = C.IvBackend(50)
    # F(x)=9*pi*x/23-2*pi-(1/2)log(8*pi/(5*x)) is increasing for x>0.
    at_six = 54*b.pi/23-2*b.pi-b.log(8*b.pi/30)/2
    if not at_six.a > 1:
        raise ArithmeticError('x>=6 exclusion failed')
    t = b.four_pi_sq/b.num(6)
    total = b.zero
    tails = []
    for builder,tailer in ((C.build_su2_terms,C.su2_tail_bounds),(C.build_su3_terms,C.su3_tail_bounds)):
        terms = builder(b,12)
        numerator = sum((dim*logdim*b.exp(-t*c2) for dim,c2,logdim in terms),b.zero)
        tail = tailer(b,t,12)['s']
        total += numerator+b.hull(b.zero.a,tail.b)
        tails.append(C._iv_pair(tail))
    if not total.b < b.num(1)/b.num(80):
        raise ArithmeticError('heat numerator bound failed')
    return {'domain':{'P_lower_inclusive':'1','MU_at_least_mZ':True,'all_three_couplings_positive':True,
                       'pure_SM_one_loop_no_thresholds':True,'transmutation_coefficient':'4'},
            'exact_elimination':exact_elimination(),
            'scale_exclusion_margin_at_x_6':C._iv_pair(at_six),
            'heat_numerator_upper_enclosure_at_alpha_2_alpha_3_1_over_6':C._iv_pair(total),
            'heat_numerator_strict_rational_upper_bound':'1/80',
            'heat_denominators_lower_bound':'1','required_pixel_entropy_lower_bound':'1/4',
            'infinite_sum_tails':tails,'no_positive_ordered_solution':True,
            'excludes_other_transmutation_or_threshold_designs':False}
