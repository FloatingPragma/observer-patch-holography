"""Exact strong-coupling enclosure of the joined original-law charge gap."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb

from ..common_source.mass_response import gap_certificate as old
from .chain import HOPPING, REPULSION


I = old.RationalInterval


@lru_cache(maxsize=32)
def dirichlet_eta(power, terms=96):
    if type(power) is not int or power < 1 or type(terms) is not int or not 1 <= terms <= 160:
        raise ValueError('positive integer eta power and bounded Euler order required')
    if power == 1: return old.log_bounds(F(2))
    total = F()
    for n in range(terms):
        difference = sum((F((-1)**j*comb(n, j), (j+1)**power) for j in range(n+1)), F())
        if not 0 <= difference <= 1: raise ArithmeticError('completely monotone Euler coefficient failed')
        total += difference/F(2**(n+1))
    # The integral of exp(-x)(1-exp(-x))^n x^(s-1)/Gamma(s)
    # proves all coefficients lie in[0,1]; the geometric tail is2^-terms.
    return I(total, total+F(1, 2**terms))


def strong_gap(hopping=HOPPING, repulsion=REPULSION, order=12):
    if type(hopping) not in (int, F) or type(repulsion) not in (int, F):
        raise ValueError('exact rational source couplings required')
    hopping, repulsion = F(hopping), F(repulsion)
    if hopping <= 0 or repulsion <= 2*hopping or type(order) is not int or not 1 <= order <= 24:
        raise ValueError('t>0,U>2t and bounded positive binomial order required')
    ratio = 2*hopping/repulsion
    result = I.enclose(repulsion-4*hopping)
    coefficient = F(1)
    for k in range(1, order+1):
        coefficient *= (F(1, 2)-(k-1))/k
        result += 8*hopping*coefficient*ratio**(2*k-1)*dirichlet_eta(2*k-1)
    next_coefficient = coefficient*(F(1, 2)-order)/(order+1)
    tail = 8*hopping*abs(next_coefficient)*ratio**(2*order+1)/(1-ratio**2)
    return I(result.lower-tail, result.upper+tail)


def certificate():
    value = strong_gap()
    if not 0 < value.lower <= value.upper < REPULSION:
        raise ArithmeticError('joined source gap or canonical-vacuum distinction failed')
    addition = (REPULSION+value)/2
    removal = (value-REPULSION)/2
    return {'thermodynamic_theorem_input': 'Lieb-Wu half-filled repulsive Hubbard-chain charge-gap theorem',
            'source_hopping': str(HOPPING), 'source_repulsion': str(REPULSION),
            'gap_definition': 'lim_L [E_L(L+1)+E_L(L-1)-2E_L(L)]',
            'exact_charge_gap_interval': value.as_json(),
            'original_H_addition_energy_interval': addition.as_json(),
            'original_H_removal_energy_interval': removal.as_json(),
            'centered_half_gap_interval': (value/2).as_json(),
            'centering_condition': 'K=H-(U/2)N is a declared reservoir/energy convention, not the original mu=0 action',
            'strong_series': 'U-4t+8t sum_(k>=1) binom(1/2,k)(2t/U)^(2k-1) eta(2k-1)',
            'binomial_terms': 12, 'Euler_terms_per_eta': 96,
            'binomial_tail_bound': '8t |binom(1/2,K+1)| r^(2K+1)/(1-r²),r=2t/U<1',
            'Euler_tail_bound': '2^-96', 'floating_quadrature_used': False,
            'common_mu_cancels_from_charge_gap': True,
            'original_H_removal_energy_is_negative': True,
            'positive_global_many_body_gap_claimed': False,
            'three_plus_one_mass_or_species_identified': False,
            'ground_state_preparation_or_detector_executed': False,
            'reference': 'https://arxiv.org/html/cond-mat/0207529'}
