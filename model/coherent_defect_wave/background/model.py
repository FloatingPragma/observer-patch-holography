"""A source-bound preparation of a pure stationary triplet background.

The preparation uses a ten-level Fourier ancilla and controlled words of the
owned link-move grammar. These are explicit optional apparatus services.
Its output is a stationary state, without a physical ground-energy claim.
"""
from collections import deque
from fractions import Fraction as F
from functools import lru_cache
import sympy as s

from ..model import JointObserver, SparseState
from ..gauss import gamma_amplitudes


def background(observer, internal=0):
    if type(internal) is not int or internal not in range(3):
        raise ValueError('this stationary background uses a triplet component')
    return SparseState({(g,p,internal,0,a):F(1) for g in observer.link_words
                        for p in range(12) for a in range(60)}, 7200)


def compiled_paths(observer):
    """One source word to every link value, ordered in five gauge pairs."""
    start = observer.initial_word
    paths, queue = {start:()}, deque([start])
    while queue:
        g = queue.popleft()
        for edge in observer.law.interior:
            h = observer.law.transport(g, edge)[0]
            if h not in paths:
                paths[h] = paths[g]+(edge,)
                queue.append(h)
    star = observer.law.vertex_stars[5]
    representatives = sorted(g for g in paths if g < g ^ star)
    order = representatives+[g ^ star for g in representatives]
    assert len(order) == 10
    return tuple((g,paths[g]) for g in order)


def apply_word(observer, g, word):
    for edge in word:
        g = observer.law.transport(g, edge)[0]
    return g


def compiled_link_event(observer, g, port, scratch, edge):
    """Reversible port copy/route/control/unroute/uncompute permutation.

    Scratch has twelve values. With initial scratch zero, this implements
    T_edge on the link and returns both the port and scratch unchanged.
    The intended preparation is in the triplet, where these port operations
    also preserve the admitted interior constraint at every suboperation.
    """
    if type(edge) is not int or edge not in observer.law.interior:
        raise ValueError('compiled link control requires an interior edge')
    if port not in range(12) or scratch not in range(12):
        raise ValueError('port or scratch outside the twelve-level register')
    scratch = (scratch+port) % 12
    if port == scratch:
        port = 5
    elif port == 5:
        port = scratch
    if port in observer.law.edges[edge]:
        g = observer.law.transport(g, edge)[0]
    if port == scratch:
        port = 5
    elif port == 5:
        port = scratch
    scratch = (scratch-port) % 12
    return g,port,scratch


def _cyclotomic_reduce(expr):
    z = s.Symbol('z')
    return s.rem(s.Poly(s.expand(expr), z), s.Poly(z**4-z**3+z**2-z+1,z)).as_expr()


def _conjugate(expr):
    z=s.Symbol('z')
    poly=s.Poly(expr,z)
    return _cyclotomic_reduce(sum(c*z**((-monomial[0]) % 10)
                                  for monomial,c in poly.terms()))


@lru_cache(maxsize=1)
def certificate():
    observer=JointObserver(); paths=compiled_paths(observer)
    star=observer.law.vertex_stars[5]
    for target,word in paths:
        assert apply_word(observer, observer.initial_word,word)==target
        assert apply_word(observer, observer.initial_word ^ star,word)==target ^ star
    # The scratch protocol is a genuine finite permutation, and uncomputes.
    for edge in observer.law.interior:
        outputs={compiled_link_event(observer,g,p,r,edge)
                 for g in observer.link_words for p in range(12) for r in range(12)}
        assert len(outputs)==10*12*12
        for g in observer.link_words:
            for p in range(12):
                assert compiled_link_event(observer,g,p,0,edge)==(
                    observer.law.transport(g,edge)[0],p,0)
    z=s.Symbol('z'); rows=[]
    for ell in range(10):
        # Input is the unnormalized two-word sum; its norm denominator is2.
        coefficients={g:s.Integer(0) for g in observer.link_words}
        for k,(target,_) in enumerate(paths):
            phase=z**((-ell*k) % 10)/10
            coefficients[target]+=phase
            coefficients[target ^ star]+=phase
        coefficients={g:_cyclotomic_reduce(c) for g,c in coefficients.items()}
        probability=_cyclotomic_reduce(sum(_conjugate(c)*c for c in coefficients.values())/2)
        assert probability == (s.Rational(1,5) if ell % 2 == 0 else 0)
        assert all(coefficients[g]==coefficients[g ^ star] for g in observer.link_words)
        rows.append({'fourier_outcome':ell,'conditional_probability':str(probability),
                     'coefficients_by_link':{str(g):str(c) for g,c in coefficients.items()}})
    assert sum(s.Rational(r['conditional_probability']) for r in rows)==1
    fixed=[]
    for internal in range(3):
        state=background(observer,internal)
        assert state.norm_squared()==1
        assert observer.step(state)==state
        assert observer.inverse_step(state)==state
        assert gamma_amplitudes(state.as_dict())==state.as_dict()
        fixed.append({'internal_component':internal,'norm_squared':'1',
                      'joint_U_eigenvalue':'1','gauss_eigenvalue':'1'})
    return {'schema':'oph.coherent_defect_wave.background.v1',
            'fixed_triplet_backgrounds':fixed,
            'source_path_words':[{'target':g,'edges':list(word)} for g,word in paths],
            'maximum_path_length':max(len(word) for _,word in paths),
            'port_scratch_permutation_checks':5*10*12*12,
            'ten_level_fourier_ancilla_outcomes':rows,
            'positive_gauss_preparation_probability':'1/2',
            'background_success_given_positive_gauss':'1/5',
            'joint_success_per_fresh_ambient_trial':'1/10',
            'fresh_trial_failure_tail':'(9/10)^k',
            'success_output':'uniform ten link words times uniform12ports times internal j<3 times parity0 times uniform60addresses',
            'preparation_services':['uniform ten-level ancilla','controlled source words',
                                    'reversible twelve-level port scratch','Fourier ancilla readout'],
            'scope':{'same_U_pure_stationary_background':True,
                     'selected_physical_vacuum':False,'unique_background':False,
                     'energy_ground_state_proved':False,'physical_mass_identified':False}}
