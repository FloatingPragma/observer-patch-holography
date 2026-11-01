"""Finite exact spin twirl, compiled selector and conditional ground words."""
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np

from model.observer.source_scalar import compiler
from model.certified_reference import certificate as cert, interval as iv
from model.ground_spectral_source import preparation as parent
from model.reference_recovery.recovery import execute_complete
from model.protected_particle_readout import execute


SELECTOR_BITS = 4
SELECTOR_WEIGHTS = tuple([5, 8, 5]*5+[0])


def quadrature_certificate():
    """Three polar nodes integrate every polynomial of degree at most five.

    x=0,+/-sqrt(3/5), with weights4/9,5/18,5/18, is exact for the
    normalized uniform measure on[-1,1]. Five azimuths kill every nonzero
    magnetic difference in[-4,4]. No infinite group average is assumed.
    """
    polar = []
    for degree in range(6):
        value = F(1) if degree == 0 else (F(5, 9)*F(3, 5)**(degree//2) if degree%2 == 0 else F())
        target = F(1, degree+1) if degree%2 == 0 else F()
        if value != target: raise ArithmeticError('polar quadrature lost polynomial exactness')
        polar.append(str(value))
    diagonal = []
    for k in range(5):
        # Integrate C(4,k)((1+x)/2)^(4-k)((1-x)/2)^k by expansion.
        value = F()
        for a in range(5-k):
            for b in range(k+1):
                degree = a+b
                moment = F(1, degree+1) if degree%2 == 0 else F()
                value += F(math.comb(4, k)*math.comb(4-k, a)*math.comb(k, b)*(-1)**b, 16)*moment
        if value != F(1, 5): raise ArithmeticError('spin-two twirl is not the full ground identity')
        diagonal.append(str(value))
    return {'ground_spin':2, 'ground_dimension':5, 'conditional_words':15,
            'selector_integer_weights':list(SELECTOR_WEIGHTS), 'selector_denominator':90,
            'polar_nodes':['-sqrt(3/5)','0','sqrt(3/5)'],
            'polar_weights':['5/18','4/9','5/18'], 'azimuths_over_pi':[str(F(2*j,5)) for j in range(5)],
            'polar_moments_through_degree_five':polar,
            'azimuth_cancellation':'For 0<|m-n|<=4 the complete five-root sum is zero.',
            'all_five_diagonal_ground_weights':diagonal,
            'exact_mixture':'P_ground/5', 'postselection_required':False}


@lru_cache(maxsize=1)
def selector_word():
    amplitudes = np.sqrt(np.asarray(SELECTOR_WEIGHTS)/90)
    # The existing binary-tree compiler is explicit at every dimension;
    # use its tree routine rather than the unrelated small-dimension QR tape.
    return tuple(compiler._compile_real_state_tree(amplitudes, SELECTOR_BITS))


def selector_certificate():
    word=selector_word(); offset=0; fixed=F(); errors=[]; nodes=0
    prefix=[0]
    for weight in SELECTOR_WEIGHTS: prefix.append(prefix[-1]+weight)
    for target in reversed(range(SELECTOR_BITS)):
        controls=SELECTOR_BITS-target-1; coefficients=[]
        for mask in range(1<<controls):
            label=['I']*SELECTOR_BITS;label[target]='Y'
            for j in range(controls):
                if mask>>j&1:label[target+j+1]='Z'
            offset,error,variables=cert.consume(word,offset,cert.pauli(label))
            if len(variables)!=1:raise ArithmeticError('one selector Walsh coefficient required')
            coefficients.append(variables[0]);fixed+=error
        maximum=F()
        for index,angle in enumerate(cert.inverse_walsh(coefficients)):
            start=index*(1<<(target+1));half=1<<target
            left=prefix[start+half]-prefix[start];right=prefix[start+2*half]-prefix[start+half]
            maximum=max(maximum,cert.two_vector_error(angle,left,right,False));nodes+=1
        errors.append(maximum)
    if offset!=len(word):raise ArithmeticError('unconsumed selector event')
    error=fixed+sum(errors,F())+cert.controls_bound(word,{})
    return {'pulses':len(word),'word_sha256':cert.word_digest(word),
            'conditional_nodes':nodes,'complete_selector_error':str(error),
            'invalid_record_15_is_retained_as_stop':True}


def polar_node(index):
    if type(index) is not int or not 0<=index<3:raise ValueError('one of three polar nodes required')
    x=iv.sqrt_fraction(F(3,5))
    return -x if index==0 else iv.Interval.point(0) if index==1 else x


def rotation_plan(record):
    if type(record) is not int or not 0<=record<15:raise ValueError('valid selector record required')
    azimuth,polar=divmod(record,3)
    x=(-math.sqrt(3/5),0.,math.sqrt(3/5))[polar]
    theta=math.acos(x);phi=2*math.pi*azimuth/5
    # In the ordered single-occupation block, XY=+sigma_y, YX=-sigma_y.
    # This is Rz(phi) Ry(theta); vacuum and doublon are fixed.
    return (('XY',-theta/4),('YX',theta/4),('ZI',phi/4),('IZ',-phi/4))


def rotation_word(record):
    return tuple(event for label,angle in rotation_plan(record)
                 for event in compiler.pauli_rotation(label,angle))


def rotation_certificate(record):
    word=rotation_word(record);offset=0;fixed=F();angles=[]
    for label in ('XY','YX'):
        offset,error,variables=cert.consume(word,offset,cert.pauli(label))
        if len(variables)!=1:raise ArithmeticError('one polar coefficient required')
        angles.append(variables[0]);fixed+=error
    if angles[0]!=-angles[1]:raise ArithmeticError('polar pair fails exact vacuum cancellation')
    azimuth,polar=divmod(record,3)
    for label,sign in (('ZI',1),('IZ',-1)):
        offset,error,variables=cert.consume(word,offset,cert.pauli(label,cert.PiMultiple(F(sign*azimuth,10))))
        if variables:raise ArithmeticError('unbound azimuth coefficient')
        fixed+=error
    if offset!=len(word):raise ArithmeticError('unconsumed ground rotation event')
    x=polar_node(polar)
    # iv.sqrt_fraction accepts rationals; monotonic endpoint enclosures
    # retain the algebraic node before the second square root.
    cosine=iv.Interval(iv.sqrt_fraction(F(iv.SCALE+x.lo,2*iv.SCALE)).lo,
                       iv.sqrt_fraction(F(iv.SCALE+x.hi,2*iv.SCALE)).hi)
    sine=iv.Interval(iv.sqrt_fraction(F(iv.SCALE-x.hi,2*iv.SCALE)).lo,
                     iv.sqrt_fraction(F(iv.SCALE-x.lo,2*iv.SCALE)).hi)
    actual_sine,actual_cosine=iv.sine_cosine(2*angles[1])
    polar_error=(actual_sine-sine).absolute_upper()+(actual_cosine-cosine).absolute_upper()
    local_error=polar_error+fixed+cert.controls_bound(word,{})
    return {'record':record,'local_word_sha256':cert.word_digest(word),
            'local_pulses':len(word),'all_eight_sites_pulses':8*len(word),
            'local_complete_operator_error':str(local_error),
            'all_eight_sites_error':str(8*local_error)}


def rotation_blocks(record):
    instrument=execute.compile_block(execute.key(rotation_word(record)),2)
    return tuple((instrument,(2*x,2*x+1)) for x in range(8))


def certificate():
    quadrature=quadrature_certificate();selector=selector_certificate()
    rotations=[rotation_certificate(j) for j in range(15)]
    old=parent.certificate()
    error=F(old['complete_preparation_state_error'])+F(selector['complete_selector_error'])+max(F(x['all_eight_sites_error']) for x in rotations)
    if error>=F(1,10**9):raise ArithmeticError('complete mixed-ground preparation accuracy failed')
    return {'quadrature':quadrature,'selector':selector,'rotations':rotations,
            'original_ground_tree_pulses':old['tree_pulses'],
            'original_inverse_frame_pulses':old['frame_pulses'],
            'complete_mixed_ground_preparation_error':str(error),
            'native_selector_and_conditional_words':True,
            'blank_initialization_and_feedforward_are_supplied_services':True,
            'ground_selected_by_autonomous_repair':False}


@lru_cache(maxsize=1)
def selector_run():
    initial=np.zeros(16);initial[0]=1
    return execute_complete(selector_word(),4,initial)
