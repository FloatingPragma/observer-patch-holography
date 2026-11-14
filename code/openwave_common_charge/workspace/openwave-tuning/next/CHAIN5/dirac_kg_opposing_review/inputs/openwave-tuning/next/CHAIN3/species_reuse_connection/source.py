"""The executed finite register effect and its explicit flavor mass slot."""
from fractions import Fraction as F
from pathlib import Path
import sys
import math
import numpy as np

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[3]
ROOT = HERE/'inputs' if (HERE/'INPUTS.json').exists() else WORKSPACE
RER = ROOT/'reverse-engineering-reality'
FIELD = ROOT/'openwave-tuning/next/CHAIN3/native_field_connection'
sys.path[:0] = [str(RER), str(FIELD)]
from model.lepton_register_round import register
from model.native_flavor_orbit import intervals as iv
import native_oscillators as field_words
from providers import native

ROUNDS = 2**20
PRECISION = 100


def cmul(a, b):
    return (iv.rounded_out(a[0]*b[0]-a[1]*b[1], 70),
            iv.rounded_out(a[0]*b[1]+a[1]*b[0], 70))


def cpow(a, power):
    result = (iv.I(1), iv.I(0))
    while power:
        if power & 1:
            result = cmul(result, a)
        power //= 2
        if power:
            a = cmul(a, a)
    return result


def midpoint(box):
    return (box.lo + box.hi)/2


def executed_register():
    integer_tick = register.full_tick(True)
    assert np.array_equal(integer_tick.T@integer_tick, 16*np.eye(432, dtype=np.int64))
    parts = register.compression(True)
    weights, disjoint = register.joint_weights()
    P = weights['repair_plus'] + weights['hold'] - 2*weights['both']
    assert disjoint and weights['total'] == 1
    # Derive all coefficients from actual integer compression arrays.
    matrices = {key: (a + math.sqrt(2)*b)/(96*math.sqrt(3))
                for key, (a, b) in parts.items()}
    t, i = register.preparation()
    prepared = np.kron((t+math.sqrt(2)*i)/math.sqrt(6),
                       (t+math.sqrt(2)*i)/math.sqrt(6))
    input_frame = np.kron(prepared[:, None], np.eye(3))
    output = integer_tick@input_frame/4
    h, ep, em = register.analyzer_parts()
    phi = (math.sqrt(2)*h+2*(ep+em))/(4*math.sqrt(3))
    analyzer = np.kron(phi[:, None], np.eye(3))
    J = np.array([register.step_of(index//36) for index in range(432)])
    assert np.linalg.norm(output.conj().T@output-np.eye(3)) < 1e-13
    assert np.linalg.norm(analyzer.conj().T@analyzer-np.eye(3)) < 1e-13
    for label, j in (('0', 0), ('+', 1), ('-', -1)):
        actual = analyzer.conj().T@(output*(J == j)[:, None])
        assert np.linalg.norm(actual-matrices[label]) < 1e-13
    return {'integer_tick': integer_tick, 'parts': parts, 'matrices': matrices,
            'P': P, 'weights': weights, 'output': output,
            'analyzer': analyzer, 'J': J}


def enclosures(data, rounds=ROUNDS):
    """Exact rational interval arithmetic, with rounded interval squaring."""
    P = data['P']
    angle = F(1, rounds)
    sine = iv.sine(iv.I(angle), terms=16)
    cosine = 1-2*iv.sine(iv.I(angle/2), terms=16)**2
    a = cpow((1-P+P*cosine, P*sine), rounds)
    root2 = iv.sqrt_rational(2, 75)
    root3 = iv.sqrt_rational(3, 75)
    boxes = {}
    for label, (A, B) in data['parts'].items():
        boxes[label] = [[(int(A[x,y])+root2*int(B[x,y]))/(96*root3)
                         for y in range(3)] for x in range(3)]
    # Fourier eigenvalues from the first column of the executed circulant,
    # using exact cube-root real/imaginary parts.
    roots = [(iv.I(1), iv.I(0)), (iv.I(F(-1,2)), root3/2),
             (iv.I(F(-1,2)), -root3/2)]
    coefficients = []
    lambdas = []
    for k in range(3):
        total = (iv.I(0), iv.I(0))
        for row in range(3):
            entry = (boxes['0'][row][0]+a[0]*(boxes['+'][row][0]+boxes['-'][row][0]),
                     a[1]*(boxes['+'][row][0]-boxes['-'][row][0]))
            value = cmul(entry, roots[(k*row)%3])
            total = (total[0]+value[0], total[1]+value[1])
        assert total[1].lo <= 0 <= total[1].hi
        assert total[0].lo > 0
        lambdas.append(iv.rounded_out(total[0], 60))
        coefficients.append(iv.rounded_out(total[0]**2, 60))
    return {'weak_amplitude': a, 'lambda': lambdas, 'mass': coefficients,
            'mass_midpoints': [midpoint(v) for v in coefficients]}


def complete_register_instrument(data, rounds):
    """All-return analyzer plus complement and every first weak failure.

    Failure j has full Kraus diag(b_J a_J**(j-1)) times the actual
    432-by-3 output frame. Its scalar formulas retain coherence between
    all J branches. No successful branch is renormalized.
    """
    P = float(data['P']); J = data['J']; B = data['output']; A = data['analyzer']
    z = np.exp(1j*J/rounds)
    one = 1-P+P*z
    fail = math.sqrt(P*(1-P))*(z-1)
    all_return = one**rounds
    returned = B*all_return[:, None]
    success = A.conj().T@returned
    complement = returned-A@success
    success_effect = success.conj().T@success
    complement_effect = complement.conj().T@complement
    # Every distinct first-failure row remains in this exact closed form.
    failure_weight = -np.expm1(rounds*np.log(np.maximum(abs(one)**2, 1e-300)))
    complete_failure_effect = B.conj().T@(B*failure_weight[:, None])
    gram = success_effect+complement_effect+complete_failure_effect
    error = float(np.linalg.norm(gram-np.eye(3), 2))
    assert error < 2e-9
    explicit = []
    if rounds <= 32:
        summed = np.zeros((3,3), complex)
        for index in range(rounds):
            kraus = B*(fail*one**index)[:, None]
            effect = kraus.conj().T@kraus
            summed += effect
            explicit.append({'first_failure_round':index+1,
                'effect_real':effect.real.tolist(), 'effect_imag':effect.imag.tolist()})
        assert np.linalg.norm(summed-complete_failure_effect) < 1e-12
    return {'rounds':rounds,'success_effect':success_effect,
            'analyzer_complement_effect':complement_effect,
            'all_first_failure_effect':complete_failure_effect,
            'complete_gram_error':error,'explicit_first_failure_effects':explicit,
            'failure_count':rounds,'failure_kraus_formula':'diag(sqrt(P(1-P))*(exp(iJ/N)-1)*(1-P+P*exp(iJ/N))^(r-1))*U0*input, r=1,...,N',
            'no_postselection':True}


def flavor_mass_word(masses, time=F(1,8), flavor=(1,2), beta=0):
    """exp(-i time Z_beta diag(m0,m1,m2,0)); pad3 is retained."""
    values = tuple(map(F, masses))+(F(0),)
    for mask in range(4):
        coefficient = -F(time)*sum((v*(-1)**((j&mask).bit_count())
                                    for j,v in enumerate(values)), F())/4
        label = ''.join('Z' if mask>>j&1 else 'I' for j in range(2))+'Z'
        yield from field_words.pauli_word(label, coefficient, tuple(flavor)+(beta,))


def flavored_macro(coexact_cap, harmonic_cap, slices, bits, masses, h=F(1,8)):
    """Actual task25 word with precisely its first scalar mass slot replaced."""
    info = field_words.layout(coexact_cap, harmonic_cap)
    flavor = (info['memory_bits'], info['memory_bits']+1)
    yield from flavor_mass_word(masses, h, flavor, info['spin'][0])
    previous = iter(field_words.macro_word(coexact_cap,harmonic_cap,slices,bits,h,F(1)))
    removed = next(previous)
    assert removed == ('Z',(info['spin'][0],),native.angle(0,0,-h))
    yield from previous


def weak_limit_bounds(P, rounds):
    # Centered Bernoulli Taylor remainder and contractive telescoping.
    delta = P*(1-P)/(2*rounds)
    # ||M+-|| <= 1, so ||delta M|| <= 2delta and both contractions have norm<=1.
    return {'all_return_amplitude_error_upper_exact':str(delta),
            'any_first_failure_probability_upper_exact':str(P*(1-P)/rounds),
            'mass_effect_operator_error_upper_exact':str(4*delta),
            'proof':'centered Bernoulli second moment, N-factor contractive telescoping; no empirical target in N',
            'state_domain':'every flavor input and arbitrary untouched reference'}
