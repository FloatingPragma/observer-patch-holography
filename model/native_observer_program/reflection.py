"""Complete event preparation and reflection from native proper-code words.

The ideal tree uses exact subtree counts. Emitted controls use certified
binary angles; actual six-level primitives retain their complete exports.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
import numpy as np
from model.observer.source_scalar import compiler
from model.native_effect_synthesis import intervals as iv, synthesis as syn

PULSE_GATE = F(2, 10**11)


def _size(events):
    if type(events) is not int or not 2 <= events <= 1260:
        raise ValueError('between two and1260 declared events required')
    bits = (events-1).bit_length()
    return bits, 1 << bits


def _walsh(values):
    """Unnormalized exact integer Walsh transform, including zero factors."""
    a = list(values)
    if not a or len(a) & (len(a)-1) or any(type(x) is not int for x in a):
        raise ValueError('complete dyadic integer angle array required')
    width = 1
    while width < len(a):
        for begin in range(0, len(a), 2*width):
            for j in range(begin, begin+width):
                x, y = a[j], a[j+width]
                a[j], a[j+width] = x+y, x-y
        width *= 2
    return a


def _coefficients(values):
    out = _walsh(values)
    assert _walsh(out) == [len(values)*x for x in values]
    return [iv.binary_exact(x, iv.AQ*len(values)) for x in out]


@lru_cache(maxsize=8)
def plan(events):
    bits, size = _size(events)
    word = []; levels = []
    weights = [int(j < events) for j in range(size)]
    for target in reversed(range(bits)):
        controls = bits-target-1
        angles = []
        for prefix in range(1 << controls):
            block = weights[prefix*(1 << (target+1)):(prefix+1)*(1 << (target+1))]
            angles.append(-iv.split_angle(sum(block[:1 << target]), sum(block[1 << target:])))
        coefficients = _coefficients(angles)
        for mask, angle in enumerate(coefficients):
            label = ['I']*bits; label[target] = 'Y'
            for control in range(controls):
                if mask >> control & 1: label[target+1+control] = 'Z'
            word.extend(compiler.pauli_rotation(label, angle))
        levels.append((target, tuple(angles)))
    phases = [0]+[iv.pi_grid()]*(size-1)
    diagonal = []
    for mask, numerator in enumerate(_walsh(phases)):
        # The2048-address identity coefficient can exceed binary64's53-bit
        # significand. Split it into exact commuting dyadic requests; do not
        # round the synthesized phase or silently replace the stored word.
        denominator=iv.AQ*size
        try:parts=[iv.binary_exact(numerator,denominator)]
        except ArithmeticError:
            low=numerator%(1<<26);high=numerator-low
            parts=[iv.binary_exact(high,denominator),iv.binary_exact(low,denominator)]
        assert sum(map(F,parts),F())==F(numerator,denominator)
        for angle in parts:
            diagonal.extend(compiler.pauli_rotation(['Z' if mask >> bit & 1 else 'I' for bit in range(bits)],angle))
    reflection = compiler.inverse(word)+diagonal+word
    compiler._validate(reflection, bits)
    return {'events':events,'bits':bits,'size':size,'preparation':word,
            'reflection':reflection,'levels':levels,'diagonal':diagonal}


def complete_bound(pulses, target_error):
    if type(pulses) is not int or pulses <= 0 or target_error < 0:
        raise ValueError('positive word length and nonnegative exact error required')
    x = pulses*PULSE_GATE
    if x >= F(1, 10): raise ValueError('declared native word exceeds the error domain')
    error = x/(1-x)+target_error
    failure = pulses*PULSE_GATE**2/(1-x)**2
    return {'stored_primitive_operator_gate':str(PULSE_GATE),
            'target_angle_operator_error_upper':str(target_error),
            'full_word_operator_error_upper':str(error),
            'all_first_processor_failure_mass_upper':str(failure),
            'complete_halfdiamond_upper':str((error*(2+error)+failure)/2),
            'complete_map_gain_upper':str(1/(1-x)**2+failure),
            'physical_hardware_precision_inferred':False}


def factor_matrix(events):
    """All columns of the declared binary-angle tree, without ideal replacement.

    Commuting Pauli factors recombine by their exact inverse Walsh identity.
    This evaluates their angle-defined block rotations. It is separate from
    direct stored six-level execution and carries its source error bound.
    """
    p = plan(events); value = np.eye(p['size'])
    for target, angles in p['levels']:
        for prefix, numerator in enumerate(angles):
            start = prefix*(1 << (target+1)); half = 1 << target
            left = slice(start,start+half); right = slice(start+half,start+2*half)
            angle = float(F(numerator,iv.AQ)); c,s = np.cos(angle),np.sin(angle)
            a,b = value[left].copy(),value[right].copy()
            value[left],value[right] = c*a+s*b,-s*a+c*b
    return value


def terminal_instrument(matrix, amplitudes, events):
    """Valid event columns plus one padding outcome, without conditioning."""
    value = np.asarray(amplitudes,dtype=complex)
    if value.ndim not in (1,2) or value.shape[0] != len(matrix) or not np.all(np.isfinite(value)):
        raise ValueError('every padded event and reference coordinate required')
    result = matrix@value
    probabilities = np.sum(abs(result.reshape(len(matrix),-1))**2,axis=1)
    return probabilities[:events],float(sum(probabilities[events:]))


def workspace_admission(value):
    """Complete optional service workspace instrument, with no conditioning."""
    a=np.asarray(value,dtype=complex)
    if a.ndim<2 or a.shape[0]<2 or not np.all(np.isfinite(a)):
        raise ValueError('complete blank and dirty workspace coordinates required')
    return {'blank_mass':float(np.vdot(a[0],a[0]).real),
            'dirty_failure_mass':float(np.vdot(a[1:],a[1:]).real),
            'normalization_performed':False}


def run_reflections():
    results = []
    for events in (60,120,1260):
        p = plan(events); bits,size = p['bits'],p['size']
        unique = {(e['kind'],e['angle']) for e in p['preparation']+p['reflection']}
        for kind,angle in unique: syn.certify_primitive(kind,angle)
        target = np.zeros(size);target[:events] = 1/np.sqrt(events)
        U = factor_matrix(events)
        phi = float(F(iv.pi_grid(),iv.AQ)); phases = np.ones(size,dtype=complex)*np.exp(1j*phi);phases[0]=1
        # All column effects use the exact low-rank ideal reflection. The
        # angle-defined factor matrix is evaluated independently above.
        expected = 2*np.outer(target,target)-np.eye(size)
        factored = (U*phases)@U.T
        prep_error = float(np.linalg.norm(U[:,0]-target))
        reflected_error = float(np.linalg.norm(factored-expected,ord='fro'))
        direct = {}
        if events <= 120:
            actual_U = compiler.execute(p['preparation'],bits)
            actual_R = compiler.execute(p['reflection'],bits)
            direct = {'actual_stored_native_complete_columns':size,
                      'preparation_column_error_diagnostic':float(np.linalg.norm(actual_U[:,0]-target)),
                      'complete_reflection_Frobenius_error_diagnostic':float(np.linalg.norm(actual_R-expected,ord='fro')),
                      'actual_factor_comparison_Frobenius_error_diagnostic':float(np.linalg.norm(actual_R-factored,ord='fro'))}
            assert direct['complete_reflection_Frobenius_error_diagnostic'] < 1e-8
        # The general reference-safe proof is operator norm, not the finite
        # Frobenius diagnostic or a renormalized selected success branch.
        assert prep_error < 2e-10 and reflected_error < 2e-8
        prep_resources = compiler.resources(p['preparation']); ref_resources = compiler.resources(p['reflection'])
        results.append({'valid_events':events,'member_role':{60:'local30-seam coherent identity/swap events',120:'local30-seam law-labelled enlargement',1260:'actual630-seam federation coherent events'}[events],'padded_events':size,'algorithmic_scratch_qubits':0,
          'complete_factor_columns_checked':size,'preparation_native_pulses':len(p['preparation']),
          'reflection_native_pulses':len(p['reflection']),'unique_full_six_level_primitives_certified':len(unique),
          'preparation_angle_defined_error_diagnostic':prep_error,
          'reflection_angle_defined_Frobenius_error_diagnostic':reflected_error,
          'preparation_resources':prep_resources,'reflection_resources':ref_resources,
          'preparation_bound':complete_bound(len(p['preparation']),F(bits,iv.AQ)),
          'reflection_bound':complete_bound(len(p['reflection']),F(2*bits+1,iv.AQ)),
          'direct_stored_execution':direct,
          'word_sha256':{key:sha256(json.dumps(p[key],sort_keys=True).encode()).hexdigest()for key in ('preparation','reflection')},
          'padding_target_action':'minus identity on every unused event address; retained as terminal failure if read',
          'dirty_workspace_policy':'No algorithmic scratch is used; any appended service workspace is admitted by P_blank with a retained I-P_blank outcome.',
          'all_processor_first_failure_outcomes_retained':True,
          'workspace_control':workspace_admission(np.array([[1,0],[0,.25j]])),
          'global_sign':'2|uniform><uniform|-I, not its negative'})
    return results
