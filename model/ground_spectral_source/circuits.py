"""Source-compiled parity-coherent charged probe and complete pointer readout.

The supplied native qubit service admits coherent N4/N5 arms. This does not
derive a physical fermion-parity phase reference from passive repair.
"""
from fractions import Fraction as F
import math

from model.observer.source_scalar import compiler
from model.certified_reference import certificate as cert
from model.protected_particle_readout import circuits as field, execute
from . import preparation as p


POINTER = 16
BITS = 17
PHASES = (F(), F(1, 2))


def majorana_plan():
    # gamma_1=c_1+c_1^dagger=Z_0 X_1. All four factors commute.
    # exp(i*pi/4*(I-Z_pointer)*(I-gamma_1)) is controlled gamma_1.
    return (('III', F(1, 4)), ('IIZ', F(-1, 4)),
            ('ZXI', F(-1, 4)), ('ZXZ', F(1, 4)))


def inverse(plan):
    return tuple((label, -angle) for label, angle in reversed(plan))


def pointer_word(phase=None):
    if phase is None:
        return [compiler.pulse('Y', -math.pi/4, [2])]+p.emit(majorana_plan())
    if type(phase) is not F or phase not in PHASES:
        raise ValueError('exact analyzer phase0 or pi/2 required')
    return (p.emit(inverse(majorana_plan()))+
            [compiler.pulse('I', float(phase)*math.pi/2, [2]),
             compiler.pulse('Z', -float(phase)*math.pi/2, [2]),
             compiler.pulse('Y', math.pi/4, [2])])


def pointer_template(phase=None):
    if phase is None:
        return [('Y', (2,), cert._pi(-1, 4))]+p.template(majorana_plan())
    if type(phase) is not F or phase not in PHASES:
        raise ValueError('exact analyzer phase0 or pi/2 required')
    return (p.template(inverse(majorana_plan()))+
            [('I', (2,), cert.PiMultiple(phase/2)),
             ('Z', (2,), cert.PiMultiple(-phase/2)), ('Y', (2,), cert._pi(1, 4))])


def pointer_blocks(phase=None):
    return ((execute.compile_block(execute.key(pointer_word(phase)), 3), (0, 1, POINTER)),)


def field_blocks(duration, steps):
    """Bind each expanded block to the original full-Fock source tape."""
    return field.field_binding(duration, steps)


def majorana(state):
    import numpy as np
    a = np.asarray(state)
    if a.shape[0] != 65536: raise ValueError('all16 field bits required')
    masks = np.arange(65536, dtype=np.int64)
    sign = 1-2*(masks & 1)
    if a.ndim == 2: sign = sign[:, None]
    return sign*a[masks ^ 2]


def exact_probe_certificate():
    weights, _ = p.signed_weights()
    occupation = sum(w for mask, w in enumerate(weights) if mask & 1)
    if occupation != 648:
        raise ArithmeticError('exact origin-up density must be1/2')
    return {'probe':'gamma_1=Z_0 X_1=c_1+c_1^dagger',
            'ground_origin_down_empty':True,
            'gamma_G_equals_normalized_down_addition':True,
            'charged_branch_particle_number':5,
            'origin_up_density':'1/2',
            'missing_Jordan_Wigner_sign_control_overlap':'0',
            'missing_Jordan_Wigner_sign_control_state_distance':'sqrt(2)',
            'relative_phase':'original E_N5-E_G; no separately inserted ground-energy scalar',
            'native_qubit_parity_coherence_service_required':True,
            'physical_parity_reference_derived':False}


def localized_defect_certificate():
    """A unitary local creation with an unchanged initial exterior state."""
    import numpy as np
    from model.interacting_vacuum import charged
    from model.field_symmetry import symmetry
    local=[symmetry._annihilator(4,j) for j in range(4)]
    gamma=local[1]+local[1].adjoint()
    identity=symmetry._identity(16)
    if not (gamma@gamma-identity).is_zero():
        raise ArithmeticError('local Majorana lost unitarity')
    for outside in local[2:]:
        if not (gamma@outside+outside@gamma).is_zero():
            raise ArithmeticError('disjoint odd CAR relation changed')
        number=outside.adjoint()@outside
        if not (gamma@number-number@gamma).is_zero():
            raise ArithmeticError('outside even number failed to commute')
    g=charged.reference_vector()
    removed=charged.car(g,1,False)
    added=charged.car(g,1,True)
    if any(np.any(x) for x in removed) or charged.inner(added,added)!=(1296,0):
        raise ArithmeticError('charged local branch is not a normalized addition')
    indices=np.flatnonzero((added[0]!=0)|(added[1]!=0))
    if any(int(x).bit_count()!=5 or not int(x)&2 for x in indices):
        raise ArithmeticError('created branch lost its number or origin-down flag')
    moments=charged.spectral_certificate()['down_addition']
    if moments['gap_mean']!='1/2' or moments['gap_second_moment']!='131/256':
        raise ArithmeticError('source-derived charged energy moments changed')
    return {'creation_is_unitary_and_origin_supported':True,
            'initial_exterior_CAR_restriction_unchanged':True,
            'exterior_proof':'Disjoint even observables commute with the local Majorana. Odd exterior expectations vanish in both the definite-parity ground and charged state; linearity covers the full exterior CAR algebra.',
            'input_number':4,'created_number':5,
            'initial_origin_down_occupation':'0 -> 1',
            'protecting_quantity':'Total particle number during original-H evolution and the stated passive N4+N5 repair waits.',
            'mean_added_energy':'1/2','added_energy_variance':'67/256',
            'energy_eigenstate_or_relativistic_rest_mass':False,
            'persistent_spatial_localization_proved':False,
            'repair_after_final_analyzer_applied':False}
