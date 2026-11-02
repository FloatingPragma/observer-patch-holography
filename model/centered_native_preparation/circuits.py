"""Explicit centered-clock correction and complete four-bin phase instrument."""
from fractions import Fraction as F
from functools import lru_cache
import math

from ..native_spectral_detector import circuits as old
from ..protected_particle_readout import execute
from ..certified_reference import certificate as cert
from . import target

M=4
TAU=F(1)
STEPS=(16,32)


def center_template(duration):
    """exp(+i duration n_field n_clock /2), on a field and clock bit."""
    return sum((cert.pauli(label,duration*sign/8)
                for label,sign in (('II',1),('ZI',-1),('IZ',-1),('ZZ',1))),[])


def ground_phase_word(duration):
    _,_,ground=target.ground();lo,hi=ground['energy_intervals'][0]
    angle=float(duration*(lo+hi)/4)
    return (old.compiler.pulse('I',angle,[0]),old.compiler.pulse('Z',-angle,[0]))


@lru_cache(maxsize=1)
def clock_blocks():
    result=[old.block(old.hadamard_template(),1,(16+j,)) for j in range(2)]
    for j,steps in enumerate(STEPS):
        duration=TAU*(1<<j);bit=16+j
        result.extend(old.controlled_field_blocks(duration,steps,bit))
        # Every exact split term preserves N, so this commutes with the
        # whole ideal splitting. Native rounding is bounded separately.
        correction=old.emit(center_template(duration))
        local=execute.compile_block(execute.key(correction),2)
        result.extend((local,(mode,bit)) for mode in range(16))
        phase=ground_phase_word(duration)
        result.append((execute.compile_block(execute.key(phase),1),(bit,)))
    result.extend(old.fourier_blocks())
    return tuple(result)


def charge_blocks():return old.charge_blocks()


def centered_clock_error_countercontrol():
    # A change of charge clock is not one common scalar energy subtraction.
    return {'omitted_N_phase_with_K_ground_reference': {'addition':'9/2','removal':'7/2'},
            'omitted_N_phase_with_H_ground_reference': {'addition':'1/2','removal':'-1/2'},
            'H_ground_reference_minus_K_ground_reference':'4'}
