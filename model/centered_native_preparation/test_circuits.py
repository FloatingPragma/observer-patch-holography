from fractions import Fraction as F
import numpy as np
import pytest
from . import circuits,tree
from ..native_spectral_detector import circuits as old
from ..protected_particle_readout import execute


@pytest.mark.parametrize('duration',[F(1),F(2)])
def test_centering_control_is_exactly_n_field_times_n_clock(duration):
    word=old.emit(circuits.center_template(duration))
    block=execute.compile_block(execute.key(word),2)
    expected=np.diag([1,1,1,np.exp(1j*float(duration)/2)])
    assert np.linalg.norm(block.success-expected)<3e-14
    assert sum(np.linalg.norm(b)**2 for group in block.failures for b in group)<1e-25


def test_missing_centering_phase_has_opposite_charged_offsets():
    control=circuits.centered_clock_error_countercontrol()
    assert control['omitted_N_phase_with_H_ground_reference']=={'addition':'1/2','removal':'-1/2'}
    assert control['omitted_N_phase_with_K_ground_reference']=={'addition':'9/2','removal':'7/2'}
    # H and N commute, but deleting the controlled N phase changes a clock
    # superposition by a relative phase. It is not a harmless global phase.
    input=np.array([0,1,0,1],complex)/np.sqrt(2)
    word=old.emit(circuits.center_template(F(1)))
    block=execute.compile_block(execute.key(word),2)
    assert np.linalg.norm(block.success@input-input)>.3
