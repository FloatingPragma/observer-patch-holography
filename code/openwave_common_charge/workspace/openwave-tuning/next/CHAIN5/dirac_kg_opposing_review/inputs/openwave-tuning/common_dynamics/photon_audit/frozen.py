"""Strict read-only adapter for the immutable quadratic spectral experiment."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import sys

sys.dont_write_bytecode = True
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PROVIDER = HERE.parent/'quadratic_field'
MANIFEST_SHA256 = '43fca237a676954e4a213b0888e2cdfee70bb9d8cefd8fbf113f5d2a60cbc81e'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_inputs():
    p = PROVIDER/'MANIFEST.json'
    assert sha(p) == MANIFEST_SHA256, 'frozen producer manifest mismatch'
    manifest = json.loads(p.read_text())
    pins = {str(p.relative_to(ROOT)):MANIFEST_SHA256}
    for group in ('files', 'source_pins'):
        for name, digest in manifest[group].items():
            assert sha(ROOT/name) == digest, name
            assert name not in pins or pins[name] == digest
            pins[name] = digest
    return pins


@lru_cache(None)
def _bundle():
    verify_inputs()
    return {name:json.loads((PROVIDER/(name+'.json')).read_text()) for name in
            ('SPECTRAL_FRAME', 'FOCK', 'NATIVE_PHOTON', 'NATIVE_PHOTON_PROGRAM')}


def spectral_certificate():
    """Return the fixed rational frame and its frozen bound, without eigh."""
    data = _bundle()
    result = {}
    for key, value in data['FOCK']['spectral_certificate'].items():
        result[key] = [F(x) for x in value] if isinstance(value, list) else F(value)
    result['denominator'] = int(result['denominator'])
    frame = data['SPECTRAL_FRAME']
    assert frame['denominator'] == result['denominator'] == 1 << 44
    result['frame_integers'] = [row.copy() for row in frame['frame_integers']]
    return result


def _native():
    # This module imports no quadratic provider module and selects no basis.
    sys.path.insert(0, str(HERE.parent/'vector_field'))
    import native_field
    return native_field


def phase_targets(time=F(1)):
    """Exact outward target intervals; this alone supplies no word at other times."""
    if type(time) is not F or not 0 < time <= 8:
        raise ValueError('exact Fraction time in (0,8] required')
    iv = _native().iv
    frequencies = [iv.sqrt_fraction(F(3, 2)*x) for x in spectral_certificate()['eigenvalues_exact']]
    vacuum = sum(frequencies, iv.Interval.point(0))/2
    energies = [vacuum]+[vacuum+x for x in frequencies]+[iv.Interval.point(0)]*65
    return [-time*x for x in energies], vacuum


def _default_time(time):
    if type(time) is not F or time != 1:
        raise ValueError('emitted frozen native word is certified at exact time1 only')


def word(setting='X', time=F(1)):
    """Replay the immutable nine-bit word. Never invoke the provider generator."""
    _default_time(time)
    if setting not in ('X', 'Y'):
        raise ValueError('configured X or Y setting required')
    data = _bundle()
    record = next(x for x in data['NATIVE_PHOTON']['experiments'] if x['setting'] == setting)
    program = [{**p, 'targets':p['targets'].copy()} for p in data['NATIVE_PHOTON_PROGRAM']['words'][setting]]
    phases, vacuum = phase_targets(time)
    midpoint = [float(F(x.lo+x.hi, 2*_native().iv.SCALE)) for x in phases]
    return program, F(record['native_word_and_pointer_bound_exact']), np.asarray(midpoint), vacuum


def photon_phase(time=F(1)):
    """Extract the common frozen eight-bit phase word and certify its Walsh phases."""
    _default_time(time)
    nf = _native()
    sys.path.insert(0, str(HERE.parent/'clock'))
    import phase_clock
    H = nf.H(0); n = len(H)
    diagonals = []
    for setting in ('X', 'Y'):
        program = word(setting, time)[0]
        clock = nf.compiler.remap(phase_clock.local_words(1, setting, time), [8])
        assert program[:n] == H and program[-len(clock):] == clock
        assert program[-len(clock)-n:-len(clock)] == H
        diagonals.append(program[n:-len(clock)-n])
    assert diagonals[0] == diagonals[1]
    targets, vacuum = phase_targets(time)
    program = diagonals[0]
    offset = 0; bound = F()
    for mask in range(256):
        label = ''.join('Z' if mask >> j & 1 else 'I' for j in range(8))
        offset, fixed, variable = nf.cert.consume(program, offset, nf.cert.pauli(label))
        assert len(variable) == 1
        target = sum((v*(-1)**((i & mask).bit_count()) for i, v in enumerate(targets)), nf.iv.Interval.point(0))/256
        bound += fixed+(nf.iv.Interval.point(variable[0])-target).absolute_upper()
    assert offset == len(program)
    bound += nf.cert.controls_bound(program, {})
    phases = [float(F(x.lo+x.hi, 2*nf.iv.SCALE)) for x in targets]
    return program, bound, np.asarray(phases), vacuum
