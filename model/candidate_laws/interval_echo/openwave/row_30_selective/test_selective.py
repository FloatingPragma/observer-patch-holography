"""Independent constraints on the physical controls and both spatial conventions."""
import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.row_30_selective import analyzer, experiment, source, staggered


def test_native_transposed_analyzers_on_complex_signal():
    setup = source.Source()
    pulses = experiment.settings(setup, source.HERE/"analyzers.json")
    field = np.random.default_rng(667).standard_normal((2, len(setup.sites)))
    signal = field[0] + 1j*field[1]
    signal /= np.linalg.norm(signal)
    for spec in pulses:
        prepared = experiment.pulse(setup, spec)
        assert abs(experiment.reverse_read(setup, spec, signal) - abs(prepared @ signal)**2) < 1e-11


def test_free_propagation_and_disabled_source_independently():
    setup = source.Source()
    specs = experiment.settings(setup, source.HERE/"controls.json")
    beat = experiment.pulse(setup, specs["beat_4"])
    native = expm_multiply(-1j * .41 * setup.h, beat)
    prediction = experiment.fly(setup, beat, np.array([.41]))[0]
    assert np.linalg.norm(native - prediction) < 1e-12
    disabled = setup.trajectory(np.zeros((source.LAYERS, len(setup.ports))))[0][-1]
    assert np.max(abs(setup.modes.T @ disabled)**2) < .34
    assert np.max(abs(abs(setup.modes.T @ disabled)**2 - setup.modes[0]**2)) < 1e-12


def test_plain_and_staggered_rest_have_nonzero_variance():
    setup = source.Source()
    signs = staggered.alternate_geometry(setup)
    _, spatial = analyzer.spatial_basis(setup)
    for state, expected in ((spatial[:, 0], 137.75), (signs*spatial[:, 0], 12.)):
        drift = setup.h @ state
        variance = drift @ drift - (state @ drift)**2
        assert abs(variance - expected) < 1e-9
