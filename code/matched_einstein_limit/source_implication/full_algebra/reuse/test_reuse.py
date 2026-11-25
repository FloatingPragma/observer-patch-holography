import json
from pathlib import Path

import numpy as np

from . import measure


def test_exact_choi_assembly():
    rng = np.random.default_rng(46)
    raw = rng.integers(-measure.UNIT//3,measure.UNIT//3,(8,4,4,2))
    r, i = measure.moment(raw)
    converted = (r.astype(float)+1j*i.astype(float))/float(measure.UNIT**2)
    ks = (raw[...,0]+1j*raw[...,1])/measure.UNIT
    assert np.max(abs(converted-sum(np.kron(k,k.conj()) for k in ks))) < 2e-15


def test_complete_instrument_changes_left_coin():
    report = measure.run()
    assert report["tick_cap"] == 12
    assert np.array(report["left_X_change_center"]).shape==(2,12)
    assert report["channel_hermitian_trace_norm_error_upper"] < 1.67e-6
    assert report["all_settings_left_X_change_upper"] < -.00051
    assert report["largest_trace_preservation_entry_residual"] < 1e-10
    saved = json.loads(Path(measure.__file__).with_name("results.json").read_text())
    np.testing.assert_allclose(report["left_X_change_center"],saved["left_X_change_center"],atol=1e-11,rtol=0)
    assert report["simulated_left_X"]["sham_plus"] == saved["simulated_left_X"]["sham_plus"]
