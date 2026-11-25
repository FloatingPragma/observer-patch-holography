import hashlib
import json
from pathlib import Path

import numpy as np

from . import clock_history as clock


def test_native_timestamp_rows_and_retained_charge_alternative():
    r=clock.run()
    saved=json.loads(Path(clock.__file__).with_suffix(".json").read_text())
    payload = json.dumps(saved, sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(payload).hexdigest() == "42618c6da61b58994d554657069494d9992d1cfb9c33372a8e4decb21f08025b"
    # Exact rational inverses recertify the original floating-point proposals.
    # All three bounds must be at least as strong as the published bound.
    numerical = ("charge_hidden_pre_cap_clock_tomography", "original_means_plus_six_bins", "per_charge_mean_tomography")
    for name in numerical:
        original, current = saved[name], r[name]
        assert current.keys() == original.keys()
        for key in original:
            if key == "smallest_singular_value_lower":
                assert current[key] >= original[key]
            elif key.endswith("norm_upper"):
                assert current[key] <= original[key]
            else:
                assert current[key] == original[key]
    assert {k: v for k, v in r.items() if k not in numerical} == {
        k: v for k, v in saved.items() if k not in numerical
    }
    assert len(set(tuple(b) for b in r["kept_clock_bins"]))==18
    assert all(t<r["cap_ticks"] for _,_,t in r["kept_clock_bins"])
    assert len(r["original_means_plus_six_bins"]["additional_bins"])==6
    assert r["witness_event"]["route_and_return_tick"]==[2,3,10]
    assert float(r["witness_event"]["probability_slope"].strip("[]").split(",")[1]) < -3e-6
    matrix=r["pooled_timestamp_matrix"]
    a=np.array(matrix["numerators"],dtype=float)/int(matrix["denominator"])
    assert np.linalg.svd(a,compute_uv=False)[-1] > r["charge_hidden_pre_cap_clock_tomography"]["smallest_singular_value_lower"]
