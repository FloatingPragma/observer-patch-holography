import json

import numpy as np

from . import centre, measure, shortened


def test_clamped_perturbations_and_commissioned_diagonal():
    null = shortened.clamped_nullspace()
    assert null["rank"]==12
    assert null["individually_invisible"]==[0,14]
    for lev in range(3):
        assert measure.free_diagonal(lev,measure.equal_degree.ports(lev))["residual"]<2e-13
        low, high=(centre.altered_matrix(lev,7,pin) for pin in (0,1))
        np.testing.assert_allclose(high-low[::-1,::-1],32*np.eye(len(low)),atol=1e-10)


def test_centre_changes_interior_transmission_and_leaves_static_work():
    record=json.loads((measure.HERE/"centre.json").read_text())
    rows=[r for b in record["boards"] for r in b if r["source"]=="sharp_intact" and r["length"]>=4]
    assert len(rows)==21 and all(abs(r["read"]["line"]-32)<.002 for r in rows)
    assert max(centre.native_columns())<1e-13
    signals=centre.endpoint_arms()
    assert max(abs(r["exact_idle_difference"]) for r in signals)<1e-12
    assert signals[2]["native_delta"]>.09
    assert record["held"]["status"]=="physical_fail"
    shortened_result=json.loads((measure.HERE/"shortened.json").read_text())
    assert shortened_result["held"]["profile_pass"]
