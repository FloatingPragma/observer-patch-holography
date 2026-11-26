from fractions import Fraction as F
import json
import math

import numpy as np
import pytest

from ..memories import ENDPOINT_GAP, SENSITIVITY, THETA, Pair
from ..native_pair import NativePair, endpoint
from .experiment import (ALPHA, BANKS_PER_LEVEL, COEFFICIENT, FRAMES, GAIN, HERE, SHAPES,
                         _probability, infer, read_count, settings, verify)


def test_fixed_volume_finite_rods_certify_the_phase_branch_before_sight_of_counts():
    e = F(1,32**6)
    along_one = 1-F(13,256*6)-e
    along_b = 1-F(9*13,256*6)-e
    assert along_b > F(9,10)
    assert 1+9*(along_one**2-1) > (F(9,10)+F(1,32**8))**2
    assert 1+9*((1+e)**2-1) < F(101,100)**2
    assert F(9,10)*F(1039,1000)**3 > 1
    assert F(1039,1000)+F(1,32**4) < F(26,25)
    assert 4*F(22,7)/15*F(101,100)*F(26,25)**4 < 1
    assert 18/(GAIN*ENDPOINT_GAP*SENSITIVITY) < 1/8


def test_executed_native_pulse_and_forward_reference_undo_give_the_amplified_read():
    code=NativePair(0,"right")
    for d in (-0.0001,0.,0.00011775475):
        m=settings(0)["coherent_repeats"]
        signal=code.gate(THETA+d)
        reference_inverse=code.gate((-THETA)%math.pi)
        amplified=np.linalg.matrix_power(signal @ reference_inverse, m-1) @ signal
        sharp=np.diag([5/8,1/8,1/8,1/8])
        state=amplified @ sharp @ amplified.conj().T
        assert abs(endpoint(state,"left")[1,1].real-_probability(m*d)) < 5e-12
        assert np.allclose(state,code.state(THETA+m*d),rtol=0,atol=5e-12)


def test_concentration_budget_and_signed_estimator_include_physical_probabilities():
    partial=0
    for j in range(8):
        cfg=settings(j)
        partial+=F(1,(j+1)*(j+2))
        assert partial==1-F(1,j+2)
        assert cfg["coherent_repeats"]*GAIN == cfg["inverse_radius"]**4
        for x in (-1/8,-0.001,0,0.052,1/8):
            N=cfg["shots_per_bank"]
            hits=round(N*_probability(x))
            _,lo,hi=read_count(hits,N-hits,cfg)
            expected=Pair(x/cfg["coherent_repeats"]).cut_entropy_change*cfg["inverse_radius"]**4
            assert lo-2e-11 <= expected <= hi+2e-11
        with pytest.raises(ValueError,match="incomplete"):
            read_count(2,0,cfg)
    assert math.isclose(ALPHA,0.01) and BANKS_PER_LEVEL==35 and len(FRAMES)==10


def test_replayed_all_frame_and_scale_estimates_distinguish_both_controls():
    output=verify()
    archive=json.loads((HERE.parent/"results.json").read_text())
    previous=float("inf")
    for row in output["levels"]:
        d=row["design"]["inverse_radius"]
        population={r["shape"]:r["trials"] for r in archive["finite_rows"] if r["inverse_ell"]==d}
        assert set(row["banks"])==set(SHAPES)
        for shape,bank in row["banks"].items():
            primary=population[shape]
            for stock,raw in bank.items():
                assert raw["ones"]+raw["zeros"]==row["design"]["shots_per_bank"]
                if stock!="blank":
                    i=list(FRAMES.values()).index(shape)
                    low,high=row["comparisons"][stock]["ten_entropy_intervals_over_ell4"][i]
                    assert low <= primary[stock]["entropy_over_ell4"] <= high
        held=row["comparisons"]["held"]
        assert held["theorem_component_bound"]<previous
        previous=held["theorem_component_bound"]
        assert row["comparisons"]["half_defect"]["stationarity_lower"]>5
        assert row["comparisons"]["blind_pulse"]["stationarity_lower"]>11
        # Population entropy has a finite-radius bias relative to the continuum action.
        errors=held["entropy_tensor_sampling_radius"]+(17/8)*held["bulk_upper"]/(2*math.pi*COEFFICIENT)
        assert abs(held["entropy_tensor_estimate"][0]-19/33) < errors
    assert previous < 0.0021


def test_tensor_reader_accepts_only_tallies_and_independent_area_and_stress_fields():
    recorded=json.loads((HERE/"results.json").read_text())["levels"][-1]
    cfg=recorded["design"]
    rows=json.loads((HERE.parent/"results.json").read_text())["finite_rows"]
    cols=("area_over_ell4","action_stress","metric_charge")
    rods={r["shape"]:{s:{k:t[k] for k in cols} for s,t in r["trials"].items()}
          for r in rows if r["inverse_ell"]==cfg["inverse_radius"]}
    counts={shape:{s:{k:v[k] for k in ("ones","zeros")} for s,v in bank.items()}
            for shape,bank in recorded["banks"].items()}
    result=infer(counts,cfg,rods)
    assert result["held"]["entropy_tensor_estimate"][0]==pytest.approx(0.57575633,abs=1e-8)
    # Retain the field, geometry, old counts and every calibration. Replace only
    # the held read by its blank control in every shape.
    for bank in counts.values():
        bank["held"]=bank["blank"]
    wrong=infer(counts,cfg,rods)
    assert wrong["held"]["entropy_tensor_estimate"]==[0.]*10
    assert wrong["held"]["stationarity_lower"]>27.
    assert wrong["held"]["area_tensor_estimate"]==result["held"]["area_tensor_estimate"]
