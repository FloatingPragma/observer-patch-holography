import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from . import observability as obs
from .. import live_source, live_survey


def test_complete_front_rank_and_native_lottery_certificate():
    saved = json.loads(Path(obs.__file__).with_suffix(".json").read_text())
    result = obs.run()
    for field in ("contacts","rank","independent_links","rank_with_retained_modular"):
        assert result[field] == saved[field]
    assert sum(r["links"] for r in result["contacts"]) == 17496
    assert result["rank"] == 6 and result["rank_with_retained_modular"] == 7
    assert result["entropy_nats_per_unit_blind_direction_lower"] == .00053
    np.testing.assert_allclose(result["singular_values_unit_HS_per_distinct_centered_operator"],
                               saved["singular_values_unit_HS_per_distinct_centered_operator"],rtol=1e-10,atol=1e-13)
    assert result["front_minor_inverse_defect_upper"] < .001
    assert result["native_lotteries"]["bank_inverse_defect_upper"] < .001
    assert result["native_lotteries"]["retained_entropy_derivative_difference_lower"] > 0
    assert result["native_lotteries"]["finite_entropy_difference_lower"] > 1.7e-15
    control = result["native_lotteries"]["finite_weight_control"]
    for sign in ("positive","negative"):
        assert sum(control[sign+"_numerators"])+control[sign+"_reference_numerator"] == 2**40
    assert control["retained_entropy_derivative_difference_lower"] == 2e-8


def test_native_preparation_replays_equal_fronts_and_an_entropy_difference():
    saved = json.loads(Path(obs.__file__).with_suffix(".json").read_text())
    recipe = saved["native_lotteries"]
    board = obs.threshold.source.formation.federation(1)
    route = obs.threshold.source.formation.region(board)
    rho,repaired = obs.threshold.source.form(board,route,-1)
    generators = live_survey.rotations(rho)
    arm_states = []
    for name in recipe["quarter_turns_about_formed_marginals"]:
        if name == "independent_repairing":
            state = repaired
        else:
            left,right = int(name[1]),int(name[3])
            u = expm(-.5j*np.pi*(left*generators["left"]+right*generators["right"]))
            state = u@rho@u.conj().T
        arm_states.append(state)
    mixtures = []
    for sign in ("positive","negative"):
        weights = np.array(recipe[sign+"_center_probabilities"])
        base = recipe[sign+"_reference_center"]
        assert min(weights) >= 0 and base >= 0 and abs(base+weights.sum()-1) < 1e-14
        mixtures.append(base*rho+np.einsum("s,sij->ij",weights,arm_states))
    diff = mixtures[0]-mixtures[1]
    for m in obs.threshold.source.formation.marginals(diff):
        assert np.linalg.norm(m) < 1e-13
    output = live_source.propagate_heads(diff)
    ref = live_source.postpump(True)["reference"]
    eig,vec = np.linalg.eigh(ref)
    slope = -np.log(eig)@np.diag(vec.conj().T@output@vec).real
    assert slope > recipe["retained_entropy_derivative_difference_lower"]
    distinct, _, _ = obs.enumerate_fronts()
    for front in distinct:
        operator = np.array(front,float)
        operator /= np.linalg.norm(operator)
        assert abs(np.trace(operator@output)) < 5e-12
    control = recipe["finite_weight_control"]
    exact_weights = (np.array(control["positive_numerators"]) - np.array(control["negative_numerators"])) / 2**40
    finite = live_source.propagate_heads(np.einsum("s,sij->ij",exact_weights,np.array(arm_states)-rho))
    slope = -np.log(eig)@np.diag(vec.conj().T@finite@vec).real
    assert slope > control["retained_entropy_derivative_difference_lower"]
    for front in distinct:
        op = np.array(front,float)/float(obs.norm_inf(front))
        assert abs(np.trace(op@finite)) < control["centered_unit_row_sum_front_difference_upper"]
