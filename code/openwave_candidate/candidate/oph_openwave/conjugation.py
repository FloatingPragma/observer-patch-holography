"""Complete local conjugation for the declared integer-load repair law.

The repaired pair balances x_a+x_b with a binary coin. Current j updates
classes (k+,k-) by (+j,-j) exactly when the declared cost gain is at most
5j(x_a-x_b-j). C(x_a,x_b,k+,k-)=(5-x_a,5-x_b,-k+,-k-) conjugates a branch
when its coin is complemented. This is class symmetry, not annihilation.
"""
from fractions import Fraction as F
from itertools import product
import json


COST = (0, 5, 8, 10, 8, 5)
STATES = tuple(product(range(6), repeat=4))
SOURCE_REFS = [
    "model/defect_integer_lift/model.py",
    "model/charge_conjugation_dynamics/symmetry.py",
    "Lean/Dynamics/ChargeConjugation.lean",
]


def conjugate(state):
    a, b, plus, minus = state
    return 5-a, 5-b, -plus % 6, -minus % 6


def step(state, coin, cost=COST):
    if len(state) != 4 or any(type(x) is not int or not 0 <= x < 6 for x in state):
        raise ValueError("four integer loads/classes in 0..5 required")
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError("binary coin required")
    if len(cost) != 6 or any(type(x) is not int for x in cost):
        raise ValueError("six integer class costs required")
    a, b, plus, minus = state
    new_a = (a+b+coin)//2
    new_b = a+b-new_a
    current = a-new_a
    proposed = ((plus+current) % 6, (minus-current) % 6)
    accepted = bool(current and cost[proposed[0]]+cost[proposed[1]]-cost[plus]-cost[minus]
                    <= 5*current*(a-b-current))
    target = (new_a, new_b, *(proposed if accepted else (plus, minus)))
    return target, (coin, current, accepted)


def branch_mismatches(cost=COST, complement_coin=True):
    failures = []
    for state in STATES:
        for coin in (0, 1):
            target, record = step(state, coin, cost)
            other, other_record = step(conjugate(state), 1-coin if complement_coin else coin, cost)
            expected_record = (1-coin, -record[1], record[2])
            if other != conjugate(target) or other_record != expected_record:
                failures.append({"state": list(state), "coin": coin,
                                 "target_mismatch": other != conjugate(target)})
    return failures


def kernel_row(state, coin_one=F(1, 2)):
    if not isinstance(coin_one, (int, F)) or not 0 <= coin_one <= 1:
        raise ValueError("exact rational coin probability in [0,1] required")
    row = {}
    for coin, probability in ((0, 1-coin_one), (1, coin_one)):
        target, _ = step(state, coin)
        row[target] = row.get(target, F(0)) + probability
    return row


def kernel_defect(coin_one):
    """Exact sup_x TV(K(Cx),C_*K(x)); targets retain all four coordinates."""
    maximum, witness = F(0), None
    for state in STATES:
        reflected = {conjugate(target): p for target, p in kernel_row(state, coin_one).items()}
        other = kernel_row(conjugate(state), coin_one)
        value = sum((abs(reflected.get(y, F(0))-other.get(y, F(0))) for y in reflected.keys() | other.keys()), F(0))/2
        if value > maximum:
            maximum, witness = value, list(state)
    return maximum, witness


def odd_class(state):
    return int(state[2] == 1)-int(state[2] == 5)


def paired_odd_mean(state, coin_one):
    return sum((p*odd_class(target)/2 for start in (state, conjugate(state))
                for target, p in kernel_row(start, coin_one).items()), F(0))


def run():
    failures = branch_mismatches()
    if failures:
        raise ArithmeticError("declared complete branch conjugation failed")
    accepted = no_op = rejected = 0
    for state in STATES:
        for coin in (0, 1):
            _, (_, current, applies) = step(state, coin)
            accepted += applies
            no_op += current == 0
            rejected += current != 0 and not applies
    fair_defect, _ = kernel_defect(F(1, 2))
    fair_pair_means = [paired_odd_mean(state, F(1, 2)) for state in STATES]
    if fair_defect or any(fair_pair_means):
        raise ArithmeticError("fair complete kernel lost class conjugation")
    cost_failures = branch_mismatches(cost=(0, 5, 8, 10, 8, 6))
    coin_failures = branch_mismatches(complement_coin=False)
    biased_defect, biased_witness = kernel_defect(F(3, 4))
    if (not any(row['target_mismatch'] for row in cost_failures)
            or not any(row['target_mismatch'] for row in coin_failures) or not biased_defect):
        raise ArithmeticError("a broken symmetry premise escaped its check")
    biased_mean = next(({"state": list(s), "odd_mean": str(paired_odd_mean(s, F(3, 4)))}
                        for s in STATES if paired_odd_mean(s, F(3, 4))), None)
    return {
        "schema": "oph.openwave.conjugation.v1", "computational_pass": True,
        "full_criterion_met": False, "criterion_rows": [8], "source_refs": SOURCE_REFS,
        "checks": {"complete_branch_conjugation": True, "fair_kernel_equivariant": True,
                   "all_1296_conjugate_preparation_pairs_have_zero_odd_mean": True,
                   "full_local_class_and_load_records_retained": True},
        "evidence": {"complete_local_branches": 2592, "accepted_twists": accepted,
                   "no_operations": no_op, "rejected_nonzero_twists": rejected,
                   "state_count": len(STATES), "fair_kernel_TV_defect_exact": str(fair_defect)},
        "mutations": {"asymmetric_cost_changes_actual_targets": True,
                      "uncomplemented_coin_changes_actual_targets": True,
                      "biased_coin_changes_kernel": True},
        "mutation_evidence": {"asymmetric_cost_detected_branches": len(cost_failures),
                      "first_asymmetric_cost_witness": cost_failures[0],
                      "coin_not_complemented_detected_branches": len(coin_failures),
                      "coin_not_complemented_changed_targets": sum(row['target_mismatch'] for row in coin_failures),
                      "biased_coin_probability": "3/4", "biased_kernel_TV_defect_exact": str(biased_defect),
                      "biased_kernel_witness": biased_witness, "biased_odd_class_mean_witness": biased_mean},
        "scope": {"source": "Declared classical integer-load repair with Z6 register and supplied cost table.",
                  "clock": "One complete local attempt; no physical duration.",
                  "ensemble_premises": "Conjugation-symmetric initial preparation, fair independent coins and equivariant schedule.",
                  "all_finite_times_statement": "The exact local intertwining supplies the premise of the cited finite-kernel Lean theorem; this script does not compile Lean.",
                  "particle_antiparticle_physical_dictionary_derived": False,
                  "annihilation_to_radiation_computed": False, "cosmological_baryogenesis_computed": False,
                  "federated_dynamics_reexecuted": False},
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
