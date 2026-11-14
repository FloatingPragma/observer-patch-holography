#!/usr/bin/env python3
"""A bounded current-record observer with an internal, fixed feedback policy.

All arithmetic is exact. The endpoint rule, tape and controller define a
finite executable state machine and its quantum-channel construction.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as Q
import hashlib
from itertools import combinations_with_replacement, product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
E = (0, 1)
F = (0, 5)
LETTERS = tuple(product((0, 1), (-1, 0, 1)))
# Fixed controller data for uniform binary total-two preparation. These values
# are issued before event 2; the independent enumerator checks calibration.
FORECAST = {(0, 0): Q(0), (1, 0): Q(5, 56),
            (0, 1): Q(0), (1, -1): Q(1, 2)}


def states(binary=False):
    """All twelve-port nonnegative integer configurations of total two.

    Lexicographic order fixes receipt indices. Binary preparation retains the
    66 configurations occupying two distinct ports; the whole carrier has 78.
    """
    configurations = []
    for a, b in combinations_with_replacement(range(12), 2):
        if binary and a == b:
            continue
        loads = [0] * 12
        loads[a] += 1
        loads[b] += 1
        configurations.append(tuple(loads))
    return tuple(sorted(configurations))


def native(loads, seam, coin):
    a, b = seam
    d = loads[a] - loads[b]
    current = (d - (2 * coin - 1) * (abs(d) % 2)) // 2
    out = list(loads)
    out[a] -= current
    out[b] += current
    return tuple(out), (coin, current)


@dataclass(frozen=True)
class State:
    loads: tuple[int, ...]
    prefix: tuple[tuple[int, int], ...] = ()

    @property
    def phase(self):
        return len(self.prefix)

    def checkpoint(self):
        """Simulator full-state snapshot; hidden loads are not observer inputs."""
        return {"phase": self.phase, "prefix": [list(r) for r in self.prefix],
                "loads": list(self.loads)}

    def controller_checkpoint(self):
        """The observer-visible controller exports no working loads."""
        return {"phase": self.phase, "prefix": [list(r) for r in self.prefix]}

    @staticmethod
    def restore_controller(checkpoint, unchanged_loads):
        """Reattach the restored controller to its unchanged physical source."""
        return State.restore(dict(checkpoint, loads=list(unchanged_loads)))

    @staticmethod
    def restore(checkpoint):
        state = State(tuple(checkpoint["loads"]),
                      tuple(tuple(r) for r in checkpoint["prefix"]))
        if (state.phase != checkpoint["phase"] or state.phase > 2
                or len(state.loads) != 12 or sum(state.loads) != 2
                or any(type(x) is not int or x < 0 for x in state.loads)
                or any(r not in LETTERS for r in state.prefix)):
            raise ValueError("Invalid complete checkpoint")
        return state


def selected_seam(state):
    """The immutable policy reads the phase and committed prefix only."""
    if state.phase == 0:
        return E
    if state.phase == 1:
        return F if state.prefix[0][1] > 0 else E
    return None


def issue_forecast(state):
    """Read a supplied controller table before the next event; never loads."""
    if state.phase != 1:
        raise ValueError("A forecast is issued only between the two events")
    return FORECAST.get(state.prefix[0])  # abstain off the declared prior support


def tick(state):
    """The same fixed rule at every invocation, including halt-on-full."""
    seam = selected_seam(state)
    if seam is None:
        return ((Q(1), state),)
    return tuple((Q(1, 2), State(out, state.prefix + (record,)))
                 for coin in (0, 1)
                 for out, record in (native(state.loads, seam, coin),))


def record_marginal(paths):
    marginal = defaultdict(Q)
    for weight, first, _, _ in paths:
        marginal[first.prefix[0]] += weight
    return dict(marginal)


def observer_step(state):
    """Issue the fixed controller's classical forecast before advancing.

    On quantum inputs its readout is a diagonal-record measurement, distinct
    from the coherent tick channel. The benchmark uses classical histories.
    """
    forecast = issue_forecast(state) if state.phase == 1 else None
    return forecast, tick(state)


def prefix_density(terms, coefficient, after_second=False):
    """Exact reduced first-letter density from equal-amplitude pure terms."""
    density = defaultdict(Q)
    for left in terms:
        for right in terms:
            if left.loads == right.loads and (
                    not after_second or left.prefix[1] == right.prefix[1]):
                density[left.prefix[0], right.prefix[0]] += coefficient
    return dict(density)


def run():
    configurations = states()
    index = {x: i for i, x in enumerate(configurations)}
    binary = [x for x in configurations if max(x) == 1]
    checks = {}

    def check(name, assertion):
        checks[name] = bool(assertion)
        if not assertion:
            raise AssertionError(name)

    check("state_space_78_and_binary_shell_66", len(configurations) == 78 and len(binary) == 66)
    # Independent endpoint characterization of the declared local repair.
    # On the selected seam compatibility means endpoint difference <= 1;
    # floor/ceiling averaging realizes it, conserving total load.
    branches = 0
    for seam in (E, F):
        a, b = seam
        for x in configurations:
            for coin in (0, 1):
                y, record = native(x, seam, coin)
                total = x[a] + x[b]
                low, high = total // 2, (total + 1) // 2
                assert (y[a], y[b]) == ((low, high) if coin == 0 else (high, low))
                assert all(y[k] == x[k] for k in range(12) if k not in seam)
                assert record[1] == x[a] - y[a]
                branches += 1
    check("all_312_native_branches_match_independent_endpoint_repair", branches == 312)

    # For fixed p, every column of V_p has two equal amplitudes 1/sqrt(2).
    # Disjoint column supports and norm one certify V_p^* V_p=I exactly,
    # including arbitrary unused record letters and all off-diagonal inputs.
    isometry_counts = {}
    for phase, prefixes in ((0, [()]), (1, [(r,) for r in LETTERS])):
        outputs = set()
        columns = 0
        for prefix in prefixes:
            for x in configurations:
                before = State(x, prefix)
                column = tick(before)
                check_outputs = {s for _, s in column}
                assert len(check_outputs) == 2 and not outputs.intersection(check_outputs)
                assert sum(w for w, _ in column) == 1
                for _, after in column:
                    assert after.phase == phase + 1
                    assert after.prefix[:phase] == prefix
                    assert sum(after.loads) == 2
                    assert sum(v*v for v in after.loads) <= sum(v*v for v in x)
                    assert all(after.loads[k] == x[k] for k in range(12) if k not in {0, 1, 5})
                outputs.update(check_outputs)
                columns += 1
        isometry_counts[str(phase)] = {"columns": columns, "nonzero_entries": len(outputs)}
    check("all_nonfull_columns_are_isometric_including_unreachable_inputs", True)
    check("load_conservation_imbalance_nonincrease_and_spectator_preservation", True)
    check("committed_classical_prefix_cylinders_preserved", True)
    halt_count = 0
    for prefix in product(LETTERS, repeat=2):
        for x in configurations:
            s = State(x, prefix)
            assert tick(s) == ((Q(1), s),)
            halt_count += 1
    check("all_2808_full_basis_states_halt", halt_count == 2808)

    # Full two-event coherent map is also an isometry on all 78 columns.
    complete_outputs = set()
    for x in configurations:
        column = [last for _, first in tick(State(x)) for _, last in tick(first)]
        assert len(set(column)) == 4 and not complete_outputs.intersection(column)
        complete_outputs.update(column)
        for last in column:
            recovered = list(last.loads)
            for event in (1, 0):
                seam = E if event == 0 else (F if last.prefix[0][1] > 0 else E)
                current = last.prefix[event][1]
                recovered[seam[0]] += current
                recovered[seam[1]] -= current
            assert tuple(recovered) == x
    check("all_312_two_event_paths_invert_and_columns_are_orthogonal", True)

    # All 264 main classical histories, with a forecast issued BEFORE tick 2.
    paths = []
    histories = []
    for x in binary:
        for first_weight, first in tick(State(x)):
            forecast, continuations = observer_step(first)
            assert forecast is not None
            restored = State.restore(json.loads(json.dumps(first.checkpoint())))
            assert tick(restored) == tick(first) and issue_forecast(restored) == forecast
            local_checkpoint = first.controller_checkpoint()
            assert set(local_checkpoint) == {"phase", "prefix"}
            restored_local = State.restore_controller(local_checkpoint, first.loads)
            assert tick(restored_local) == tick(first)
            for second_weight, last in continuations:
                weight = first_weight * second_weight / len(binary)
                paths.append((weight, first, last, forecast))
                histories.append({"initial_state": index[x], "first_record": list(first.prefix[0]),
                                  "issued_before_event_2": str(forecast),
                                  "second_seam": list(selected_seam(first)),
                                  "second_record": list(last.prefix[1]),
                                  "final_state": index[last.loads], "probability": str(weight)})
    check("all_main_paths_retained_without_filtering", len(paths) == 264 and sum(p[0] for p in paths) == 1)
    check("complete_checkpoint_preserves_conditional_future_and_forecast", True)
    check("controller_checkpoint_contains_no_loads_and_restores_with_unchanged_source", True)
    marginal = record_marginal(paths)
    successes = defaultdict(Q)
    for weight, first, last, _ in paths:
        successes[first.prefix[0]] += weight * (last.prefix[1][1] > 0)
    check("preissued_forecasts_calibrated_on_every_positive_record",
          all(FORECAST[r] == successes[r] / marginal[r] for r in marginal))
    event_rate = sum(successes.values())
    risk = sum(w * (forecast - (last.prefix[1][1] > 0)) ** 2
               for w, _, last, forecast in paths)
    prior_risk = event_rate * (1 - event_rate)
    shuffled_risk = sum(w * mr * (FORECAST[r] - (last.prefix[1][1] > 0)) ** 2
                        for w, _, last, _ in paths for r, mr in marginal.items())
    action_mass = defaultdict(Q)
    action_success = defaultdict(Q)
    for w, first, last, _ in paths:
        seam = selected_seam(first)
        action_mass[seam] += w
        action_success[seam] += w * (last.prefix[1][1] > 0)
    action_prior_risk = sum(w * (action_success[selected_seam(first)]
                               / action_mass[selected_seam(first)]
                               - (last.prefix[1][1] > 0)) ** 2
                            for w, first, last, _ in paths)
    action_shuffle_risk = sum(
        w * mr / action_mass[selected_seam(first)]
        * (FORECAST[r] - (last.prefix[1][1] > 0)) ** 2
        for w, first, last, _ in paths for r, mr in marginal.items()
        if (F if r[1] > 0 else E) == selected_seam(first))
    check("positive_predictive_gain_against_prior_and_matched_independent_shuffle",
          risk == Q(395, 7392) and prior_risk == Q(305, 4356)
          and shuffled_risk == Q(21125, 243936) and risk < prior_risk < shuffled_risk)
    check("record_predicts_beyond_knowledge_of_supplied_action",
          risk < action_prior_risk < action_shuffle_risk
          and action_shuffle_risk - risk == 2 * (action_prior_risk - risk))
    occupation = sum(w * last.loads[0] for w, _, last, _ in paths)
    true_f_mass = sum(w for w, first, _, _ in paths if selected_seam(first) == F)

    # Controls alter the selector connection, while actual native records remain.
    # Shuffled selector reads an independent letter with the same marginal as R1.
    control = {}
    for mode in ("record_disconnected", "independent_matched_record"):
        mass = Q(0)
        mean = Q(0)
        f_mass = Q(0)
        branches = 0
        for x in binary:
            for p1, first in tick(State(x)):
                selectors = [(None, Q(1))] if mode == "record_disconnected" else list(marginal.items())
                for replacement, pr in selectors:
                    seam = F if replacement is not None and replacement[1] > 0 else E
                    for coin in (0, 1):
                        y, _ = native(first.loads, seam, coin)
                        w = p1 * pr / (2 * len(binary))
                        mass += w
                        mean += w * y[0]
                        f_mass += w * (seam == F)
                        branches += 1
        control[mode] = {"mass": str(mass), "final_occupation_0": str(mean),
                         "second_seam_f_mass": str(f_mass), "branches": branches}
        assert mass == 1 and mean == Q(1, 6)
    check("record_feedback_changes_the_occupation_audit", occupation == Q(35, 264))
    check("matched_selector_has_identical_seam_frequency",
          true_f_mass == Q(5, 66) == Q(control["independent_matched_record"]["second_seam_f_mass"]))
    check("all_control_branches_retained", True)

    # A reachable coherence explicitly refutes full-quantum-prefix protection.
    a = tuple(int(k in (0, 2)) for k in range(12))
    b = tuple(int(k in (1, 2)) for k in range(12))
    first_terms = [first for x in (a, b) for _, first in tick(State(x))]
    second_terms = [last for first in first_terms for _, last in tick(first)]
    before = prefix_density(first_terms, Q(1, 4))
    after = prefix_density(second_terms, Q(1, 8), after_second=True)
    witness = ((0, 1), (0, 0))
    check("reachable_quantum_prefix_coherence_changes",
          before[witness] == Q(1, 4) and after[witness] == Q(1, 8))
    check("same_coherent_input_keeps_classical_record_diagonal",
          all(before.get((r, r), 0) == after.get((r, r), 0) for r in LETTERS))

    return {"schema": "oph.native_observer.v1", "all_pass": all(checks.values()),
            "checks": checks, "state_space": [list(x) for x in configurations],
            "architecture": {"candidate": "N1 current C, not N1-E contrast E or N4",
                             "record_letters": [list(r) for r in LETTERS],
                             "source_dimension": 78, "record_dimension": 6,
                             "capacity": 2, "total_hilbert_dimension": 3354,
                             "isometry_columns": isometry_counts, "halt_basis_states": halt_count,
                             "protected_algebra": "committed classical first-prefix cylinders at phase >= 1",
                             "pointer_parameter": "two local patch events; not physical duration"},
            "forecast": {"target": "second native current j > 0", "issued_at": "phase 1 before event 2",
                         "table": [{"record": list(r), "record_probability": str(marginal[r]),
                                    "forecast": str(FORECAST[r])} for r in sorted(marginal)],
                         "off_prior_letters": "abstain", "event_rate": str(event_rate),
                         "brier_risk": str(risk), "prior_only_brier_risk": str(prior_risk),
                         "independent_shuffle_brier_risk": str(shuffled_risk),
                         "action_known_prior_brier_risk": str(action_prior_risk),
                         "within_action_shuffle_brier_risk": str(action_shuffle_risk),
                         "gain_over_action_known_prior": str(action_prior_risk - risk),
                         "gain_over_within_action_shuffle": str(action_shuffle_risk - risk),
                         "gain_over_prior": str(prior_risk - risk),
                         "gain_over_shuffle": str(shuffled_risk - risk)},
            "feedback_audit": {"actual_final_occupation_0": str(occupation),
                               "actual_second_seam_f_mass": str(true_f_mass),
                               "difference_from_controls": str(occupation - Q(1, 6)),
                               "controls": control,
                               "scope": "load expectation is an external mathematical audit, not a native meter"},
            "quantum_prefix_counterexample": {"initial_equal_superposition_states": [index[a], index[b]],
                                              "record_matrix_unit": [list(r) for r in witness],
                                              "before": str(before[witness]), "after": str(after[witness])},
            "main_histories": histories,
            "scope": "Finite observer functions with explicit tape, prior, forecast table, selector and halt. Constructed software with no learning, physical hardware execution, calibrated physical clock or full observer algebra identification claimed."}


def canonical(result):
    return json.dumps(result, indent=2, sort_keys=True) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path, help="exclusively create a new deterministic result")
    group.add_argument("--check", type=Path, help="compare against exact canonical result bytes")
    args = parser.parse_args()
    result = run()
    encoded = canonical(result)
    if args.write:
        with args.write.open("x") as stream:
            stream.write(encoded)
    if args.check and args.check.read_text() != encoded:
        raise SystemExit("Result mismatch")
    print(json.dumps({"all_pass": result["all_pass"], "checks": len(result["checks"]),
                      "result_sha256": hashlib.sha256(encoded.encode()).hexdigest()}))
