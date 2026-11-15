"""Same-clock periodic C10 field refinement with all integer fluxes retained."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import math
import sys

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE = ROOT / "openwave-tuning"
sys.path.insert(0, str(BASE / "next/C10"))
from joint_source import JointNativeSource
from charged_reservoir import transport, FiberPreparation
sys.path.insert(0, str(BASE / "next/ROW1_NATIVE/routes/quantum_spatial"))
from quantum_spatial import lift_record, D as SPATIAL_D

MASS, KAPPA, G = .23, .37, .61
BETA = Fraction(361, 50)  # six half-edge hops, plus two face shifts
T = Fraction(37, 100)
RATE = Fraction(4, 5)
EDGES = ((0, 3), (1, 3), (1, 4), (2, 4), (2, 5), (0, 5))
SEAMS = ((0, 1), (1, 2), (2, 0))
DELTA = np.array([1, -1] * 3, dtype=int)
CHORD = 3


def centered(x, modulus):
    return (x + modulus // 2) % modulus - modulus // 2


def integer_record(mask, cycle):
    """Exact source tree (chord 3), including integer B=3-popcount(mask)."""
    b = 3 - int(mask).bit_count()
    res = [int(z < 3) - ((mask >> z) & 1) - int(z == 0) * b for z in range(6)]
    f = [0] * 6
    f[CHORD] = cycle
    a, v = EDGES[CHORD]
    res[a] -= cycle
    res[v] += cycle
    parent = {3: (0, 0), 5: (0, 5), 1: (3, 1), 2: (5, 4), 4: (1, 2)}
    for z in (4, 2, 1, 5, 3):
        p, edge = parent[z]
        a, _ = EDGES[edge]
        f[edge] = res[z] if a == z else -res[z]
        res[p] += res[z]
        res[z] = 0
    assert res[0] == 0
    return tuple(f)


def gauss(mask, f):
    out = np.array([(mask >> z) & 1 for z in range(6)])
    out[0] += 3 - int(mask).bit_count()
    for edge, (a, b) in enumerate(EDGES):
        out[a] += f[edge]
        out[b] -= f[edge]
    return tuple(out)


class Family:
    """Periodic source at modulus M, or integer Dirichlet diagnostic |k|<=K."""
    def __init__(self, modulus=None, cutoff=None):
        if (modulus is None) == (cutoff is None):
            raise ValueError("exactly one finite carrier required")
        if modulus is not None and (type(modulus) is not int or modulus < 6):
            raise ValueError("integer modulus at least six required")
        self.modulus, self.cutoff = modulus, cutoff
        self.cycles = tuple(range(-(modulus // 2), (modulus - 1) // 2 + 1)) if modulus else tuple(range(-cutoff, cutoff + 1))
        self.states = tuple((mask, k) for mask in range(64) for k in self.cycles)
        self.index = {s: i for i, s in enumerate(self.states)}
        self.n = len(self.states)
        self.records = np.array([integer_record(mask, k) for mask, k in self.states])
        self.display = np.array([[centered(int(x), modulus) for x in f] for f in self.records]) if modulus else self.records.copy()
        electric = .5 * KAPPA * np.sum(self.display**2, axis=1)
        mass = np.array([MASS * sum((1 if z < 3 else -1) * ((mask >> z) & 1) for z in range(6)) for mask, _ in self.states])
        self.electric = electric
        rows, cols, values = [], [], []
        hop_rows, hop_cols, hop_values = [], [], []
        for col, (mask, k) in enumerate(self.states):
            for h, (port, transit) in enumerate(EDGES):
                for src, dst, step in ((port, transit, 1), (transit, port, -1)):
                    moved = transport(mask, src, dst, 6)
                    if moved is None:
                        continue
                    target = self.resolve(moved[0], k + (step if h == CHORD else 0))
                    if target is not None:
                        hop_rows.append(target); hop_cols.append(col)
                        hop_values.append((1 if h % 2 == 0 else -1) * moved[1])
            for direction in (-1, 1):
                target = self.resolve(mask, k + direction)
                if target is not None:
                    rows.append(target); cols.append(col); values.append(-G)
        shape = (self.n, self.n)
        self.hop = sparse.csr_matrix((hop_values, (hop_rows, hop_cols)), shape=shape)
        self.magnetic = sparse.csr_matrix((values, (rows, cols)), shape=shape) + 2 * G * sparse.eye(self.n)
        self.field = sparse.diags(electric) + self.magnetic
        self.h = (self.hop + sparse.diags(mass) + self.field).tocsr()

    def resolve(self, mask, cycle):
        if self.modulus:
            cycle = centered(cycle, self.modulus)
        return self.index.get((mask, cycle))

    def lift(self, columns, source_states):
        out = np.zeros((self.n, columns.shape[1]), complex)
        for i, state in enumerate(source_states):
            out[self.index[state]] += columns[i]
        return out

    def evolve(self, columns, time=float(T)):
        return expm_multiply(-1j * time * self.h, columns)

    def instruments(self):
        """Complete MIX instrument, with every seam/law/coin/value label."""
        if self.modulus is None:
            raise ValueError("use a larger periodic member for terminal instrument diagnostics")
        groups = {}
        for col, (mask, k) in enumerate(self.states):
            for e, (a, b) in enumerate(SEAMS):
                contrast = ((mask >> a) & 1) - ((mask >> b) & 1)
                for law in ("C", "E"):
                    for coin in (0, 1):
                        current = (contrast - (2 * coin - 1) * abs(contrast)) // 2
                        target_mask, phase = mask, 1
                        if current:
                            src, dst = (a, b) if current == 1 else (b, a)
                            target_mask, phase = transport(mask, src, dst, 6)
                        target_k = k + (-current if 2 * e + 1 == CHORD else 0)
                        row = self.resolve(target_mask, target_k)
                        key = (e, law, coin, current if law == "C" else contrast)
                        groups.setdefault(key, ([], [], []))
                        groups[key][0].append(row); groups[key][1].append(col)
                        groups[key][2].append(phase / math.sqrt(12))
        return [(key, sparse.csr_matrix((v, (r, c)), shape=(self.n, self.n)))
                for key, (r, c, v) in sorted(groups.items())]


def tail_bound(x, first):
    """Rational majorant of 2 sum_{n>=first} x^n/n!, capped at TV=1."""
    if first <= 0 or x >= first + 1:
        return Fraction(1)
    value = 2 * x**first / math.factorial(first) / (1 - x / (first + 1))
    return min(Fraction(1), value)


def bound(modulus, preparation_radius, events=None, terminal_events=0):
    safe = (modulus - 1) // 2
    if events is None:
        # Poisson mixture: convolution of quantum Dyson order and event count.
        return tail_bound((BETA + RATE) * T, safe - preparation_radius - terminal_events + 1)
    return tail_bound(BETA * T, safe - preparation_radius - events - terminal_events + 1)


def outward(value):
    if value in (0, 1):
        return float(value)
    return float(np.nextafter(float(value), np.inf))


def receipt_bound(value):
    return {"exact_rational_upper": str(value), "outward_float_upper": outward(value)}


def admit_actual_input(hamiltonian, preparation, duration, model):
    """Admission regenerates source data, rather than trusting source labels."""
    rho, _, _ = model.reference()
    expected = sum(k @ rho @ k.T for _, k in model.source.reservoir_instrument())
    if (duration != T or hamiltonian.shape != (384, 384)
            or not np.array_equal(hamiltonian, model.hamiltonian)
            or preparation.shape != expected.shape
            or not np.allclose(preparation, expected, atol=1e-13, rtol=0)):
        raise ValueError("actual C10 generator, complete preparation, and common duration required")


def normdiff(left, left_family, right, right_family):
    keys = sorted(set(left_family.states) | set(right_family.states))
    idx = {s: i for i, s in enumerate(keys)}
    diff = np.zeros((len(keys), left.shape[1]), complex)
    for s, vector in zip(left_family.states, left):
        diff[idx[s]] += vector
    for s, vector in zip(right_family.states, right):
        diff[idx[s]] -= vector
    return float(np.linalg.norm(diff))


def execute():
    actual = JointNativeSource()
    prep_tree = FiberPreparation(actual.source.construction, 6)
    assert prep_tree.chords == (CHORD,)
    family6 = Family(6)
    order = [family6.index[(mask, centered(f[CHORD], 6))] for mask, f in actual.states]
    actual_h = actual.hamiltonian
    h_error = float(np.max(abs(family6.h.toarray()[np.ix_(order, order)] - actual_h)))
    # Exact preparation prescription is the actual global ground projector,
    # then all four source-reservoir instrument outcomes. Numerical vectors
    # implement its rank-one instance; the proof does not use rounded weights.
    rho, ground, multiplicity = actual.reference()
    assert multiplicity == 1
    vals, vecs = np.linalg.eigh(rho)
    psi = vecs[:, -1]
    original_columns = np.column_stack([k @ psi for _, k in actual.source.reservoir_instrument()])
    source_states = [(mask, centered(f[CHORD], 6)) for mask, f in actual.states]
    rho_check = original_columns @ original_columns.conj().T
    expected_rho = sum(k @ rho @ k.T for _, k in actual.source.reservoir_instrument())
    admit_actual_input(actual_h, expected_rho, T, actual)
    clock_spec = json.loads((HERE.parent / "clock/CLOCK_CONTRACT.json").read_text())
    with np.load(HERE.parent / "clock/common_clock_states.npz") as clock_states:
        clock_preparation_error = float(np.max(abs(expected_rho - clock_states["input_state"])))
        for meta, k in actual.source.reservoir_instrument():
            key = f"reservoir_{meta['reservoir_initial']}_{meta['reservoir_final']}_input"
            clock_preparation_error = max(clock_preparation_error, float(np.max(abs(k @ rho @ k.T - clock_states[key]))))
    h_hash = hashlib.sha256(np.asarray(actual_h, dtype="<f8").tobytes()).hexdigest()
    clock_match = (h_hash == clock_spec["source_H_array_sha256_le_float64"]
                   and float(T) == clock_spec["numerical_member"]["duration"]
                   and float(RATE) == clock_spec["numerical_member"]["native_attempt_rate"]
                   and 2 * float(RATE) == clock_spec["numerical_member"]["uniformization_rate"]
                   and clock_spec["numerical_member"]["law"] == "MIX")
    preparation_error = float(np.max(abs(rho_check - expected_rho)))
    integer_f = np.array([integer_record(mask, k) for mask, k in source_states])
    radii = np.max(abs(integer_f), axis=1)
    radius = int(radii.max())  # all 384 labels, independent of numerical zeros
    weights = np.sum(abs(original_columns)**2, axis=1)
    support = [{"radius": r, "outside_probability_diagnostic": float(weights[radii > r].sum()),
                "marked_clock_good_part_bound": receipt_bound(bound(6, r)),
                "full_preparation_comparison_formula": "min(1, delta_good + 2 sqrt(p_outside)); no conditioning or branch removal"}
               for r in range(radius + 1)]
    spatial_f = np.array([lift_record(list(f), 3 - mask.bit_count())["lifted_records"] for mask, f in actual.states])
    chart_cycle = integer_f - spatial_f
    chart_degree_error = int(np.max(abs(chart_cycle @ SPATIAL_D.T)))
    chart_residue_error = int(np.max(abs(chart_cycle % 6)))
    actual_ops = {tuple(meta[key] for key in ("seam", "law", "coin", "record_value")): k / math.sqrt(3)
                  for seam in actual.field.operators_by_law["MIX"] for meta, k in seam}
    operators6 = family6.instruments()
    instrument_error = max(float(np.max(abs(k.toarray()[np.ix_(order, order)] - actual_ops[key]))) for key, k in operators6)
    completeness = (sum(k.T @ k for _, k in operators6) - sparse.eye(family6.n)).tocoo()
    complete_error = float(np.max(abs(completeness.data))) if completeness.nnz else 0.
    # Integer chart mismatch is a measured obstruction on the actual full carrier.
    integer_energy = .5 * KAPPA * np.sum(integer_f**2, axis=1)
    periodic_energy = actual.field.electric_values
    chart_difference = integer_energy - periodic_energy
    wrap_mass = float(weights[abs(chart_difference) > 1e-12].sum())
    face_targets = np.argmax(actual.field.W, axis=0)
    coordinate_defect = int(np.max(abs(actual.field.representatives[face_targets] - actual.field.representatives - DELTA)))
    uniform_fiber_leak = float(np.max(periodic_energy.reshape(64, 6).std(axis=1)))
    reference = Family(cutoff=26)
    initial_ref = reference.lift(original_columns, source_states)
    final_ref = reference.evolve(initial_ref)
    ref_bound = tail_bound(BETA * T, 26 - 3 + 1)
    runs, arrays = [], {"actual_preparation_columns": original_columns,
                       "actual_source_masks": np.array([s[0] for s in actual.states]),
                       "actual_source_modular_records": np.array([s[1] for s in actual.states])}
    one_event = []
    for modulus in (6, 12, 18, 30, 42):
        model = family6 if modulus == 6 else Family(modulus)
        initial = model.lift(original_columns, source_states)
        final = model.evolve(initial)
        ops = model.instruments()
        rows = []
        for key, k in ops:
            branch = k @ final
            rows.append({"seam": key[0], "law": key[1], "coin": key[2], "record_value": key[3],
                         "probability": float(np.linalg.norm(branch)**2)})
        initial_field = float(np.vdot(initial, model.field @ initial).real)
        final_field = float(np.vdot(final, model.field @ final).real)
        mask_prob = np.array([np.sum(abs(final[i * modulus:(i + 1) * modulus])**2) for i in range(64)])
        runs.append({"modulus": modulus, "dimension": model.n,
                     "same_preparation_label_probabilities": [float(np.linalg.norm(v)**2) for v in initial.T],
                     "norm_diagnostic": float(np.linalg.norm(final)**2),
                     "passive_purified_vector_error_vs_integer_cutoff26": normdiff(final, model, final_ref, reference),
                     "passive_TV_theorem_bound": receipt_bound(bound(modulus, radius, events=0)),
                     "complete_marked_clock_TV_theorem_bound": receipt_bound(bound(modulus, radius)),
                     "marked_clock_plus_terminal_MIX_TV_bound": receipt_bound(bound(modulus, radius, terminal_events=1)),
                     "field_energy_change": final_field - initial_field,
                     "complete_terminal_occupation_probabilities": mask_prob.tolist(),
                     "complete_terminal_MIX_records": rows,
                     "terminal_probability_sum": sum(row["probability"] for row in rows)})
        if modulus in (6, 42):
            arrays[f"M{modulus}_passive_final_columns"] = final
            # One complete native MIX event at the fixed absolute time .17.
            # Every one of its 30 marks and four reservoir labels is retained.
            mid = model.evolve(initial, .17)
            marked = np.column_stack([k @ mid for _, k in ops])
            continued = model.evolve(marked, .20)
            one_event.append({"modulus": modulus, "attempt_time": .17,
                              "continuation_duration": .20, "number_of_native_marks": len(ops),
                              "column_order": "native mark outer, four reservoir labels inner, both in source iteration order",
                              "native_marks": [{"seam": key[0], "law": key[1], "coin": key[2], "record_value": key[3]} for key, _ in ops],
                              "retained_branch_probabilities": [float(np.linalg.norm(continued[:, j])**2) for j in range(continued.shape[1])],
                              "probability_sum": float(np.linalg.norm(continued)**2),
                              "same_clock_one_event_TV_bound": receipt_bound(bound(modulus, radius, events=1))})
            arrays[f"M{modulus}_one_event_final_columns"] = continued
    # Entire C10 source contains components with a nonzero commutator:
    # field and matter exchange energy in the same autonomous evolution.
    commutator = actual.field.field @ actual.field.source - actual.field.source @ actual.field.field
    # Removing the electric term leaves the original hopping fiber invariant,
    # but it changes this actual source and its prepared response.
    bad_h = family6.h - sparse.diags(family6.electric)
    initial6 = family6.lift(original_columns, source_states)
    bad_final = expm_multiply(-1j * float(T) * bad_h, initial6)
    good_final = arrays["M6_passive_final_columns"]
    mutation_distance = float(np.linalg.norm(bad_final - good_final))
    rejected = {}
    for name, h, state, duration in (
            ("electric_term_removed", actual_h - actual.field.electric, expected_rho, T),
            ("reference_substituted_for_reservoir_output", actual_h, rho, T),
            ("clock_duration_changed", actual_h, expected_rho, T * 2)):
        try:
            admit_actual_input(h, state, duration, actual)
            rejected[name] = False
        except ValueError:
            rejected[name] = True
    checks = {
        "actual384_H_recovered": h_error < 1e-14,
        "actual_complete_MIX_recovered": instrument_error < 1e-14,
        "complete_MIX_isometry": complete_error < 1e-14,
        "actual_global_reference_reservoir_preparation_recovered": preparation_error < 1e-13,
        "actual_common_clock_H_preparation_and_parameters_recovered": clock_match and clock_preparation_error < 1e-12,
        "exact_integer_Gauss_all_retained_labels": all(gauss(mask, integer_record(mask, k)) == (1, 1, 1, 0, 0, 0) for mask, k in reference.states),
        "actual_spatial_and_boundary_degrees_preserved": chart_degree_error == 0 and chart_residue_error == 0,
        "all_terminal_probabilities_retained": all(abs(row["terminal_probability_sum"] - 1) < 1e-11 for row in runs),
        "complete_one_event_histories_retained": all(abs(row["probability_sum"] - 1) < 1e-11 and len(row["retained_branch_probabilities"]) == 120 for row in one_event),
        "interacting_field_matter_commutator_nonzero": float(np.linalg.norm(commutator)) > 1,
        "actual_M6_integer_chart_energy_identification_rejected": wrap_mass > 0,
        "drop_electric_generator_mutation_rejected": mutation_distance > 1e-3,
        "source_preparation_clock_mutations_rejected": all(rejected.values()),
        "full_carrier_coordinate_commutator_transfer_rejected": coordinate_defect == 6,
        "uniform_fiber_invariance_transfer_rejected": uniform_fiber_leak > 1,
        "all_passive_errors_below_proved_bounds": all(row["passive_purified_vector_error_vs_integer_cutoff26"] <= (2. if row["passive_TV_theorem_bound"]["exact_rational_upper"] == "1" else row["passive_TV_theorem_bound"]["outward_float_upper"]) + outward(ref_bound) + 1e-10 for row in runs),
        "complete_marked_clock_refinement_nonvacuous": float(bound(42, radius)) < 1e-4,
    }
    result = {"schema": "oph.c10.dynamical_integer_field.v1", "checks": checks,
              "parameters": {"mass": MASS, "kappa": KAPPA, "face_coefficient": G,
                             "duration": str(T), "native_attempt_rate": str(RATE), "uniformization_rate": "8/5",
                             "bounded_hopping_norm_upper": str(BETA), "initial_integer_record_radius_all384": radius},
              "baseline": {"dimension": actual.n, "H_max_difference": h_error,
                           "H_sha256_le_float64": h_hash,
                           "instrument_max_difference": instrument_error, "instrument_completeness_error": complete_error,
                           "preparation_max_difference": preparation_error, "global_reference_energy": ground,
                           "common_clock_preparation_archive_max_difference": clock_preparation_error,
                           "common_clock_H_duration_and_rates_match": clock_match},
              "preparation_outcomes": [{**meta, "probability": float(np.linalg.norm(original_columns[:, j])**2)}
                                       for j, (meta, _) in enumerate(actual.source.reservoir_instrument())],
              "actual_preparation_chart": {"exact_lift": "source tree chord3; k=centered(Fmod6[3]); exact integer B=3-N",
                                           "spatial_chord5_chart_difference_max": int(abs(chart_cycle).max()),
                                           "spatial_chart_degree_difference": chart_degree_error,
                                           "spatial_chart_residue_difference": chart_residue_error,
                                           "outside_support_diagnostics": support,
                                           "periodic_vs_integer_electric_energy_max_difference": float(abs(chart_difference).max()),
                                           "periodic_vs_integer_electric_energy_preparation_mean_difference": float(weights @ chart_difference),
                                           "nonzero_energy_chart_defect_preparation_probability": wrap_mass,
                                           "M6_full_preparation_complete_clock_bound": receipt_bound(bound(6, radius))},
              "integer_diagnostic_cutoff": {"cycle_cutoff": 26, "dimension": reference.n,
                                            "passive_error_to_infinite_TV_bound": receipt_bound(ref_bound),
                                            "boundary_handling": "Hermitian Dirichlet cut only in diagnostic; omitted infinite amplitude bounded, never renormalized"},
              "members": runs,
              "complete_one_event_executions": one_event,
              "obstructions": {"field_matter_commutator_Frobenius": float(np.linalg.norm(commutator)),
                               "drop_electric_mutation_vector_distance": mutation_distance,
                               "admission_mutations_rejected": rejected,
                               "periodic_coordinate_commutator_wrap_defect_operator_norm": coordinate_defect,
                               "uniform_fiber_electric_leak_operator_norm": uniform_fiber_leak,
                               "spatial_cycle_rank": 1,
                               "integer_energy_minimization_used": False},
              "scope": "Actual C10 M6 source and full marked-clock preparation extend by a declared modulus family to an interacting integer-flux field on the same finite triangle. All four reservoir labels, native no-current outcomes, integer boundary B, and Poisson count tails retained. No oscillator, spatial Maxwell, physical clock, or particle-force identification.",
              "official_rows_changed": []}
    return result, arrays


if __name__ == "__main__":
    result, arrays = execute()
    np.savez_compressed(HERE / "dynamical_field_states.npz", **arrays)
    (HERE / "RESULT.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"checks": result["checks"], "radius": result["parameters"]["initial_integer_record_radius_all384"],
                      "M42_clock_bound": result["members"][-1]["complete_marked_clock_TV_theorem_bound"]}, indent=2))
