"""Marked Poisson clock on the complete interacting C10 carrier.

Hamiltonian durations and event rate are supplied apparatus parameters.
Every native attempt, including zero current, costs one native commit pulse.
No native pulse compiler for exp(-itH) is asserted.
"""
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
C10 = ROOT / "openwave-tuning/next/C10"
sys.path.insert(0, str(C10))
from joint_source import JointNativeSource


def applied_current(meta):
    value, coin = meta["record_value"], meta["coin"]
    return value if meta["law"] == "C" else (value - (2 * coin - 1) * abs(value)) // 2


class CommonClock:
    """Exact instrument and sparse number-sector semigroup adapters.

    Number-sector compression is used only for the declared block-diagonal
    preparation. trajectory() and native() act on the entire 384-state space.
    """

    def __init__(self, law="MIX"):
        if law not in ("C", "E", "MIX"):
            raise ValueError("complete C, E or MIX law required")
        self.model = JointNativeSource()
        self.law = law
        self.h = self.model.hamiltonian
        self.n = self.model.n
        self.ops = [
            ({**meta, "applied_current": applied_current(meta),
              "native_attempts": 1, "native_commit_pulses": 1,
              "boundary_rejected": False}, sparse.csr_matrix(k / math.sqrt(3)))
            for seam in self.model.field.operators_by_law[law]
            for meta, k in seam
        ]
        self.sectors = {}
        for q, (ids, energies, vectors) in self.model.spectral.items():
            h = sparse.csr_matrix(self.h[np.ix_(ids, ids)])
            eye = sparse.eye(len(ids), format="csr")
            # Column-major vectorization: vec(A X B)=(B.T tensor A)vec(X).
            coherent = -1j * (sparse.kron(eye, h) - sparse.kron(h.T, eye))
            jump = sparse.csr_matrix((len(ids)**2, len(ids)**2), dtype=complex)
            move = jump.copy()
            for meta, k in self.ops:
                local = k[ids][:, ids]
                term = sparse.kron(local.conjugate(), local, format="csr")
                jump += term
                if meta["applied_current"]:
                    move += term
            self.sectors[q] = (ids, energies, vectors, coherent.tocsr(), jump, move)

    def preparation(self):
        """Recompute all four reservoir branches from the actual global reference."""
        reference, _, _ = self.model.reference()
        branches = [(meta, k @ reference @ k.T)
                    for meta, k in self.model.source.reservoir_instrument()]
        return sum(state for _, state in branches).astype(complex), branches

    def evolve(self, rho, duration, h_scale=1.0):
        """Exact full-space unitary service, including cross-number coherences."""
        if not math.isfinite(duration) or duration < 0:
            raise ValueError("finite nonnegative supplied duration required")
        if not math.isfinite(h_scale) or h_scale < 0:
            raise ValueError("finite nonnegative Hamiltonian scale required")
        u = np.zeros((self.n, self.n), complex)
        for ids, energies, vectors, *_ in self.sectors.values():
            u[np.ix_(ids, ids)] = (vectors * np.exp(-1j * h_scale * duration * energies)) @ vectors.T
        return u @ rho @ u.conj().T

    def native(self, rho, detailed=False):
        out = np.zeros((self.n, self.n), complex)
        rows = []
        for meta, k in self.ops:
            branch = (k @ rho) @ k.T
            out += branch
            if detailed:
                rows.append({**meta, "probability": float(np.trace(branch).real),
                             "unnormalized_final_energy": float(np.trace(self.h @ branch).real)})
        return (out, rows) if detailed else out

    def require_block_diagonal(self, rho):
        off = np.array(rho, copy=True)
        for ids, *_ in self.sectors.values():
            off[np.ix_(ids, ids)] = 0
        if np.max(abs(off)) > 1e-13:
            raise ValueError("semigroup compression requires number-block-diagonal input")

    def semigroup(self, rho, duration, rate, *, attempt_mark=1.0,
                  move_mark=1.0, idle_rate=0.0, idle_mark=1.0, h_scale=1.0):
        """E[z^N x^M y^I rho_record], with every no-current attempt included.

        N is native attempt count, M the nonzero-current count, and I the
        independent uniformization idle count. For marks one this is CPTP.
        Marks in [0,1] give positive, subnormalized generating densities.
        """
        if not all(math.isfinite(v) and v >= 0 for v in (duration, rate, idle_rate, h_scale)):
            raise ValueError("nonnegative finite clock parameters required")
        if not all(math.isfinite(v) and 0 <= v <= 1 for v in (attempt_mark, move_mark, idle_mark)):
            raise ValueError("positive generating marks in [0,1] required")
        self.require_block_diagonal(rho)
        out = np.zeros_like(rho, dtype=complex)
        for ids, _, _, coherent, jump, move in self.sectors.values():
            state = rho[np.ix_(ids, ids)]
            if not np.any(state):
                continue
            eye = sparse.eye(len(ids)**2, format="csr")
            gen = (h_scale * coherent + rate * attempt_mark * (jump + (move_mark - 1) * move)
                   + (-rate + idle_rate * (idle_mark - 1)) * eye)
            vec = expm_multiply(duration * gen, state.reshape(-1, order="F"))
            out[np.ix_(ids, ids)] = vec.reshape(state.shape, order="F")
        return out

    def count_split(self, rho, duration, rate):
        """Retain exact N=0, N=1 and N>=2 density maps with no renormalization.

        The final class is a complete infinite-history complement, not a
        runtime abort or an assertion that counts beyond one do not occur.
        """
        self.require_block_diagonal(rho)
        one = np.zeros_like(rho, dtype=complex)
        for ids, _, _, coherent, jump, _ in self.sectors.values():
            state = rho[np.ix_(ids, ids)]
            if not np.any(state):
                continue
            size = len(ids)**2
            gen0 = coherent - rate * sparse.eye(size, format="csr")
            gen = sparse.bmat([[gen0, None], [rate * jump, gen0]], format="csr")
            vec = np.r_[state.reshape(-1, order="F"), np.zeros(size, complex)]
            result = expm_multiply(duration * gen, vec)
            one[np.ix_(ids, ids)] = result[size:].reshape(state.shape, order="F")
        zero = math.exp(-rate * duration) * self.evolve(rho, duration)
        total = self.semigroup(rho, duration, rate)
        return zero, one, total - zero - one

    def trajectory(self, rho, duration, rate, uniform_rate, events):
        """Execute one specified marked path; keep all alternatives at each node.

        events are (absolute_time, 'idle' or native operator index). The
        path density is with respect to ordered event times. This is one
        path witness, not a replacement for the complete semigroup output.
        """
        if not (0 < rate <= uniform_rate and math.isfinite(uniform_rate)):
            raise ValueError("0 < native rate <= finite uniformization rate required")
        state = np.asarray(rho, complex).copy()
        if abs(np.trace(state) - 1) > 1e-12:
            raise ValueError("normalized input required")
        records, segments = [], []
        time, native_count, idle_count, moving_count = 0.0, 0, 0, 0
        selected_weight = 1.0
        work = 0.0
        for at, mark in events:
            if not (time < at < duration):
                raise ValueError("strictly ordered interior event times required")
            delta = at - time
            state = self.evolve(state, delta)
            segments.append(delta)
            _, native_rows = self.native(state, detailed=True)
            alternatives = [{**row, "mark": i, "probability": rate / uniform_rate * row["probability"]}
                            for i, row in enumerate(native_rows)]
            alternatives.append({"mark": "idle", "probability": 1 - rate / uniform_rate,
                                 "native_attempts": 0, "native_commit_pulses": 0})
            before = float(np.trace(self.h @ state).real)
            if mark == "idle":
                probability = 1 - rate / uniform_rate
                idle_count += 1
            else:
                if type(mark) is not int or not 0 <= mark < len(self.ops):
                    raise ValueError("native mark index outside complete instrument")
                meta, k = self.ops[mark]
                branch = (k @ state) @ k.T
                native_probability = float(np.trace(branch).real)
                probability = rate / uniform_rate * native_probability
                if native_probability <= 0:
                    raise ValueError("chosen zero-probability path has no conditional state")
                state = branch / native_probability
                native_count += 1
                moving_count += int(meta["applied_current"] != 0)
            if probability <= 0:
                raise ValueError("chosen path has zero density")
            after = float(np.trace(self.h @ state).real)
            work += after - before
            selected_weight *= probability
            records.append({"at": at, "selected_mark": mark,
                            "selected_conditional_probability": probability,
                            "complete_alternatives": alternatives,
                            "conditional_energy_change": after - before})
            time = at
        segments.append(duration - time)
        state = self.evolve(state, duration - time)
        return state, {
            "events": records, "coherent_segment_durations": segments,
            "supplied_coherent_duration": sum(segments),
            "native_attempts": native_count, "native_zero_current_attempts": native_count - moving_count,
            "native_moving_attempts": moving_count, "uniformization_idles": idle_count,
            "native_commit_pulses": native_count, "coherent_service_calls": len(segments),
            "coherent_service_native_pulse_compilation": None,
            "conditional_energy_change_sum": work,
            "path_density_wrt_ordered_times": math.exp(-uniform_rate * duration) * uniform_rate**len(events) * selected_weight,
            "complete_branching": all(abs(sum(x["probability"] for x in r["complete_alternatives"]) - 1) < 1e-12 for r in records),
        }


def source_pins():
    paths = [Path(__file__), C10 / "joint_source.py", C10 / "modular_field_adapter.py",
             C10 / "charged_reservoir.py", C10 / "joint_source_states.npz",
             ROOT / "openwave-tuning/next/C9/fermionic_shared_records.py",
             ROOT / "openwave-tuning/next/C9/fermionic_uniform_reference.py",
             ROOT / "openwave-tuning/validation/bridges/providers/fock.py",
             ROOT / "openwave-tuning/candidate/oph_openwave/construction.py",
             ROOT / "openwave-tuning/candidate/oph_openwave/fiber_preparation.py",
             ROOT / "openwave-tuning/candidate/oph_openwave/modular_native.py",
             ROOT / "openwave-tuning/candidate/oph_openwave/quantum_records.py",
             ROOT / "openwave-tuning/candidate/oph_openwave/uniform_native.py",
             ROOT / "openwave-tuning/candidate/oph_openwave/geometry.py"]
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def run():
    clock = CommonClock()
    rho, preparation = clock.preparation()
    archive = np.load(C10 / "joint_source_states.npz")
    old = sum(archive[f"global_allN_reservoir_{a}_{b}"] for a in (0, 1) for b in (0, 1))
    preparation_error = float(np.max(abs(rho - old)))
    duration, rate, uniform_rate = .37, .8, 1.6
    final = clock.semigroup(rho, duration, rate)
    zero, one, tail = clock.count_split(rho, duration, rate)
    tilted = clock.semigroup(rho, duration, rate, attempt_mark=.43, idle_rate=uniform_rate-rate, idle_mark=.61)
    expected_tilt_trace = math.exp(duration * (rate * (.43 - 1) + (uniform_rate - rate) * (.61 - 1)))
    idle_invariant = clock.semigroup(rho, duration, rate, idle_rate=uniform_rate-rate)
    no_move = clock.semigroup(rho, duration, rate, move_mark=0)
    coherent = clock.evolve(rho, duration)
    _, initial_records = clock.native(rho, detailed=True)
    _, coherent_records = clock.native(coherent, detailed=True)
    no_attempt_readout_difference = max(abs(a["probability"] - b["probability"])
                                      for a, b in zip(initial_records, coherent_records))
    # A concrete path retains one native zero-current attempt, one controller
    # idle and one nonzero native current. Select populated marks explicitly.
    mid = clock.evolve(rho, .07)
    _, rows = clock.native(mid, detailed=True)
    noop = max((i for i, r in enumerate(rows) if r["applied_current"] == 0), key=lambda i: rows[i]["probability"])
    _, k = clock.ops[noop]
    mid = (k @ mid) @ k.T
    mid /= np.trace(mid)
    mid = clock.evolve(mid, .21 - .07)
    _, rows = clock.native(mid, detailed=True)
    move = max((i for i, r in enumerate(rows) if r["applied_current"] != 0), key=lambda i: rows[i]["probability"])
    path_state, path = clock.trajectory(rho, duration, rate, uniform_rate, [( .07, noop), (.13, "idle"), (.21, move)])
    slowdown = .5
    rate_only = clock.semigroup(rho, duration, slowdown * rate)
    common_scaled = clock.semigroup(rho, slowdown * duration, rate)
    both_scaled = clock.semigroup(rho, duration, slowdown * rate, h_scale=slowdown)
    raw_pulse_error = abs(math.exp(duration * rate * (.43 - 1)) - math.exp(duration * rate * .5 * (.43 - 1)))
    unnormalized_failure = float(np.trace(tail).real)
    checks = {
        "same_actual384H_and_complete_charged_preparation": preparation_error < 1e-12 and np.array_equal(archive["hamiltonian"], clock.h),
        "complete_native_pulse_clock_including_zero_current": path["native_commit_pulses"] == 2 and path["native_zero_current_attempts"] == 1,
        "all_trajectory_alternatives_retained": path["complete_branching"],
        "path_duration_exact_sum_and_service_cost": abs(path["supplied_coherent_duration"] - duration) < 1e-15 and path["coherent_service_calls"] == 4,
        "trace_preserving_joined_semigroup": abs(np.trace(final) - 1) < 1e-12,
        "uniformization_idle_rate_invariance": np.max(abs(final - idle_invariant)) < 1e-12,
        "independent_Poisson_attempt_and_idle_generating_counts": abs(np.trace(tilted) - expected_tilt_trace) < 1e-12,
        "complete_count_split_with_nonzero_tail": abs(np.trace(zero) - math.exp(-rate*duration)) < 1e-12 and abs(np.trace(one) - rate*duration*math.exp(-rate*duration)) < 1e-12 and unnormalized_failure > .01,
        "positive_count_complement": np.linalg.eigvalsh(tail).min() > -1e-12,
        "no_attempt_branch_has_coherent_native_readout_change": no_attempt_readout_difference > 1e-5,
        "rate_only_idle_intervention_not_common_time_rescaling": np.linalg.norm(rate_only-common_scaled) > 1e-3,
        "simultaneous_H_and_rate_scaling_is_exact": np.max(abs(both_scaled-common_scaled)) < 1e-12,
        "native_move_count_cannot_replace_attempt_count": float(np.trace(no_move).real) - float(np.trace(zero).real) > .05,
    }
    arrays = {"input_state": rho, "final_state": final, "count0_state": zero,
              "count1_state": one, "count_ge2_state": tail, "no_moving_attempt_state": no_move,
              "path_final_state": path_state, "coherent_no_attempt_conditional_state": coherent,
              "rate_only_slowdown_state": rate_only, "common_scaled_state": common_scaled}
    preparation_rows = []
    complete_preparation_final = np.zeros_like(rho)
    for meta, state in preparation:
        key = f"reservoir_{meta['reservoir_initial']}_{meta['reservoir_final']}"
        arrays[key+"_input"] = state
        outputs = clock.count_split(state, duration, rate)
        for suffix, output in zip(("_count0", "_count1", "_count_ge2"), outputs):
            arrays[key+suffix] = output
        complete_preparation_final += sum(outputs)
        preparation_rows.append({**meta, "probability": float(np.trace(state).real),
                                 "input_state_archive_key": key+"_input",
                                 "complete_final_count_state_archive_keys": [key+s for s in ("_count0", "_count1", "_count_ge2")]})
    checks["all_reservoir_record_and_count_joint_states_retained"] = np.max(abs(complete_preparation_final-final)) < 1e-12
    np.savez_compressed(HERE / "common_clock_states.npz", **arrays)
    receipt = {
        "schema": "oph.common_c10.marked_clock.v1", "checks": {k: bool(v) for k, v in checks.items()},
        "parameters": {"duration": duration, "native_attempt_rate": rate, "uniformization_rate": uniform_rate,
                       "law": "MIX", "dimension": clock.n, "Hamiltonian_inputs": {"mass": .23, "kappa": .37, "face_coefficient": .61},
                       "preparation": "Recomputed global interacting ground reference followed by complete fair-0/1 charged reservoir instrument; all four branches retained",
                       "source_pulse_cost": "One native commit pulse per native attempt, including every zero-current branch",
                       "coherent_service_cost": "One exp(-i H delta) service per marked waiting interval; native compilation cost unspecified"},
        "preparation_branches": preparation_rows,
        "operational_relation": {
            "generator": "L rho=-i[H,rho]+lambda*(sum_r K_r rho K_r* -rho)",
            "marked_generator": "L(z,x,y)=-i ad_H +lambda*z*(J_noop+x J_move)-lambda Id +(Lambda-lambda)*(y-1)Id",
            "attempt_law": "N~Poisson(lambda*T); every no-current native event contributes one",
            "idle_law": "I~Poisson((Lambda-lambda)*T), independent of N",
            "duration": "sum of all coherent waiting durations is T for every retained history",
            "resource_law": "native commit pulses=N; coherent interval services=N+I+1 under the declared segmentation",
            "expected_native_commit_pulses": rate*duration, "expected_coherent_services": uniform_rate*duration+1,
            "physical_clock_equality": False,
        },
        "count_probabilities": {"N=0": float(np.trace(zero).real), "N=1": float(np.trace(one).real),
                                "N>=2": unnormalized_failure, "zero_moving_attempts": float(np.trace(no_move).real)},
        "path_witness": path,
        "mutations": {"N0_coherent_readout_change": no_attempt_readout_difference,
                      "rate_only_vs_common_time_Frobenius_difference": float(np.linalg.norm(rate_only-common_scaled)),
                      "false_half_rate_count_PGF_difference": raw_pulse_error,
                      "tail_discarded_probability": unnormalized_failure},
        "errors": {"preparation_vs_pinned_archive": preparation_error,
                   "count_PGF": float(abs(np.trace(tilted)-expected_tilt_trace)),
                   "joint_rate_H_scaling": float(np.max(abs(both_scaled-common_scaled))),
                   "count_complement_min_eigenvalue": float(np.linalg.eigvalsh(tail).min())},
        "scope": "Exact finite hybrid architecture on the same C10 preparation, H and complete native instrument. The Poisson rate, coherent duration and reservoir preparation are explicit inputs. Native pulse cost counts actual commit attempts; coherent-H pulse compilation and physical elapsed-time calibration are not supplied. No modular, cubic, passive-rest-clock, field-limit or particle dictionary is transferred by this clock adapter.",
        "source_sha256": source_pins(), "state_archive": "common_clock_states.npz",
        "state_archive_sha256": hashlib.sha256((HERE / "common_clock_states.npz").read_bytes()).hexdigest(),
        "official_rows_changed": [],
    }
    (HERE / "COMMON_CLOCK.json").write_text(json.dumps(receipt, indent=2) + "\n")
    contract = {
        "schema": "oph.source_state_clock_contract.v1",
        "clock_id": "C10_complete_marked_Poisson_v1",
        "source_id": "C10_interacting_allnumber384_Z6",
        "source_dimension": clock.n,
        "source_H_array_sha256_le_float64": hashlib.sha256(np.asarray(clock.h, dtype="<f8").tobytes()).hexdigest(),
        "preparation_array_sha256_le_complex128": hashlib.sha256(np.asarray(rho, dtype="<c16").tobytes()).hexdigest(),
        "basis": "C10 JointNativeSource.states order; (occupation mask, six modular half-records)",
        "preparation": receipt["parameters"]["preparation"],
        "preparation_labels": preparation_rows,
        "dynamics": receipt["operational_relation"],
        "parameter_units": {"T": "supplied inverse-C10-H unit", "lambda": "native attempts per supplied T unit",
                            "Lambda": "marked opportunities per supplied T unit", "N": "integer native attempts",
                            "native_commit_pulses": "integer C/E commit operations, including zero-current outcomes"},
        "numerical_member": receipt["parameters"],
        "api": {
            "constructor": "common_clock.CommonClock(law='MIX')",
            "preparation": "preparation() -> (rho, list[(reservoir_metadata, subnormalized_state)])",
            "unconditional_or_tilted": "semigroup(rho,T,lambda,attempt_mark=z,move_mark=x,idle_rate=Lambda-lambda,idle_mark=y) -> subnormalized full384 density",
            "count_classes": "count_split(rho,T,lambda) -> (rho_N0,rho_N1,rho_Nge2)",
            "marked_history": "trajectory(rho,T,lambda,Lambda,[(absolute_time,native_index_or_idle),...]) -> (conditional_full384_density,complete_ledger)",
            "terminal_native_readout": "native(rho,detailed=True) -> (poststate,complete_actual_C_E_records)",
            "terminal_number_energy": "model.energy_number_readout(rho) -> complete actual-H grouped probabilities",
        },
        "state_admission": {"semigroup_and_count_split": "number-block-diagonal densities; cross-number coherences rejected",
                            "full_history_and_coherent_service": "full384 density, retaining cross-number coherences"},
        "transfer_boundaries": {"physical_time_equality": False, "native_compiler_for_coherent_H": False,
                               "passive192CAR_clock_equivalence": False, "routed_flux_clock_equivalence": False,
                               "BR16_log_native_transfer_equivalence": False, "continuum_or_particle_award": False},
        "source_sha256": receipt["source_sha256"],
        "receipt": "COMMON_CLOCK.json", "receipt_sha256": hashlib.sha256((HERE/"COMMON_CLOCK.json").read_bytes()).hexdigest(),
        "state_archive": receipt["state_archive"], "state_archive_sha256": receipt["state_archive_sha256"],
        "official_rows_changed": [],
    }
    (HERE / "CLOCK_CONTRACT.json").write_text(json.dumps(contract, indent=2)+"\n")
    return receipt


if __name__ == "__main__":
    result = run()
    print(json.dumps({"checks": result["checks"], "counts": result["count_probabilities"], "mutations": result["mutations"]}, indent=2))
    if not all(result["checks"].values()):
        raise SystemExit("clock check failed")
