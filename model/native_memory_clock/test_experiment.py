"""Adversarial custody, complete-outcome and declared-apparatus checks.

Temporary freezes and receipts belong only to isolated test fixtures. These
tests do not create, overwrite or amend the package's scientific evidence.
"""

from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from math import ceil, log, sqrt
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from . import experiment as e
from . import clock


def _json_bytes(value):
    return (json.dumps(value, separators=(",", ":"), allow_nan=False) + "\n").encode("utf8")


class ExperimentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # A bounded exploratory calculation, never a record-mode execution.
        cls.science = json.loads(json.dumps(e.compute(False), allow_nan=False))

    def _isolated_package(self):
        temporary = tempfile.TemporaryDirectory(prefix="oph-memory-clock-tests-")
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        here = root / "model" / "native_memory_clock"
        here.mkdir(parents=True)
        (root / "model" / "__init__.py").write_text("", encoding="utf8")
        for name in e.BINDINGS:
            target = root / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(e.ROOT / name, target)
        freeze = {"schema": "native-memory-clock-freeze.v1",
                  "files": {name: e.digest(root / name) for name in e.BINDINGS}}
        (here / "freeze.json").write_bytes(_json_bytes(freeze))
        return root, here, freeze

    def _cli(self, root, *arguments):
        return subprocess.run(
            [sys.executable, "-B", "-O", "-m", "model.native_memory_clock.experiment", *arguments],
            cwd=root, text=True, capture_output=True, check=False, timeout=60,
        )

    def test_recursive_comparison_rejects_types_fields_outcomes_and_shots(self):
        baseline = {"valid": True, "count": 2, "residual": 0.0,
                    "sampling": {"performed": False, "bound": None},
                    "outcomes": [[0, 1, [], [], 0.5], [1, 0, [0], [1], 0.5]],
                    "independent_shot_outcome_indices": [0, 1, 1]}
        e.compare(baseline, deepcopy(baseline))
        mutations = []
        for key, value in (("valid", 1), ("count", 2.0), ("residual", float("nan")),
                           ("residual", float("inf")), ("residual", 1e-5)):
            changed = deepcopy(baseline)
            changed[key] = value
            mutations.append(changed)
        missing = deepcopy(baseline)
        del missing["sampling"]
        mutations.append(missing)
        added = deepcopy(baseline)
        added["unreviewed_claim"] = "accepted"
        mutations.append(added)
        changed = deepcopy(baseline)
        changed["independent_shot_outcome_indices"][0] = 1
        mutations.append(changed)
        changed = deepcopy(baseline)
        changed["outcomes"][1][2] = [1]
        mutations.append(changed)
        changed = deepcopy(baseline)
        changed["outcomes"].pop()
        mutations.append(changed)
        for changed in mutations:
            with self.subTest(mutation=changed), self.assertRaises(ValueError):
                e.compare(baseline, changed)
        close = deepcopy(baseline)
        close["residual"] = e.TOL / 4
        e.compare(baseline, close)

    def test_aggregate_probability_corruption_cannot_hide_in_scalar_tolerance(self):
        count = 4096
        table = [[1, index, [0], [0], 1 / count] for index in range(count)]
        expected = {"cases": [{"records": [{"outcomes": table}]}]}
        for balanced in (False, True):
            changed = deepcopy(expected)
            for index, row in enumerate(changed["cases"][0]["records"][0]["outcomes"]):
                row[-1] += (-1 if balanced and index >= count // 2 else 1) * 1e-10
            with self.assertRaisesRegex(ValueError, "incomplete outcome law|aggregate probability discrepancy"):
                e.compare(expected, changed)
        for field in ("x_record_law", "uniform_clock_same_read_time_law",
                      "interior_measured_return_clock_law"):
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "invalid probability"):
                e.compare({field: [0.0, 1.0]}, {field: [-1e-12, 1.0 + 1e-12]})
        for field in ("endpoint_source_law", "complete_classical_records"):
            with self.subTest(field=field), self.assertRaises(ValueError):
                e.compare({field: [[0, 0.5], [1, 0.5]]}, {field: [[0, 0.5], [1, 0.50001]]})

    def test_every_required_source_is_bound_and_byte_changes_fail(self):
        root, here, freeze = self._isolated_package()
        self.assertIn("model/native_memory_clock/PROOF.md", e.BINDINGS)
        self.assertIn("model/native_memory_clock/protocol.json", e.BINDINGS)
        self.assertIn("model/native_memory_clock/provenance.json", e.BINDINGS)
        with patch.object(e, "ROOT", root), patch.object(e, "HERE", here):
            self.assertEqual(e.frozen(), hashlib.sha256((here / "freeze.json").read_bytes()).hexdigest())
            for name in e.BINDINGS:
                path = root / name
                original = path.read_bytes()
                path.write_bytes(original + b"\n")
                try:
                    with self.subTest(name=name), self.assertRaisesRegex(ValueError, "frozen source changed"):
                        e.frozen()
                finally:
                    path.write_bytes(original)
            for mutation in ("missing", "additional", "schema"):
                changed = deepcopy(freeze)
                if mutation == "missing":
                    changed["files"].pop(next(iter(changed["files"])))
                elif mutation == "additional":
                    changed["files"]["../../unbound.py"] = "0" * 64
                else:
                    changed["schema"] = "unrecognized"
                (here / "freeze.json").write_bytes(_json_bytes(changed))
                with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                    e.frozen()

    def test_actual_cli_rejects_overwrite_and_bad_freeze_under_optimization(self):
        root, here, freeze = self._isolated_package()
        # Non-overwrite is checked before execution and preserves even failed evidence.
        for name, options in (("results.json", ("--record",)),
                              ("extended_results.json.gz", ("--record", "--extended"))):
            target = here / name
            raw = b"retained failed receipt\n"
            target.write_bytes(raw)
            result = self._cli(root, *options)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn("refusing to overwrite retained evidence", result.stderr)
            self.assertEqual(target.read_bytes(), raw)
        changed = deepcopy(freeze)
        changed["files"].pop("model/native_memory_clock/native.py")
        (here / "freeze.json").write_bytes(_json_bytes(changed))
        result = self._cli(root)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("entire declared source", result.stderr)
        (here / "freeze.json").write_bytes(_json_bytes(freeze))
        source = here / "native.py"
        source.write_bytes(source.read_bytes() + b"\n# fixture mutation\n")
        result = self._cli(root)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("frozen source changed", result.stderr)

    def test_actual_cli_rejects_scientific_receipt_mutation_read_only(self):
        root, here, _ = self._isolated_package()
        science = deepcopy(self.science)
        science["cases"][0]["records"][0]["outcomes"][0][-1] += 0.001
        payload = {"freeze_sha256": e.digest(here / "freeze.json"), "science": science,
                   "execution": {"scope": "temporary adversarial fixture"}}
        target = here / "results.json"
        raw = _json_bytes(payload)
        target.write_bytes(raw)
        result = self._cli(root)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("incomplete outcome law", result.stderr)
        self.assertEqual(target.read_bytes(), raw)
        self.assertNotIn('"verified": true', result.stdout)

    def test_protocol_matches_executed_cases_and_sampling_apparatus(self):
        protocol = json.loads((e.HERE / "protocol.json").read_text(encoding="utf8"))
        choices, gates = protocol["choices"], protocol["gates"]
        cases = self.science["cases"]
        self.assertEqual(sorted(set(row["horizon"] for row in cases)), choices["headline_horizons"])
        self.assertEqual(sorted(set(tuple(row["background"]) for row in cases)),
                         sorted(tuple(row) for row in choices["backgrounds"]))
        self.assertEqual(len(cases), 4)
        self.assertEqual(set(row["law"] for row in cases), {"native", "uphill"})
        self.assertEqual([row["horizon"] for row in self.science["complete_spectrum"]], choices["spectral_horizons"])
        self.assertEqual(gates["joint_probability_l1"], e.TOL)
        self.assertEqual(gates["norm_and_free_energy"], e.TOL)
        self.assertEqual(gates["complete_clock_spectrum_max_error"], e.TOL)
        self.assertEqual(self.science["physical_rows_validated"], [])
        self.assertEqual(self.science["empirical_mass_inputs"], [])
        total_indices = 0
        for case in cases:
            self.assertEqual([row["time_over_pi_T"] for row in case["records"]], choices["source_time_over_pi_T"])
            for record, seed in zip(case["records"], choices["seeds"]):
                sampling = record["sampling"]
                indices = record["independent_shot_outcome_indices"]
                count = choices["shots_per_terminal_law"]
                self.assertTrue(sampling["performed"])
                self.assertEqual(sampling["count"], count)
                self.assertEqual(sampling["seed"], seed)
                self.assertEqual(len(indices), count)
                self.assertTrue(all(type(index) is int and 0 <= index < len(record["outcomes"]) for index in indices))
                ideal_radius = sqrt(log(2 * len(record["outcomes"]) * 12 / 0.01) / (2 * count))
                self.assertAlmostEqual(sampling["hoeffding_radius"], ideal_radius, places=14)
                self.assertGreaterEqual(sampling["cdf_grid_bias_bound"], 0)
                self.assertGreaterEqual(sampling["simultaneous_per_outcome_error_bound"], ideal_radius)
                self.assertLessEqual(sampling["largest_sample_frequency_error"],
                                     sampling["simultaneous_per_outcome_error_bound"])
                self.assertLess(abs(sampling["cdf_endpoint_defect"]), e.TOL)
                total_indices += len(indices)
                e.check_probability_vectors([row[-1] for row in record["outcomes"]],
                                            [row[-1] for row in record["outcomes"]], "complete fixture")
        self.assertEqual(total_indices, 3072)

    def test_zero_shots_are_unsampled_and_do_not_report_zero_uncertainty(self):
        result = e.finite_run(clock.Construction(1), 0)
        for row in result["records"]:
            sampling = row["sampling"]
            self.assertFalse(sampling["performed"])
            self.assertEqual(sampling["count"], 0)
            self.assertIsNone(sampling["hoeffding_radius"])
            self.assertIsNone(sampling["simultaneous_per_outcome_error_bound"])
            self.assertIsNone(sampling["largest_sample_frequency_error"])
            self.assertEqual(row["independent_shot_outcome_indices"], [])
            self.assertAlmostEqual(sum(item[-1] for item in row["outcomes"]), 1, places=12)

    def test_sampler_bias_encloses_exact_finite_grid_distribution(self):
        grid = 2 ** 53
        for probabilities in ([1 / 3] * 3,
                              [2 ** -54, 0.25, 0.75 - 2 ** -54],
                              [0.1] * 9 + [0.1 - 1e-12]):
            table = [[0, i, [], [], probability] for i, probability in enumerate(probabilities)]
            budget = e.sampling_budget(table, [], 6100)
            certified = Fraction(budget["cdf_grid_bias_upper_rational"])
            previous = 0
            for probability, cumulative in zip(probabilities, np.cumsum(probabilities)):
                # Count all k/2^53 < the actual binary CDF boundary, exactly.
                below = max(0, min(grid, ceil(Fraction.from_float(float(cumulative)) * grid)))
                actual_probability = Fraction(below - previous, grid)
                self.assertLessEqual(abs(actual_probability - Fraction.from_float(probability)), certified)
                previous = below
            missing = Fraction(grid - previous, grid)
            self.assertLessEqual(missing, Fraction.from_float(budget["missing_probability_upper_bound"]))
            self.assertGreaterEqual(Fraction.from_float(budget["cdf_grid_bias_bound"]), certified)

    def test_binary_operator_budget_distinguishes_exact_uniform_profile(self):
        uniform = e.binary_coupling_bound(clock.Construction(3, profile="uniform"))
        self.assertEqual(Fraction(uniform["operator_norm_error_upper_rational"]), 0)
        transfer = e.binary_coupling_bound(clock.Construction(3))
        self.assertGreater(Fraction(transfer["operator_norm_error_upper_rational"]), 0)
        self.assertTrue(transfer["numerical_roundoff_is_not_certified_by_coupling_bound"])

    def test_coherent_decoder_is_same_apparatus_and_needed_for_interference(self):
        controls = e.coherent_control()
        self.assertEqual([row["phase_twist"] for row in controls], [False, True])
        self.assertAlmostEqual(controls[0]["x_record_law"][1], 0, places=12)
        self.assertAlmostEqual(controls[1]["x_record_law"][1], 1 / 12, places=12)
        for row in controls:
            self.assertAlmostEqual(sum(row["x_record_law"]), 1, places=12)
            self.assertLess(row["classical_l1_from_native"], e.TOL)
            self.assertAlmostEqual(sum(outcome[-1] for outcome in row["complete_classical_records"]), 1, places=12)
        with patch.object(e, "step", lambda c, b, direction: (b, 1)):
            with self.assertRaisesRegex(ValueError, "native coherence lost|coherent alternative invisible"):
                e.coherent_control()

    def test_dephasing_has_nonzero_work_for_a_coherent_clock_preparation(self):
        c = clock.Construction(2)
        states, h, _ = clock.prepared_span(c, {1: 1.0})
        first = clock.Basis(0, (1, 0, 0), (0, 0), (0, 0))
        second = clock.step(c, first, 1)[0]
        state = np.zeros(len(states), complex)
        state[states.index(first)] = 1 / sqrt(2)
        state[states.index(second)] = 1 / sqrt(2)
        before = float(np.vdot(state, h @ state).real)
        after = 0.0
        for n in range(3):
            branch = state * np.array([b.clock == n for b in states])
            after += float(np.vdot(branch, h @ branch).real)
        self.assertAlmostEqual(before, c.coupling(0), places=14)
        self.assertGreater(before, 0.7)
        self.assertAlmostEqual(after, 0.0, places=14)
        self.assertAlmostEqual(after - before, -c.coupling(0), places=14)

    def test_exact_native_history_counter_retains_all_addresses_and_no_operations(self):
        for n in (0, 1, 2, 3):
            for background in e.BACKGROUNDS:
                records = e.native_records(1, n, background, "native")
                self.assertEqual(len(records), 12 ** n)
                self.assertEqual(sum(records.values()), Fraction(1))
                self.assertEqual({len(key[2]) for key in records}, {n})
                self.assertEqual({len(key[3]) for key in records}, {n})
        first = e.native_records(1, 1, (0, 0, 0, 0), "native")
        self.assertTrue(any(key[1] == 1 for key in first))
        self.assertTrue(any(key[1] != 1 for key in first))


if __name__ == "__main__":
    unittest.main()
