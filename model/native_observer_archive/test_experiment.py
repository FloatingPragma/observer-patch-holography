"""Adversarial replay and complete factored-law checks.

Exploratory calculations and adversarial receipts use isolated temporary
fixtures. No test records, overwrites, or repairs package evidence.
"""

from copy import deepcopy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from . import experiment as e
from .packet import Rail


def _json_bytes(value, allow_nan=False):
    return (json.dumps(value, separators=(",", ":"), allow_nan=allow_nan) + "\n").encode("utf8")


def _fraction(value):
    return Fraction(int(value["numerator"]), int(value["denominator"]))


def _replace(value, path, replacement):
    for key in path[:-1]:
        value = value[key]
    value[path[-1]] = replacement


def _rational_paths(value, path=()):
    if type(value) is dict:
        if set(value) == {"numerator", "denominator"}:
            yield path
        else:
            for key, child in value.items():
                yield from _rational_paths(child, path + (key,))
    elif type(value) is list:
        for index, child in enumerate(value):
            yield from _rational_paths(child, path + (index,))


class ExperimentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.science = json.loads(json.dumps(e.compute(False), allow_nan=False))
        cls.protocol = json.loads((e.HERE / "protocol.json").read_text(encoding="utf8"))
        e.validate_science(cls.science)

    def _fixture(self):
        temporary = tempfile.TemporaryDirectory(prefix="oph-observer-archive-tests-")
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        here = root / "model" / "native_observer_archive"
        here.mkdir(parents=True)
        for name in (*e.BINDINGS, "model/__init__.py", "tools/pin_scope.py", "code/pin_scope.py"):
            destination = root / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(e.ROOT / name, destination)
        freeze = {"schema": "native-observer-archive-freeze.v2",
                  "files": {name: e.digest(root / name) for name in e.BINDINGS}}
        (here / "freeze.json").write_bytes(_json_bytes(freeze))
        return root, here, freeze

    def _receipt(self, here, science, **extra):
        payload = {"freeze_sha256": e.digest(here / "freeze.json"), "science": science,
                   "execution": {"scope": "temporary adversarial fixture"}, **extra}
        raw = gzip.compress(_json_bytes(payload, allow_nan=True), mtime=0)
        (here / "results.json.gz").write_bytes(raw)
        return raw

    def _cli(self, root, *arguments):
        return subprocess.run([sys.executable, "-B", "-O", "-m",
                               "model.native_observer_archive.experiment", *arguments],
                              cwd=root, text=True, capture_output=True, check=False, timeout=60)

    def test_protocol_exhausts_declared_cases_modes_and_all_branches(self):
        expected_cases = {(tuple(background), source, variant)
                          for background in self.protocol["backgrounds"]
                          for source in self.protocol["sources"]
                          for variant in self.protocol["variants"]}
        actual_cases = {(tuple(row["background"]), row["source"], row["variant"])
                        for row in self.science["cases"]}
        self.assertEqual(actual_cases, expected_cases)
        self.assertEqual(len(self.science["cases"]), len(expected_cases))
        self.assertEqual(len(actual_cases), 30)
        clocks = self.science["clock_laws"]
        self.assertEqual({(row["direction"], row["time"]) for row in clocks},
                         {(direction, tau) for direction in self.protocol["packet_directions"]
                          for tau in self.protocol["read_times"]})
        self.assertEqual(len(clocks), 10)
        self.assertEqual(len(self.science["spectrum"]["values"]), self.protocol["rail"]["sites"])
        for case in self.science["cases"]:
            self.assertEqual(_fraction(case["address_weight"]), Fraction(1, 12))
            records = case["complete_prefix_records"]
            self.assertEqual(len(records), 12)
            self.assertTrue(all(len(branch) == 5 for branch in records))
            for stage in range(5):
                # Retained address registers distinguish all twelve branches,
                # including no-operations and unfinished observer records.
                self.assertEqual(len({json.dumps(branch[stage]) for branch in records}), 12)
            self.assertEqual(len(case["public_laws"]), len(clocks))
        self.assertEqual(self.science["physical_rows_validated"], [])
        self.assertEqual(self.science["empirical_calibration_inputs"], [])

    def test_reconstruct_complete_factored_law_and_public_marginals(self):
        rail = Rail(**self.protocol["rail"])
        stages = np.array([rail.stage(n) for n in range(rail.sites)])
        clocks = {(row["direction"], row["time"]): row for row in self.science["clock_laws"]}
        for case in self.science["cases"]:
            records = case["complete_prefix_records"]
            for public in case["public_laws"]:
                clock = clocks[public["direction"], public["time"]]
                probabilities = np.array(clock["probabilities"])
                self.assertEqual(probabilities.size, rail.sites)
                complete = np.repeat(probabilities / 12, 12)
                self.assertEqual(complete.size, 12 * rail.sites)
                self.assertTrue(np.all(complete >= 0))
                self.assertAlmostEqual(float(complete.sum()), 1.0, places=13)
                stage_law = np.bincount(stages, weights=probabilities, minlength=5)
                np.testing.assert_allclose(stage_law, clock["stage_probabilities"], atol=2e-14, rtol=0)
                reconstructed = np.zeros(16)
                for branch in records:
                    labels = np.array([4 * (2 * row[0] + row[1]) + 2 * row[2] + row[3]
                                       for row in branch])
                    reconstructed += np.bincount(labels[stages], weights=probabilities / 12,
                                                 minlength=16)
                np.testing.assert_allclose(reconstructed, public["probabilities"], atol=2e-14, rtol=0)
                ready = sum(reconstructed[4 * a + b] for a in (2, 3) for b in (2, 3))
                self.assertAlmostEqual(ready, public["both_valid"], places=13)
                self.assertEqual(public["valid_disagreement"], reconstructed[11] + reconstructed[14])
                if case["variant"] in ("no_communication", "no_publish"):
                    self.assertEqual(ready, 0)
        # Distinct branch labels make complete-record L1 equal clock-factor L1.
        left = np.array(clocks[1, 256]["probabilities"])
        right = np.array(clocks[-1, 256]["probabilities"])
        self.assertAlmostEqual(float(np.abs(np.repeat(left / 12, 12) - np.repeat(right / 12, 12)).sum()),
                               float(np.abs(left - right).sum()), places=13)

    def test_discrete_missing_mode_branch_domain_and_scalar_mutations_fail(self):
        mutations = [
            (("cases",), self.science["cases"][:-1]),
            (("cases", 0, "complete_prefix_records"), self.science["cases"][0]["complete_prefix_records"][:-1]),
            (("cases", 0, "complete_prefix_records", 0), self.science["cases"][0]["complete_prefix_records"][0][:-1]),
            (("spectrum", "values"), self.science["spectrum"]["values"][:-1]),
            (("clock_laws", 0, "probabilities"), self.science["clock_laws"][0]["probabilities"][:-1]),
            (("cases", 0, "complete_prefix_records", 0, 0, 0), True),
            (("cases", 0, "complete_prefix_records", 0, 0, 0), 2),
            (("cases", 0, "source"), True),
            (("cases", 0, "background"), [0, 0, 0, 1]),
            (("cases", 0, "public_laws", 0, "both_valid"), 1.0),
            (("cases", 0, "public_laws", 0, "probabilities"), [1.0] + [0.0] * 14 + [True]),
            (("clock_laws", 0, "probabilities", 0), float("nan")),
            (("clock_laws", 0, "probabilities", 0), float("inf")),
            (("clock_laws", 0, "probabilities", 0), -1e-30),
            (("literal_operator_check", "vector_l2_error"), float("nan")),
            (("physical_rows_validated",), [1]),
        ]
        for path, replacement in mutations:
            changed = deepcopy(self.science)
            _replace(changed, path, replacement)
            with self.subTest(path=path), self.assertRaises(ValueError):
                e.compare(self.science, changed)
        for mutate in (lambda value: value.pop("resources"),
                       lambda value: value.update({"source_equivalence_claimed": True})):
            changed = deepcopy(self.science)
            mutate(changed)
            with self.assertRaisesRegex(ValueError, "keys"):
                e.compare(self.science, changed)

    def test_aggregate_l1_cannot_hide_in_per_scalar_tolerance(self):
        for name, count in (("probabilities", 4096), ("stage_probabilities", 5)):
            original = {name: [1 / count] * count}
            for balanced in (False, True):
                changed = deepcopy(original)
                for i in range(count):
                    changed[name][i] += (-1 if balanced and i >= count // 2 else 1) * e.TOL / 2
                self.assertTrue(all(abs(a - b) <= e.TOL for a, b in zip(original[name], changed[name])))
                with self.subTest(name=name, balanced=balanced), self.assertRaisesRegex(ValueError, "normalization|aggregate"):
                    e.compare(original, changed)

    def test_every_rational_bound_and_weight_is_compared_exactly(self):
        paths = list(_rational_paths(self.science))
        self.assertGreater(len(paths), 30)
        for path in paths:
            value = self.science
            for key in path:
                value = value[key]
            changed = deepcopy(self.science)
            # Even a change far below floating replay tolerance must fail.
            _replace(changed, path + ("numerator",), str(int(value["numerator"]) - 1))
            with self.subTest(path=path), self.assertRaises(ValueError):
                e.compare(self.science, changed)
        certificate = self.science["interval_certificate"]
        self.assertGreater(_fraction(certificate["boundary_upper"]), 0)
        self.assertLess(_fraction(certificate["failure_upper"]), Fraction(1, 1000))
        self.assertGreater(_fraction(self.science["binary_packet_vector_error_upper"]), 0)

    def test_scientific_scalar_probability_and_error_domains(self):
        for field, baseline in (("both_valid", 0.0), ("valid_disagreement", 0.0),
                                ("vector_l2_error", 0.0)):
            with self.subTest(field=field), self.assertRaises(ValueError):
                e.compare({field: baseline}, {field: -1e-12})
        for field in ("both_valid", "valid_disagreement"):
            with self.subTest(field=field), self.assertRaises(ValueError):
                e.compare({field: 1.0}, {field: 1.0 + 1e-12})

    def test_stored_gates_are_independent_of_replay_tolerance(self):
        residual = self.science["literal_operator_check"]["vector_l2_error"]
        over_gate = e.TOL + residual / 2
        changed = deepcopy(self.science)
        changed["literal_operator_check"]["vector_l2_error"] = over_gate
        self.assertGreater(over_gate, e.TOL)
        # This is close enough to replay, but it is not a passing experiment.
        e.compare(self.science, changed)
        with self.assertRaisesRegex(ValueError, "stored literal propagation residual gate"):
            e.validate_science(changed)
        no_communication = next(i for i, row in enumerate(self.science["cases"])
                                if row["variant"] == "no_communication")
        mutations = [
            (("clock_laws", 0, "raw_norm"), 1.0 + e.TOL, "stored packet norm residual"),
            (("literal_operator_check", "norm_squared"), 1.0 + e.TOL, "stored literal norm"),
            (("literal_operator_check", "initial_energy"), 2 * e.TOL, "energy/meter"),
            (("literal_operator_check", "final_energy"), 2 * e.TOL, "energy/meter"),
            (("literal_operator_check", "pointer_measurement_energy"), 2 * e.TOL, "energy/meter"),
            (("clock_laws", 0, "mean_hamiltonian"), 2 * e.TOL, "packet energy"),
            (("cases", no_communication, "public_laws", 0, "both_valid"), 1e-10, "readiness disagrees"),
            (("cases", no_communication, "public_laws", 0, "both_valid"), 1e-16, "invalid blank pointers"),
            (("cases", 0, "public_laws", 0, "valid_disagreement"), 1e-16, "pointer disagreement"),
            (("interval_certificate", "normalized_binary_packet_failure_upper"),
             {"numerator": "1", "denominator": "1000"}, "interval budget"),
            (("controls", "completed_source_tv_native_feedback_off"),
             {"numerator": "0", "denominator": "1"}, "same-apparatus"),
        ]
        for path, replacement, message in mutations:
            changed = deepcopy(self.science)
            _replace(changed, path, replacement)
            with self.subTest(path=path), self.assertRaisesRegex(ValueError, message):
                e.validate_science(changed)
        # Consistently falsifying both summary and law cannot weaken an exact
        # no-communication/no-publication gate through the summation allowance.
        for variant in ("no_communication", "no_publish"):
            changed = deepcopy(self.science)
            case = next(row for row in changed["cases"] if row["variant"] == variant)
            row = case["public_laws"][0]
            row["probabilities"][0] -= 1e-16
            row["probabilities"][10] += 1e-16
            row["both_valid"] = 1e-16
            with self.assertRaisesRegex(ValueError, "invalid blank pointers"):
                e.validate_science(changed)

    def test_public_laws_retain_structural_zeros_and_record_factorization(self):
        # B.valid with A.invalid does not enter the both-valid summary. It is
        # nevertheless forbidden, even below the numerical replay allowance.
        for variant, outcome, message in (
            ("no_communication", 2, "disabled receiver"),
            ("no_publish", 2, "disabled receiver"),
            ("no_publish", 8, "disabled publisher"),
            ("native", 2, "outside its complete prefix records"),
        ):
            for mass in (5e-11, 1e-16):
                changed = deepcopy(self.science)
                case = next(row for row in changed["cases"] if row["variant"] == variant)
                law = case["public_laws"][0]["probabilities"]
                law[0] -= mass
                law[outcome] += mass
                e.compare(self.science, changed)
                with self.subTest(variant=variant, mass=mass), self.assertRaisesRegex(ValueError, message):
                    e.validate_science(changed)
        # An allowed record cannot acquire an independent probability either:
        # the clock factors and complete prefix records determine its weight.
        changed = deepcopy(self.science)
        case = next(row for row in changed["cases"] if row["variant"] == "native")
        law = case["public_laws"][0]["probabilities"]
        law[0] -= 5e-11
        law[8] += 5e-11
        e.compare(self.science, changed)
        with self.assertRaisesRegex(ValueError, "clock and record factors"):
            e.validate_science(changed)

    def test_invalid_regenerated_science_cannot_be_recorded(self):
        root, here, _ = self._fixture()
        changed = deepcopy(self.science)
        changed["literal_operator_check"]["vector_l2_error"] = 2 * e.TOL
        with patch.object(e, "ROOT", root), patch.object(e, "HERE", here), \
             patch.object(e, "compute", return_value=changed), \
             patch.object(sys, "argv", ["experiment", "--record"]):
            with self.assertRaisesRegex(ValueError, "stored literal propagation residual gate"):
                e.main()
        self.assertFalse((here / "results.json.gz").exists())

    def test_freeze_requires_exact_members_and_rejects_each_source_change(self):
        root, here, freeze = self._fixture()
        for name in ("model/native_memory_clock/__init__.py", "model/native_memory_clock/native.py"):
            self.assertIn(name, e.BINDINGS)
            self.assertEqual((root / name).read_bytes(), (e.ROOT / name).read_bytes())
        with patch.object(e, "ROOT", root), patch.object(e, "HERE", here):
            self.assertEqual(e.frozen(), e.digest(here / "freeze.json"))
            for name in e.BINDINGS:
                path = root / name
                original = path.read_bytes()
                path.write_bytes(original + b"\n")
                try:
                    with self.subTest(name=name), self.assertRaisesRegex(ValueError, "frozen source changed"):
                        e.frozen()
                finally:
                    path.write_bytes(original)
            mutations = []
            changed = deepcopy(freeze)
            changed["files"].pop("model/native_memory_clock/native.py")
            mutations.append(changed)
            changed = deepcopy(freeze)
            changed["files"]["unbound.py"] = "0" * 64
            mutations.append(changed)
            changed = deepcopy(freeze)
            changed["files"][next(iter(changed["files"]))] = True
            mutations.append(changed)
            mutations.extend(({**freeze, "schema": "unknown"},
                              {**freeze, "schema": "native-observer-archive-freeze.v1"},
                              {**freeze, "ignored_claim": True}))
            for changed in mutations:
                (here / "freeze.json").write_bytes(_json_bytes(changed))
                with self.assertRaises(ValueError):
                    e.frozen()

    def test_optimized_cli_preserves_failed_evidence_and_checks_binding_first(self):
        root, here, freeze = self._fixture()
        for name, arguments in (("results.json.gz", ("--record",)),
                                ("extended_results.json.gz", ("--record", "--extended"))):
            target = here / name
            original = b"retained invalid evidence\n"
            target.write_bytes(original)
            result = self._cli(root, *arguments)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("refusing to overwrite retained evidence", result.stderr)
            self.assertEqual(target.read_bytes(), original)
        changed = deepcopy(freeze)
        changed["files"].pop("model/native_memory_clock/native.py")
        (here / "freeze.json").write_bytes(_json_bytes(changed))
        result = self._cli(root)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("complete source binding", result.stderr)
        (here / "freeze.json").write_bytes(_json_bytes(freeze))
        source = root / "model" / "native_memory_clock" / "native.py"
        source.write_bytes(source.read_bytes() + b"\n# temporary source mutation\n")
        result = self._cli(root)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("frozen source changed", result.stderr)

    def test_optimized_cli_rejects_changed_executable_package_initializers(self):
        root, here, _ = self._fixture()
        self._receipt(here, self.science)
        for name in ("model/native_memory_clock/__init__.py",):
            path = root / name
            original = path.read_bytes()
            path.write_bytes(original + b"\n_import_executed_before_replay = True\n")
            try:
                result = self._cli(root)
                with self.subTest(name=name):
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("frozen source changed: " + name, result.stderr)
            finally:
                path.write_bytes(original)

    def test_optimized_cli_replays_complete_fixture_and_rejects_scientific_corruption(self):
        root, here, _ = self._fixture()
        baseline = self._receipt(here, self.science)
        result = self._cli(root)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('"verified": true', result.stdout)
        self.assertEqual((here / "results.json.gz").read_bytes(), baseline)
        mutations = [
            (("cases", 0, "complete_prefix_records"), self.science["cases"][0]["complete_prefix_records"][:-1]),
            (("spectrum", "values"), self.science["spectrum"]["values"][:-1]),
            (("interval_certificate", "failure_upper", "numerator"), "0"),
            (("cases", 0, "complete_prefix_records", 0, 0, 0), True),
            (("cases", 0, "public_laws", 0, "both_valid"), 1.0),
            (("clock_laws", 0, "probabilities", 0), float("nan")),
            (("literal_operator_check", "vector_l2_error"),
             e.TOL + self.science["literal_operator_check"]["vector_l2_error"] / 2),
            (("clock_laws", 0, "raw_norm"), 1.0 + e.TOL),
            (("literal_operator_check", "norm_squared"), 1.0 + e.TOL),
            (("cases", 2, "public_laws", 0, "both_valid"), 1e-10),
        ]
        for path, replacement in mutations:
            changed = deepcopy(self.science)
            _replace(changed, path, replacement)
            raw = self._receipt(here, changed)
            result = self._cli(root)
            with self.subTest(path=path):
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertNotIn('"verified": true', result.stdout)
                self.assertEqual((here / "results.json.gz").read_bytes(), raw)
        self._receipt(here, self.science, unbound_scientific_claim=True)
        result = self._cli(root)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("receipt schema or source pin", result.stderr)


if __name__ == "__main__":
    unittest.main()
