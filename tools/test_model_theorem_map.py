"""Synthetic reference/coverage/freshness tests; no scientific producer runs."""

from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
import io
import json
import subprocess
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import build_model_theorem_map as builder


class ModelTheoremMapTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "proof.lean").write_text("theorem exact_join : True := trivial\n")
        self.premises = {"rows": [{"id": "PR-01"}, {"id": "PR-85"}]}
        self.claims = {"claims": [
            {"claim_id": "C-SUPPLIED", "assumptions": ["declared_action", "shared"],
             "imported_results": ["spectral_theorem"],
             "premise_dependencies": {"classification": "explicit_non_consumer",
                                      "consumed": [], "open": [], "boundary": []}},
            {"claim_id": "C-BRANCH", "assumptions": ["shared"],
             "imported_results": ["spectral_theorem", "external_scale"],
             "null_decision_assumptions": ["coverage", "shared"],
             "threshold_kinematic_assumptions": ["dispersion"],
             "premise_dependencies": {"classification": "explicit_edges",
                                      "consumed": ["PR-01"], "open": ["PR-01"],
                                      "boundary": ["PR-85"]}},
        ]}
        self.data = {
            "schema": builder.SCHEMA, "as_of": "2026-10-11",
            "coverage_note": "Navigation only; physical identification is excluded.",
            "contract_coverage": {
                "schema": builder.CONTRACT_SCHEMA,
                "claim_ids": ["C-SUPPLIED", "C-BRANCH"], "lean_files": [],
                "domains": {"native": {field: "declared " + field for field in builder.DOMAIN_FACETS}},
                "bindings": [], "connections": [],
            },
            "groups": [{
                "id": "RE-01", "title": "Supplied scalar | branch",
                "links": sorted(builder.MODEL_IDS), "scope": "active",
                "status": "conditional", "premises": ["PR-01", "PR-85"],
                "constructed": "The stated mathematical implication is proved.",
                "remaining": "Supply an implementation on the same domain.",
                "closure": "Bound the implementation error without circular inference.",
                "sources": [{"path": "proof.lean", "label": "Exact implication",
                             "anchor": "theorem exact_join"}],
                "claims": ["C-SUPPLIED"], "assumptions": ["declared_action", "coverage"],
            }],
        }

    def errors(self, data=None):
        return builder.validate(self.data if data is None else data,
                                self.premises, self.claims, self.root)

    def test_unstaged_lean_source_is_inventoried_but_ignored_build_tree_is_not(self):
        subprocess.run(['git','init','-q',str(self.root)],check=True)
        (self.root/'Lean').mkdir()
        (self.root/'Lean'/'Fresh.lean').write_text('theorem fresh : True := trivial\n')
        (self.root/'.gitignore').write_text('Lean/generated/\n')
        (self.root/'Lean'/'generated').mkdir()
        (self.root/'Lean'/'generated'/'Ignored.lean').write_text('theorem ignored : True := trivial\n')
        self.assertEqual(builder.lean_inventory(self.root),{'Lean/Fresh.lean':1})

    def reject(self, mutate, expected):
        data = deepcopy(self.data)
        mutate(data)
        self.assertIn(expected, "\n".join(self.errors(data)))

    def test_complete_fixture_passes_without_changing_claim_edges(self):
        before = deepcopy(self.claims)
        self.assertEqual(self.errors(), [])
        rendered = builder.render(self.data, builder.inventory(self.data, self.premises, self.claims, self.root))
        self.assertEqual(self.claims, before)
        self.assertIn("explicit_non_consumer", rendered)
        self.assertIn("no proof-discharge edges", rendered)
        self.assertIn("not a global completeness", rendered.replace("\n", " "))
        self.assertIn("Supplied scalar \\| branch", rendered)

    def test_missing_and_unknown_canonical_ids_fail_separately(self):
        self.reject(lambda d: d["groups"][0]["premises"].pop(), "omitted canonical premises ['PR-85']")
        self.reject(lambda d: d["groups"][0]["premises"].append("PR-404"), "unknown references ['PR-404']")

    def test_missing_and_unknown_model_criteria(self):
        self.reject(lambda d: d["groups"][0]["links"].remove("M12"), "missing model criteria ['M12']")
        self.reject(lambda d: d["groups"][0]["links"].append("M13"), "unknown references ['M13']")

    def test_duplicate_group_id_fails(self):
        self.reject(lambda d: d["groups"].append(deepcopy(d["groups"][0])), "duplicate group id RE-01")

    def test_shared_premises_between_distinct_groups_are_allowed(self):
        other = deepcopy(self.data["groups"][0])
        other["id"] = "RE-02"
        other["scope"] = "downstream"
        self.data["groups"].append(other)
        self.assertEqual(self.errors(), [])
        self.assertEqual(builder.inventory(self.data, self.premises, self.claims, self.root)["covered_premises"], 2)

    def test_duplicate_references_within_group_fail(self):
        for field in ("premises", "links", "claims", "assumptions"):
            with self.subTest(field=field):
                self.reject(lambda d: d["groups"][0][field].append(d["groups"][0][field][0]),
                            f"{field}: duplicate entries")

    def test_source_anchor_is_exact_and_case_sensitive(self):
        self.reject(lambda d: d["groups"][0]["sources"][0].update(anchor="Theorem exact_join"),
                    "exact substring absent")
        self.reject(lambda d: d["groups"][0]["sources"][0].update(anchor=" "), "nonempty exact substring")

    def test_source_fragment_must_name_an_explicit_anchor(self):
        (self.root / "page.md").write_text('# Page\n\n<a id="gravity-rungs"></a>\n\n## Ladder\n')
        self.data["groups"][0]["sources"].append({"path": "page.md#gravity-rungs", "label": "Ladder"})
        self.assertEqual(self.errors(), [])
        rendered = builder.render(self.data, builder.inventory(self.data, self.premises, self.claims, self.root))
        self.assertIn("- [Ladder](../page.md#gravity-rungs).", rendered)
        self.reject(lambda d: d["groups"][0]["sources"][1].update(path="page.md#ladder"),
                    "anchor #ladder absent from page.md")
        self.reject(lambda d: d["groups"][0]["sources"][1].update(path="page.md#Gravity_Rungs"),
                    "fragment must be lowercase words joined by hyphens")
        self.reject(lambda d: d["groups"][0]["sources"][1].update(path="page.md#"),
                    "fragment must be lowercase words joined by hyphens")
        self.reject(lambda d: d["groups"][0]["sources"][1].update(path="absent.md#gravity-rungs"),
                    "missing file")
        self.reject(lambda d: d["groups"][0]["sources"][1].update(path="../page.md#gravity-rungs"),
                    "inside the repository")

    def test_source_without_anchor_is_allowed_but_missing_file_is_not(self):
        del self.data["groups"][0]["sources"][0]["anchor"]
        self.assertEqual(self.errors(), [])
        self.reject(lambda d: d["groups"][0]["sources"][0].update(path="absent.lean"), "missing file")

    def test_escaping_paths_and_external_symlink_are_rejected(self):
        for path in ("../proof.lean", "/proof.lean", "dir\\proof.lean"):
            with self.subTest(path=path):
                self.reject(lambda d: d["groups"][0]["sources"][0].update(path=path), "inside the repository")
        with tempfile.TemporaryDirectory() as outside:
            external = Path(outside) / "source.lean"
            external.write_text("theorem exact_join")
            (self.root / "link.lean").symlink_to(external)
            self.reject(lambda d: d["groups"][0]["sources"][0].update(path="link.lean"), "outside repository")

    def test_empty_required_text_and_unknown_keys_fail(self):
        for key in ("title", "constructed", "remaining", "closure"):
            with self.subTest(key=key):
                self.reject(lambda d: d["groups"][0].update({key: "\n "}), f"{key}: expected nonempty string")
        self.reject(lambda d: d.update(coverage_note=""), "coverage_note")
        self.reject(lambda d: d["groups"][0].update(proved=True), "unknown keys ['proved']")
        self.reject(lambda d: d.update(discharge_edges=[]), "unknown keys ['discharge_edges']")

    def test_status_scope_and_date_validation(self):
        self.reject(lambda d: d["groups"][0].update(status="closed"), "closed is not supported")
        self.reject(lambda d: d["groups"][0].update(scope="all_proofs"), "scope: expected")
        self.reject(lambda d: d.update(as_of="2026-13-11"), "expected YYYY-MM-DD")
        self.reject(lambda d: d.update(as_of="20261011"), "expected YYYY-MM-DD")

    def test_optional_references_must_resolve(self):
        self.reject(lambda d: d["groups"][0]["claims"].append("invented"), "claims: unknown references")
        self.reject(lambda d: d["groups"][0]["assumptions"].append("invented"), "assumptions: unknown references")

    def test_extra_assumption_fields_are_independently_inventoried(self):
        inv = builder.inventory(self.data, self.premises, self.claims, self.root)
        self.assertEqual(inv["ordinary_assumptions"], 2)
        self.assertEqual(inv["imported_names"], 2)
        self.assertEqual(inv["role_counts"], {"consumed": 1, "open": 1, "boundary": 1})
        self.assertEqual(inv["extra_assumption_fields"]["null_decision_assumptions"],
                         {"claims_with_field": 1, "occurrences": 2, "unique_names": 2, "also_ordinary": 1})
        self.assertEqual(inv["extra_assumption_fields"]["threshold_kinematic_assumptions"]["unique_names"], 1)
        text = builder.render(self.data, inv)
        self.assertIn("`threshold_kinematic_assumptions` | 1 | 1 | 1 | 0", text)
        self.assertIn("cannot be added to obtain a certified total", text)

    def test_registry_shape_and_unknown_role_premise_fail_closed(self):
        self.claims["claims"][1]["premise_dependencies"]["open"].append("PR-404")
        self.assertIn("unknown premises ['PR-404']", "\n".join(self.errors()))
        self.claims["claims"][1]["null_decision_assumptions"] = "not a list"
        self.assertIn("null_decision_assumptions: expected a list", "\n".join(self.errors()))

    def test_duplicate_json_keys_are_not_silently_overwritten(self):
        path = self.root / "duplicate.json"
        path.write_text('{"schema": "first", "schema": "second"}')
        with self.assertRaisesRegex(ValueError, "duplicate JSON key: schema"):
            builder.read_json(path)

    def _write_fixture(self, register=None):
        if register is None:
            register = dict(self.premises, **{builder.SECTION: self.data})
        for relative, value in ((builder.REGISTER, register), (builder.CLAIMS, self.claims)):
            target = self.root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(json.dumps(value))

    def test_map_is_the_model_supply_section_of_the_premise_register(self):
        self.assertEqual(builder.REGISTER, Path("tracking/premise_register.json"))
        self.assertEqual(builder.SECTION, "model_supply")
        live = builder.read_json(builder.ROOT / builder.REGISTER)
        self.assertEqual(live[builder.SECTION]["schema"], builder.SCHEMA)
        self.assertEqual(builder.validate(live[builder.SECTION], live,
                                          builder.read_json(builder.ROOT / builder.CLAIMS),
                                          builder.ROOT), [])

    def test_missing_section_fails_closed_without_writing(self):
        self._write_fixture(register=self.premises)
        output = self.root / builder.OUTPUT
        stderr = io.StringIO()
        with patch.object(builder, "ROOT", self.root), redirect_stderr(stderr):
            self.assertEqual(builder.main([]), 1)
        self.assertIn("missing model_supply section", stderr.getvalue())
        self.assertFalse(output.exists())

    def test_cli_generation_and_read_only_stale_check(self):
        self._write_fixture()
        output = self.root / builder.OUTPUT
        with patch.object(builder, "ROOT", self.root), redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            self.assertEqual(builder.main(["--check"]), 1)
            self.assertFalse(output.exists())
            self.assertEqual(builder.main([]), 0)
            before = output.read_bytes()
            mtime = output.stat().st_mtime_ns
            self.assertEqual(builder.main(["--check"]), 0)
            self.assertEqual(output.stat().st_mtime_ns, mtime)
            output.write_text("stale\n")
            self.assertEqual(builder.main(["--check"]), 1)
            self.assertEqual(output.read_text(), "stale\n")
            self.assertEqual(builder.main([]), 0)
            self.assertEqual(output.read_bytes(), before)

    def test_invalid_map_cannot_overwrite_preexisting_output(self):
        self.data["groups"][0]["status"] = "closed"
        self._write_fixture()
        output = self.root / builder.OUTPUT
        output.parent.mkdir(parents=True)
        output.write_text("preserve me\n")
        with patch.object(builder, "ROOT", self.root), redirect_stderr(io.StringIO()):
            self.assertEqual(builder.main([]), 1)
        self.assertEqual(output.read_text(), "preserve me\n")

    def joined_fixture(self):
        coverage = self.data["contract_coverage"]
        coverage["domains"]["observer"] = deepcopy(coverage["domains"]["native"])
        upstream = {
            "id": "native-map", "domain": "native", "claims": ["C-SUPPLIED"],
            "establishes": ["complete instrument"], "requires": [],
            "assumptions": ["declared_action", "shared"], "imports": ["spectral_theorem"],
            "additional_assumptions": {},
            "scope": "One finite source clause.", "limitations": "No physical selection.",
            "sources": [{"path": "proof.lean", "anchor": "theorem exact_join"}],
        }
        downstream = deepcopy(upstream)
        downstream.update(id="observer-map", domain="observer", establishes=["complete history"],
                          requires=["complete instrument"])
        coverage["bindings"] = [upstream, downstream]
        coverage["connections"] = [{
            "from": "native-map", "to": "observer-map", "fact": "complete instrument",
            "coordinate_maps": {}, "sources": [{"path": "proof.lean", "anchor": "theorem exact_join"}],
        }]
        self.assertEqual(self.errors(), [])
        return coverage

    def test_every_claim_including_nonconsumers_requires_explicit_inventory(self):
        self.reject(lambda d: d["contract_coverage"]["claim_ids"].remove("C-BRANCH"),
                    "exact claim inventory/order changed")
        self.claims["claims"].append(dict(self.claims["claims"][0], claim_id="NEW-FINITE"))
        self.assertIn("exact claim inventory/order changed", "\n".join(self.errors()))

    def test_new_lean_source_cannot_disappear_behind_registry_coverage(self):
        lean = self.root / "Lean"
        lean.mkdir()
        (lean / "New.lean").write_text("/- theorem commented : True := trivial -/\n"
                                      "private theorem hidden : True := trivial\n"
                                      "theorem actual : True := trivial\n")
        self.assertEqual(builder.lean_inventory(self.root), {"Lean/New.lean": 1})
        self.assertIn("exact Lean source inventory changed", "\n".join(self.errors()))
        self.data["contract_coverage"]["lean_files"] = ["Lean/New.lean"]
        self.assertEqual(self.errors(), [])

    def test_lean_inventory_has_one_case_sensitive_order_on_all_platforms(self):
        expected = [
            "Lean/Dynamics/Exact.lean",
            "Lean/Dynamics.lean",
            "Lean/Locality.lean",
            "Lean/Main.lean",
            "Lean/ObservableNormalForms/ObservableNormalForms/Exact.lean",
            "Lean/ObservableNormalForms/lakefile.lean",
            "Lean/lakefile.lean",
        ]
        # Windows Path comparison folds case and orders lakefile before
        # Locality; both checkout and fixture discovery must preserve the
        # canonical case-sensitive component order instead. Directory/file
        # prefixes also retain the established POSIX Path ordering.
        for name in reversed(expected):
            target = self.root / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text("theorem actual : True := trivial\n")
        git_result = builder.subprocess.CompletedProcess(
            args=[], returncode=0,
            stdout=b"\0".join(name.encode() for name in reversed(expected))+b"\0")
        for mode in ("git", "fixture"):
            with self.subTest(mode=mode):
                options = ({"return_value": git_result} if mode == "git"
                           else {"side_effect": OSError("git unavailable")})
                with patch.object(builder.subprocess, "run", **options):
                    actual = builder.lean_inventory(self.root)
                self.assertEqual(list(actual), expected)
                self.assertEqual(list(actual.values()), [1]*len(expected))

    def test_binding_cannot_erase_canonical_inputs_or_claim_discharge(self):
        self.joined_fixture()
        for field in ("assumptions", "imports"):
            with self.subTest(field=field):
                self.reject(lambda d: d["contract_coverage"]["bindings"][0][field].pop(),
                            "every canonical input must remain, without discharge")
        self.reject(lambda d: d["contract_coverage"]["bindings"][0].update(discharged=True),
                    "unknown keys ['discharged']")

    def test_malformed_domain_or_fact_fields_return_errors_without_crashing(self):
        self.joined_fixture()
        self.reject(lambda d: d["contract_coverage"]["bindings"][0].update(domain=[]),
                    "unknown model domain")
        self.reject(lambda d: d["contract_coverage"]["connections"][0].update({"from": []}),
                    "unknown binding endpoint")
        self.reject(lambda d: d["contract_coverage"]["bindings"][0].update(establishes=None),
                    "expected a list")
        self.reject(lambda d: d["contract_coverage"]["bindings"][1].update(requires=[{}]),
                    "expected a list")

    def test_extra_assumption_fields_cannot_be_silently_folded_or_lost(self):
        coverage = self.joined_fixture()
        binding = coverage["bindings"][0]
        binding.update(claims=["C-BRANCH"], assumptions=["shared"],
                       imports=["external_scale", "spectral_theorem"])
        self.assertIn("all separate canonical assumption fields must remain", "\n".join(self.errors()))
        binding["additional_assumptions"] = {
            "null_decision_assumptions": ["coverage", "shared"],
            "threshold_kinematic_assumptions": ["dispersion"]}
        self.assertEqual(self.errors(), [])

    def test_each_changed_coordinate_needs_its_own_sourced_mapping(self):
        coverage = self.joined_fixture()
        for facet in builder.DOMAIN_FACETS:
            with self.subTest(facet=facet):
                self.reject(lambda d: d["contract_coverage"]["domains"]["observer"].update({facet: "different"}),
                            "explicit maps required for exactly changed coordinates")
        coverage["domains"]["observer"]["law"] = "different law"
        coverage["connections"][0]["coordinate_maps"] = {"law": "Explicit restriction, not silent substitution."}
        self.assertEqual(self.errors(), [])
        self.reject(lambda d: d["contract_coverage"]["connections"][0]["sources"][0].update(anchor="invented"),
                    "exact source anchor absent")
        self.reject(lambda d: d["contract_coverage"]["connections"][0]["coordinate_maps"].update(clock="Unneeded"),
                    "exactly changed coordinates")

    def test_missing_mismatched_and_duplicate_fact_deliveries_fail(self):
        self.joined_fixture()
        self.reject(lambda d: d["contract_coverage"]["connections"].clear(), "missing connection")
        self.reject(lambda d: d["contract_coverage"]["connections"][0].update(fact="marginal only"),
                    "produced and required fact must agree exactly")
        self.reject(lambda d: d["contract_coverage"]["connections"].append(
            deepcopy(d["contract_coverage"]["connections"][0])), "duplicate input delivery")

    def test_circular_tower_support_cannot_establish_itself(self):
        coverage = self.joined_fixture()
        coverage["bindings"][0]["requires"] = ["complete history"]
        coverage["connections"].append({"from": "observer-map", "to": "native-map",
            "fact": "complete history", "coordinate_maps": {},
            "sources": [{"path": "proof.lean", "anchor": "theorem exact_join"}]})
        self.assertIn("circular connection cannot establish its own input", "\n".join(self.errors()))

    def test_compact_view_retains_complete_machine_inventory_without_duplicate_rows(self):
        self.joined_fixture()
        before = deepcopy(self.data)
        text = builder.render_contracts(self.data, builder.inventory(
            self.data, self.premises, self.claims, self.root))
        self.assertEqual(self.data, before)
        joined = "\n".join(text)
        self.assertIn("model_supply.contract_coverage", joined)
        self.assertIn("Claims with open common-model mapping | 1", joined)
        self.assertIn("native-map", joined)
        self.assertNotIn("Every canonical contract", joined)
        self.assertNotIn("Every Lean source module", joined)
        self.assertIn("C-BRANCH", self.data["contract_coverage"]["claim_ids"])

    def test_rendered_open_mapping_does_not_invent_a_physical_witness(self):
        self.joined_fixture()
        text = builder.render(self.data, builder.inventory(self.data, self.premises, self.claims, self.root))
        self.assertIn("Claims with open common-model mapping | 1", text)
        self.assertIn("Open: no same-model instantiation recorded", text)
        self.assertIn("strings and source\nanchors are not proofs", text)
        self.assertIn("No entry silently discharges a canonical premise", text)
        self.assertIn("explicit_non_consumer", text)


if __name__ == "__main__":
    unittest.main()
