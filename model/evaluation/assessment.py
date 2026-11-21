"""Requirement definitions, candidate evidence and quantitative readouts for the public comparison."""
from __future__ import annotations

import json
import math
from pathlib import Path
from urllib.parse import urlsplit

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
STATUS = {"demonstrated": "Demonstrated", "partial": "Partial", "open": "Open", "counterexample": "Counterexample"}
FILES = ("theorem_criteria.json", "observer_criteria.json", "openwave_criteria.json")
READOUTS = {
    "interval_echo": ("Traveller/home elapsed ratio on alternating 2- and 6-service handoffs",
                      "model/candidate_laws/interval_echo/RESULTS.json",
                      "cohorts/panels/alternating_2_6/three_bits/ratio",
                      "cohorts/panels/alternating_2_6/three_bits/standard_error"),
    "contrast_memory_spherical_v1": ("Largest probability of losing pair adjacency across nine sampled waits on 20 cells at κ = 8",
                                     "model/candidate_laws/contrast_memory_spherical_v1/RESULTS.json",
                                     "particle_stability/law/20 cells/by_kappa/8.0/max_no_pair", None),
    "scalar_channel": ("Cage clock frequency in its supplied reference units, equal-step schedule",
                       "model/experiments/cage/results.json", "clock/steps_of_equal_length/line/frequency",
                       "clock/steps_of_equal_length/line/standard_error"),
}
SHORT = {"interval_echo": "Interval Echo", "contrast_memory_spherical_v1": "Contrast memory", "scalar_channel": "Scalar channel"}
DECISIONS = {
    "interval_echo": "Full federation generator, a transported and jointly read patch clock, and a current instrument on the bound pair. The 19-port ring comparator contradicts the common rest-energy and inertia assignment for its free clusters.",
    "contrast_memory_spherical_v1": "Reference for the shared seam dynamics, with an all-size binding theorem and local reader controls. The 19-port ring resolves different rest-energy and kinetic ratios for free clusters.",
    "scalar_channel": "Independent test of a stable cage clock and finite mode dispersion. Its registered generator covers one cell. The larger cage construction has its own apparatus.",
}


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def check_source(source):
    if not isinstance(source, str) or not source.strip():
        raise ValueError(f"Missing evidence path: {source!r}")
    parsed = urlsplit(source)
    if parsed.scheme == "https" and parsed.netloc:
        return
    path = ROOT / source.partition("#")[0]
    if parsed.scheme or Path(source).is_absolute() or not path.resolve().is_relative_to(ROOT) or not path.is_file():
        raise ValueError(f"Unknown repository evidence: {source}")


def readout(name):
    if name not in READOUTS:
        return None
    label, source, pointer, uncertainty = READOUTS[name]
    receipt = read_json(ROOT / source)

    def value(pointer):
        item = receipt
        for key in pointer.split("/"):
            item = item[key]
        if isinstance(item, bool) or not isinstance(item, (float, int)) or not math.isfinite(item):
            raise ValueError(f"Expected a finite measurement at {source}#/{pointer}")
        return item

    return {"label": label, "value": value(pointer), "standard_error": value(uncertainty) if uncertainty else None,
            "source": source, "pointer": "/" + pointer, "error_pointer": "/" + uncertainty if uncertainty else None}


def species_group(candidate):
    catalogue = read_json(HERE / "particles.json")
    group = {k: v for k, v in catalogue.items() if k != "families"}
    assignments = getattr(candidate, "PARTICLE_IDENTIFICATIONS", {})
    rows = []
    for family in catalogue["families"]:
        for key, name, distinction in family["types"]:
            record = assignments.get(key, {"status": "open", "note": "No species assignment is registered.", "evidence": []})
            rows.append({"id": "SM-" + key, "title": name, "domain": family["name"],
                         "requires": family["criterion"] + " " + distinction, "sources": catalogue["sources"],
                         "results": {candidate.NAME: record}})
    if set(assignments) - {r['id'][3:] for r in rows}:
        raise ValueError(f"{candidate.NAME}: unknown particle type")
    return {**group, "criteria": rows}


def definitions():
    return [group for file in FILES for group in read_json(HERE / file)["groups"]]


def for_candidate(candidate, register_ids, known_laws, shared=None):
    groups = [*(definitions() if shared is None else shared), species_group(candidate)]
    seen_groups, seen_ids = set(), set()
    assessments = []
    for group in groups:
        group_id = group["id"]
        if group_id in seen_groups or not group["criteria"]:
            raise ValueError(f"Duplicate or empty requirement group {group_id}")
        seen_groups.add(group_id)
        for source in group.get("sources", []):
            check_source(source)
        rows = []
        for definition in group["criteria"]:
            key = definition["id"]
            if key in seen_ids or not definition["requires"].strip() or not definition["sources"]:
                raise ValueError(f"Duplicate or incomplete requirement {key}")
            seen_ids.add(key)
            declared = definition["results"]
            if set(declared) - set(known_laws) or set(definition.get("premises", [])) - register_ids:
                raise ValueError(f"{key}: unregistered candidate or premise")
            record = declared.get(candidate.NAME, {"status": "open", "note": "Awaiting candidate-specific assessment.", "evidence": []})
            if (record["status"] not in STATUS or not record["note"].strip()
                    or (record["status"] != "open" and not record["evidence"])):
                raise ValueError(f"{key}: unsupported verdict for {candidate.NAME}")
            for source in [*definition["sources"], *record["evidence"]]:
                check_source(source)
            rows.append({**{k: v for k, v in definition.items() if k != "results"}, **record})
        assessments.append({**{k: v for k, v in group.items() if k != "criteria"}, "criteria": rows})
    return assessments
