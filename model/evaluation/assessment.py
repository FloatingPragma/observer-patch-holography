"""Requirement definitions, candidate evidence and quantitative readouts for the public comparison."""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
STATUS = {"demonstrated": "Demonstrated", "partial": "Partial", "open": "Open", "counterexample": "Counterexample"}
FILES = ("theorem_criteria.json", "observer_criteria.json", "openwave_criteria.json", "gravitational_particle.json")
READOUTS = {
    "null_hopping": ("Calculated first-read probability that no seam reports the released pair after one exchange-time unit on one cell",
                     "model/candidate_laws/null_hopping/results.json", "kappa_0/no_pair", None),
    "private_record_walk": ("Home-pair probability with eight sharp comparison signs at exposure 37",
                            "model/candidate_laws/private_record_walk/results.json",
                            "appointments/37/sharp_home", None),
    "interval_echo": ("Traveller/home ratio for the designed closing filter on alternating 2- and 6-window handoffs",
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


def check_source(source, *, root=None, local=False):
    if not isinstance(source, str) or not source.strip():
        raise ValueError(f"Missing evidence path: {source!r}")
    parsed = urlsplit(source)
    if not local and parsed.scheme == "https" and parsed.hostname and parsed.path and not parsed.username:
        return
    root = ROOT if root is None else root
    path = root / parsed.path
    if (parsed.scheme or parsed.netloc or parsed.query or Path(parsed.path).is_absolute()
            or not path.resolve().is_relative_to(root) or not path.is_file()):
        raise ValueError(f"Unknown repository evidence: {source}")
    if re.match(r'^L\d', parsed.fragment):
        lines = re.fullmatch(r'L([0-9]+)(?:-L([0-9]+))?', parsed.fragment)
        if not lines:
            raise ValueError(f"Invalid line reference: {source}")
        first, last = int(lines[1]), int(lines[2] or lines[1])
        with path.open('rb') as stream:
            if not 1 <= first <= last <= sum(1 for _ in stream):
                raise ValueError(f"Line reference outside evidence: {source}")
    elif parsed.fragment and path.suffix.lower() == '.md':
        headings = re.findall(r'^#{1,6}\s+(.+?)\s*$', path.read_text(encoding='utf-8'), re.M)
        anchors = {re.sub(r'[^\w -]', '', title.lower()).replace(' ', '-') for title in headings}
        if unquote(parsed.fragment) not in anchors:
            raise ValueError(f"Unknown Markdown evidence heading: {source}")


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
    group = {k: v for k, v in catalogue.items() if k not in {"families", "observations"}}
    observations = catalogue.get("observations", {}).get(candidate.NAME, {})
    assignments = getattr(candidate, "PARTICLE_IDENTIFICATIONS", {})
    rows = []
    for family in catalogue["families"]:
        for key, name, distinction in family["types"]:
            record = assignments.get(key, {"status": "open", "note": "No species assignment is registered.", "evidence": []})
            if key in observations:
                observed = observations[key]
                record = {**record, "note": record["note"] + " " + observed["note"],
                          "evidence": list(dict.fromkeys([*record["evidence"], *observed["evidence"]]))}
                if 'comparison' in observed:
                    previous = record.get('comparison', {})
                    revised = observed['comparison']
                    if any(previous.get(field) != revised.get(field) for field in ('source', 'law')):
                        raise ValueError(f"{key}: comparison revision must identify the original law and source")
                    record['comparison'] = revised
            rows.append({"id": "SM-" + key, "title": name, "domain": family["name"],
                         "requires": family["criterion"] + " " + distinction, "sources": catalogue["sources"],
                         "results": {candidate.NAME: record}})
    if (set(assignments) | set(observations)) - {r['id'][3:] for r in rows}:
        raise ValueError(f"{candidate.NAME}: unknown particle type")
    return {**group, "criteria": rows}


def comparison_studies(record):
    additional = record.get("additional_comparisons", [])
    if not isinstance(additional, list):
        raise ValueError("Additional comparisons must be a list of studies")
    primary = [record["comparison"]] if "comparison" in record else []
    return primary + additional


def definitions():
    return [group for file in FILES for group in read_json(HERE / file)["groups"]]


EVIDENCE_KINDS = {'record': 'Observer records', 'mixed': 'Records on supplied settings',
                  'state': 'Operator or state calculation', 'construction': 'Construction identity'}


def check_observer_record(key, support):
    if support['kind'] == 'record' and any(
        not isinstance(support.get(field), str) or not support[field].strip()
        for field in ('observer', 'instrument', 'clock', 'outcomes')
    ):
        raise ValueError(f"{key}: identify the observer, instrument, clock and compared outcomes")


def validate_support(key, record, definition, premise_status):
    evidence = record['evidence']
    if not isinstance(evidence, list) or len(evidence) != len(set(evidence)):
        raise ValueError(f"{key}: evidence must list distinct sources")
    filing = record.get('filing', {})
    for source in evidence:
        # The upstream decision or review can document filing status; the candidate's physical evidence is local.
        if source.startswith('https://') and source not in (*definition['sources'], filing.get('review')):
            raise ValueError(f"{key}: external evidence must be a listed criterion source or filing review")
        check_source(source)
    if record['status'] != 'open' and not any(not s.startswith('https://') for s in evidence):
        raise ValueError(f"{key}: candidate assessment needs repository evidence")
    if 'full_credit' in record:
        if record['status'] != 'demonstrated':
            raise ValueError(f"{key}: full-credit certificate supplied for a lower grade")
        support = record['full_credit']
        if (not isinstance(support, dict) or support.get('kind') not in EVIDENCE_KINDS
                or not isinstance(support.get('scope'), str) or not support['scope'].strip()
                or support.get('source') not in evidence or support['source'].startswith('https://')):
            raise ValueError(f"{key}: full credit needs a classified, scoped local certificate or receipt")
        path = Path(support['source'].partition('#')[0])
        if path.suffix.lower() not in {'.py', '.lean', '.json', '.gz', '.npz', '.csv', '.npy', '.jl'}:
            raise ValueError(f"{key}: full credit needs an executable certificate or recorded data")
        check_observer_record(key, support)
        if re.match(r'^(OW-|SM-|OBS-|G[1-6]$)', key) and support['kind'] != 'record':
            raise ValueError(f"{key}: full physical identification or observer credit needs observer records")
        pending = [p for p in definition.get('premises', []) if premise_status.get(p) != 'fulfilled']
        if pending:
            raise ValueError(f"{key}: full credit has undischarged candidate premises {pending}")


def for_candidate(candidate, premise_status, known_laws, shared=None):
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
            if set(declared) - set(known_laws) or set(definition.get("premises", [])) - set(premise_status):
                raise ValueError(f"{key}: unregistered candidate or premise")
            record = declared.get(candidate.NAME, {"status": "open", "note": "Awaiting candidate-specific assessment.", "evidence": []})
            if (record["status"] not in STATUS or not record["note"].strip()
                    or (record["status"] != "open" and not record["evidence"])):
                raise ValueError(f"{key}: unsupported verdict for {candidate.NAME}")
            if "demonstrated_measurement" in record:
                measurement = record["demonstrated_measurement"]
                if (not isinstance(measurement, dict) or not isinstance(measurement.get("label"), str)
                        or measurement.get('kind') not in EVIDENCE_KINDS
                        or not measurement["label"].strip() or measurement.get("source") not in record["evidence"]):
                    raise ValueError(f"{key}: demonstrated measurement needs a label and a listed evidence source")
                check_observer_record(key, measurement)
            for comparison in comparison_studies(record):
                if (not isinstance(comparison, dict) or
                        any(not isinstance(comparison.get(field), str) or not comparison[field].strip()
                            for field in ("law", "finding", "source"))):
                    raise ValueError(f"{key}: a comparison needs its law, finding and source")
                check_source(comparison["source"])
            validate_support(key, record, definition, premise_status)
            for source in definition["sources"]:
                check_source(source)
            rows.append({**{k: v for k, v in definition.items() if k != "results"}, **record})
        assessments.append({**{k: v for k, v in group.items() if k != "criteria"}, "criteria": rows})
    return assessments
