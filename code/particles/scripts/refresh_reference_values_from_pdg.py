#!/usr/bin/env python3
"""Refresh the pinned particle reference set from the official PDG API.

This keeps `/particles` self-contained: we pin a compact machine-readable
reference file locally, then build status tables from that file without
needing network access on every run.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import urllib.request
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "data" / "particle_reference_values.json"
PDG_API_BASE = "https://pdgapi.lbl.gov/summaries"
PDG_INFO_URL = "https://pdg.lbl.gov/api"


def _unit_scale_to_gev(unit: str) -> Optional[float]:
    normalized = unit.strip()
    if normalized == "GeV":
        return 1.0
    if normalized == "MeV":
        return 1.0e-3
    if normalized == "eV":
        return 1.0e-9
    return None


def _fetch_json(url: str) -> Dict[str, Any]:
    with urllib.request.urlopen(url, timeout=30) as response:
        return json.load(response)


def _fetch_summary(summary_id: str) -> Dict[str, Any]:
    url = f"{PDG_API_BASE}/{summary_id}"
    payload = _fetch_json(url)
    values = payload.get("pdg_values") or []
    if not values:
        raise RuntimeError(f"PDG summary {summary_id} returned no pdg_values")
    value = values[0]
    unit = value.get("unit")
    scale = _unit_scale_to_gev(unit) if unit else None
    value_gev = float(value["value"]) * scale if scale is not None and value.get("value") is not None else None
    err_plus = value.get("error_positive")
    err_minus = value.get("error_negative")
    err_plus_gev = float(err_plus) * scale if scale is not None and err_plus is not None else None
    err_minus_gev = float(err_minus) * scale if scale is not None and err_minus is not None else None
    return {
        "source": {
            "label": "PDG 2025 API",
            "summary_id": summary_id,
            "url": url,
            "info_url": PDG_INFO_URL,
        },
        "description": payload.get("description"),
        "reference_kind": "upper_limit" if value.get("is_upper_limit") else "value",
        "display": (
            f"{value.get('value_text')} {unit}".strip()
            if value.get("value_text") and unit
            else value.get("value_text")
        ),
        "value_gev": value_gev,
        "error_plus_gev": err_plus_gev,
        "error_minus_gev": err_minus_gev,
        "raw_value": value.get("value"),
        "raw_unit": unit,
        "value_text": value.get("value_text"),
        "comment": value.get("comment"),
        "type": value.get("type"),
        "is_limit": bool(value.get("is_limit", False)),
        "is_upper_limit": bool(value.get("is_upper_limit", False)),
        "confidence_level": value.get("confidence_level"),
    }


def _manual_reference(
    *,
    label: str,
    reference_kind: str,
    display: str,
    notes: str,
    url: str,
    value_gev: Optional[float] = None,
) -> Dict[str, Any]:
    return {
        "source": {
            "label": label,
            "url": url,
            "info_url": PDG_INFO_URL,
        },
        "description": None,
        "reference_kind": reference_kind,
        "display": display,
        "value_gev": value_gev,
        "error_plus_gev": None,
        "error_minus_gev": None,
        "raw_value": None,
        "raw_unit": None,
        "value_text": None,
        "comment": notes,
        "type": None,
        "is_limit": False,
        "is_upper_limit": False,
        "confidence_level": None,
    }


def build_reference_payload() -> Dict[str, Any]:
    fetched: Dict[str, Dict[str, Any]] = {
        "photon": _fetch_summary("S000M"),
        "w_boson": _fetch_summary("S043M"),
        "z_boson": _fetch_summary("S044M"),
        "higgs": _fetch_summary("S126M"),
        "electron": _fetch_summary("S003M"),
        "muon": _fetch_summary("S004M"),
        "tau": _fetch_summary("S035M"),
        "up_quark": _fetch_summary("Q002M"),
        "down_quark": _fetch_summary("Q001M"),
        "strange_quark": _fetch_summary("Q003M"),
        "charm_quark": _fetch_summary("Q004M"),
        "bottom_quark": _fetch_summary("Q005M"),
        "top_quark": _fetch_summary("Q007TP4"),
        "top_quark_direct_aux": _fetch_summary("Q007TP"),
        "proton": _fetch_summary("S016M"),
        "neutron": _fetch_summary("S017M"),
        "neutral_pion": _fetch_summary("S009M"),
        "rho_770_0": _fetch_summary("M009M0"),
    }

    manual: Dict[str, Dict[str, Any]] = {
        "gluon": _manual_reference(
            label="PDG 2025 particle listings context",
            reference_kind="no_direct_free_particle_mass_measurement",
            display="no direct free-particle mass measurement",
            notes="Free gluons are confined; there is no direct measured gluon rest mass entry comparable to other particle masses.",
            url="https://pdg.lbl.gov/2025/listings/particle_properties.html",
        ),
        "graviton": _manual_reference(
            label="GW dispersion observational context",
            reference_kind="upper_limit",
            display="<1e-32 GeV",
            notes="No direct graviton rest-mass measurement exists. Gravitational-wave propagation constrains any hard graviton mass to be extremely small; OPH predicts an exact structural zero on the dynamical-metric branch.",
            url="https://floatingpragma.io",
            value_gev=1.0e-32,
        ),
        "electron_neutrino": _manual_reference(
            label="PDG 2025 neutrino properties",
            reference_kind="not_directly_measured",
            display="not directly measured",
            notes="Individual flavor neutrino masses are not directly measured as standalone masses; PDG quotes oscillation data, effective masses, and cosmological bounds instead.",
            url="https://pdg.lbl.gov/2025/reviews/contents_sports.html",
        ),
        "muon_neutrino": _manual_reference(
            label="PDG 2025 neutrino properties",
            reference_kind="not_directly_measured",
            display="not directly measured",
            notes="Individual flavor neutrino masses are not directly measured as standalone masses; PDG quotes oscillation data, effective masses, and cosmological bounds instead.",
            url="https://pdg.lbl.gov/2025/reviews/contents_sports.html",
        ),
        "tau_neutrino": _manual_reference(
            label="PDG 2025 neutrino properties",
            reference_kind="not_directly_measured",
            display="not directly measured",
            notes="Individual flavor neutrino masses are not directly measured as standalone masses; PDG quotes oscillation data, effective masses, and cosmological bounds instead.",
            url="https://pdg.lbl.gov/2025/reviews/contents_sports.html",
        ),
    }

    return {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source": {
            "label": "Particle Data Group",
            "edition": "2025",
            "api_info_url": PDG_INFO_URL,
        },
        "entries": {
            **fetched,
            **manual,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Refresh the pinned PDG reference set used by /particles.")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Output JSON path.")
    args = parser.parse_args()

    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    payload = build_reference_payload()
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
