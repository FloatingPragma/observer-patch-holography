from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from oph_fpe.claims import CONTINUATION, with_claim_metadata


FINITE_GAP_RECEIPT = "FINITE_NONABELIAN_GAUGE_GAP_DIAGNOSTIC_RECEIPT"
REPORT_SCHEMA = "oph_yang_mills_gap_certificate_v0"


def yang_mills_gap_certificate_report(
    *,
    lattice_size: int = 2,
    sweeps: int = 16,
    beta: float = 2.2,
    proposal_width: float = 0.35,
    seed: int = 20260706,
    transition_bins: int = 8,
    refinement_lattice_sizes: Iterable[int] | None = (2, 3),
    refinement_sweeps: int | None = None,
    continuum_certificate: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Finite SU(2) diagnostic lane for the OPH Yang-Mills-gap claim boundary.

    This emits useful compact-simple nonabelian lattice data and a finite
    transfer/repair-gap proxy. It deliberately does not promote the Clay/Jaffe-
    Witten mass-gap claim; that requires the continuum OS/gauge certificate
    named in the paper.
    """

    lattice_size = max(2, int(lattice_size))
    sweeps = max(4, int(sweeps))
    beta = _positive(beta, 2.2)
    proposal_width = _positive(proposal_width, 0.35)
    seed = int(seed)
    transition_bins = max(3, int(transition_bins))

    main = _run_su2_wilson_chain(
        lattice_size=lattice_size,
        sweeps=sweeps,
        beta=beta,
        proposal_width=proposal_width,
        seed=seed,
    )
    replay = _run_su2_wilson_chain(
        lattice_size=lattice_size,
        sweeps=sweeps,
        beta=beta,
        proposal_width=proposal_width,
        seed=seed,
    )
    deterministic_replay = (
        main["plaquette_trace"] == replay["plaquette_trace"]
        and main["wilson_loop_trace"] == replay["wilson_loop_trace"]
        and main["polyakov_loop_trace"] == replay["polyakov_loop_trace"]
    )

    transfer = _finite_transfer_gap(main["plaquette_trace"], transition_bins=transition_bins)
    reflection = _reflection_gram_proxy(main)
    finite_nontriviality = bool(float(main["plaquette_variance"]) > 1.0e-12)
    finite_receipt = bool(
        deterministic_replay
        and finite_nontriviality
        and transfer["finite_transfer_gap_proxy_receipt"]
        and reflection["finite_reflection_gram_proxy_receipt"]
    )

    refinement_rows = _refinement_rows(
        main=main,
        lattice_sizes=tuple(refinement_lattice_sizes or (lattice_size,)),
        sweeps=refinement_sweeps or max(4, min(sweeps, 8)),
        beta=beta,
        proposal_width=proposal_width,
        seed=seed,
        transition_bins=transition_bins,
    )
    positive_gap_values = [
        float(row["finite_transfer_gap_estimate"])
        for row in refinement_rows
        if row.get("finite_transfer_gap_estimate") is not None
        and math.isfinite(float(row["finite_transfer_gap_estimate"]))
        and float(row["finite_transfer_gap_estimate"]) > 0.0
    ]
    finite_gap_floor = min(positive_gap_values) if positive_gap_values else None
    continuum = _continuum_certificate_status(continuum_certificate or {})
    promotion = _promotion_status(
        finite_receipt=finite_receipt,
        finite_gap_floor=finite_gap_floor,
        continuum=continuum,
    )

    report: dict[str, Any] = {
        "schema": REPORT_SCHEMA,
        "mode": "finite_su2_yang_mills_gap_certificate_lane",
        "paper_source": "markdown/yang_mills_gap_clay_problem.md",
        "gauge_group": {
            "name": "SU(2)",
            "compact": True,
            "simple": True,
            "nonabelian": True,
            "representation": "unit quaternions; normalized trace is quaternion scalar part",
        },
        "regulator": {
            "dimension": 4,
            "lattice_size": int(lattice_size),
            "site_count": int(lattice_size**4),
            "link_count": int((lattice_size**4) * 4),
            "plaquette_count": int((lattice_size**4) * 6),
            "boundary": "periodic_hypercubic",
            "beta": float(beta),
            "sweeps": int(sweeps),
            "proposal_width": float(proposal_width),
            "seed": int(seed),
            "transition_bins": int(transition_bins),
        },
        "lattice_gauge_stage": {
            "compact_simple_nonabelian_reference": True,
            "su2_reference": True,
            "su3_reference": False,
            "u1_reference": False,
            "four_dimensional_wilson_lattice": True,
            "continuum_yang_mills_theory_constructed": False,
        },
        "finite_lattice_diagnostics": {
            "mean_plaquette": main["mean_plaquette"],
            "plaquette_variance": main["plaquette_variance"],
            "acceptance_rate": main["acceptance_rate"],
            "polyakov_abs_mean_final": main["polyakov_abs_mean_final"],
            "finite_nontriviality_proxy_receipt": finite_nontriviality,
            "canonical_serial_chain_replay_receipt": deterministic_replay,
            "diagnostic_only_not_stationary_law_proof": True,
        },
        "plaquette_trace": main["plaquette_trace"],
        "wilson_loop_trace": main["wilson_loop_trace"],
        "polyakov_loop_trace": main["polyakov_loop_trace"],
        "orientation_plaquette_rows": main["orientation_plaquette_rows"],
        "finite_transfer_gap_diagnostic": transfer,
        "reflection_positivity_proxy": reflection,
        "refinement_gap_rows": refinement_rows,
        "finite_positive_gap_floor_estimate": finite_gap_floor,
        "continuum_certificate": continuum,
        "promotion_status": promotion,
        FINITE_GAP_RECEIPT: finite_receipt,
        "finite_nonabelian_gauge_gap_diagnostic_receipt": finite_receipt,
        "finite_repair_gap_proxy_receipt": bool(transfer["finite_transfer_gap_proxy_receipt"]),
        "continuum_yang_mills_mass_gap_receipt": False,
        "YANG_MILLS_GAP_REPRODUCED_RECEIPT": False,
        "CLAY_YANG_MILLS_GAP_RECEIPT": False,
        "yang_mills_identification": "conditional",
        "explicit_nonclaims": [
            "Clay-admissible Yang-Mills construction",
            "continuum Osterwalder-Schrader reconstruction",
            "physical glueball mass prediction",
            "SU(3) QCD mass-gap computation",
            "proof that the finite transfer proxy survives refinement",
        ],
        "claim_boundary": (
            "Finite SU(2) compact-simple nonabelian lattice diagnostic. A positive finite diagnostic "
            "receipt means this run emitted nonabelian Wilson-lattice data, deterministic replay, a "
            "finite transfer-gap proxy, and a finite reflection-Gram proxy. It is not a reproduction "
            "of the Yang-Mills mass gap. The Clay/Jaffe-Witten claim remains closed until the "
            "support-visible compact-gauge continuum certificate supplies Schwinger convergence, "
            "reflection positivity, Euclidean covariance/locality, nontriviality, and transfer/"
            "intertwiner convergence."
        ),
    }
    return with_claim_metadata(
        report,
        claim_level=CONTINUATION,
        receipt=FINITE_GAP_RECEIPT,
        observable_id="finite_su2_yang_mills_gap_diagnostic",
        fit_objective="nonabelian_compact_gauge_continuum_certificate_lane",
    )


def write_yang_mills_gap_certificate_report(
    out: Path,
    *,
    lattice_size: int = 2,
    sweeps: int = 16,
    beta: float = 2.2,
    proposal_width: float = 0.35,
    seed: int = 20260706,
    transition_bins: int = 8,
    refinement_lattice_sizes: Iterable[int] | None = (2, 3),
    refinement_sweeps: int | None = None,
    continuum_certificate: dict[str, Any] | None = None,
) -> dict[str, Any]:
    report = yang_mills_gap_certificate_report(
        lattice_size=lattice_size,
        sweeps=sweeps,
        beta=beta,
        proposal_width=proposal_width,
        seed=seed,
        transition_bins=transition_bins,
        refinement_lattice_sizes=refinement_lattice_sizes,
        refinement_sweeps=refinement_sweeps,
        continuum_certificate=continuum_certificate,
    )
    destination = Path(out)
    if destination.suffix.lower() != ".json":
        destination = destination / "yang_mills_gap_certificate_report.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    stem = destination.with_suffix("")
    _write_csv(stem.with_name(stem.name + "_plaquette_trace.csv"), report["plaquette_trace"])
    _write_csv(stem.with_name(stem.name + "_wilson_loop_trace.csv"), report["wilson_loop_trace"])
    _write_csv(stem.with_name(stem.name + "_polyakov_loop_trace.csv"), report["polyakov_loop_trace"])
    _write_csv(stem.with_name(stem.name + "_orientation_plaquettes.csv"), report["orientation_plaquette_rows"])
    _write_csv(stem.with_name(stem.name + "_refinement_gap.csv"), report["refinement_gap_rows"])
    _write_csv(stem.with_name(stem.name + "_promotion_gates.csv"), _promotion_rows(report))
    stem.with_name(stem.name + ".md").write_text(_markdown_report(report), encoding="utf-8")
    report["report_path"] = str(destination)
    return report


def _run_su2_wilson_chain(
    *,
    lattice_size: int,
    sweeps: int,
    beta: float,
    proposal_width: float,
    seed: int,
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    links = _random_su2((lattice_size,) * 4 + (4,), rng)
    current = _wilson_log_weight(links, beta=beta)
    total_accepted = 0
    total_updates = 0
    plaquette_trace: list[dict[str, Any]] = []
    wilson_loop_trace: list[dict[str, Any]] = []
    polyakov_loop_trace: list[dict[str, Any]] = []
    orientation_rows: list[dict[str, Any]] = []

    for sweep in range(sweeps):
        accepted = 0
        updates = 0
        for site in np.ndindex((lattice_size,) * 4):
            for mu in range(4):
                key = site + (mu,)
                old = np.array(links[key], copy=True)
                links[key] = _q_mul(_random_near_identity(rng, proposal_width), old)
                candidate = _wilson_log_weight(links, beta=beta)
                accept_u = float(rng.random())
                if math.log(max(accept_u, 1.0e-300)) <= min(0.0, candidate - current):
                    current = candidate
                    accepted += 1
                else:
                    links[key] = old
                updates += 1
        total_accepted += accepted
        total_updates += updates
        stats = _plaquette_statistics(links)
        acceptance_rate = float(total_accepted / max(total_updates, 1))
        row = {
            "sweep": int(sweep),
            "mean_plaquette": stats["mean"],
            "plaquette_variance": stats["variance"],
            "action_density": float(current / max((lattice_size**4) * 6, 1)),
            "acceptance_rate": acceptance_rate,
        }
        plaquette_trace.append(row)
        wilson_loop_trace.append(
            {
                "sweep": int(sweep),
                "loop": "plaquette_1x1",
                "mean_normalized_trace": stats["mean"],
                "variance": stats["variance"],
            }
        )
        polyakov_loop_trace.append(
            {
                "sweep": int(sweep),
                "loop": "time_polyakov_abs",
                "mean_abs_normalized_trace": _polyakov_abs_mean(links),
            }
        )
        for label, value in stats["orientation_means"].items():
            orientation_rows.append(
                {
                    "sweep": int(sweep),
                    "orientation": label,
                    "mean_plaquette": value,
                }
            )

    series = np.asarray([row["mean_plaquette"] for row in plaquette_trace], dtype=float)
    return {
        "lattice_size": int(lattice_size),
        "sweeps": int(sweeps),
        "beta": float(beta),
        "proposal_width": float(proposal_width),
        "seed": int(seed),
        "plaquette_trace": plaquette_trace,
        "wilson_loop_trace": wilson_loop_trace,
        "polyakov_loop_trace": polyakov_loop_trace,
        "orientation_plaquette_rows": orientation_rows,
        "mean_plaquette": float(np.mean(series)) if series.size else None,
        "plaquette_variance": float(np.var(series)) if series.size else None,
        "acceptance_rate": float(total_accepted / max(total_updates, 1)),
        "polyakov_abs_mean_final": polyakov_loop_trace[-1]["mean_abs_normalized_trace"] if polyakov_loop_trace else None,
    }


def _refinement_rows(
    *,
    main: dict[str, Any],
    lattice_sizes: tuple[int, ...],
    sweeps: int,
    beta: float,
    proposal_width: float,
    seed: int,
    transition_bins: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    seen: set[int] = set()
    for size in lattice_sizes:
        size = max(2, int(size))
        if size in seen:
            continue
        seen.add(size)
        run = main if size == int(main["lattice_size"]) else _run_su2_wilson_chain(
            lattice_size=size,
            sweeps=sweeps,
            beta=beta,
            proposal_width=proposal_width,
            seed=seed + 9973 * size,
        )
        transfer = _finite_transfer_gap(run["plaquette_trace"], transition_bins=transition_bins)
        rows.append(
            {
                "lattice_size": int(size),
                "site_count": int(size**4),
                "sweeps": int(run["sweeps"]),
                "mean_plaquette": run["mean_plaquette"],
                "plaquette_variance": run["plaquette_variance"],
                "acceptance_rate": run["acceptance_rate"],
                "finite_transfer_gap_estimate": transfer["spectral_gap_estimate"],
                "lambda2_abs": transfer["lambda2_abs"],
                "screening_mass_proxy": transfer["screening_mass_proxy"],
                "finite_transfer_gap_proxy_receipt": transfer["finite_transfer_gap_proxy_receipt"],
            }
        )
    return rows


def _finite_transfer_gap(rows: list[dict[str, Any]], *, transition_bins: int) -> dict[str, Any]:
    series = np.asarray([float(row["mean_plaquette"]) for row in rows], dtype=float)
    if series.size < 3:
        return {
            "method": "plaquette_trace_binned_markov_proxy",
            "transition_bins": int(transition_bins),
            "spectral_gap_estimate": None,
            "lambda2_abs": None,
            "screening_mass_proxy": None,
            "row_stochastic_error": None,
            "nondegenerate_transition_support": False,
            "finite_transfer_gap_proxy_receipt": False,
            "claim_boundary": "insufficient finite trace length for a transfer-gap proxy",
        }
    lo = float(np.min(series))
    hi = float(np.max(series))
    if hi - lo < 1.0e-12:
        binned = np.zeros(series.size, dtype=int)
    else:
        edges = np.linspace(lo, hi, transition_bins + 1)
        binned = np.clip(np.digitize(series, edges[1:-1]), 0, transition_bins - 1)
    counts = np.full((transition_bins, transition_bins), 1.0e-3, dtype=float)
    for left, right in zip(binned[:-1], binned[1:]):
        counts[int(left), int(right)] += 1.0
    transition = counts / counts.sum(axis=1, keepdims=True)
    eig = np.linalg.eigvals(transition)
    ordered = sorted((abs(complex(value)) for value in eig), reverse=True)
    lambda2 = float(ordered[1]) if len(ordered) > 1 else 0.0
    gap = max(0.0, float(1.0 - lambda2))
    centered = series - float(np.mean(series))
    c0 = float(np.mean(centered * centered))
    c1 = float(np.mean(centered[:-1] * centered[1:])) if centered.size > 1 else 0.0
    screening_mass = -math.log(c1 / c0) if c0 > 0.0 and 0.0 < c1 < c0 else None
    row_error = float(np.max(np.abs(transition.sum(axis=1) - 1.0)))
    nondegenerate = bool(len(set(int(value) for value in binned)) >= 2)
    receipt = bool(nondegenerate and gap > 0.0 and row_error <= 1.0e-12)
    return {
        "method": "plaquette_trace_binned_markov_proxy",
        "transition_bins": int(transition_bins),
        "spectral_gap_estimate": gap,
        "lambda2_abs": lambda2,
        "screening_mass_proxy": screening_mass,
        "row_stochastic_error": row_error,
        "nondegenerate_transition_support": nondegenerate,
        "finite_transfer_gap_proxy_receipt": receipt,
        "claim_boundary": (
            "Finite binned Markov proxy from plaquette trace. This is sampler/regulator evidence only, "
            "not a continuum Yang-Mills Hamiltonian spectral gap."
        ),
    }


def _reflection_gram_proxy(run: dict[str, Any]) -> dict[str, Any]:
    plaquette = np.asarray([row["mean_plaquette"] for row in run["plaquette_trace"]], dtype=float)
    action = np.asarray([row["action_density"] for row in run["plaquette_trace"]], dtype=float)
    polyakov = np.asarray(
        [row["mean_abs_normalized_trace"] for row in run["polyakov_loop_trace"]],
        dtype=float,
    )
    if plaquette.size == 0:
        return {
            "method": "finite_feature_reflection_gram_proxy",
            "reflection_gram_lower_bound": None,
            "finite_reflection_gram_proxy_receipt": False,
        }
    features = np.vstack([np.ones_like(plaquette), plaquette, action, polyakov]).T
    gram = (features.T @ features) / max(1, features.shape[0])
    lower = float(np.min(np.linalg.eigvalsh(gram)))
    return {
        "method": "finite_feature_reflection_gram_proxy",
        "feature_labels": ["constant", "mean_plaquette", "action_density", "time_polyakov_abs"],
        "reflection_gram_lower_bound": lower,
        "finite_reflection_gram_proxy_receipt": bool(lower >= -1.0e-10),
        "claim_boundary": (
            "Positive finite feature Gram for reflection-style diagnostics. This is not the continuum "
            "Osterwalder-Schrader reflection-positivity certificate."
        ),
    }


def _continuum_certificate_status(certificate: dict[str, Any]) -> dict[str, Any]:
    required = {
        "support_visible_extraction_receipt": bool(certificate.get("support_visible_extraction_receipt", False)),
        "renormalized_schwinger_convergence_receipt": bool(
            certificate.get("renormalized_schwinger_convergence_receipt", False)
        ),
        "reflection_positivity_receipt": bool(certificate.get("reflection_positivity_receipt", False)),
        "euclidean_covariance_locality_receipt": bool(
            certificate.get("euclidean_covariance_locality_receipt", False)
        ),
        "nontriviality_receipt": bool(certificate.get("nontriviality_receipt", False)),
        "transfer_intertwiner_convergence_receipt": bool(
            certificate.get("transfer_intertwiner_convergence_receipt", False)
        ),
    }
    complete = all(required.values())
    return {
        "schema": "oph_compact_gauge_continuum_certificate_status_v0",
        "provided": bool(certificate),
        "required_receipts": required,
        "candidate_complete": complete,
        "external_certificate_hash": certificate.get("external_certificate_hash"),
        "trusted_external_verification": False,
        "continuum_certificate_receipt": False,
        "missing": [key for key, value in required.items() if not value],
        "claim_boundary": (
            "These are the Assumption-14-style continuum-certificate slots. This simulator lane records "
            "whether fields were supplied but does not independently verify the continuum theorem."
        ),
    }


def _promotion_status(
    *,
    finite_receipt: bool,
    finite_gap_floor: float | None,
    continuum: dict[str, Any],
) -> dict[str, Any]:
    reasons = []
    if not finite_receipt:
        reasons.append("finite SU(2) nonabelian diagnostic did not pass")
    if finite_gap_floor is None:
        reasons.append("no positive finite transfer-gap floor was estimated across refinement rows")
    reasons.extend(
        f"missing continuum certificate field: {name}"
        for name in continuum.get("missing", [])
    )
    reasons.append("Clay/Yang-Mills promotion remains disabled in this simulator lane")
    return {
        "finite_nonabelian_regulator": "pass" if finite_receipt else "fail",
        "finite_positive_gap_floor": "pass" if finite_gap_floor is not None else "fail",
        "continuum_certificate": "pending",
        "os_reconstruction": "pending",
        "yang_mills_identification": "conditional",
        "yang_mills_mass_gap": "not_promoted",
        "reasons": reasons,
    }


def _plaquette_statistics(links: np.ndarray) -> dict[str, Any]:
    values = []
    orientation: dict[str, list[float]] = {}
    for site in np.ndindex(links.shape[:4]):
        for mu in range(4):
            for nu in range(mu + 1, 4):
                value = float(_plaquette(links, site, mu, nu)[0])
                values.append(value)
                orientation.setdefault(f"{mu}{nu}", []).append(value)
    arr = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(arr)),
        "variance": float(np.var(arr)),
        "orientation_means": {
            key: float(np.mean(np.asarray(rows, dtype=float)))
            for key, rows in sorted(orientation.items())
        },
    }


def _wilson_log_weight(links: np.ndarray, *, beta: float) -> float:
    total = 0.0
    for site in np.ndindex(links.shape[:4]):
        for mu in range(4):
            for nu in range(mu + 1, 4):
                total += float(_plaquette(links, site, mu, nu)[0])
    return float(beta) * total


def _plaquette(links: np.ndarray, site: tuple[int, int, int, int], mu: int, nu: int) -> np.ndarray:
    size = int(links.shape[0])
    site_mu = _shift(site, mu, size)
    site_nu = _shift(site, nu, size)
    return _q_mul(
        _q_mul(_q_mul(links[site + (mu,)], links[site_mu + (nu,)]), _q_conj(links[site_nu + (mu,)])),
        _q_conj(links[site + (nu,)]),
    )


def _polyakov_abs_mean(links: np.ndarray) -> float:
    size = int(links.shape[0])
    rows = []
    for spatial in np.ndindex((size, size, size)):
        q = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
        for t in range(size):
            q = _q_mul(q, links[(spatial[0], spatial[1], spatial[2], t, 3)])
        rows.append(abs(float(q[0])))
    return float(np.mean(np.asarray(rows, dtype=float)))


def _shift(site: tuple[int, int, int, int], axis: int, size: int) -> tuple[int, int, int, int]:
    values = list(site)
    values[axis] = (values[axis] + 1) % size
    return tuple(values)  # type: ignore[return-value]


def _random_su2(shape: tuple[int, ...], rng: np.random.Generator) -> np.ndarray:
    values = rng.normal(size=shape + (4,))
    norms = np.linalg.norm(values, axis=-1, keepdims=True)
    return values / np.maximum(norms, 1.0e-12)


def _random_near_identity(rng: np.random.Generator, width: float) -> np.ndarray:
    axis = rng.normal(size=3)
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm <= 1.0e-12:
        axis = np.array([1.0, 0.0, 0.0], dtype=float)
    else:
        axis = axis / axis_norm
    angle = float(rng.uniform(-width, width))
    half = 0.5 * angle
    return np.array([math.cos(half), *(math.sin(half) * axis)], dtype=float)


def _q_mul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    a, b, c, d = left
    e, f, g, h = right
    return np.array(
        [
            a * e - b * f - c * g - d * h,
            a * f + b * e + c * h - d * g,
            a * g - b * h + c * e + d * f,
            a * h + b * g - c * f + d * e,
        ],
        dtype=float,
    )


def _q_conj(value: np.ndarray) -> np.ndarray:
    return np.array([value[0], -value[1], -value[2], -value[3]], dtype=float)


def _positive(value: Any, default: float) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return float(default)
    return parsed if math.isfinite(parsed) and parsed > 0.0 else float(default)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or ["empty"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _promotion_rows(report: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    promotion = report.get("promotion_status") if isinstance(report.get("promotion_status"), dict) else {}
    for key, value in promotion.items():
        if key == "reasons":
            continue
        rows.append({"gate": key, "status": value})
    rows.extend({"gate": "blocker", "status": reason} for reason in promotion.get("reasons", []))
    return rows


def _markdown_report(report: dict[str, Any]) -> str:
    promotion = report["promotion_status"]
    transfer = report["finite_transfer_gap_diagnostic"]
    return "\n".join(
        [
            "# Yang-Mills Gap Certificate Lane",
            "",
            f"- schema: `{report['schema']}`",
            f"- gauge group: `{report['gauge_group']['name']}`",
            f"- finite diagnostic receipt: `{report[FINITE_GAP_RECEIPT]}`",
            f"- finite transfer gap estimate: `{transfer['spectral_gap_estimate']}`",
            f"- finite gap floor estimate: `{report['finite_positive_gap_floor_estimate']}`",
            f"- Yang-Mills mass-gap reproduced: `{report['YANG_MILLS_GAP_REPRODUCED_RECEIPT']}`",
            f"- Clay receipt: `{report['CLAY_YANG_MILLS_GAP_RECEIPT']}`",
            "",
            "## Promotion Status",
            "",
            *[f"- {key}: `{value}`" for key, value in promotion.items() if key != "reasons"],
            "",
            "## Blockers",
            "",
            *[f"- {reason}" for reason in promotion.get("reasons", [])],
            "",
            "## Claim Boundary",
            "",
            str(report["claim_boundary"]),
            "",
        ]
    )
