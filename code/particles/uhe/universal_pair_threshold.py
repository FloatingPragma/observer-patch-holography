#!/usr/bin/env python3
"""Conditional leading-order photon/electron pair-threshold receipt.

The calculation assumes a preferred frame, rotationally isotropic leading
dispersion, ordinary additive energy-momentum conservation, a negligible LIV
correction for the soft photon, and one common negative dimension-six
coefficient for the hard photon, electron, and positron.  It is kinematics,
not an interaction, opacity, flux, or exclusion calculation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


HBAR_C_EV_M = 1.973269804e-7
PLANCK_LENGTH_M = 1.616255e-35
PLANCK_ENERGY_EV = 1.220890e28
ELECTRON_MASS_EV = 0.51099895069e6
PRODUCER_PATH = Path(__file__).resolve()
REPOSITORY_ROOT = PRODUCER_PATH.parents[3]
LEAN_SOURCE_PATH = (
    REPOSITORY_ROOT / "Lean" / "Screen" / "SeamCurrentPhotonLeptonThreshold.lean"
)


def _require_positive_finite(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def delta_abs_from_a(a_evinv: float) -> float:
    """Return ``d=a^2/20`` in the retained quartic dispersion.

    The exact carrier expansion is written as
    ``E^2=p^2+m^2-d*p^4+...``.  The threshold convention is
    ``E^2=p^2+m^2-d*E^4`` through first order in ``d`` and at leading
    ultrarelativistic order.  Replacing ``p^4`` by ``E^4`` here discards
    mixed mass/dispersion and higher-dispersion terms; it is not an exact
    identity between the two finite dispersions.
    """

    _require_positive_finite("a", a_evinv)
    return a_evinv**2 / 20.0


def universal_threshold_at_share(
    hard_energy_ev: float,
    mass_ev: float,
    delta_abs_ev_minus2: float,
    u: float,
) -> float:
    """Leading soft-photon threshold at u=x(1-x)."""

    _require_positive_finite("hard energy", hard_energy_ev)
    _require_positive_finite("mass", mass_ev)
    _require_positive_finite("|delta|", delta_abs_ev_minus2)
    if not math.isfinite(u) or not 0 < u <= 0.25:
        raise ValueError("u must be finite and lie in (0,1/4]")
    return (
        mass_ev**2 / (4.0 * hard_energy_ev * u)
        + 0.75 * delta_abs_ev_minus2 * hard_energy_ev**3 * u
    )


def transition_energy(mass_ev: float, delta_abs_ev_minus2: float) -> float:
    """Energy where the minimizing share leaves the x=1/2 boundary."""

    _require_positive_finite("mass", mass_ev)
    _require_positive_finite("|delta|", delta_abs_ev_minus2)
    return 2.0 * math.sqrt(mass_ev) * (3.0 * delta_abs_ev_minus2) ** -0.25


def minimizing_u(hard_energy_ev: float, mass_ev: float, delta_abs_ev_minus2: float) -> float:
    _require_positive_finite("hard energy", hard_energy_ev)
    _require_positive_finite("mass", mass_ev)
    _require_positive_finite("|delta|", delta_abs_ev_minus2)
    unconstrained = mass_ev / (
        math.sqrt(3.0 * delta_abs_ev_minus2) * hard_energy_ev**2
    )
    return min(0.25, unconstrained)


def universal_threshold_envelope(
    hard_energy_ev: float, mass_ev: float, delta_abs_ev_minus2: float
) -> float:
    """Minimum over the physical domain ``0 < u <= 1/4``.

    Below the transition this is the equal-share endpoint, not the
    unconstrained AM--GM lower bound.  At and above the transition the
    physical minimum attains that AM--GM bound.
    """

    return universal_threshold_at_share(
        hard_energy_ev,
        mass_ev,
        delta_abs_ev_minus2,
        minimizing_u(hard_energy_ev, mass_ev, delta_abs_ev_minus2),
    )


def unconstrained_am_gm_lower_bound(
    hard_energy_ev: float, mass_ev: float, delta_abs_ev_minus2: float
) -> float:
    """The ``2*sqrt(A*B)`` bound over all positive ``u``.

    It is a physical-domain minimum only when its equality witness satisfies
    ``u <= 1/4``, equivalently at or above ``transition_energy``.
    """

    _require_positive_finite("hard energy", hard_energy_ev)
    _require_positive_finite("mass", mass_ev)
    _require_positive_finite("|delta|", delta_abs_ev_minus2)
    return (
        math.sqrt(3.0)
        / 2.0
        * mass_ev
        * hard_energy_ev
        * math.sqrt(delta_abs_ev_minus2)
    )


def photon_only_threshold(
    hard_energy_ev: float, mass_ev: float, delta_abs_ev_minus2: float
) -> float:
    """Leading photon-only threshold with Lorentz-invariant charged leptons."""

    _require_positive_finite("hard energy", hard_energy_ev)
    _require_positive_finite("mass", mass_ev)
    _require_positive_finite("|delta|", delta_abs_ev_minus2)
    return mass_ev**2 / hard_energy_ev + 0.25 * delta_abs_ev_minus2 * hard_energy_ev**3


def photon_only_window(
    soft_energy_ev: float, mass_ev: float, delta_abs_ev_minus2: float
) -> dict[str, float | str | None]:
    """Return the finite leading photon-only window, if one exists."""

    _require_positive_finite("soft energy", soft_energy_ev)
    k_min = (4.0 * mass_ev**2 / (3.0 * delta_abs_ev_minus2)) ** 0.25
    epsilon_min = photon_only_threshold(k_min, mass_ev, delta_abs_ev_minus2)
    if soft_energy_ev < epsilon_min:
        return {
            "status": "no_leading_kinematic_window",
            "lower_root_ev": None,
            "upper_root_ev": None,
            "epsilon_min_ev": epsilon_min,
            "k_at_epsilon_min_ev": k_min,
        }
    if math.isclose(soft_energy_ev, epsilon_min, rel_tol=2e-14):
        return {
            "status": "double_root",
            "lower_root_ev": k_min,
            "upper_root_ev": k_min,
            "epsilon_min_ev": epsilon_min,
            "k_at_epsilon_min_ev": k_min,
        }

    residual = lambda k: photon_only_threshold(  # noqa: E731
        k, mass_ev, delta_abs_ev_minus2
    ) - soft_energy_ev
    lower = _bisect_root(residual, k_min * 1e-18, k_min)
    upper = _bisect_root(residual, k_min, k_min * 1e18)
    return {
        "status": "finite_leading_kinematic_window",
        "lower_root_ev": lower,
        "upper_root_ev": upper,
        "epsilon_min_ev": epsilon_min,
        "k_at_epsilon_min_ev": k_min,
    }


def global_minimum(mass_ev: float, delta_abs_ev_minus2: float) -> tuple[float, float]:
    """Global minimum (hard energy, soft energy) of the leading envelope."""

    _require_positive_finite("mass", mass_ev)
    _require_positive_finite("|delta|", delta_abs_ev_minus2)
    k_min = (16.0 * mass_ev**2 / (9.0 * delta_abs_ev_minus2)) ** 0.25
    epsilon_min = 4.0 * mass_ev**2 / (3.0 * k_min)
    return k_min, epsilon_min


def _bisect_root(function, lo: float, hi: float) -> float:
    flo = function(lo)
    fhi = function(hi)
    if flo == 0:
        return lo
    if fhi == 0:
        return hi
    if flo * fhi > 0:
        raise ValueError("root is not bracketed")
    for _ in range(160):
        mid = math.sqrt(lo * hi)
        fmid = function(mid)
        if fmid == 0 or hi / lo - 1.0 < 2e-14:
            return mid
        if flo * fmid > 0:
            lo, flo = mid, fmid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def kinematic_window(
    soft_energy_ev: float, mass_ev: float, delta_abs_ev_minus2: float
) -> dict[str, float | str | None]:
    """Return the finite hard-energy window allowed by leading kinematics."""

    _require_positive_finite("soft energy", soft_energy_ev)
    _require_positive_finite("mass", mass_ev)
    _require_positive_finite("|delta|", delta_abs_ev_minus2)
    k_min, epsilon_min = global_minimum(mass_ev, delta_abs_ev_minus2)
    k_transition = transition_energy(mass_ev, delta_abs_ev_minus2)
    epsilon_transition = universal_threshold_envelope(
        k_transition, mass_ev, delta_abs_ev_minus2
    )
    if soft_energy_ev < epsilon_min:
        return {
            "status": "no_leading_kinematic_window",
            "lower_root_ev": None,
            "upper_root_ev": None,
            "epsilon_min_ev": epsilon_min,
            "k_at_epsilon_min_ev": k_min,
        }
    if math.isclose(soft_energy_ev, epsilon_min, rel_tol=2e-14):
        return {
            "status": "double_root",
            "lower_root_ev": k_min,
            "upper_root_ev": k_min,
            "epsilon_min_ev": epsilon_min,
            "k_at_epsilon_min_ev": k_min,
        }

    residual = lambda k: universal_threshold_envelope(  # noqa: E731
        k, mass_ev, delta_abs_ev_minus2
    ) - soft_energy_ev
    lower = _bisect_root(residual, k_min * 1e-18, k_min)
    if soft_energy_ev <= epsilon_transition:
        upper = _bisect_root(residual, k_min, k_transition)
    else:
        # On the asymmetric branch epsilon=(sqrt(3)/2)m k sqrt(d).
        upper = 2.0 * soft_energy_ev / (
            math.sqrt(3.0) * mass_ev * math.sqrt(delta_abs_ev_minus2)
        )
    u_upper = minimizing_u(upper, mass_ev, delta_abs_ev_minus2)
    x_small = (1.0 - math.sqrt(max(0.0, 1.0 - 4.0 * u_upper))) / 2.0
    return {
        "status": "finite_leading_kinematic_window",
        "lower_root_ev": lower,
        "upper_root_ev": upper,
        "epsilon_min_ev": epsilon_min,
        "k_at_epsilon_min_ev": k_min,
        "transition_energy_ev": k_transition,
        "epsilon_at_transition_ev": epsilon_transition,
        "u_at_upper_root": u_upper,
        "smaller_outgoing_share_at_upper_root": x_small,
        "smaller_outgoing_energy_ev": x_small * upper,
    }


def a_from_quartic_scale(quartic_scale_ev: float) -> tuple[float, float]:
    """Translate ``d=1/scale^2=a^2/20`` into the carrier coordinate."""

    _require_positive_finite("quartic scale", quartic_scale_ev)
    a_m = math.sqrt(20.0) / quartic_scale_ev * HBAR_C_EV_M
    return a_m, a_m / PLANCK_LENGTH_M


def build_receipt() -> dict[str, object]:
    target_photons = {
        "representative_radio_4e-8_eV": 4.0e-8,
        "CMB_mean_6.34e-4_eV": 6.34e-4,
        "selected_CMB_Wien_tail_3e-3_eV": 3.0e-3,
    }
    windows: dict[str, object] = {}
    for a_over_lp in (1.0, 10.0, 100.0, 1000.0, 10000.0):
        a_evinv = a_over_lp * PLANCK_LENGTH_M / HBAR_C_EV_M
        delta_abs = delta_abs_from_a(a_evinv)
        entry: dict[str, object] = {
            "a_over_planck_length": a_over_lp,
            "a_evinv": a_evinv,
            "delta_abs_ev_minus2": delta_abs,
            "transition_energy_ev": transition_energy(ELECTRON_MASS_EV, delta_abs),
        }
        for label, epsilon in target_photons.items():
            window = kinematic_window(epsilon, ELECTRON_MASS_EV, delta_abs)
            if window["upper_root_ev"] is not None:
                window["a_times_upper_root_dimensionless"] = (
                    a_evinv * float(window["upper_root_ev"])
                )
            entry[label] = window
        windows[f"a={a_over_lp:g}_lP"] = entry

    planck_a_evinv = PLANCK_LENGTH_M / HBAR_C_EV_M
    planck_delta_abs = delta_abs_from_a(planck_a_evinv)
    photon_only_comparison = {
        label: photon_only_window(epsilon, ELECTRON_MASS_EV, planck_delta_abs)
        for label, epsilon in target_photons.items()
    }

    # Use the sign-specific subluminal value in the LHAASO paper's conclusion,
    # rather than rounding its 6.9e11 GeV result upward via the one-significant-
    # figure abstract value 6e-8 E_Pl.
    lhaaso_eqg2_ev = 6.9e11 * 1.0e9
    fermi_eqg2_ev = 1.3e11 * 1.0e9
    crab_mlv_ev = 2.0e16 * 1.0e9
    lhaaso_a_m, lhaaso_a_lp = a_from_quartic_scale(lhaaso_eqg2_ev)
    fermi_a_m, fermi_a_lp = a_from_quartic_scale(fermi_eqg2_ev)
    crab_a_m, crab_a_lp = a_from_quartic_scale(crab_mlv_ev)
    auger_a_m = math.sqrt(20.0e-58) * HBAR_C_EV_M

    return {
        "schema": "oph.universal_pair_threshold_leading.v1",
        "epistemic_status": {
            "retrospective": True,
            "conditional_kinematic_theorem": True,
            "physical_scale_selected_by_oph": False,
            "photon_and_lepton_identification_derived": False,
            "preferred_frame_and_additive_conservation_derived": False,
            "interaction_vertex_or_cross_section_derived": False,
            "opacity_or_flux_prediction": False,
            "auger_exclusion_recomputed_for_universal_sector": False,
            "frozen_prediction_score": False,
        },
        "assumptions": [
            "preferred frame and rotationally isotropic leading dispersion",
            "ordinary additive energy-momentum conservation",
            "soft-photon LIV neglected",
            "common delta=-a^2/20 for hard photon, electron, and positron",
            "ultrarelativistic leading-order expansion",
        ],
        "analytic_result": {
            "threshold_at_u": "m^2/(4 k u) + (3 d k^3/4) u",
            "u_domain": "0 < u=x(1-x) <= 1/4",
            "interior_minimizer": "u*=m/(sqrt(3d) k^2)",
            "unconstrained_am_gm_lower_bound": "epsilon=(sqrt(3)/2) m k sqrt(d)",
            "physical_constrained_minimum": (
                "equal-share value below 3 d k^4=16 m^2; AM-GM value at and "
                "above that transition"
            ),
            "group_velocity_low_ak": "v_g=1-(3/2)d p^2+O(d^2 p^4)=1-3a^2p^2/40+...",
            "quartic_notation_bridge": (
                "p^4 in the carrier dispersion may be replaced by E^4 only "
                "at the retained leading ultrarelativistic/first-d order"
            ),
        },
        "source_binding": {
            "producer_path": str(PRODUCER_PATH.relative_to(REPOSITORY_ROOT)),
            "producer_sha256": _sha256(PRODUCER_PATH),
            "lean_source_path": str(LEAN_SOURCE_PATH.relative_to(REPOSITORY_ROOT)),
            "lean_source_sha256": _sha256(LEAN_SOURCE_PATH),
        },
        "conditional_external_bounds": {
            "LHAASO_GRB221009A_quadratic_photon_TOF": {
                "input_sub_luminal_EQG2_lower_bound_ev": lhaaso_eqg2_ev,
                "a_max_m": lhaaso_a_m,
                "a_max_over_lP": lhaaso_a_lp,
                "source": "https://doi.org/10.1103/PhysRevLett.133.071501",
                "scope": (
                    "95% subluminal 6.9e11 GeV result; leading isotropic "
                    "photon/clock attachment; source timing and the paper's "
                    "afterglow/EBL model remain analysis assumptions"
                ),
            },
            "Fermi_four_GRB_quadratic_photon_TOF": {
                "input_sub_luminal_EQG2_lower_bound_ev": fermi_eqg2_ev,
                "a_max_m": fermi_a_m,
                "a_max_over_lP": fermi_a_lp,
                "source": "https://arxiv.org/abs/1305.3463",
                "scope": "subluminal 95% result under the no-source-intrinsic-dispersion scenario",
            },
            "Crab_quadratic_electron_synchrotron": {
                "input_MLV_lower_bound_ev": crab_mlv_ev,
                "a_max_m": crab_a_m,
                "a_max_over_lP": crab_a_lp,
                "source": "https://arxiv.org/abs/1207.0670",
                "scope": (
                    "95% M_LV>2e16 GeV result for the paper's quartic "
                    "electron dispersion, conditional on its CPT/parity, "
                    "Hamiltonian, and Crab synchrotron model assumptions"
                ),
            },
            "Auger_photon_only_translation_not_universal": {
                "input_delta_gamma_2_lower_bound_ev_minus2": -1.0e-58,
                "a_max_m": auger_a_m,
                "a_max_over_lP": auger_a_m / PLANCK_LENGTH_M,
                "source": "https://arxiv.org/abs/2112.06773",
                "scope": (
                    "numeric translation of delta_gamma,2>-1e-58 eV^-2 in "
                    "the paper's alternative UHECR scenario with a "
                    "subdominant proton component and photon-only LIV; not "
                    "transferable to a common photon/electron coefficient "
                    "without recomputing propagation"
                ),
            },
        },
        "representative_target_windows": windows,
        "photon_only_lorentz_invariant_lepton_comparison_at_a=1_lP": (
            photon_only_comparison
        ),
        "nonclaims": [
            "A representative photon energy is not a radiation spectrum.",
            "A kinematic window is not an optical depth or flux.",
            "No Auger bound on the universal branch follows without cross sections, backgrounds, cascades, sources, composition, and detector response.",
            "The low-ak coefficient map is not the complete reciprocal-cosine dispersion at arbitrary momentum.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rendered = json.dumps(build_receipt(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
