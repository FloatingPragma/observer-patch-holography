"""Complete local pairing kernel on the frozen compact defect source.

The scalar site coordinates are independent real interventions. A scalar
Hamiltonian, kinetic metric and selected collective mode are not assumed.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import linalg as la

from model.compact_defect_pair.spectral import source
from model.compact_defect_pair.topology import spectral_bounds

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SETTINGS = {"radius": 3, "mass": 64, "half_separation": 1}
TOLERANCE = 2e-10


@lru_cache(maxsize=1)
def build():
    """Diagonalize every spatial row, retaining the full four-mode band."""
    raw = source(**SETTINGS)
    h, c = raw["h"].toarray(), raw["C"].toarray()
    energies, vectors = la.eigh(h)
    mask = abs(energies) < 2
    if np.count_nonzero(mask) != 4:
        raise ArithmeticError("the frozen source must have four band modes")
    w, q = vectors[:, mask], vectors[:, ~mask]
    pair = w.conj().T @ c @ w.conj()
    sites = w.reshape(len(raw["points"]), 2, 4)
    density = np.einsum("xsi,xsj->xij", sites.conj(), sites)
    kernel = density @ pair
    checks = {
        "band_invariance": la.norm(h @ w - w * energies[mask]),
        "particle_hole_band": la.norm(c @ w.conj() - w @ pair),
        "band_metric": la.norm(w.conj().T @ w - np.eye(4)),
        "complete_site_density": la.norm(density.sum(axis=0) - np.eye(4)),
        "uniform_pair_unitarity": la.norm(pair.conj().T @ pair - np.eye(4)),
        "every_site_pair_antisymmetry": float(np.max(abs(kernel + kernel.transpose(0, 2, 1)))),
    }
    if max(checks.values()) > TOLERANCE:
        raise ArithmeticError("complete spatial kernel identities failed")
    return dict(raw=raw, h=h, c=c, energies=energies, vectors=vectors,
                mask=mask, W=w, Q=q, pair=pair, density=density,
                kernel=kernel, checks=checks)


def profile_response(profile, *, value=None):
    value = build() if value is None else value
    f = np.asarray(profile)
    if f.shape != (len(value["raw"]["points"]),) or not np.isrealobj(f) or not np.isfinite(f).all():
        raise ValueError("one finite real scalar coordinate per source site required")
    f = f.astype(float)
    density, pair = value["density"], value["pair"]
    compressed = np.einsum("x,xij->ij", f, density)
    vertex = compressed @ pair
    second = np.einsum("x,xij->ij", f*f, density)
    leakage_gram = pair.conj().T @ (second - compressed @ compressed) @ pair
    leakage_gram = (leakage_gram + leakage_gram.conj().T) / 2
    direct = np.repeat(f, 2)[:, None] * (value["c"] @ value["W"].conj())
    omitted = direct - value["W"] @ vertex
    if la.norm(omitted.conj().T @ omitted - leakage_gram) > TOLERANCE:
        raise ArithmeticError("complete complement leakage identity failed")
    singulars = la.svdvals(vertex)
    leakage = float(la.norm(omitted, 2))
    return dict(vertex=vertex, density=compressed, leakage_gram=leakage_gram,
                singulars=singulars, leakage=leakage,
                euclidean_profile_norm=float(la.norm(f)),
                maximum_site_amplitude=float(np.max(abs(f))))


def diagnostic_profiles(value=None):
    value = build() if value is None else value
    points = np.asarray(value["raw"]["points"])
    distance = np.max(abs(points - np.array([0, 0, 1])), axis=1)
    return {"uniform": np.ones(len(points)),
            "positive_core": (distance == 0).astype(float),
            "positive_core_patch": (distance <= 1).astype(float)}


def core_residue(value=None):
    """Canonical norms from the full exterior resolvent at every band energy."""
    value = build() if value is None else value
    core = np.asarray(value["raw"]["core_indices"])
    exterior = np.setdiff1d(np.arange(value["h"].shape[0]), core)
    b = value["h"][np.ix_(exterior, core)]
    d = value["h"][np.ix_(exterior, exterior)]
    de, du = la.eigh(d)
    rotated = du.conj().T @ b
    rows = []
    for j, energy in enumerate(value["energies"][value["mask"]]):
        inverse_b = du @ (rotated / (de - energy)[:, None])
        metric = np.eye(4) + inverse_b.conj().T @ inverse_b
        u = value["W"][core, j]
        wave = np.zeros(value["h"].shape[0], complex)
        wave[core], wave[exterior] = u, -inverse_b @ u
        error = la.norm(wave - value["W"][:, j])
        normalization = float(np.real(u.conj() @ metric @ u))
        if error > TOLERANCE or abs(normalization - 1) > TOLERANCE:
            raise ArithmeticError("full exterior reconstruction or residue normalization failed")
        rows.append({"energy": float(energy), "raw_core_norm_squared": float(np.vdot(u, u).real),
                     "residue_metric_eigenvalues": la.eigvalsh(metric).tolist(),
                     "residue_normalized_norm_squared": normalization,
                     "full_wave_reconstruction_error": float(error)})
    return rows


def texture_controls():
    """The same fixed patch and geometry at three declared texture amplitudes."""
    rows = []
    for mass in (32, 64, 128):
        raw = source(radius=3, mass=mass, half_separation=1)
        energies, vectors = la.eigh(raw["h"].toarray())
        w = vectors[:, abs(energies) < 2]
        if w.shape[1] != 4:
            raise ArithmeticError("complete comparison band must have four modes")
        points = np.asarray(raw["points"])
        f = (np.max(abs(points-np.array([0, 0, 1])), axis=1) <= 1).astype(float)
        vertex = w.conj().T @ (np.repeat(f, 2)[:, None] * (raw["C"] @ w.conj()))
        singulars = la.svdvals(vertex)
        # Literal first hopping out of the negative core into the patch.
        negative = np.flatnonzero(np.all(points == (0, 0, -1), axis=1))[0]
        columns = raw["h"][:, 2*negative:2*negative+2].toarray()
        first_hop = np.repeat(f, 2)[:, None] * columns
        first_gram = first_hop.conj().T @ first_hop
        if la.norm(first_gram - np.eye(2)/4) > 1e-12:
            raise ArithmeticError("the complete negative-core first-hop Gram differs")
        rows.append({"texture_amplitude": mass, "complete_one_body_rows": len(energies),
                     "singular_values": singulars.tolist(),
                     "weak_strength_times_amplitude_squared": float(singulars[-1]*mass*mass),
                     "first_hop_squared_overlap_coefficient_exact": "1/4",
                     "leading_weak_strength": 1/(4*mass*mass)})
    return rows


def full_intervention(profile, amplitude=0.125, *, value=None):
    """Compare all eight low Nambu rows with their exact band compression.

    The source is B(t)=diag(h,-h*)+t[[0,fC],[(fC)*,0]]. The internal
    unit invariant Higgs tensor factors out; no family matrix enters.
    """
    value = build() if value is None else value
    row = profile_response(profile, value=value)
    f = np.asarray(profile, float)
    t = float(amplitude)
    if not np.isfinite(t):
        raise ValueError("finite scalar intervention required")
    h, c = value["h"], value["c"]
    delta = np.repeat(f, 2)[:, None] * c
    full = np.block([[h, t * delta], [t * delta.conj().T, -h.conj()]])
    a = np.diag(value["energies"][value["mask"]])
    small = np.block([[a, t*row["vertex"]], [t*row["vertex"].conj().T, -a]])
    actual = la.eigh(full, eigvals_only=True, subset_by_value=[-2., 2.])
    compressed = la.eigvalsh(small)
    if len(actual) != 8:
        raise ArithmeticError("the complete low Nambu band changed rank")
    # Exact parent estimates: r=(7/4)/(64-4)=7/240 and complement >=62.
    radius = 7/240
    coupling = abs(t) * row["maximum_site_amplitude"]
    denominator = 62 - radius - 2*coupling
    if denominator <= 0:
        raise ValueError("intervention exceeds the declared separated-band domain")
    schur = t*t*row["leakage"]**2/denominator
    residue = t*t*row["leakage"]**2/denominator**2
    error = float(np.max(abs(actual-compressed)))
    if error > schur + TOLERANCE:
        raise ArithmeticError("full spectrum exceeds its Schur comparison allowance")
    return {"amplitude": t, "full_Nambu_rows": len(full), "retained_low_rows": len(actual),
            "complete_low_energies": actual.tolist(), "compressed_low_energies": compressed.tolist(),
            "maximum_spectral_difference": error,
            "Schur_self_energy_norm_bound_diagnostic": schur,
            "residue_metric_distance_from_identity_bound_diagnostic": residue,
            "analytic_complement_denominator": denominator,
            "bound_uses_numerical_leakage_not_outward_interval": True}


def source_bindings():
    paths = ["model/compact_defect_pair/spectral.py", "model/compact_defect_pair/topology.py",
             "model/compact_defect_pair/PROOF.md", "model/defect_mode_overlap/overlap.py",
             "model/defect_mode_overlap/PROOF.md"]
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}


def calculate():
    value = build()
    profiles = diagnostic_profiles(value)
    diagnostics = {}
    for name, f in profiles.items():
        row = profile_response(f, value=value)
        normalized = profile_response(f/la.norm(f), value=value)
        diagnostics[name] = {"site_count": int(np.count_nonzero(f)),
            "singular_values": row["singulars"].tolist(),
            "complement_leakage_norm": row["leakage"],
            "euclidean_profile_norm": row["euclidean_profile_norm"],
            "unit_euclidean_profile_singular_values": normalized["singulars"].tolist(),
            "profile_normalization_is_physical_scalar_residue": False}
    patch = profiles["positive_core_patch"]
    patch_row = profile_response(patch, value=value)
    diagnostics["positive_core_patch"]["largest_to_smallest_ratio"] = float(
        patch_row["singulars"][0]/patch_row["singulars"][-1])
    site_rows = []
    for point, density, vertex in zip(value["raw"]["points"], value["density"], value["kernel"], strict=True):
        gram = value["pair"].conj().T @ (density-density@density) @ value["pair"]
        site_rows.append({"site": list(point), "singular_values": la.svdvals(vertex).tolist(),
            "complement_leakage_norm": float(np.sqrt(max(0, la.eigvalsh(gram)[-1])))})
    return {"schema": "oph.defect_mode_overlap.v1", "source": SETTINGS,
            "one_body_rows": len(value["h"]), "sites": len(site_rows), "complete_band_rank": 4,
            "band_energies": value["energies"][value["mask"]].tolist(),
            "parent_exact_bounds": spectral_bounds(64, 1, 3),
            "identities_maximum_numerical_error": max(value["checks"].values()),
            "source_bindings": source_bindings(), "all_site_diagnostics": site_rows,
            "profiles": diagnostics, "texture_amplitude_controls": texture_controls(),
            "complete_exterior_residue": core_residue(value),
            "full_interventions": {name: full_intervention(profiles[name]/la.norm(profiles[name]), value=value)
                for name in ("uniform", "positive_core_patch")},
            "scope": {"all_spatial_rows_and_scalar_site_insertions_retained": True,
                "full_local_leakage_kernel_reconstructible": True,
                "empirical_inputs": [], "inserted_family_coefficients": False,
                "selected_dynamical_scalar_mode": False, "scalar_kinetic_residue_supplied": False,
                "three_generations_derived": False, "physical_Yukawa_prediction": False,
                "local_Gauss_field_transfer_proved": False}}


def write():
    result = calculate()
    (HERE/"results.json").write_text(json.dumps(result, indent=2)+"\n")
    value = build()
    np.savez_compressed(HERE/"site_kernel.npz", points=np.asarray(value["raw"]["points"]),
        band_energies=value["energies"][value["mask"]], W=value["W"],
        pair=value["pair"], density=value["density"], kernel=value["kernel"])
    return result


def validate_archive(packet, *, value=None):
    value = build() if value is None else value
    shapes = {"points": (343, 3), "band_energies": (4,), "W": (686, 4),
              "pair": (4, 4), "density": (343, 4, 4), "kernel": (343, 4, 4)}
    if set(packet) != set(shapes):
        raise ArithmeticError("complete scalar archive key inventory differs")
    for key, shape in shapes.items():
        array = packet[key]
        if not isinstance(array, np.ndarray) or array.shape != shape or not np.isfinite(array).all():
            raise ArithmeticError("invalid complete scalar archive array: "+key)
    if not np.array_equal(packet["points"], np.asarray(value["raw"]["points"])):
        raise ArithmeticError("complete source site alphabet differs")
    if not np.allclose(packet["band_energies"], value["energies"][value["mask"]], atol=TOLERANCE, rtol=0):
        raise ArithmeticError("complete band energies differ")
    # Stored coordinates may choose a different unitary frame for the same band.
    change = value["W"].conj().T @ packet["W"]
    if (la.norm(packet["W"] - value["W"] @ change) > TOLERANCE or
            la.norm(change.conj().T @ change - np.eye(4)) > TOLERANCE):
        raise ArithmeticError("stored band is not the complete source band")
    expected = {"pair": change.conj().T @ value["pair"] @ change.conj(),
                "density": change.conj().T @ value["density"] @ change,
                "kernel": change.conj().T @ value["kernel"] @ change.conj()}
    for key, actual in expected.items():
        if not np.allclose(packet[key], actual, atol=TOLERANCE, rtol=0):
            raise ArithmeticError("complete scalar site kernel differs")


def check():
    saved = json.loads((HERE/"results.json").read_text())
    fresh = calculate()
    if saved["source_bindings"] != fresh["source_bindings"]:
        raise ArithmeticError("source bindings differ")
    # Numerical spectra are portable diagnostics, compared at their stated tolerance.
    def compare(a, b):
        if isinstance(a, dict):
            return set(a) == set(b) and all(compare(a[k], b[k]) for k in a)
        if isinstance(a, list):
            return len(a) == len(b) and all(compare(x, y) for x, y in zip(a, b, strict=True))
        if isinstance(a, float):
            return abs(a-b) <= TOLERANCE*max(1, abs(a), abs(b))
        return a == b
    if not compare(saved, fresh):
        raise ArithmeticError("replayed overlap receipt differs")
    value = build()
    with np.load(HERE/"site_kernel.npz") as packet:
        validate_archive(packet, value=value)
    return fresh
