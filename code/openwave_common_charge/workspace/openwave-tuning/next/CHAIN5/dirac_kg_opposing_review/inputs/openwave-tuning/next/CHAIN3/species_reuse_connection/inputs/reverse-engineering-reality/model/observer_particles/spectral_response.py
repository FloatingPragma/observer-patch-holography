"""Charge-resolved spectral weights of a supplied finite Hermitian generator.

For canonical finite field coordinates and fixed interpolator columns B,
G(z)=B^dagger (z I-H)^(-1) B. Exact joint energy/charge projectors give
positive residue matrices and sum to B^dagger B. This module diagnoses that
finite algebra; it does not infer H, a vacuum, CAR, a physical clock, a
regulator, or a particle species from detector records.

Nearby numerical eigenvalues are retained as unresolved clusters with their
full frequency lists. A cluster's total weight is not an assertion that all
its frequencies coincide. Dark modes are retained, and no weight is
renormalized after restricting the interpolators.
"""
from __future__ import annotations

from numbers import Real

import numpy as np


def _positive(value, name):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
            or not np.isfinite(value) or value <= 0):
        raise ValueError(f"{name} must be finite and positive")
    return float(value)


def _hermitian(value, name, atol):
    matrix = np.asarray(value, dtype=complex)
    if (matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]
            or matrix.shape[0] == 0 or not np.isfinite(matrix).all()):
        raise ValueError(f"{name} must be a nonempty finite square matrix")
    error = float(np.linalg.norm(matrix-matrix.conj().T, "fro"))
    if error > atol*max(1., float(np.linalg.norm(matrix, "fro"))):
        raise ValueError(f"{name} must be Hermitian; damping is outside this spectral contract")
    return (matrix+matrix.conj().T)/2, error


def _clusters(values, resolution):
    """Bound each cluster's full diameter; do not chain neighboring gaps."""
    start = 0
    for stop in range(1, len(values)+1):
        if stop == len(values) or values[stop]-values[start] > resolution:
            yield slice(start, stop)
            start = stop


def charge_resolved_spectrum(hamiltonian, charge, interpolators, *,
                             atol=1e-11, energy_resolution=1e-10,
                             charge_resolution=1e-10):
    """Retain all joint spectral subspaces and their interpolator weights.

    Interpolators is an n-by-k matrix B, with no normalization required.
    Charge must commute with H within the declared numerical admission.
    Input Hermitian-part and charge-block projection errors are returned
    explicitly. The matrices and residuals use binary64 arithmetic, not
    outward-rounded bounds. Resolutions are absolute in the supplied
    frequency and charge coordinates.
    """
    atol = _positive(atol, "atol")
    energy_resolution = _positive(energy_resolution, "energy_resolution")
    charge_resolution = _positive(charge_resolution, "charge_resolution")
    h, h_error = _hermitian(hamiltonian, "Hamiltonian", atol)
    q, q_error = _hermitian(charge, "charge", atol)
    b = np.asarray(interpolators, dtype=complex)
    if (q.shape != h.shape or b.ndim != 2 or b.shape[0] != len(h)
            or b.shape[1] == 0 or not np.isfinite(b).all()):
        raise ValueError("charge and finite interpolator columns must share the full generator domain")
    commutator = float(np.linalg.norm(h@q-q@h, "fro"))
    if commutator > atol*max(1., float(np.linalg.norm(h, "fro")*np.linalg.norm(q, "fro"))):
        raise ValueError("charge is not conserved by the supplied generator")
    charges, charge_basis = np.linalg.eigh(q)
    rows = []
    reconstructed = np.zeros_like(h)
    total = np.zeros((b.shape[1], b.shape[1]), complex)
    gram = b.conj().T@b
    visible_threshold = atol*max(1., float(np.linalg.norm(gram, "fro")))
    for charge_slice in _clusters(charges, charge_resolution):
        basis = charge_basis[:, charge_slice]
        block = basis.conj().T@h@basis
        reconstructed += basis@block@basis.conj().T
        frequencies, energy_basis = np.linalg.eigh((block+block.conj().T)/2)
        amplitudes = (basis@energy_basis).conj().T@b
        for energy_slice in _clusters(frequencies, energy_resolution):
            selected = amplitudes[energy_slice]
            residue = selected.conj().T@selected
            total += residue
            eigenweights = np.linalg.eigvalsh(residue)
            rows.append({
                "charge_interval": [float(charges[charge_slice][0]), float(charges[charge_slice][-1])],
                "frequency_interval": [float(frequencies[energy_slice][0]), float(frequencies[energy_slice][-1])],
                "frequencies": frequencies[energy_slice].copy(),
                "cluster_dimension": energy_slice.stop-energy_slice.start,
                "exact_energy_degeneracy_asserted": False,
                "residue": residue,
                "residue_trace": float(np.trace(residue).real),
                "numerical_residue_rank": int(np.count_nonzero(eigenweights > visible_threshold)),
                "dark_at_numerical_admission": bool(np.max(eigenweights) <= visible_threshold),
            })
    projection_error = float(np.linalg.norm(h-reconstructed, "fro"))
    if projection_error > atol*max(1., float(np.linalg.norm(h, "fro"))):
        raise ValueError("charge-sector projection would discard a resolved part of the generator")
    return {
        "dimension": len(h), "interpolator_count": b.shape[1], "clusters": rows,
        "interpolator_gram": gram, "residue_sum": total,
        "residue_sum_error_Frobenius": float(np.linalg.norm(total-gram, "fro")),
        "Hamiltonian_Hermitian_part_error_Frobenius": h_error/2,
        "charge_Hermitian_part_error_Frobenius": q_error/2,
        "charge_commutator_error_Frobenius": commutator,
        "charge_block_projection_error_Frobenius": projection_error,
        "atol": atol, "energy_resolution": energy_resolution,
        "charge_resolution": charge_resolution,
    }


def spectrum_summary(result):
    """Serialize the computed weights without dropping dark/negative modes."""
    def encode(value):
        if isinstance(value, np.ndarray):
            if np.iscomplexobj(value):
                return {"real": value.real.tolist(), "imag": value.imag.tolist()}
            return value.tolist()
        if isinstance(value, dict):
            return {key: encode(item) for key, item in value.items()}
        if isinstance(value, list):
            return [encode(item) for item in value]
        return value
    output = encode(result)
    output["scope"] = {
        "input": "supplied finite Hermitian generator, conserved charge and fixed interpolators in canonical field coordinates",
        "all_generator_modes_retained": True,
        "residues_normalized_after_projection": False,
        "nearby_frequencies_identified_as_one_exact_pole": False,
        "exact_residue_ranks_certified": False,
        "native_detector_record_reconstruction_performed": False,
        "physical_particle_species_or_mass_identified": False,
        "outward_numerical_certificate": False,
    }
    return output
