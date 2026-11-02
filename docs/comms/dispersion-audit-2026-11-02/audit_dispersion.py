"""Independent, data-free numerical and symbolic audit of the IPI note.

Run from any directory with the research environment's numpy, sympy, mpmath.
Does not import OPH's geometry or certificate implementation, read experimental
comparison data, or modify frozen artifacts. Numerical sampling is a control,
not a substitute for the exact identities described in the note.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
import subprocess

import mpmath as mp
import numpy as np
import sympy as s

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


def geometry():
    phi = (1 + np.sqrt(5.0)) / 2
    vertices = []
    for x, y in itertools.product((-1, 1), repeat=2):
        vertices.extend(((0, x, y * phi), (x, y * phi, 0), (y * phi, 0, x)))
    vertices = np.array(vertices, dtype=float)
    vertices /= np.linalg.norm(vertices, axis=1)[:, None]
    distances = np.linalg.norm(vertices[:, None] - vertices[None, :], axis=2)
    edge_length = distances[distances > 1e-8].min()
    edges = [(i, j) for i in range(12) for j in range(i + 1, 12)
             if abs(distances[i, j] - edge_length) < 1e-12]
    edge_set = set(edges)
    faces = [(i, j, k) for i, j, k in itertools.combinations(range(12), 3)
             if all(pair in edge_set for pair in ((i, j), (i, k), (j, k)))]
    ev = np.array([vertices[j] - vertices[i] for i, j in edges])
    fv = np.array([vertices[i] + vertices[j] + vertices[k] for i, j, k in faces])
    ev /= np.linalg.norm(ev, axis=1)[:, None]
    fv /= np.linalg.norm(fv, axis=1)[:, None]
    assert (len(vertices), len(ev), len(fv)) == (12, 30, 20)
    return vertices, ev, fv


def main():
    rng = np.random.default_rng(20261102)
    v, e, f = geometry()
    adjacency = s.Matrix(12, 12, lambda i,j: int(i != j and abs(float(v[i]@v[j])-1/np.sqrt(5)) < 1e-12))
    ell = s.symbols("ell")
    laplacian = 5*s.eye(12)-adjacency
    characteristic = laplacian.charpoly(ell).as_expr()
    assert s.expand(characteristic-ell*(ell-6)**5*(ell**2-10*ell+20)**3) == 0
    n = rng.normal(size=(512, 3))
    n /= np.linalg.norm(n, axis=1)[:, None]

    def i6(x):
        t = np.asarray(x) @ v.T
        return 25 / 132 * ((231*t**6 - 315*t**4 + 105*t**2 - 5)/16).sum(axis=-1)

    iv = i6(n)
    moment_residual = 0.0
    for directions, beta in ((v, 64/175), (e, -2/7), (f, -64/189)):
        dot = n @ directions.T
        count = len(directions)
        for order, expected in ((2, count/3), (4, count/5),
                                (6, count/7 + beta*iv)):
            moment_residual = max(moment_residual, float(np.max(abs((dot**order).sum(axis=1) - expected))))
    assert moment_residual < 2e-12
    assert np.max(abs(i6(v) - 1)) < 1e-12
    assert np.max(abs(i6(e) + 5/16)) < 1e-12
    assert np.max(abs(i6(f) + 5/9)) < 1e-12

    # Exact independent cosine/square-root/log algebra; a=1 for this control.
    k, A, B, C, gamma = s.symbols("k A B C gamma", real=True)
    omega = s.series(s.sqrt(1+A*k**2+(B+C)*k**4), k, 0, 6).removeO()*k
    radial_speed = s.diff(omega, k)
    assert s.expand(radial_speed).coeff(k, 2) == 3*A/2
    assert s.expand(radial_speed).coeff(k, 4) == 5*(B+C)/2-5*A**2/8
    native_log = s.series(-s.log(1-gamma*(k**2+A*k**4))/gamma, k, 0, 6)
    assert s.expand(native_log.removeO()).coeff(k, 4) == A + gamma/2

    # Positive multi-radius controls. These do not assume the claimed band.
    ratios = []
    for _ in range(512):
        radii = np.exp(rng.uniform(-2, 2, 5))
        weights = np.exp(rng.uniform(-2, 2, 5))
        seeds = rng.normal(size=(5, 3))
        seeds /= np.linalg.norm(seeds, axis=1)[:, None]
        mu2, mu4, mu6 = [np.sum(weights*radii**p) for p in (2, 4, 6)]
        isotropic = 10/21*mu2*mu6/mu4**2
        anisotropy = 16/75*np.sum(weights*radii**6*i6(seeds))/mu6
        assert isotropic >= 10/21-1e-13
        assert -16/135-1e-13 <= anisotropy <= 16/75+1e-13
        ratios.append((float(isotropic), float(anisotropy)))

    # A valid positive vertex/face mixture has no rank-six anisotropy.
    zero_mix = s.Rational(16, 75)*(25*12 + 27*20*s.Rational(-5, 9))/(25*12+27*20)
    assert zero_mix == 0
    # The eight-order lock fails for different radii, within the valid class.
    def mean(power):
        return (12+20*2**power*s.Rational(-5, 9))/(12+20*2**power)
    lock_residual = s.Rational(64, 125)*mean(8) - s.Rational(12, 5)*s.Rational(16, 75)*mean(6)
    assert lock_residual == -s.Rational(57344, 10360225)
    # Stability does not imply positive shell weights. This A5-invariant signed
    # stencil is nonnegative term by term: (1-cos x)(4/5-cos x/5)>=0.
    x = s.symbols("x", real=True)
    signed = (1-s.cos(x))-s.Rational(1, 10)*(1-s.cos(2*x))
    assert s.trigsimp(signed-(1-s.cos(x))*(s.Rational(4, 5)-s.cos(x)/5)) == 0
    signed_C4 = -s.Rational(1,20)*(1-s.Rational(16,10))/(1-s.Rational(4,10))
    assert signed_C4 == s.Rational(1,20)

    # Long-wavelength full edge symbol vs sixth-order polynomial, evaluated
    # at high precision to avoid cancellation. Geometry is reconstructed above.
    mp.mp.dps = 70
    phi = (1+mp.sqrt(5))/2
    raw = []
    for b, c in itertools.product((-1, 1), repeat=2):
        raw.extend(((0,b,c*phi),(b,c*phi,0),(c*phi,0,b)))
    vm = [mp.matrix(p)/mp.sqrt(1+phi**2) for p in raw]
    em = []
    for i, j in itertools.combinations(range(12), 2):
        d = vm[j]-vm[i]
        if abs(mp.norm(d)-2/mp.sqrt(1+phi**2)) < mp.mpf("1e-60"):
            em.append(d/mp.norm(d))
    max_fraction = mp.mpf(0)
    for direction in n[:32]:
        nm = mp.matrix([mp.mpf(float(t)) for t in direction]); nm /= mp.norm(nm)
        dots = [mp.fdot(p,nm) for p in vm]
        im = mp.mpf(25)/132*sum((231*t**6-315*t**4+105*t**2-5)/16 for t in dots)
        for q in map(mp.mpf, ("0.01", "0.1", "0.5", "1")):
            exact = sum(1-mp.cos(q*mp.fdot(p,nm)) for p in em)/5
            p6 = q**2-q**4/20+(mp.mpf(1)/840-im/12600)*q**6
            bound = mp.mpf(7)/388800*q**8
            assert abs(exact-p6) <= bound
            max_fraction = max(max_fraction, abs(exact-p6)/bound)

    registry = json.loads((REPO/"claims/frozen_prediction_register.json").read_text())
    custody_checks = {}
    for row in ("FZ-11", "FZ-12"):
        entry = registry["external_custody_contracts"][row]
        folder = REPO/"evidence/custody"/entry["custody_path"]
        for name, expected in entry["artifact_sha256"].items():
            actual = hashlib.sha256((folder/name).read_bytes()).hexdigest()
            assert actual == expected, f"immutable hash drift: {row}/{name}"
            custody_checks[f"{row}/{name}"] = actual

    input_paths = [
        "docs/AXIOM_REFERENCE.md",
        "claims/frozen_prediction_register.json",
        "flagship/from_observer_consensus_to_standard_physics.tex",
        "code/a5_fingerprint/carrier_class_dispersion_certificate.py",
        "code/a5_fingerprint/runtime/carrier_class_dispersion_receipt.json",
        "code/a5_fingerprint/fz12_full_symbol_remainder_certificate.py",
        "code/a5_fingerprint/fz12_observation_map_certificate.py",
        "code/a5_fingerprint/fz12_free_photon_hamiltonian.py",
        "model/PHYSICAL_AUDIT.md", "model/native_effective_fields/construction.py",
        "model/native_effective_fields/README.md",
        "model/physical_space_time/construction.py",
        "model/physical_space_time/README.md",
        "model/physical_space_time/OCCUPATION_FIELD.md",
        "model/common_source/chiral_candidate/dispersion_identifiability.py",
        "model/energy_identifiability/README.md",
        "model/full_fock_field/README.md",
        "Lean/Screen/DispersionArmingInterface.lean",
    ]
    for directory in ("code/a5_fingerprint", "model/native_effective_fields",
                      "model/physical_space_time", "model/energy_identifiability",
                      "model/common_source/chiral_candidate"):
        input_paths.extend(str(path.relative_to(REPO)) for path in (REPO/directory).glob("*.py"))
    input_paths = sorted(set(input_paths))
    # SI constants are exact defining constants; Planck length is illustrative.
    c0, h, elementary_charge = 299792458., 6.62607015e-34, 1.602176634e-19
    hbar_c_eVm = h/(2*np.pi)*c0/elementary_charge
    wavelength = 1064e-9
    wave_number = 2*np.pi/wavelength
    q_planck = 1.616255e-35*wave_number
    output = {
        "date": "2026-11-02",
        "scope": "Data-free IPI note audit; numerical controls supplement exact repo checks.",
        "source_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
        "working_tree_notice": "Audited working-tree bytes, not a claim that all inputs equal HEAD.",
        "input_sha256": {p: hashlib.sha256((REPO/p).read_bytes()).hexdigest() for p in input_paths},
        "original_custody_hashes_verified": custody_checks,
        "external_bitcoin_anchoring_reverified": False,
        "independent_geometry": {"vertices":12,"edges":30,"faces":20,"directions_tested":512,
            "max_moment_abs_error":moment_residual, "sampled_extrema": [float(iv.min()),float(iv.max())],
            "laplacian_characteristic_polynomial_exact":str(s.factor(characteristic)),
            "smallest_positive_laplacian_eigenvalue":"5-sqrt(5), multiplicity 3"},
        "positive_mixture_controls": {"count":512,"min_isotropic_ratio":min(t[0] for t in ratios),
            "min_anisotropic_ratio":min(t[1] for t in ratios), "max_anisotropic_ratio":max(t[1] for t in ratios)},
        "zero_anisotropy_positive_mixture": str(zero_mix),
        "multi_radius_cross_order_lock_residual": str(lock_residual),
        "stable_signed_stencil_C4_over_a2": str(signed_C4),
        "edge_remainder_128_samples_max_bound_fraction": str(max_fraction),
        "frequency_series": str(omega), "radial_group_speed_series": str(radial_speed),
        "log_generator_quartic": str(A+gamma/2),
        "units": {"hbar_c_eV_m":hbar_c_eVm,"1064nm_energy_eV":h*c0/wavelength/elementary_charge,
                  "1064nm_planck_q":q_planck,
                  "1064nm_edge_planck_phase_anisotropy_peak_to_peak":q_planck**4/16200},
        "verification_log": "verification.txt",
        "conclusion": "Algebra reproduced on its domain. No physical arming or OPH-wide exclusion theorem established."
    }
    target = HERE/"audit_results.json"
    target.write_text(json.dumps(output, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS", "output":str(target),
                      "moment_error":moment_residual,"custody_files_verified":len(custody_checks),
                      "edge_remainder_fraction":float(max_fraction)}))


if __name__ == "__main__":
    main()
