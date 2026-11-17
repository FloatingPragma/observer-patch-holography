"""Exact finite family instrument and declared unconditional collision clock."""
from functools import lru_cache
import mpmath as mp
import numpy as np
import sympy as s
from model.coherent_seam_repair.program import Observer
from . import primitives as fc


def sn(matrix):
    return np.array(matrix.evalf(18).tolist(), dtype=complex)


def spectrum(masses):
    vals = sorted(mp.mpf(x) for x in masses)
    return {"m2/m1": str(vals[1] / vals[0]), "m3/m1": str(vals[2] / vals[0]),
            "Q": str(sum(vals) / sum(mp.sqrt(x) for x in vals) ** 2)}


def native_event_check():
    observer = Observer()
    local = fc.repair()
    checked = 0
    for address, (a, b) in enumerate(observer.addresses):
        idle = next(port for port in range(12) if port not in (a,b))
        ports = (a,b,idle)
        for load in range(3):
            for record in range(2):
                actual = observer.repair(observer.basis_state(ports[load],record,address)).amplitudes
                expected = {}
                for target in range(3):
                    for rout in range(2):
                        coefficient = local[fc.index(target,rout,0),fc.index(load,record,0)]
                        if coefficient:
                            expected[observer.index(ports[target],rout,address)] = coefficient
                assert {i:v for i,v in enumerate(actual) if v} == expected
                checked += 1
    return {"actual_Observer_repair_columns": checked, "exact_mismatches": 0,
            "two_address_reflection_is_original_60_address_reflection": False}


@lru_cache(maxsize=1)
def calculate():
    mp.mp.dps = 70
    native = native_event_check()
    W = fc.repair() * fc.address_reflection()
    V = fc.oriented_record()
    B = V * W
    assert B.H * B == s.eye(12)
    prep = fc.state(s.sqrt(s.Rational(1, 3)), 0, s.sqrt(s.Rational(2, 3)))
    after = B * prep
    phi = fc.analyzer()
    assert (phi.H * phi)[0] == 1
    M = fc.compression(phi, after)
    C = fc.circulant()
    assert (M - s.sqrt(2) * C / 3).applyfunc(s.simplify) == s.zeros(3)
    assert M * fc.shift() == fc.shift() * M
    assert fc.shift()**3 == s.eye(3)
    adjoint=M.H.subs(s.conjugate(fc.Z),1/fc.Z)
    assert (adjoint-M).applyfunc(s.simplify)==s.zeros(3)
    assert s.simplify(s.trace(adjoint*M))==s.Rational(2,9)
    # Exact two-tick identity for every phase and family character. Each
    # character reduces the controlled shift to a scalar zeta.
    zeta = s.Symbol("zeta", nonzero=True)
    D = s.diag(*[zeta ** ((1 if i % 2 == 0 else -1) if (i // 2) % 2 else 0)
                 for i in range(12)])
    c2, e2 = D*B*D*B*prep, B*B*prep
    flags = fc.flags()
    k2 = ((flags["C_idle"].H*c2)[0]*(flags["E_agree"].H*e2)[0]
          + (flags["C_rep+"].H*c2)[0]*(flags["E_idle"].H*e2)[0]
          + (flags["C_rep-"].H*c2)[0]*(flags["E_idle"].H*e2)[0])/s.sqrt(3)
    assert s.simplify(k2-s.sqrt(3)/9) == 0
    spectral_variable = s.Symbol("eigenvalue")
    poly = (D*B).charpoly(spectral_variable).as_expr()
    expected_poly = ((spectral_variable-1)**2*(spectral_variable+1)**2
                     *(spectral_variable**2+1)**2*(spectral_variable**2+zeta)
                     *(zeta*spectral_variable**2-1)/zeta)
    assert s.factor(poly-expected_poly) == 0

    # Character-block form of the actual 432-state U. It reduces to three
    # 144-state blocks, because U commutes with the neutral family shift.
    theta = mp.mpf(2) / 9
    z = complex(mp.exp(1j * theta))
    Bn = sn(B)
    base = np.kron(Bn, Bn)
    p = sn(s.kronecker_product(prep, prep)).ravel()
    ph = sn(phi).ravel()
    blocks = []
    for k in range(3):
        char = np.exp(2j * np.pi * k / 3)
        diagonal = []
        for i in range(12):
            rec, orient = (i // 2) % 2, i % 2
            step = (1 if orient == 0 else -1) if rec else 0
            diagonal.extend([(z * char) ** step] * 12)
        uk = np.asarray(diagonal)[:, None] * base
        assert np.max(np.abs(uk.conj().T @ uk - np.eye(144))) < 2e-14
        blocks.append(uk)
    lambdas = [1 / mp.sqrt(6) + 1 / mp.sqrt(3) * mp.cos(theta + 2 * mp.pi * k / 3)
               for k in range(3)]
    mu = np.array([np.vdot(ph, uk @ p) for uk in blocks])
    target = np.array([float(s.sqrt(2) / 3) * float(x) for x in lambdas])
    assert np.max(np.abs(mu - target)) < 2e-15
    multi = []
    for n in range(1, 5):
        vals = [np.vdot(ph, np.linalg.matrix_power(uk, n) @ p) for uk in blocks]
        multi.append({"ticks": n, "compression_eigenvalues": [[v.real, v.imag] for v in vals],
                      "distance_from_M_power": float(np.max(np.abs(np.array(vals) - mu ** n))),
                      "squared_amplitude_ratios": {key:float(value) for key,value in
                                                   spectrum([abs(v) ** 2 for v in vals]).items()}})
    assert abs(complex(*multi[1]["compression_eigenvalues"][0])
               - complex(*multi[1]["compression_eigenvalues"][1])) < 2e-15

    # A single selected matrix element is an effect, with a complete binary
    # outcome instrument {M, sqrt(I-M^dag M)} available by dilation.
    energy = np.abs(mu) ** 2
    success = {"family_character_k": list(range(3)), "probability": energy.tolist(),
               "no_operation_probability": (1 - energy).tolist(),
               "tracial_success_probability": float(sum(energy) / 3)}

    # Declared compiler: P_k = U_k^dag |phi><phi| U_k is a projector.
    # Each fresh ancillary object pair undergoes exp(-i dt P_k), is retained
    # in the complete output, and is traced only for the family marginal.
    # Its compression is A=M^dag M. No success branch is selected.
    us = [uk.conj().T @ ph for uk in blocks]
    assert all(abs(np.vdot(v, v) - 1) < 2e-14 for v in us)
    p_projected = [v * np.vdot(v, p) for v in us]
    assert np.max(np.abs([np.vdot(p, q) for q in p_projected] - energy)) < 2e-15
    t = 1.0
    exact_schur = np.exp(-1j * t * (energy[:, None] - energy[None, :]))
    collision_rows = []
    for n in (1, 2, 10, 100, 1000, 10000):
        delta = t / n
        states = np.asarray([p + np.expm1(-1j * delta) * q for q in p_projected])
        # r[j,k] = <ancilla_k|ancilla_j>; all ancilla outputs retained.
        r = states @ states.conj().T
        rn = r ** n
        assert np.max(np.abs(np.diag(r) - 1)) < 2e-14
        assert np.min(np.linalg.eigvalsh(r)) > -1e-13
        err = float(np.max(np.abs(rn - exact_schur)))
        # This entrywise norm is a diagnostic. The analytic diamond bound
        # proved in PROOF.md is 4*t^2/n and applies to arbitrary references.
        collision_rows.append({"N": n, "max_coherence_multiplier_error": err,
                               "diamond_error_upper_bound": min(2.0, 4 * t * t / n),
                               "min_coherence_modulus": float(np.min(np.abs(rn)))})
        assert err <= 4 * t * t / n + 1e-10

    # Exact polynomial test at z=1: U has rational entries and the complete
    # instrument compression is checked in symbolic radical arithmetic.
    exact0 = {"B_unitary": True, "analyzer_norm": "1", "compression": str(M),
              "M_equals_sqrt2_over3_C": True, "R_order": 3,
              "trace_M_dagger_M": "2/9", "tracial_success_probability": "2/27",
              "two_tick_compression_for_all_phases": "sqrt(3)/9 * I3",
              "family_characters": ["1", "exp(2*pi*i/3)", "exp(4*pi*i/3)"],
              "physical_generations_from_characters": False,
              "read_probabilities": {"R+": "1/12", "H": "1/6", "independent_XOR": "2/9"},
              "clock_coefficient_inferred_from_read_probability": False}
    out = {"domain": "declared FR family instrument and explicit collision extension",
           "native_event": native,
           "exact_instrument": exact0, "one_tick_squared_amplitude_spectrum": spectrum([x*x for x in lambdas]),
           "C_eigenvalues": [str(x) for x in lambdas], "multi_tick": multi,
           "actual_FR_Floquet": {"character_sector_factorization": "U_k=(D(z*omega^k)B) tensor B",
                                 "C_factor_characteristic_polynomial": str(s.factor(poly)),
                                 "E_factor_characteristic_polynomial": "(eigenvalue^4-1)^3",
                                 "FR_report_arccos_readout_is_eigenphase_of_this_U": False},
           "neutral_background_comparison": {"Y=y*C_rest_parameters": spectrum(lambdas),
                                              "Y=y*C_squared_rest_parameters": spectrum([x*x for x in lambdas]),
                                              "these_are_different_Yukawa_operators": True,
                                              "v_convention": "v=v_EW/sqrt(2) for H=(0,v)"},
           "complete_binary_instrument": success,
           "collision_clock": {"Hamiltonian": "A=M^dag M=(2/9) C(theta)^2",
                               "ancilla_generator": "P=U^dag(|Phi><Phi| tensor I)U",
                               "input": "fresh copies of the same prepared object pair",
                               "discarded_outcomes": [], "time": t, "rows": collision_rows,
                               "fixed_FR_tick_has_this_generator": False,
                               "primary_Wilson_source_transfer": False},
           "empirical_inputs": [], "phase_status": "declared theta=2/9; no fit in this calculation",
           "architecture_ancestry": "balanced circulant is an inherited mass-informed architecture; no blind-prediction claim",
           "register_and_Fock_are_distinct": {
               "one_family_lepton_modes":6,"Fock_K_tensor_common_register_dimension":3*2**6,
               "Fock_K_tensor_three_family_one_particle_space_dimension":2**18,
               "register_identified_with_independently_occupiable_generations":False}}
    return out
