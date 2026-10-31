"""Static charged quasiparticles of the actual calibrated Weyl/Higgs action.

The original 180 Nambu directions are retained.  At rest six unpaired
neutrino annihilator modes fix the zero-energy vacuum convention.  The
other 84 positive directions and their particle-hole partners diagonalize
the supplied quadratic action.  No target mass enters the eigensolver.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from ..common_source import calibrated_matter as matter
from ..common_source.chiral_candidate import action, kinetic


MODE_COUNT = 90
NAMBU_COUNT = 180
ABSOLUTE_OPERATOR_GUARD = 2e-10
ORTHOGONALITY_GUARD = 2e-11
GROUP_TOLERANCE = 2e-9


def _numeric(matrix):
    return np.asarray(matrix.evalf(17) if hasattr(matrix, "evalf") else matrix,
                      dtype=complex)


def _immutable(value):
    value = np.array(value, copy=True)
    value.flags.writeable = False
    return value


def particle_hole(columns):
    """The antiunitary Nambu involution C(u,v)=(conj(v),conj(u))."""
    value = np.asarray(columns)
    if value.ndim not in (1, 2) or value.shape[0] != NAMBU_COUNT:
        raise ValueError("180-row Nambu vector or columns required")
    return np.concatenate((value[MODE_COUNT:].conj(), value[:MODE_COUNT].conj()), axis=0)


def _phase_columns(value):
    result = value.copy()
    for column in range(result.shape[1]):
        pivot = int(np.argmax(abs(result[:, column])))
        phase = result[pivot, column]
        if abs(phase):
            result[:, column] *= phase.conjugate()/abs(phase)
    return result


def _clusters(values):
    groups = []
    for index, energy in enumerate(values):
        if not groups or abs(energy-values[groups[-1][0]]) > GROUP_TOLERANCE:
            groups.append([index])
        else:
            groups[-1].append(index)
    return groups


@lru_cache(maxsize=1)
def quasiparticles():
    """Return the common static action, canonical chart and complete probes.

    ``energies`` and ``charge6`` index the 90 quasiparticle colors.
    ``positive_basis`` maps them into (c,c-dagger) Nambu columns.
    ``full_basis`` appends every particle-hole partner; no negative or zero
    direction is omitted.  Degenerate eigenvectors are an arbitrary chart;
    ``projectors`` are the invariant energy/charge blocks.
    """
    delta = _numeric(matter.pairing_kernel())
    if delta.shape != (MODE_COUNT, MODE_COUNT):
        raise ArithmeticError("the calibrated field changed its 90-mode domain")
    zero = np.zeros_like(delta)
    bdg = np.block([[zero, delta], [delta.conj().T, zero]])
    original_charge = _numeric(action.spin_electric_charge())
    charge6_original = np.rint(6*np.diag(original_charge).real).astype(int)
    if np.linalg.norm(original_charge-np.diag(charge6_original/6)) > ORTHOGONALITY_GUARD:
        raise ArithmeticError("the original electric-charge chart is not exact sixth-integer diagonal")
    charge6_nambu = np.concatenate((charge6_original, -charge6_original))
    charge_nambu = np.diag(charge6_nambu/6)
    if (np.linalg.norm(delta+delta.T) > ABSOLUTE_OPERATOR_GUARD
            or np.linalg.norm(bdg@charge_nambu-charge_nambu@bdg) > ABSOLUTE_OPERATOR_GUARD):
        raise ArithmeticError("the source pairing is not antisymmetric and electrically invariant")
    columns, energies, charges = [], [], []
    zero_modes = np.flatnonzero(charge6_original == 0)
    if len(zero_modes) != 6 or np.any(delta[zero_modes]) or np.any(delta[:, zero_modes]):
        raise ArithmeticError("the declared six massless neutrino annihilator modes changed")
    for charge6 in sorted(set(charge6_nambu)):
        support = np.flatnonzero(charge6_nambu == charge6)
        if charge6 == 0:
            if len(support) != 12 or np.any(bdg[support]) or np.any(bdg[:, support]):
                raise ArithmeticError("the entire neutral Nambu block must remain zero")
            block_columns = np.eye(NAMBU_COUNT, dtype=complex)[:, zero_modes]
            values = np.zeros(len(zero_modes))
        else:
            values, vectors = np.linalg.eigh(bdg[np.ix_(support, support)])
            half = len(values)//2
            if (len(values) % 2 or values[half-1] >= -GROUP_TOLERANCE
                    or values[half] <= GROUP_TOLERANCE):
                raise ArithmeticError("charged block lost its separated positive/negative spectrum")
            if np.max(abs(values[:half]+values[half:][::-1])) > ABSOLUTE_OPERATOR_GUARD:
                raise ArithmeticError("charged block lost its opposite-energy pairing")
            values = values[half:]
            block_columns = np.zeros((NAMBU_COUNT, half), complex)
            block_columns[support] = _phase_columns(vectors[:, half:])
        columns.append(block_columns)
        energies.extend(values)
        charges.extend([charge6]*len(values))
    positive = np.concatenate(columns, axis=1)
    energies = np.asarray(energies)
    charges = np.asarray(charges, dtype=int)
    if positive.shape != (NAMBU_COUNT, MODE_COUNT):
        raise ArithmeticError("the full positive/zero half has not been retained")
    complete = np.concatenate((positive, particle_hole(positive)), axis=1)
    diagonal = np.diag(np.concatenate((energies, -energies)))
    if np.linalg.norm(complete.conj().T@complete-np.eye(NAMBU_COUNT)) > ORTHOGONALITY_GUARD:
        raise ArithmeticError("the quasiparticle chart is not a canonical Bogoliubov basis")
    if np.linalg.norm(bdg@complete-complete@diagonal) > ABSOLUTE_OPERATOR_GUARD:
        raise ArithmeticError("the full action is not reconstructed by all quasiparticle directions")
    electric = positive.conj().T@charge_nambu@positive
    if np.linalg.norm(electric-np.diag(charges/6)) > ORTHOGONALITY_GUARD:
        raise ArithmeticError("quasiparticle electric charges do not match the full probe")
    blocks = []
    for charge6 in sorted(set(charges)):
        support = np.flatnonzero(charges == charge6)
        for group in _clusters(energies[support]):
            indices = support[group]
            vectors = positive[:, indices]
            blocks.append({"charge6": int(charge6), "energy": float(np.mean(energies[indices])),
                           "multiplicity": len(indices), "indices": tuple(map(int, indices)),
                           "projector": _immutable(vectors@vectors.conj().T)})
    lepton_blocks = [block for block in blocks if block["charge6"] == -6]
    if len(lepton_blocks) != 3 or [b["multiplicity"] for b in lepton_blocks] != [2, 2, 2]:
        raise ArithmeticError("the charged-lepton quasiparticle sectors changed")
    return {"energies": _immutable(energies), "charge6": _immutable(charges),
            "positive_basis": _immutable(positive), "full_basis": _immutable(complete),
            "bdg": _immutable(bdg), "pairing": _immutable(delta),
            "electric_charge": _immutable(electric),
            "original_nambu_charge": _immutable(charge_nambu),
            "projectors": tuple(blocks), "zero_annihilator_indices": tuple(map(int, zero_modes)),
            "ramsey_indices": (lepton_blocks[0]["indices"][0], lepton_blocks[1]["indices"][0]),
            "vacuum_energy": -float(sum(energies))/2,
            "operator_error_guard": ABSOLUTE_OPERATOR_GUARD,
            "orthogonality_guard": ORTHOGONALITY_GUARD}


def transform_one_body_probe(operator):
    """Retain all quasiparticle-number-preserving and pair-changing blocks.

    A Hermitian ordinary one-body probe O has Nambu representation
    diag(O,-O.T), plus its usual scalar offset in the quadratic operator.
    Broken gauge/Higgs probes generally do not preserve quasiparticle N=2.
    """
    operator = _numeric(operator)
    if (operator.shape != (MODE_COUNT, MODE_COUNT) or not np.isfinite(operator).all()
            or np.linalg.norm(operator-operator.conj().T) > ORTHOGONALITY_GUARD):
        raise ValueError("a finite Hermitian 90-mode one-body probe required")
    zero = np.zeros_like(operator)
    original = np.block([[operator, zero], [zero, -operator.T]])
    basis = quasiparticles()["full_basis"]
    transformed = basis.conj().T@original@basis
    return {"full": _immutable(transformed),
            "number_preserving": _immutable(transformed[:MODE_COUNT, :MODE_COUNT]),
            "pair_creation": _immutable(transformed[:MODE_COUNT, MODE_COUNT:]),
            "original_trace_offset": float(np.trace(operator).real)/2,
            "quasiparticle_vacuum_offset": float((np.trace(operator)
                -np.trace(transformed[:MODE_COUNT, :MODE_COUNT])).real)/2,
            "pair_block_norm": float(np.linalg.norm(transformed[:MODE_COUNT, MODE_COUNT:]))}


def moving_symbol(momentum):
    """The same calibrated mass action with supplied c=1 Weyl kinematics.

    This finite continuum symbol does not identify native source positions
    with momentum.  A physical source realization at nonzero p requires
    both p and -p CAR modes; it is not the 90-mode rest process.
    """
    return kinetic.bdg_symbol(momentum, matter.mass_kernel(), speed=1.)


def validate_chart(positive, energies, charge6):
    """Independently admit a complete chart, including degenerate rotations.

    This checks the CAR, the entire action and the electric observable. It
    does not compare eigenvectors with those of a preferred diagonalizer.
    All six neutral annihilators and their conjugates must be present.
    """
    positive, energies, charge6 = map(np.asarray, (positive, energies, charge6))
    if (positive.shape != (NAMBU_COUNT, MODE_COUNT) or energies.shape != (MODE_COUNT,)
            or charge6.shape != (MODE_COUNT,) or not np.isfinite(positive).all()
            or not np.isfinite(energies).all() or not np.isfinite(charge6).all()
            or np.any(energies < 0) or np.any(charge6 != np.rint(charge6))):
        raise ValueError("a complete finite 90-mode chart with nonnegative energies required")
    reference = quasiparticles()
    basis = np.column_stack((positive, particle_hole(positive)))
    if np.linalg.norm(basis.conj().T@basis-np.eye(NAMBU_COUNT)) > ORTHOGONALITY_GUARD:
        raise ValueError("chart violates CAR completeness")
    spectrum = np.r_[energies, -energies]
    if np.linalg.norm(reference['bdg']@basis-basis*spectrum) > ABSOLUTE_OPERATOR_GUARD:
        raise ValueError("chart does not reconstruct the supplied action")
    if np.linalg.norm(reference['original_nambu_charge']@positive-positive*(charge6/6)) > ORTHOGONALITY_GUARD:
        raise ValueError("chart does not reconstruct the supplied electric observable")
    return basis


@lru_cache(maxsize=1)
def current_probe_data():
    """All twelve full Nambu currents and the exact Higgs stabilizer.

    The supplied primitive generator basis mixes broken and unbroken
    directions. Its nine-dimensional unbroken span is determined from the
    actual fundamental generators annihilating the chosen neutral Higgs
    direction, without thresholding a numerical current spectrum.
    """
    import sympy as s

    fundamental = action._representation()[1]
    higgs = s.eye(5)[:, 4]
    columns = []
    for generator in fundamental:
        column = generator*higgs
        columns.append(column.applyfunc(s.re).col_join(column.applyfunc(s.im)))
    higgs_map = s.Matrix.hstack(*columns)
    kernel = higgs_map.nullspace()
    if higgs_map.rank() != 3 or len(kernel) != 9:
        raise ArithmeticError("the exact neutral-Higgs stabilizer changed")
    probes = tuple(transform_one_body_probe(1j*np.kron(np.eye(2), _numeric(g)))
                   for g in action.generators())
    anomalous = np.column_stack([np.r_[p['pair_creation'].real.ravel(),
                                        p['pair_creation'].imag.ravel()] for p in probes])
    singular = np.linalg.svd(anomalous, compute_uv=False)
    if singular[2] < 1 or singular[3] > ORTHOGONALITY_GUARD:
        raise ArithmeticError("the full current anomalous map does not have rank three")
    coefficients = np.column_stack([_numeric(v).real[:, 0] for v in kernel])
    if np.linalg.norm(anomalous@coefficients) > ORTHOGONALITY_GUARD:
        raise ArithmeticError("the exact Higgs stabilizer does not preserve quasiparticle number")
    return {'probes': probes, 'exact_stabilizer_coefficients': tuple(tuple(map(str, v)) for v in kernel),
            'stabilizer_coefficients': _immutable(coefficients),
            'anomalous_singular_values': _immutable(singular),
            'exact_higgs_orbit_rank': 3, 'exact_stabilizer_dimension': 9}
