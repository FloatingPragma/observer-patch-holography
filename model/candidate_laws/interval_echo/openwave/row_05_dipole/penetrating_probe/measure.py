"""One-opening reflection null and two-opening current-sensitive scattering."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy import linalg

from model.executable_model import cell, federation
from model.experiments.material_inertia import neutral_defect
from model.experiments.material_inertia.pair_spectrum.cell_field_probe.measure import pointer_decoder
from ..experiment import gates_of, source


HERE = Path(__file__).parent
ARCHIVE = HERE / "results.json"
TIMES = (0, 8, 16, 32)
REFLECTION = (9, 1, 2, 10, 11, 5, 6, 8, 7, 0, 3, 4)
TOLERANCE = 2e-9


def reflection_certificate():
    import sympy as sp

    primitive = cell()
    permutation = np.array(REFLECTION)
    adjacency = primitive.adjacency().toarray().astype(int)
    if (not np.array_equal(permutation[permutation], np.arange(12))
            or not np.array_equal(adjacency[np.ix_(permutation, permutation)], adjacency)
            or permutation[1] != 1 or set(permutation[[3, 10]]) != {3, 10}):
        raise ArithmeticError("The proposed one-opening reflection is not a native graph symmetry")
    equations = sp.Matrix(np.vstack((np.eye(12, dtype=int)+adjacency, np.eye(12, dtype=int)[[1, 3, 10]])).tolist())
    basis = sp.Matrix.hstack(*equations.nullspace())
    if basis.cols != 2:
        raise ArithmeticError("The stationary current source has a different compact dimension")
    reflected = basis[list(REFLECTION), :]
    action = (basis.T*basis).inv()*basis.T*reflected
    if (basis*action != reflected or action.trace() != 0 or action.det() != -1
            or action*action != sp.eye(2)):
        raise ArithmeticError("The compact source reflection is not orientation reversing")
    board = federation(0)
    if gates_of(board, 0) != (1, 3, 10):
        raise ArithmeticError("The source gateways differ from the declared apparatus")
    plus, _, _, _ = source((1, 3, 10))
    phase = np.vdot(plus.conj(), plus[permutation])
    error = float(np.linalg.norm(plus[permutation]-phase*plus.conj()))
    if error > 2e-12 or abs(abs(phase)-1) > 2e-12:
        raise ArithmeticError("The reflection does not exchange the actual source preparations")
    return {"port_permutation": list(REFLECTION), "exact_graph_automorphism": True,
            "exact_compact_action": [[str(action[i, j]) for j in range(2)] for i in range(2)],
            "compact_trace": 0, "compact_determinant": -1,
            "conjugate_source_phase": [float(phase.real), float(phase.imag)],
            "source_mapping_error": error}, plus, phase


def region(board, x, y, kinetic, h, owners):
    ports = np.flatnonzero(np.isin(np.arange(board.n)//12, owners))
    selected = np.flatnonzero(np.isin(x, ports) & np.isin(y, ports))
    return x[selected], y[selected], kinetic[selected][:, selected].toarray(), h[selected][:, selected].toarray()


def initial(x, y, plus, two_openings):
    probe = {82: 1/np.sqrt(2), 48: 1j/np.sqrt(2)} if two_openings else {82: 1.}
    columns = np.zeros((len(x), 2), complex)
    for sign, vacancy in enumerate((plus, plus.conj())):
        for q, amplitude in probe.items():
            keep = (x == q) & (y < 12)
            columns[keep, sign] = amplitude*vacancy[y[keep]]
    if not np.allclose((abs(columns)**2).sum(axis=0), 1, rtol=0, atol=2e-13):
        raise ArithmeticError("Probe/source product preparation lost a branch")
    return columns


def terminal(x, y, state, mode, detector):
    weights = abs(state)**2
    if abs(weights.sum()-1) > 2e-12:
        raise ArithmeticError("Scattering flight does not normalize")
    excess = float(weights[x == detector].sum())
    vacancy = float(weights[y == detector].sum())
    host = float(weights[(x != detector) & (y != detector)].sum())
    exterior = x >= 12
    inside = float(weights[~exterior].sum())
    # The source analyzer is only offered after the excess is outside its
    # owner. Then every compact vacancy address is available and this is an
    # ordinary orthogonal source-mode projector, without hard-core overlap.
    amplitudes = {}
    opposite_weight, other_weight = 0.0, 0.0
    for q in np.unique(x[exterior]):
        keep = x == q
        row_mode = np.zeros(keep.sum(), complex)
        local = y[keep] < 12
        row_mode[local] = mode[y[keep][local]]
        amplitude = np.vdot(row_mode, state[keep])
        opposite_amplitude = np.vdot(row_mode.conj(), state[keep])
        amplitudes[int(q)] = amplitude
        opposite_weight += float(abs(opposite_amplitude)**2)
        other_weight += float(np.linalg.norm(state[keep]-amplitude*row_mode-opposite_amplitude*row_mode.conj())**2)
    retained = float(sum(abs(value)**2 for value in amplitudes.values()))
    joint = float(abs(amplitudes.get(detector, 0))**2)
    source_miss = inside+opposite_weight+other_weight
    success = {"port_load_zero_host_excess": [vacancy, host, excess],
               "probe_exterior_inside": [float(weights[exterior].sum()), inside],
               "source_original_opposite_other_then_probe_inside": [retained, opposite_weight, other_weight, inside],
               "correct_source_and_exterior_then_complement": [retained, source_miss],
               "correct_source_and_excess_at_detector_then_complement": [joint, source_miss+float(sum(abs(value)**2 for q, value in amplitudes.items() if q != detector))]}
    for values in success.values():
        if min(values) < 0 or abs(sum(values)-1) > 2e-12:
            raise ArithmeticError("A complete scattering detector has an invalid outcome")
    return success


def run():
    reflection, plus, phase = reflection_certificate()
    native_pointer = pointer_decoder()
    board = federation(0)
    xfull, yfull, _, kinetic, h, _ = neutral_defect.generators(board)
    studies = []
    for owners, detector in (((0, 6), 82), ((0, 4, 6), 48)):
        two = len(owners) == 3
        x, y, zero, native = region(board, xfull, yfull, kinetic, h, owners)
        columns = initial(x, y, plus, two)
        lookup = {(int(a), int(b)): j for j, (a, b) in enumerate(zip(x, y, strict=True))}
        reflected_x = np.where(x < 12, np.array(REFLECTION)[np.minimum(x, 11)], x)
        reflected_y = np.where(y < 12, np.array(REFLECTION)[np.minimum(y, 11)], y)
        pair_map = np.array([lookup[int(a), int(b)] for a, b in zip(reflected_x, reflected_y, strict=True)])
        settings = []
        for strength, matrix in ((8, native), (0, zero)):
            exact_commutes = np.array_equal(matrix[np.ix_(pair_map, pair_map)], matrix)
            if not two and not exact_commutes:
                raise ArithmeticError("The complete one-opening native operator breaks the exact reflection")
            energies, modes = linalg.eigh(matrix)
            residual = float(linalg.norm(matrix@modes-modes*energies))
            orthogonality = float(linalg.norm(modes.T@modes-np.eye(len(modes))))
            allowance = 2*(max(TIMES)*residual+orthogonality+1e-11)
            if allowance > 1e-9:
                raise ArithmeticError("The finite spectral propagation exceeds its numerical allowance")
            coefficients = modes.T@columns
            readings = []
            for exposure in TIMES:
                wave = modes@(np.exp(-1j*(energies-energies[0])*exposure)[:, None]*coefficients)
                arms = [terminal(x, y, wave[:, j], mode, detector)
                        for j, mode in enumerate((plus, plus.conj()))]
                source_relation = float(np.linalg.norm(wave[:, 1]-phase.conjugate()*wave[pair_map, 0]))
                if not two and source_relation > 2e-11:
                    raise ArithmeticError("The one-opening propagated states violate source reflection")
                difference = arms[0]["port_load_zero_host_excess"][2]-arms[1]["port_load_zero_host_excess"][2]
                readings.append({"exposure": exposure, "source_plus_minus": arms,
                                 "excess_probability_difference": difference,
                                 "one_opening_state_symmetry_error": source_relation if not two else None,
                                 "equal_source_mixture": {key: [(a+b)/2 for a, b in zip(arms[0][key], arms[1][key], strict=True)] for key in arms[0]}})
            settings.append({"kappa": strength, "exact_reflection_commutes": exact_commutes,
                             "eigensystem_residual_frobenius": residual,
                             "orthogonality_error_frobenius": orthogonality,
                             "numerical_probability_allowance": allowance,
                             "readings": readings})
        studies.append({"owners": list(owners), "ports": 12*len(owners), "dimension": len(x),
                        "detector_port": detector, "settings": settings})
    primary = studies[1]["settings"][0]["readings"][-1]
    gate = bool(abs(primary["excess_probability_difference"])-studies[1]["settings"][0]["numerical_probability_allowance"] > .01)
    return {"schema": "oph.interval_echo.penetrating_current_probe.v1",
            "protocol_sha256": sha256((HERE/"PROTOCOL.md").read_bytes()).hexdigest(),
            "measure_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "reflection_certificate": reflection, "native_local_decoder": native_pointer,
            "scope": "Principal restrictions with complete native diagonals. Supplied number apertures, coherent sources, phase and exposure appointments. Ideal finite-record probabilities, no sampled trials, imposed flux, static dipole B or electron g claim.",
            "source_return_scope": "Supplied ideal coherent exterior-number flag and local compact-mode analyzers; native gate compilation is outside this study. Correct original compact vacancy mode AND exterior excess, with all other outcomes retained. No unsupported bare source projector when the excess occupies its owner.",
            "studies": studies, "two_opening_reversal_gate": {"minimum_absolute_difference": .01, "passed": gate}}


def compare(actual, archived):
    if type(actual) is not type(archived):
        raise ArithmeticError("Receipt type differs")
    if isinstance(actual, dict):
        if actual.keys() != archived.keys():
            raise ArithmeticError("Receipt fields differ")
        for key in actual:
            compare(actual[key], archived[key])
    elif isinstance(actual, list):
        if len(actual) != len(archived):
            raise ArithmeticError("Receipt outcome partition differs")
        for left, right in zip(actual, archived, strict=True):
            compare(left, right)
    elif isinstance(actual, float):
        if not np.isfinite(archived) or abs(actual-archived) > TOLERANCE:
            raise ArithmeticError("Receipt probability differs")
    elif actual != archived:
        raise ArithmeticError("Receipt metadata differs")


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record", type=Path, help="explicit exclusive first-record creation")
    args = cli.parse_args()
    result = run()
    payload = json.dumps(result, indent=2, allow_nan=False)+"\n"
    if args.record:
        with args.record.open("x", encoding="utf-8") as destination:
            destination.write(payload)
    else:
        compare(result, json.loads(ARCHIVE.read_text()))
    if not result["two_opening_reversal_gate"]["passed"]:
        raise ArithmeticError("The two-opening source-reversal gate failed")
    for study in result["studies"]:
        for setting in study["settings"]:
            print(study["owners"], setting["kappa"], setting["readings"][-1]["excess_probability_difference"])


if __name__ == "__main__":
    main()
