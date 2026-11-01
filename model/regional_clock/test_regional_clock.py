from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import actuator as a, independent as ind
from .verify import retained_weights
from ..local_repair_family import modular_clock, learned_flow
from ..native_program import interpreter

TRUE = (25, 25)+(8,)*10
R, S = a.REGIONS


@pytest.fixture(scope='module')
def learned():
    return retained_weights()


@pytest.fixture(scope='module')
def arrays(learned):
    return {ports: a.executed(learned, ports) for ports in a.REGIONS}


def test_all_number_sectors_and_exterior_units(learned):
    assert tuple(modular_clock.law.basis()) == ind.STATES
    for ports in a.REGIONS:
        certificate = ind.algebra_certificate(learned, ports)
        assert certificate['regional_matrix_units_by_number'] == [1, 4, 1]
        assert certificate['exterior_matrix_units'] == 1+100+2025
        rho = ind.reduced_reference(TRUE, ports)
        assert sum(rho.values()) == 1 and len(rho) == 4
        assert all(x > 0 for x in rho.values())
    # Independent integer census of the actual stationary law.
    assert sum(math.prod(w**n for w, n in zip(TRUE, x, strict=True)) for x in ind.STATES) == 7505


def test_global_restriction_is_insufficient():
    global_ratios = tuple(math.prod(F(w)**n for w, n in zip(TRUE, x, strict=True)) for x in ind.STATES)
    rho = ind.reduced_reference(TRUE, R)
    for _, u, v, support in ind.units(R):
        assert all(global_ratios[i]/global_ratios[j] == rho[u]/rho[v] for i, j in support)
    with pytest.raises(ArithmeticError, match='exterior'):
        ind.algebra_certificate(TRUE, R, global_ratios)


def test_wrong_reference_law_and_reversed_phase():
    with pytest.raises(ArithmeticError, match='reference'):
        ind.algebra_certificate(TRUE, R, ind.target((1,)*12, R))
    with pytest.raises(ArithmeticError, match='reference'):
        ind.algebra_certificate(TRUE, R, tuple(1/x for x in ind.target(TRUE, R)))
    fair = modular_clock.law.Parameters()
    biased = modular_clock.parameters()
    # The same target is not the stationary clock for the original fair law.
    assert modular_clock.law.coin_probability(fair, R, 1) != modular_clock.law.coin_probability(biased, R, 1)


def test_omitted_sectors_do_not_define_a_complete_channel():
    one = np.diag([int(sum(x[p] for p in R) == 1) for x in ind.STATES])
    assert np.trace(one) == 20
    assert not np.array_equal(one.conj().T@one, np.eye(66))
    assert sorted(sum(sum(x[p] for p in R) == n for x in ind.STATES) for n in range(3)) == [1, 20, 45]


def test_independent_angle_and_literal_admission(learned):
    for ports in a.REGIONS:
        tape = a.word(learned, ports)
        cert = ind.angle_certificate(tape, learned, ports)
        assert len(tape) == 2042 and cert['ideal_CNOT_count'] == 510
        assert F(cert['analytic_operator_error_upper']) < F(1, 10**11)
        admitted = a.admitted(learned, ports)
        assert admitted['declared_resources']['proper_code_load_events'] == 2042
        assert admitted['declared_resources']['proper_code_export_events'] == 2042
        # Remove a real phase, retain the exact circuit: target check must bite.
        changed = [dict(x) for x in tape]
        index = max(range(len(changed)), key=lambda i: abs(changed[i]['angle']) if changed[i]['kind'] == 'Z' else -1)
        changed[index]['angle'] += .01
        with pytest.raises(ArithmeticError, match='analytic cap'):
            ind.angle_certificate(changed, learned, ports)
    assert a.prepare_request(learned, R)['status'] == 'ADMITTED'
    stopped = a.prepare_request(learned, R, F(2))
    assert stopped['status'] == 'COMPILER_STOP' and stopped['pulses_executed'] == 0


def test_all512_word_matches_independent_target_and_padding(learned, arrays):
    for ports, raw in arrays.items():
        expected = np.ones(512, complex)
        for address, q in zip(ind.CODES, ind.target(learned, ports), strict=True):
            expected[address] = np.exp(-1j*math.log(float(q)))
        assert np.linalg.norm(raw-np.diag(expected), 'fro') < 1e-9
        outside = tuple(x for x in range(512) if x not in ind.CODES)
        assert len(outside) == 446
        assert np.linalg.norm(raw[np.ix_(outside, outside)]-np.eye(446), 'fro') < 1e-9


def test_integrated_coefficient_mutation_exceeds_complete_source_gate(monkeypatch):
    compiler = modular_clock.common.full.compiler
    original = compiler.integrated_controls
    def changed(kind, angle):
        coefficients = original(kind, angle).copy()
        if kind == 'H':
            coefficients[0] += .01
        return coefficients
    monkeypatch.setattr(compiler, 'integrated_controls', changed)
    with pytest.raises(ArithmeticError, match='analytic cap'):
        ind.angle_certificate(a.word(TRUE, R), TRUE, R)


def exterior_reference_marginal(density, ports, r=2):
    """Independent full occupation trace, retaining all exterior numbers/R."""
    outside = tuple(p for p in range(12) if p not in ports)
    keys = sorted({tuple(x[p] for p in outside) for x in ind.STATES})
    index = {x: i for i, x in enumerate(keys)}
    result = np.zeros((len(keys)*r, len(keys)*r), complex)
    for i, x in enumerate(ind.STATES):
        for j, y in enumerate(ind.STATES):
            if tuple(x[p] for p in ports) == tuple(y[p] for p in ports):
                u, v = index[tuple(x[p] for p in outside)], index[tuple(y[p] for p in outside)]
                result[u*r:(u+1)*r, v*r:(v+1)*r] += density[i*r:(i+1)*r, j*r:(j+1)*r]
    return result


def test_entangled_exterior_and_reference_complete_export(learned, arrays):
    rng = np.random.default_rng(1138)
    psi = rng.normal(size=132)+1j*rng.normal(size=132)
    psi /= np.linalg.norm(psi)
    rho = np.outer(psi, psi.conj())
    for ports, raw in arrays.items():
        phases = np.array([np.exp(-1j*math.log(float(x))) for x in ind.target(learned, ports)])
        ideal_psi = np.repeat(phases, 2)*psi
        ideal = np.outer(ideal_psi, ideal_psi.conj())
        assert np.linalg.norm(exterior_reference_marginal(ideal, ports)-exterior_reference_marginal(rho, ports)) < 1e-13
        output = a.address_export(raw, rho, 2)
        assert set(output) == {('SUCCESS',), ('ADDRESS_ABORT',)}
        assert np.linalg.norm(output[('SUCCESS',)]-ideal, 'nuc')/2 < 1e-9
        assert abs(sum(np.trace(x) for x in output.values())-1) < 1e-9
        assert np.linalg.norm(output[('ADDRESS_ABORT',)]) < 1e-18


def test_success_conditioning_and_omitted_compiler_failures(monkeypatch):
    compiler = modular_clock.common.full.compiler
    # A hostile response rotates valid address bit 0 into processor complement
    # row 2. Its complete export is TP; its success alone is trace decreasing.
    leaked = np.eye(6, dtype=complex)
    c, s = math.cos(.4), math.sin(.4)
    leaked[0, 0] = leaked[2, 2] = c
    leaked[2, 0], leaked[0, 2] = s, -s
    monkeypatch.setattr(compiler, 'native_pulse', lambda *args: leaked.copy())
    tape = [compiler.pulse('Z', 0., [0])]
    branches = dict(a.first_failure_kraus(tape))
    assert len(branches) == 5  # four complement rows AND success
    effect = sum(k.conj().T@k for k in branches.values())
    assert np.linalg.norm(effect-np.eye(512)) < 1e-13
    success = branches[('SUCCESS',)]
    assert np.linalg.norm(success.conj().T@success-np.eye(512)) > .1
    # R entangled with addressed bit; discarded/normalized success changes R.
    psi = np.zeros((512, 2), complex)
    psi[0, 0] = psi[1, 1] = 1/math.sqrt(2)
    outputs = {key: k@psi for key, k in branches.items()}
    ref = sum(v.conj().T@v for v in outputs.values())
    assert np.linalg.norm(ref-np.eye(2)/2) < 1e-13
    survivor = outputs[('SUCCESS',)]
    posterior = survivor.conj().T@survivor/np.linalg.norm(survivor)**2
    assert np.linalg.norm(posterior-np.eye(2)/2) > .05
    complete = dict(a.complete_kraus(tape))
    assert len(complete) == 4+1+446
    assert np.linalg.norm(sum(k.conj().T@k for k in complete.values())-np.eye(66)) < 1e-13
    assert np.linalg.norm(complete[('SUCCESS',)].conj().T@complete[('SUCCESS',)]-np.eye(66)) > .1
    # The existing interpreter also retains a malformed-word rejection.
    rejected = interpreter.submit_word({}, np.eye(2)/2)
    assert rejected['status'] == 'REJECTED_PROGRAM'


def test_overlap_join_and_distinct_expectation_obligation(learned, arrays):
    union, overlap = tuple(sorted(set(R) | set(S))), tuple(sorted(set(R) & set(S)))
    qr, qs, qu, qi = (ind.target(learned, p) for p in (R, S, union, overlap))
    assert all(x*y == u*i for x, y, u, i in zip(qr, qs, qu, qi, strict=True))
    assert np.linalg.norm(arrays[R]@arrays[S]-arrays[S]@arrays[R], 'fro') < 1e-9
    f = tuple(F(x[7]) for x in ind.STATES)
    expect = lambda p, x: ind.expectation_diagonal(TRUE, p, x)
    assert expect(overlap, expect(R, f)) == expect(overlap, f)
    assert expect(R, expect(overlap, f)) == expect(overlap, f)
    assert expect(R, expect(S, f)) != expect(overlap, f)
    assert expect(R, expect(S, f)) != expect(S, expect(R, f))


def test_regional_native_covariance_does_not_imply_stationarity(learned):
    for ports in a.REGIONS:
        q = ind.target(learned, ports)
        for name in ('C', 'E'):
            for edge in modular_clock.law.SEAMS:
                for _, _, pairs in modular_clock.law.raw_maps(name, edge):
                    assert len({q[j]/q[i] for i, j in pairs}) <= 1
    assert ind.target(learned, R) != ind.target(TRUE, R)


def test_finite_prepared_closeness_is_not_exact_regional_invariance():
    # Faithful state arbitrarily close to the stationary reference. Perturb
    # one spectator context only; the same regional matrix unit now has two
    # modular frequencies. This blocks inference of exact invariance from an
    # error bound; it makes no unexecuted assertion about the actual marginal.
    q = list(modular_clock.law.reference(modular_clock.parameters()))
    i = next(j for j, x in enumerate(ind.STATES) if x[1] and x[3])
    j = next(j for j, x in enumerate(ind.STATES) if x[2] and x[3])
    eps = min(q)/1000000
    q[i] += eps
    q[j] -= eps
    assert all(x > 0 for x in q) and sum(q) == 1
    with pytest.raises(ArithmeticError, match='reference'):
        ind.algebra_certificate(TRUE, R, q)


@pytest.mark.parametrize('ports', [(7, 1), (1, 1), (True, 7), (1, 12)])
def test_invalid_region_is_rejected_before_execution(ports):
    assert a.prepare_request(TRUE, ports)['status'] == 'COMPILER_STOP'
