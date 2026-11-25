from fractions import Fraction
import json
from math import factorial

import numpy as np

from ....receipt_replay import pinned_receipt, verify_recorded_band
from .. import station_region
from . import joint, total


def density_entropy(state):
    values = np.linalg.eigvalsh(state)
    assert np.all(values > 0)
    return -np.dot(values, np.log(values))


def terminal_blocks(flight, route, head):
    selected = np.take(flight, route, axis=0).reshape(16, 10, 4).transpose(1, 0, 2)
    remaining = np.delete(flight, route, axis=0).transpose(0,2,1,3).reshape(-1,4,4)

    def reduced(maps):
        return (maps @ head @ maps.conj().transpose(0,2,1)).sum(axis=0)/10

    return reduced(selected), reduced(remaining)


def test_closed_flight_and_retained_head_entropy():
    full = total.run()
    logged = pinned_receipt(total.HERE/"total.json", "87d126bf09616813e55dda706f631558261b3f1e6ab25ad31c43d24e6fb8cb5f")
    result = joint.run()
    archive = pinned_receipt(total.HERE/"joint.json", "80b1e5ec8805737b7ad158934affd895d37fa486b667d4c598e0dc67c1196d4a")
    verify_recorded_band(logged["source_slope"], full["source_slope"], significant_digits=7)
    verify_recorded_band(archive["reference_slope"], result["reference_slope"], significant_digits=7)
    assert full["control"] == result["control"] == "fixed_all_number_weights"
    assert result["retained_ports"] == [0, 1, 38, 39] and result["retained_heads"] == 2
    assert result["charge_erased"] and result["fair_charge_weights"] == ["1/2","1/2"]
    pi = station_region.parse_interval
    assert all(pi(x).a > lo for x,lo in zip(full["reference_conditional_eigenvalue_lower"],(.20254,.21036),strict=True))
    assert pi(full["tangent_HS_norm_lower"]).a > .9979
    assert pi(full["entropy_change_at_0p01_upper"]).b < -.00074047
    assert pi(result["entropy_change_at_0p01_upper"]).b < -.00016035
    assert pi(result["quadratic_entropy_curvature_lower"]).a > 3.2204
    assert pi(result["centered_tangent_error_gain"]).b < .366
    assert max(pi(x["maximum_log_error"]).b for x in result["conditionals"]) < 1.55e-9
    assert pi(result["reference_slope"]).a - 1.36457e-6 > 5.63e-6
    region = station_region
    source = region.measure
    board = source.bridge.native.iso.pairing.formation.federation(1)
    route = source.bridge.native.iso.pairing.formation.region(board)
    rings = source.bridge.beam.screens(board, route)
    assert len(set(sum(rings, []))) == 10
    start = np.zeros((4*board.n,40))
    start[source.bridge.beam.indices(sum(rings, []))] = np.eye(40)
    reference, _, _ = source.source()
    tangent, _ = region.controls()["fixed_all_number_weights"]
    full_slope = full_change = joint_slope = joint_change = 0
    for odd, weight in enumerate((.8657,.1343)):
        rho = np.array(reference[odd], complex)*weight
        delta = np.array(tangent[odd], complex)
        spectrum, vectors = np.linalg.eigh(rho)
        full_slope -= np.dot(np.log(spectrum), np.diag(vectors.conj().T@delta@vectors).real)
        full_change += density_entropy(rho+.01*delta)-density_entropy(rho)
        homes, changes = [], []
        for q in (-1, 1):
            h = source.bridge.native.generator(board, route, q, bool(odd))
            # Independent small-step Taylor solver, with an analytic local
            # remainder and a direct Gram check on all forty launch columns.
            assert abs(h).sum(axis=1).max() <= 43
            remainder = Fraction(43,32)**23/factorial(23)
            assert (1+remainder)**64-1 < Fraction(1,10**16)
            op = -1j*h/32
            flight = start.astype(complex)
            for _ in range(64):
                term = flight
                updated = flight.copy()
                for order in range(1, 23):
                    term = (op@term)/order
                    updated += term
                flight = updated
            np.testing.assert_allclose(flight.conj().T@flight, np.eye(40), rtol=0, atol=2e-13)
            flight = flight.reshape(board.n, 4, 10, 4)
            homes.append(terminal_blocks(flight, route, rho))
            changes.append(terminal_blocks(flight, route, delta))

        def pooled(values):
            return [values[0][0]/2,values[1][0]/2,(values[0][1]+values[1][1])/2]

        for block, diff in zip(pooled(homes),pooled(changes),strict=True):
            ev, v = np.linalg.eigh(block)
            joint_slope -= np.dot(np.log(ev), np.diag(v.conj().T@diff@v).real)
            joint_change += density_entropy(block+.01*diff)-density_entropy(block)
    assert abs(joint_slope-float(pi(result["nominal_slope"]).mid)) < 1e-10
    assert abs(joint_slope-float(pi(archive["nominal_slope"]).mid)) < 1e-10
    assert full_slope in pi(full["source_slope"])
    assert full_change < float(pi(full["entropy_change_at_0p01_upper"]).a)
    assert joint_slope in pi(result["reference_slope"])
    assert joint_change < float(pi(result["entropy_change_at_0p01_upper"]).a)


def test_conditional_memory_gain():
    region = station_region
    joint_report = json.loads((total.HERE/"joint.json").read_text())
    loads = json.loads((region.measure.HERE/"station_region.json").read_text())
    only = loads["regions"]["stations"]["fair"]["fixed_all_number_weights"]
    pi = region.parse_interval
    slope = pi(joint_report["reference_slope"]) - pi(only["slope"])
    assert slope.a > 5.63e-6
    delta_upper = pi(joint_report["entropy_change_at_0p01_upper"]) - pi(only["finite_change_0_to_0p01"])
    assert delta_upper.b < -1.6030e-4
