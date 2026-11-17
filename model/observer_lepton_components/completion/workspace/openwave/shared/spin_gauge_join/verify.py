#!/usr/bin/env python3
"""Exact same-field carrier join and a complete finite current-path reader.

Run from any directory. No upstream/status files are written.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import replace
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RER = ROOT / 'reverse-engineering-reality'
sys.path.insert(0, str(RER))
from model.native_electric_gauss import source as g
from model.defect_field_bridge import full_domain as car
from model.common_source import slow_geometry

TRIANGLE = (0, 5, 11)
EDGES = ((0, 5), (5, 11), (0, 11))
SPINS = tuple(product((0, 1), repeat=2))
SIGMA = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], dtype=complex)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encoded(s, spin):
    """Positive-tag creation first, then negative compensator creation."""
    other = next(x for x in s.pair if x != s.tagged_port)
    mask = (1 << (2*s.tagged_port+spin[0])) | (1 << (2*other+spin[1]))
    return mask, 1 if s.tagged_port < other else -1


def car_swap(mask, seam, selector):
    a, b = seam
    na, nb = (((mask >> (2*x)) & 3).bit_count() for x in seam)
    # A coherent identity/swap selector is not the original endpoint coin.
    coin = 1-selector if na == 1 and nb == 0 else selector if na == 0 and nb == 1 else 0
    _, target, sign, current = car.branch_action(
        mask, 'C', coin, sites=12, flavors=2, seam=seam)
    return target, sign, current


def joined(s, spin, seam, selector=1, limit=1):
    """Actual existing Gauss action plus the field's canonical CAR sign."""
    hop = g.apply_native_current(s, seam, selector, limit)
    mask, _ = encoded(s, spin)
    target, sign, current = car_swap(mask, seam, selector)
    if hop.overflow:
        target, sign, current = mask, 1, 0
    assert target == encoded(hop.state, spin)[0]
    assert current == (0 if hop.overflow else hop.proposed_number_current)
    return hop, sign, target


def state_key(s, spin):
    return s.pair, s.tagged_port, spin, s.charge_number, s.memory, s.flux


def local_audit():
    count = 0
    coin_columns = 0
    negative = 0
    overflows = 0
    controls = {'omit_car_sign': 0, 'omit_flux_update': 0,
                'use_proposed_current_at_overflow': 0, 'flip_transported_spin': 0}
    witnesses = {}
    branch_counts = []
    for seam in g.SEAMS:
        e = g.EDGE_INDEX[seam]
        memories = range(6) if seam == (0, 5) else (0,)
        coherent_images = [set(), set()]
        raw_images = {(law, coin, raw, overflow): set()
                      for law in ('C', 'E') for coin in (0, 1)
                      for raw in range(16) for overflow in (False, True)}
        for q, pair, tag_index, spin, memory, value in product(
                (-1, 0, 1), g.PAIRS, (0, 1), SPINS, memories, (-1, 0, 1)):
            flux = [0]*30
            flux[e] = value
            s = g.ElectricState(pair, pair[tag_index], q, memory, tuple(flux))
            mask, epsilon = encoded(s, spin)
            hops = {}
            for selector in (0, 1):
                hop, sign, target = joined(s, spin, seam, selector)
                hops[selector] = hop, sign, target
                t = hop.state
                assert epsilon*sign == encoded(t, spin)[1]
                inverse, inverse_sign, _ = joined(t, spin, seam, selector)
                assert inverse.state == s and sign*inverse_sign == 1
                out = state_key(t, spin)
                assert out not in coherent_images[selector]
                coherent_images[selector].add(out)
                assert g.gauss(t) == g.gauss(s)
                assert g.electric_energy(t)-g.electric_energy(s) == g.KAPPA*(
                    F(hop.realized_signed_current**2, 2)-value*hop.realized_signed_current)
                assert spin == out[2]  # independent noncommuting spin fiber is transported intact
                negative += sign < 0
                overflows += hop.overflow
                count += 1
                if sign < 0:
                    controls['omit_car_sign'] += 1
                    witnesses.setdefault('omit_car_sign', {'pair': pair, 'tag': s.tagged_port,
                        'spin': spin, 'seam': seam, 'q': q, 'memory': memory, 'local_flux': value})
                if t.pair != s.pair and q:
                    assert g.gauss(replace(t, flux=s.flux)) != g.gauss(s)
                    controls['omit_flux_update'] += 1
                    moved_tag = int(s.tagged_port not in seam)
                    changed_spin = list(spin)
                    changed_spin[moved_tag] = 1-changed_spin[moved_tag]
                    assert abs(sum(1-2*x for x in changed_spin)-sum(1-2*x for x in spin)) == 2
                    controls['flip_transported_spin'] += 1
                if hop.overflow:
                    assert hop.proposed_signed_current != hop.realized_signed_current
                    controls['use_proposed_current_at_overflow'] += 1
                    witnesses.setdefault('overflow', {'pair': pair, 'tag': s.tagged_port,
                        'seam': seam, 'q': q, 'local_flux': value,
                        'proposed': hop.proposed_signed_current, 'realized': hop.realized_signed_current})
            for law in ('C', 'E'):
                total = F(0)
                for coin in (0, 1):
                    raw, proposed_mask, proposed_sign, current = car.branch_action(
                        mask, law, coin, sites=12, flavors=2, seam=seam)
                    hop, sign, target = hops[abs(current)]
                    if not hop.overflow:
                        assert target == proposed_mask and sign == proposed_sign
                    key = law, coin, raw, hop.overflow
                    out = state_key(hop.state, spin)
                    assert out not in raw_images[key]
                    raw_images[key].add(out)
                    total += F(1, 2)
                    coin_columns += 1
                assert total == 1
        branch_counts.append({'seam': list(seam), 'memories_enumerated': list(memories),
            'complete_raw_slots_per_law': 64,
            'nonempty_raw_slots': {law: sum(bool(rows) for key, rows in raw_images.items()
                                         if key[0] == law) for law in ('C', 'E')},
            'coherent_columns_per_selector': len(coherent_images[0])})
    assert all(controls.values())
    return {'coherent_basis_columns_checked': count,
        'native_coin_basis_columns_checked': coin_columns,
        'negative_car_sign_columns': negative, 'overflow_columns': overflows,
        'all_exact_intertwining_gauss_spin_involution_checks_pass': True,
        'all_raw_effect_sums_exactly_identity': True,
        'all_branches_partial_injections': True,
        'controls_failing_columns': controls, 'control_witnesses': witnesses,
        'per_seam': branch_counts,
        'memory_extension': 'All six values executed on seam(0,5). On every seam the shift is independent of input memory, so equality on one value tensors the exact Z6 translation on all six.',
        'complement': 'Flagged identity on the nonadmitted full Fock/tag sector; no outcome discarded. Orthogonal admission projector commutes with total field spin.',
        'reference_scope': 'Exact signed-column intertwiners fix every matrix unit and arbitrary untouched reference, not only diagonal states.'}


def spin_loop(vector, spin, inverse=False, limit=2):
    out = {}
    edges = tuple(reversed(EDGES)) if inverse else EDGES
    for s, amplitude in vector.items():
        sign = 1
        t = s
        for seam in edges:
            hop, factor, _ = joined(t, spin, seam, limit=limit)
            sign *= factor
            t = hop.state
        g.add_state(out, t, g.scale(sign, amplitude))
    return out


def spin_reference(vector, spin, phase, limit=2):
    moved = spin_loop(vector, spin, limit=limit)
    outputs = []
    for sign in (1, -1):
        out = {s: g.scale(F(1, 2), amplitude) for s, amplitude in vector.items()}
        for s, amplitude in moved.items():
            g.add_state(out, s, g.scale(F(sign, 2), g.mul(g.ROOTS[phase], amplitude)))
        outputs.append(out)
    assert sum(g.probability(out) for out in outputs) == g.probability(vector)
    return outputs


def loop_audit():
    source_receipt = g.loop_experiment()
    records = []
    boundary_checks = 0
    for spin in SPINS:
        initial = g.character_preparation()
        first = spin_reference(initial, spin, 1)
        assert first == list(g.loop_reference_measurement(initial, 1, 2))
        for phase in (0, 1):
            for a, first_vector in enumerate(first):
                for b, vector in enumerate(spin_reference(first_vector, spin, phase)):
                    assert all(not any(g.gauss(s)) for s in vector)
                    records.append({'spin_at_positive_negative': spin,
                        'second_phase_sixths': phase, 'first_pointer': a, 'second_pointer': b,
                        'probability': str(g.probability(vector)),
                        'joint_electric_energy': str(g.energy_mean(vector)),
                        'terminal_amplitudes_Q_zeta_over_sqrt6': [
                            {'car_mask': encoded(s, spin)[0], 'tag': s.tagged_port,
                             'memory': s.memory, 'flux': list(s.flux),
                             'amplitude': list(map(str, amp))}
                            for s, amp in sorted(vector.items(), key=lambda row:(row[0].flux,row[0].memory))]})
        boundary = initial
        for _ in range(4):
            moved = spin_loop(boundary, spin, limit=1)
            assert spin_loop(moved, spin, inverse=True, limit=1) == boundary
            boundary = moved
            boundary_checks += 1
    return {'all_four_spin_inputs_exact': True,
        'arbitrary_spin_coherence': 'The spin fiber is the identity in each path amplitude; equality therefore holds on all sixteen spin matrix units and arbitrary references.',
        'complete_joint_pointer_records': records,
        'boundary_forward_inverse_checks': boundary_checks,
        'pointer_law_phase0': ['5/16', '3/16', '3/16', '5/16'],
        'pointer_law_phase1': ['3/8', '1/8', '1/8', '3/8'],
        'nonselective_electric_energy': [source_receipt[k] for k in
            ('initial_electric_energy', 'after_first_nonselective_energy', 'after_second_nonselective_energy')],
        'energy_status': 'Existing supplied kappa=1/8. No new magnetic generator, energy normalization or physical time.'}


def finite_slice(q):
    states = []
    for position, local_flux, memory in product(TRIANGLE, product((-1, 0, 1), repeat=3), range(6)):
        flux = list(g.pair_string((0,1), q, memory).flux)
        for seam, value in zip(EDGES, local_flux, strict=True):
            flux[g.EDGE_INDEX[seam]] = value
        states.append(g.ElectricState(tuple(sorted((position, 1))), position, q, memory, tuple(flux)))
    return states


def path_audit(q):
    states = finite_slice(q)
    index = {s:i for i,s in enumerate(states)}
    n = len(states)
    targets = np.zeros((3,n), dtype=np.int64)
    currents = np.zeros((3,n), dtype=np.int64)
    signs = np.zeros((3,n), dtype=np.int64)
    overflow = np.zeros((3,n), dtype=bool)
    operators = []
    for e, seam in enumerate(EDGES):
        matrix = np.zeros((n,n), dtype=complex)
        for a,s in enumerate(states):
            hop, sign, _ = joined(s,(0,1),seam,limit=1)
            b = index[hop.state]
            current = 0 if hop.overflow else hop.proposed_number_current
            targets[e,a], currents[e,a], signs[e,a], overflow[e,a] = b,current,sign,hop.overflow
            matrix[b,a] = 1j*current*sign
        assert np.array_equal(matrix,matrix.conj().T)
        assert len(set(targets[e])) == n
        operators.append(matrix)
    frame = np.array([[float(a)+float(b)*np.sqrt(5) for a,b in row]
                      for row in slow_geometry.raw_frame()])
    frame /= np.linalg.norm(frame[0])
    steps = np.array([frame[b]-frame[a] for a,b in EDGES])
    x = np.einsum('ei,eab->iab',steps,np.array(operators))
    d = sum(np.kron(SIGMA[i], x[i]) for i in range(3))
    scalar = sum(a@a for a in x)
    # INDEPENDENT path enumeration. Each entry records two actual calls of
    # the native coherent swap/current operation. No [X_i,X_j] is used here.
    area = np.zeros((3,n,n),dtype=complex)
    scalar_paths = np.zeros((n,n),dtype=complex)
    path_records = []
    failures = []
    for first in range(3):
        for second in range(3):
            for a in range(n):
                middle = targets[first,a]
                b = targets[second,middle]
                j1,j2 = currents[first,a],currents[second,middle]
                factor = -int(j1*j2*signs[first,a]*signs[second,middle])
                # Ordered matrix product is J_second J_first.
                scalar_paths[b,a] += np.dot(steps[second],steps[first])*factor
                area[:,b,a] += np.cross(steps[second],steps[first])*factor
                path_records.append((first,second,a,middle,b,j1,j2,factor))
                failures.append((j1==0,j1!=0 and j2==0,j1!=0 and j2!=0))
    cross = 1j*sum(np.kron(SIGMA[i],area[i]) for i in range(3))
    errors = {'scalar_two_path_reconstruction': float(np.max(np.abs(scalar-scalar_paths))),
        'clifford_square_from_independent_path_reader': float(np.max(np.abs(d@d-np.kron(np.eye(2),scalar_paths)-cross))),
        'path_cross_hermiticity': float(np.max(np.abs(cross-cross.conj().T)))}
    assert max(errors.values()) < 1e-10
    # A complete unitary completion of each current contraction. Moving
    # domain gives J; hold/overflow routes to a retained second output.
    completion_checks = 0
    for current_operator in operators:
        effect = current_operator.conj().T@current_operator
        assert np.array_equal(effect,np.diag(np.diag(effect)))
        residual = np.eye(n)-effect
        block = np.block([[current_operator,residual],[residual,-current_operator]])
        assert np.array_equal(block.conj().T@block,np.eye(2*n))
        completion_checks += 1
    # Pair the complete current-query words with a reference. Neither the
    # path zeros nor the failure outputs are conditioned out. The path cross
    # entries below are the independently measured signed amplitude sum.
    word_pointer = []
    basis_input = index[g.pair_string((0,1),q,0)]
    for e in range(3):
        for f in range(e+1,3):
            v = np.zeros(n,complex);v[basis_input]=1
            a = operators[f]@operators[e]@v
            b = operators[e]@operators[f]@v
            for quarter_phase in (0,1):
                phase = (1j)**quarter_phase
                plus = np.vdot((a+phase*b)/2,(a+phase*b)/2).real
                minus = np.vdot((a-phase*b)/2,(a-phase*b)/2).real
                failure = 1-(np.vdot(a,a).real+np.vdot(b,b).real)/2
                assert min(plus,minus,failure) >= -1e-14
                assert abs(plus+minus+failure-1)<1e-14
                word_pointer.append({'first_seam':e,'second_seam':f,
                    'quarter_phase':quarter_phase,'plus':float(plus),
                    'minus':float(minus),'all_query_failure':float(failure)})
    arr = {'targets':targets,'realized_number_currents':currents,
           'car_signs':signs,'overflow':overflow,
           'path_records':np.array(path_records,dtype=np.int64),
           'complete_path_failure_stages':np.array(failures,dtype=bool),
           'frame_steps':steps,
           'path_area_entries':area}
    return {'charge_number':q,'source_dimension':n,'spin_current_dimension':2*n,
        'spectator_spin': 'Independent full two-state factor tensors identity; all four spin labels checked in carrier/loop audits.',
        'basis': [{'pair':list(s.pair),'tag':s.tagged_port,'memory':s.memory,
                  'local_flux':[s.flux[g.EDGE_INDEX[e]] for e in EDGES]} for s in states],
        'ordered_two_step_records':len(path_records),
        'query_completion_unitarity_checks':completion_checks,
        'all_ordered_path_failure_stages_complete':bool(np.all(np.sum(failures,axis=1)==1)),
        'moving_columns_per_seam':[int(np.count_nonzero(c)) for c in currents],
        'overflow_columns_per_seam':[int(np.count_nonzero(c)) for c in overflow],
        'numeric_errors':errors,
        'path_cross_frobenius_norm':float(np.linalg.norm(cross)),
        'path_cross_max_entry':float(np.max(np.abs(cross))),
        'complete_path_pointer_control':word_pointer,
        'kinetic_status':'D is a declared coherent spin-current reader. No source clock makes D or D squared an energy in this packet.',
        'curvature_status':'Oriented two-step path amplitude, including native geometry, memory and boundary. Not yet a Maxwell B field; q=0 control is retained.'}, arr


CITED = [
    'Lean/EventAlgebra/PauliMinimalCoupling.lean',
    'Lean/Screen/SeamCurrentCarrierQuotient.lean',
    'Lean/Screen/SeamCurrentHomogeneousAction.lean',
    'Lean/Screen/SeamCurrentDirichletGenerator.lean',
    'Lean/Screen/SeamCurrentPhysicalMetricAttachment.lean',
    'model/native_electric_gauss/PROOF.md',
    'model/native_observer_program/PROOF.md',
    'model/native_circulation/README.md',
    'model/observer/holonomy_reconstruction.py',
    'code/a5_closure/sl2f5_port_spin_bridge_certificate.py',
]


def pins():
    paths = set()
    for module in tuple(sys.modules.values()):
        filename = getattr(module,'__file__',None)
        if filename:
            path = Path(filename).resolve()
            if path.is_relative_to(RER) and path.is_file(): paths.add(path)
    for name in CITED: paths.add(RER/name)
    paths.update((RER/'evidence/native_geometric_source_20260925/spec.json',
                  HERE/'DESIGN.md',HERE/'TASK.md',HERE/'verify.py',HERE/'preregistration.json',
                  HERE.parent/'magnetic_source_join/MANIFEST.json',
                  HERE.parent/'magnetic_source_join/PROOF.md',
                  HERE.parent/'magnetic_source_join/NEXT_INTERFACE.md'))
    return {str(p.relative_to(ROOT)):digest(p) for p in sorted(paths)}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true')
    args=parser.parse_args()
    prereg=json.loads((HERE/'preregistration.json').read_text())
    assert all(digest(HERE/name)==value for name,value in prereg['files'].items())
    slow_geometry.registered_binding()
    before=pins()
    started=time.time()
    local=local_audit();print('Exact carrier and raw instruments passed.',flush=True)
    loop=loop_audit();print('Exact same-field coherent loop passed.',flush=True)
    responses=[];arrays={}
    for q in (0,1):
        response,arr=path_audit(q);responses.append(response)
        arrays.update({f'q{q}_{key}':value for key,value in arr.items()})
        print(f'Complete q={q} current/path reader passed.',flush=True)
    after=pins();assert before==after, 'Source inputs changed during run.'
    result={'schema':'oph.spin_gauss_join.v1','local_exact':local,'same_field_loop':loop,
        'finite_current_path_readers':responses,
        'nature_measurements_used':[],
        'new_native_M6_pulse_execution':False,
        'stationary_sourced_dipole_g_larmor': 'Not established by this packet.',
        'preregistration_sha256':digest(HERE/'preregistration.json')}
    if args.record:
        (HERE/'results.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
        (HERE/'SOURCE_PINS.json').write_text(json.dumps(after,sort_keys=True,indent=2)+'\n')
        np.savez_compressed(HERE/'path_records.npz',**arrays)
    else:
        assert result==json.loads((HERE/'results.json').read_text())
        assert after==json.loads((HERE/'SOURCE_PINS.json').read_text())
        old=np.load(HERE/'path_records.npz')
        assert set(old.files)==set(arrays)
        for key,value in arrays.items():assert np.array_equal(value,old[key]),key
    print(json.dumps({'all_pass':True,'seconds':round(time.time()-started,3),
        'source_inputs_pinned':len(after),'coherent_columns':local['coherent_basis_columns_checked'],
        'coin_columns':local['native_coin_basis_columns_checked'],
        'path_cross_norms':[r['path_cross_frobenius_norm'] for r in responses]},indent=2))


if __name__=='__main__':main()
