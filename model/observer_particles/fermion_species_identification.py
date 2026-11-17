#!/usr/bin/env python3
"""Exact finite source symmetry and conditional fermion-species discriminators.

No empirical masses or SM assignment enter the source band or lift readout.
The source-native E instrument is executed separately from the declared
two-port hopping/statistics comparison. Neither is substituted for the other.
"""
from __future__ import annotations
from evidence_contract import public_receipt
from collections import Counter
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PUBLIC = ROOT
sys.path.insert(0, str(PUBLIC))
sys.path.insert(0, str(PUBLIC / 'code/a5_closure'))
import family_band_attachment_certificate as bands
import noncentral_seam_reduction_certificate as lifts
import super_tannakian_matter_lift_certificate as matter
from model.observer import native
from model.observer.source_discrimination import circuit
from model.observer.source_scalar import compiler


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exact(value):
    return [str(value.a), str(value.b)]


def pointer_law(character, dimension):
    mean = character / bands.F5(dimension)
    plus = (bands.ONE + mean) / bands.F5(2)
    minus = bands.ONE - plus
    assert not (-plus).is_positive() and not (-minus).is_positive()
    return {'plus': exact(plus), 'minus': exact(minus)}


def decode_character(law, dimension):
    plus = bands.F5(*map(F, law['plus']))
    minus = bands.F5(*map(F, law['minus']))
    assert plus + minus == bands.ONE
    return bands.F5(dimension)*(plus-minus)


def permutation_order(g):
    h = tuple(range(len(g)))
    for k in range(1, 121):
        h = tuple(g[j] for j in h)
        if h == tuple(range(len(g))):
            return k
    raise AssertionError('unclosed permutation')


def source_band_readout():
    carrier, rotations, pin = bands.load_carrier()
    adjacency = bands.lift(bands.adjacency_int(carrier))
    projectors = bands.spectral_projectors(adjacency)
    resolution = bands.verify_spectral_resolution(adjacency, projectors)
    bands.verify_equivariance(projectors, rotations)
    orbits = bands.verify_pair_orbits(carrier, rotations)
    characters = {}
    records = []
    for name, projector in projectors.items():
        dimension = int(bands.mat_trace(projector).a)
        values = []
        for g in rotations:
            # tr(P U_g), independent of the choice of basis inside the band.
            chi = sum((projector[i][g[i]] for i in range(12)), bands.ZERO)
            law = pointer_law(chi, dimension)
            decoded = decode_character(law, dimension)
            assert decoded == chi
            values.append(decoded)
            records.append({'band': name, 'rotation': list(g), 'law': law})
        characters[name] = values
    gram = {}
    for left, x in characters.items():
        gram[left] = {}
        for right, y in characters.items():
            product = sum((a*b for a,b in zip(x,y)), bands.ZERO) / bands.F5(60)
            assert product == bands.F5(int(left == right))
            gram[left][right] = exact(product)
    # The class is distinguished by the actual character as well as order:
    # A5 has two order-five classes with different triplet characters.
    table = Counter((permutation_order(g), tuple(exact(characters['3'][k])))
                    for k,g in enumerate(rotations))
    classes = []
    for (order, chi3), size in sorted(table.items()):
        i = next(i for i,g in enumerate(rotations)
                 if permutation_order(g)==order and tuple(exact(characters['3'][i]))==chi3)
        classes.append({'order': order, 'size': size,
                        'characters': {name: exact(c[i]) for name,c in characters.items()}})
    # Three copies of an invariant scalar line can have the same one-pole
    # time response as the rank-three band, but rotations distinguish them.
    mismatch = next(i for i,g in enumerate(rotations) if permutation_order(g)==2)
    chi = characters['3'][mismatch]
    assert chi == bands.F5(-1)
    covariance_checks = sum(len(rotations) * 12 * 12 for _ in projectors)
    return {
        'source_carrier_pin': pin, 'resolution': resolution,
        'character_field': '[a,b] means a+b*sqrt(5)',
        'classes': classes, 'character_gram': gram,
        'pair_orbits': orbits, 'equivariance_entry_checks': covariance_checks,
        'readout_settings': records,
        'readout_scope': 'Exact controlled-rotation probabilities with prepared normalized band projector; no native preparation or rotation pointer word newly compiled.',
        'same_time_response_control': {
            'source_band_dimension': 3, 'single_pole': '5-sqrt(5)',
            'alternative': 'three invariant copies with that identical supplied single pole',
            'rank3_source_order2_pointer': pointer_law(chi, 3),
            'three_scalar_copies_order2_pointer': pointer_law(bands.F5(3), 3),
            'source_rotation_commutant_dimension': 1,
            'three_scalar_copies_commutant_dimension': 9,
            'conclusion': 'Time-pole multiplicity does not identify generations. The actual selected triplet is an irreducible spatial source-rotation representation; physical family attachment remains separate.'},
        'interpretation': 'Integer-spin-type finite A5 representations of port coefficients; no continuous spatial-spin or physical species attachment inferred.'
    }, rotations


def parse_q5(text):
    # Parse only the two exact rational formats emitted by the pinned artifact.
    if '*sqrt(5)' not in text:
        return (F(text), F(0))
    a, b = text.split(' + ')
    return (F(a), F(b.removesuffix('*sqrt(5)')))


def source_spin_lift(rotations):
    artifact = json.loads((PUBLIC / 'code/a5_closure/manifests/spin_statistics_semantic_artifact.json').read_text())
    generators, _ = lifts.pinned_lift_quaternions(artifact)
    group = lifts.build_binary_icosahedral(generators, artifact)
    frame = [tuple(parse_q5(v) for v in artifact['port_vertex_frame'][p])
             for p in artifact['carrier_binding']['port_order']]
    frame_index = {v: i for i,v in enumerate(frame)}
    projection = []
    for q in group.elements:
        conjugate = (q[0], *[(-a,-b) for a,b in q[1:]])
        out = []
        for v in frame:
            image = lifts.q5_quaternion_multiply(
                lifts.q5_quaternion_multiply(q, ((F(0),F(0)), *v)), conjugate)
            assert image[0] == (F(0),F(0))
            out.append(frame_index[image[1:]])
        projection.append(tuple(out))
    assert set(projection) == set(rotations)
    assert set(Counter(projection).values()) == {2}
    for i in range(120):
        for j in range(120):
            assert projection[group.table[i][j]] == tuple(projection[i][k] for k in projection[j])
    spin_characters = [bands.F5(2*q[0][0], 2*q[0][1]) for q in group.elements]
    norm = sum((x*x for x in spin_characters), bands.ZERO) / bands.F5(120)
    assert norm == bands.ONE
    z = lifts.minus_one_index(group)
    assert projection[z] == tuple(range(12)) and spin_characters[z] == bands.F5(-2)
    involutions = [i for i,g in enumerate(projection) if permutation_order(g)==2]
    assert len(involutions)==30 and all(group.table[i][i]==z for i in involutions)
    records = [{'quaternion': [[str(a),str(b)] for a,b in q],
                'deck_rotation': list(projection[i]), 'law': pointer_law(spin_characters[i],2)}
               for i,q in enumerate(group.elements)]
    return {
        'group_order': group.order, 'order_profile': group.order_profile(),
        'projection_products_checked': 14400, 'cover_fiber_size': 2,
        'spinor_character_norm_squared': exact(norm),
        'deck_involution_lifts_checked': len(involutions),
        'all_deck_involution_lifts_square_to_minus_one': True,
        'central_pointer': {'scalar_port_representation': pointer_law(bands.ONE,1),
                            'fundamental_lift_spinor': pointer_law(bands.F5(-2),2)},
        'readout_settings': records,
        'positive_result': 'The actual carrier admits a non-split binary-icosahedral transport cover and its irreducible two-dimensional unit-quaternion representation.',
        'scope': 'Central sign is a relative phase observable only with controlled reference access. This executes exact finite representation/readout algebra, not a physical matter spinor preparation, continuum spin-statistics theorem, or spin assignment to every source pole. Other half-integer-spin representations share central sign.'}


def grading_attachment_control():
    # Call the actual public implementer and exterior-power code: a concrete
    # test of what its central lift does, not a convention inferred from labels.
    identity3 = [[bands.F5(int(i==j)) for j in range(3)] for i in range(3)]
    minus_identity2 = [[matter.C5(bands.F5(-int(i==j)), bands.ZERO)
                        for j in range(2)] for i in range(2)]
    pi = matter.carrier_implementer(identity3, minus_identity2)
    fock = matter.FockSpace()
    lifted = fock.exterior_lift(pi)
    parity = fock.parity()
    rows = []
    for i,s in enumerate(fock.subsets):
        if len(s) not in (2,4):
            continue
        assert all(lifted[j][i].is_zero() for j in range(32) if j!=i)
        sign = int(lifted[i][i].re.a)
        assert parity[i][i].re == bands.ONE
        rows.append({'auxiliary_subset': list(s), 'weak_occupancy': sum(j>=3 for j in s),
                     'actual_lift_centre': sign, 'auxiliary_car_parity': 1})
    signs = Counter(row['actual_lift_centre'] for row in rows)
    assert signs == {1:7, -1:8}
    return {
        'dimension': len(rows), 'states': rows,
        'actual_centre_eigenvalue_multiplicities': dict(signs),
        'auxiliary_car_even_on_all15': True,
        'uniform_physical_fermion_parity_emitted_by_this_operator': False,
        'scope': 'The internal rank15 gauge module may consistently be typed as15 external Weyl fields. Its auxiliary even exterior degree and diag(R3,U2) transport are not themselves the common physical Spin centre. No contradiction to a separately attached Spin factor is claimed.'}


def native_exchange_control():
    rows = []
    maximum_error = 0.0
    for coin in (0,1):
        word = circuit.repair_word('E',coin)
        actual = compiler.execute(word,6,np.eye(64,dtype=complex)[:,:6])
        values = circuit.record_values('E',coin)
        inverse = {value:code for code,value in values.items()}
        for initial, local in enumerate(circuit.SEAM_STATES):
            load = local + ((2-sum(local)),) + (0,)*9
            output, (_, current) = native.native(load,(0,1),coin)
            record = local[0]-local[1]
            target = circuit.SEAM_STATES.index(output[:2])+8*inverse[record]
            expected = np.eye(64)[:,target]
            error = float(np.linalg.norm(actual[:,initial]-expected))
            maximum_error = max(maximum_error,error)
            law = abs(actual[:,initial])**2
            rows.append({'coin':coin,'coin_probability':'1/2','input_local_occupation':list(local),
                         'raw_law_64_outputs':law.tolist(),
                         'raw_record_precontrast':record,'current':current,
                         'target_output_address':target,'target_local_occupation':list(output[:2]),
                         'state_error_L2':error, 'raw_probability_sum':float(law.sum())})
    assert maximum_error < 3e-11
    exports = {}
    for d in (2,4):
        ks = compiler.transfer_kraus(d,True)
        completeness = sum(k.conj().T@k for k in ks)
        assert np.array_equal(completeness,np.eye(6))
        exports[str(d)] = {'all_branches':len(ks),'failure_branches':6-d,
                          'completeness_exact':True,'ideal_valid_input_failure_probability':0}
    return {'law':'E','scheduler':'one fixed seam, fair classical coin retained',
            'native_steps':1,'local_columns_executed':12,'source_pulses_per_coin':5373,
            'records':rows,'maximum_state_error_L2':maximum_error,
            'all_reserved_padding_record_outcomes_retained':True,'complete_export':exports,
            'result':'Each of20,11,02 repairs to11; identical coincidence1 does not distinguish CAR from this dissipative balancing rule.',
            'scope':'Actual source words under the existing protected-code/load/export apparatus. The fixed-seam repair is not the supplied unitary hopping factor used by the statistics discriminator.'}


def declared_exchange_readout():
    path = PUBLIC / 'code/spin_exchange/build.py'
    spec = importlib.util.spec_from_file_location('species_exchange_build',path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=module
    spec.loader.exec_module(module)
    receipt = module.build()
    # Independent tensor-space propagation, different implementation from
    # creation polynomials. All three occupation bins retained.
    u = np.array([[3,-4j],[-4j,3]],complex)/5
    unitary = np.kron(u,u)
    controls = {}
    for name, state in [('fermion',np.array([0,1,-1,0])/np.sqrt(2)),
                        ('boson',np.array([0,1,1,0])/np.sqrt(2)),
                        ('distinguishable',np.array([0,1,0,0]))]:
        out = unitary@state
        laws = {'20':abs(out[0])**2,'11':abs(out[1])**2+abs(out[2])**2,'02':abs(out[3])**2}
        errors = [abs(laws[k]-float(F(v))) for k,v in receipt['results'][name]['probabilities'].items()]
        assert max(errors)<1e-14
        controls[name]={'exact_probabilities':receipt['results'][name]['probabilities'],
                        'independent_tensor_law_error':float(max(errors))}
    return {'declared_splitter':receipt['splitter'], 'statistics_models':controls,
            'identical_one_body_density':receipt['results']['fermion']['one_body_density'],
            'same_one_particle_bilinears':True,
            'coincidence_gap':receipt['gap']['fermion_minus_boson_coincidence'],
            'partial_distinguishability':receipt['partial_distinguishability'],
            'scope':'Replayed ideal probabilities for a supplied identical-channel preparation and hopping action; no source selection of statistics, native exchange experiment, empirical outcome, or inherited quantum-link Gauss.'}


def main():
    band,rotations = source_band_readout()
    result = {'schema':'oph.finite.fermion_species_identification.v1',
              'status':'pass_scoped_source_symmetry_and_conditional_statistics_readouts',
              'public_head':'bb469bd80432fa24672bbb605afc1d171e46fa9d',
              'source_band_readout':band,'source_spin_lift':source_spin_lift(rotations),
              'grading_attachment_control':grading_attachment_control(),
              'native_exchange_control':native_exchange_control(),
              'declared_exchange_readout':declared_exchange_readout(),
              'SM_particle_identification_from_native_poles':False,
              'next_same_response_attachment':'A physical charged interpolator whose source time response, gauge representation and common Spin action are measured together, then an identical-channel two-excitation exchange readout. The existing conditional internal module and geometric spin cover supply concrete inputs, not that joint emission.'}
    dependencies = set()
    for module in tuple(sys.modules.values()):
        file = getattr(module,'__file__',None)
        if file and Path(file).is_file() and Path(file).is_relative_to(PUBLIC):
            dependencies.add(Path(file))
    for relative in ['code/a5_closure/manifests/echosahedral_federation_reference.json',
                     'code/a5_closure/manifests/spin_statistics_semantic_artifact.json',
                     'code/a5_closure/manifests/super_tannakian_matter_reference.json',
                     'code/a5_closure/manifests/family_band_attachment_reference.json',
                     'code/sm_fermion_current/quantum_link_receipt.json']:
        dependencies.add(PUBLIC/relative)
    for relative in compiler.source_summary()['source_sha256']:
        dependencies.add(PUBLIC/relative)
    result['source_sha256'] = {str(p.relative_to(ROOT)):sha(p) for p in sorted(dependencies)}
    result['producer_sha256']=sha(__file__)
    target=HERE/'results/fermion_species_identification.json'
    target.write_text(json.dumps(public_receipt(result),indent=2,sort_keys=True)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status':result['status'],'output':str(target),
                      'binary_icosahedral_order':result['source_spin_lift']['group_order'],
                      'native_state_error':result['native_exchange_control']['maximum_state_error_L2'],
                      'actual_internal_centre':result['grading_attachment_control']['actual_centre_eigenvalue_multiplicities']}))


if __name__=='__main__':
    main()
