"""Separate full-pair geometry, instrument and budget verification."""
from fractions import Fraction as F
import hashlib, json, sys, time, argparse, math
from pathlib import Path
import numpy as np
sys.set_int_max_str_digits(0)
HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / '2026-11-16_codex_aligned_pair_spatial'
SEAL = 'db48ac7ee3c80f75a63ea991907f984c54135282316a8ef023878fd014535ba8'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read(p):
    return strict_json(p)

def tv(x, y):
    return float(np.abs(x - y).sum() / 2)

def exact_R(d):
    a, b, c = (int(d[:, 0] @ d[:, 0]), int(d[:, 0] @ d[:, 1]), int(d[:, 1] @ d[:, 1]))
    det = a * c - b * b
    return [[F(c * int(d[i, 0]) * int(d[j, 0]) - b * (int(d[i, 0]) * int(d[j, 1]) + int(d[i, 1]) * int(d[j, 0])) + a * int(d[i, 1]) * int(d[j, 1]), det) for j in range(12)] for i in range(12)]

def exact_points(rotate_second):
    pairs = [(-1, 0), (1, 0), (0, 0), (0, 1), (0, -1)]
    minus, plus, zero, phi, nphi = pairs
    local = [(minus, phi, zero), (plus, phi, zero), (minus, nphi, zero), (plus, nphi, zero), (zero, minus, phi), (zero, plus, phi), (zero, minus, nphi), (zero, plus, nphi), (phi, zero, minus), (phi, zero, plus), (nphi, zero, minus), (nphi, zero, plus)]
    out = []
    for cell in range(20):
        for a, b, c in local:
            if rotate_second and cell == 1:
                b, c = (tuple((-x for x in c)), b)
            out.append((a, (b[0] + 10 * cell, b[1]), c))
    return out

def permutation(points):
    lookup = {q: i for i, q in enumerate(points)}
    if not len(lookup) == 240:
        raise ValueError('Required check at ancestor line 32: len(lookup) == 240')
    out = [lookup[tuple((-v for v in a)), b, tuple((-v for v in c))] for a, b, c in points]
    if not all((out[out[k]] == k for k in range(240))):
        raise ValueError('Required check at ancestor line 34: all((out[out[k]] == k for k in range(240)))')
    return out + list(range(240, 256))

def superoperator_pump(Q, seed, dephase=False):
    escape = np.eye(12) - Q.conj().T @ Q
    M = np.kron(Q, Q.conj()) + np.outer(np.eye(12).reshape(144) / 12, escape.T.reshape(144))
    if dephase:
        M *= np.eye(12).reshape(144)[:, None]
    v = np.zeros(144, complex)
    v[13 * seed] = 1
    return (np.linalg.matrix_power(M, 512) @ v).reshape(12, 12)

def local_flag_state(columns, rho, E):
    A = columns[:, :, :12].reshape(512, 12)
    block = (A @ rho @ A.conj().T).reshape(2, 256, 2, 256)
    return np.einsum('ab,biaj->ij', E, block)

def all_addresses(rho, effects):
    result = []
    for word in range(256):
        E = np.ones((1, 1), complex)
        for bit in range(7, -1, -1):
            E = np.kron(E, effects[word >> bit & 1])
        result.append(np.sum(E.T * rho).real)
    return np.asarray(result)

def reconstruct():
    start = time.perf_counter()
    I = SOURCE / 'inputs'
    mfile = SOURCE / 'MANIFEST.sha256'
    if not sha(mfile) == SEAL:
        raise ValueError('Required check at ancestor line 56: sha(mfile) == SEAL')
    count = 0
    for line in mfile.read_text().splitlines():
        digest, rel = line.split('  ', 1)
        if not sha(SOURCE / rel) == digest:
            raise ValueError('Required check at ancestor line 59: sha(SOURCE / rel) == digest')
        count += 1
    receipt = read(SOURCE / 'RECEIPT.json')
    saved = np.load(SOURCE / 'PAIR_MAPS.npz')
    require(receipt['all_gates_pass'] is True and len(receipt['gates']) == 12 and all((v is True for v in receipt['gates'].values())), 'complete saved verdict')
    common = read(I / 'common_RECEIPT.json')
    native = np.load(I / 'common_MAPS.npz')
    first = np.load(I / 'first_maps.npz')
    second = np.load(I / 'second_maps.npz')
    graph = read(I / 'graph.json')['seams']
    L = np.zeros((240, 240), int)
    for a, b in graph:
        L[a, a] += 1
        L[b, b] += 1
        L[a, b] -= 1
        L[b, a] -= 1
    basis = [first['dark_integer_basis'], second['dark_integer_basis']]
    if not all((np.array_equal(L @ d, 6 * d) for d in basis)):
        raise ValueError('Required check at ancestor line 67: all((np.array_equal(L @ d, 6 * d) for d in basis))')
    r = [exact_R(basis[0][:12]), exact_R(basis[1][12:24])]
    points = exact_points(True)
    p = permutation(points)
    oldp = permutation(exact_points(False))
    if not receipt['coords_exact_a_plus_b_phi'] == json.loads(json.dumps(points)):
        raise ValueError("Required check at ancestor line 70: receipt['coords_exact_a_plus_b_phi'] == json.loads(json.dumps(points))")
    if not np.array_equal(p, saved['global_Ry_address_permutation']):
        raise ValueError("Required check at ancestor line 71: np.array_equal(p, saved['global_Ry_address_permutation'])")
    p0 = p[:12]
    p1 = [v - 12 for v in p[12:24]]
    for a in range(12):
        for b in range(12):
            for c in range(12):
                for d in range(12):
                    if not r[0][a][c] * r[1][b][d] == r[0][p0[a]][p0[c]] * r[1][p1[b]][p1[d]]:
                        raise ValueError('Required check at ancestor line 77: r[0][a][c] * r[1][b][d] == r[0][p0[a]][p0[c]] * r[1][p1[b]][p1[d]]')
    oldtv = sum((abs(r[0][a][a] * r[1][b][b] - r[0][oldp[a]][oldp[a]] * r[1][oldp[12 + b] - 12][oldp[12 + b] - 12]) for a in range(12) for b in range(12))) / 8
    if not oldtv == F(1, 3) == F(receipt['unaligned_joint_ideal_TV']):
        raise ValueError("Required check at ancestor line 79: oldtv == F(1, 3) == F(receipt['unaligned_joint_ideal_TV'])")
    if not (p[0] == 1 and p[1] == 0 and (np.linalg.det(np.array(receipt['global_rotation_matrix'])) == 1)):
        raise ValueError("Required check at ancestor line 80: p[0] == 1 and p[1] == 0 and (np.linalg.det(np.array(receipt['global_rotation_matrix'])) == 1)")
    rho0 = superoperator_pump(first['certified_Q_proposal'], 0)
    rho1 = superoperator_pump(second['certified_Q_proposal'], 2)
    control = superoperator_pump(first['certified_Q_proposal'], 0, True)
    state_error = max(np.max(abs(rho0 - native['first_final_density'])), np.max(abs(rho1 - native['second_final_density'])), np.max(abs(control - native['dephasing_final_density'])))
    if not state_error < 1e-12:
        raise ValueError('Required check at ancestor line 85: state_error < 1e-12')
    if not np.max(abs(np.kron(rho0, rho1) - saved['prepared_full144_density'])) < 1e-12:
        raise ValueError("Required check at ancestor line 86: np.max(abs(np.kron(rho0, rho1) - saved['prepared_full144_density'])) < 1e-12")
    pointer = np.load(I / 'pointer_witness.npz')['pointer_effects_C']

    def pure(k):
        a = np.zeros((12, 12), complex)
        a[k, k] = 1
        return a
    require(saved['selected_full144_density'].shape == (144, 144), 'selected full density shape')
    require(np.max(abs(saved['selected_full144_density'] - np.kron(np.array(r[0], float) / 2, np.array(r[1], float) / 2))) < 1e-14, 'selected full density')
    cases = {'native': (rho0, rho1), 'dephasing': (control, rho1), 'polar_plant': (pure(0), pure(2)), 'rotated_polar_plant': (pure(p0[0]), pure(p1[2]))}
    law_error = 0
    recomputed = {}
    for name, states in cases.items():
        require(np.max(abs(saved[name + '_full144_density'] - np.kron(*states))) < 1e-12, 'full prepared case ' + name)
        local = []
        flags = []
        dirty = []
        scratch = []
        for i, state in enumerate(states):
            blocks = native[f'loader{i}_clean_columns']
            outs = [local_flag_state(blocks, state, E) for E in pointer]
            flags.append([float(np.trace(z).real) for z in outs])
            local.append(all_addresses(outs[1], pointer))
            dirty.append(float(np.trace(native[f'loader{i}_dirty_gram'][:12, :12] @ state).real))
            scratch.append([float(np.trace(g[:12, :12] @ state).real) for g in native[f'loader{i}_scratch_word_grams']])
        full = np.multiply.outer(*local)
        reference = saved[name + '_all65536_position_centres']
        err = float(np.max(abs(full - reference)))
        law_error = max(law_error, err)
        if not err < 1e-12:
            raise ValueError('Required check at ancestor line 106: err < 1e-12')
        row = receipt['laws'][name]
        rejects = [flags[0][a] * flags[1][b] for a, b in [(0, 0), (0, 1), (1, 0)]]
        require(abs(dirty[0] - row['dirty_first_centre']) < 1e-12 and abs(sum(flags[0]) * dirty[1] - row['dirty_second_after_first_clean_centre']) < 1e-12, 'dirty complement ' + name)
        require(abs(float(full.sum()) - row['complete_position_mass_centre']) < 1e-12, 'position mass ' + name)
        require(abs(float(np.trace(full[:240, :240])) - row['collision_address_centre']) < 1e-12, 'collision ' + name)
        require(abs(float(full.sum() - full[:240, :240].sum()) - row['invalid_address_centre']) < 1e-12, 'invalid ' + name)
        residual = 1 - full.sum() - sum(rejects) - dirty[0] - sum(flags[0]) * dirty[1]
        if not abs(residual - row['unresolved_processor_plus_rounding_residual']) < 1e-12:
            raise ValueError("Required check at ancestor line 114: abs(residual - row['unresolved_processor_plus_rounding_residual']) < 1e-12")
        if not np.max(abs(np.asarray(scratch) - row['scratch_words_comparison_centres'])) < 1e-12:
            raise ValueError("Required check at ancestor line 115: np.max(abs(np.asarray(scratch) - row['scratch_words_comparison_centres'])) < 1e-12")
        if not np.max(abs(np.asarray(rejects) - row['region_reject_centres'])) < 1e-12:
            raise ValueError("Required check at ancestor line 116: np.max(abs(np.asarray(rejects) - row['region_reject_centres'])) < 1e-12")
        fulltv = tv(full, full[np.ix_(p, p)])
        marg = full.sum(axis=1)
        if not (abs(fulltv - row['joint_halfturn_TV_centre']) < 1e-12 and fulltv + 1e-12 >= tv(marg, marg[p])):
            raise ValueError("Required check at ancestor line 118: abs(fulltv - row['joint_halfturn_TV_centre']) < 1e-12 and fulltv + 1e-12 >= tv(marg, marg[p])")
        recomputed[name] = fulltv
    unequal = 0.75 * saved['polar_plant_all65536_position_centres'] + 0.25 * saved['rotated_polar_plant_all65536_position_centres']
    require(np.max(abs(unequal - saved['unequal_plant_all65536_position_centres'])) < 1e-12, 'unequal law')
    require(recomputed['polar_plant'] > 0.99 and tv(unequal, unequal[np.ix_(p, p)]) > 0.49, 'directional reader sensitivity')
    bit = F(35, 2 ** 31) + F(200000000001, 20000000000000000000000)
    receipts = [read(I / 'first_RECEIPT.json')['scientific'], read(I / 'second_RECEIPT.json')['scientific']]
    flow = F(read(I / 'SHARP_FLOW_CERTIFICATE.json')['512_complete_U4_halfdiamond_upper'])
    raw = [flow + F(x['512epoch_apparatus_error_upper']) for x in receipts]
    gains = [x / (1 - 2 * x) for x in raw]
    bads = [F(x['512epoch_ideal_bad_weight_upper']) for x in receipts]
    joint_bad = 1 - (1 - bads[0]) * (1 - bads[1])
    if not F(receipt['ideal_pair_bad_weight_upper']) == joint_bad:
        raise ValueError("Required check at ancestor line 128: F(receipt['ideal_pair_bad_weight_upper']) == joint_bad")
    loading = F(common['complete_broad_loader_error_upper'])
    route = F(read(I / 'route_RECEIPT.json')['per_trial_timeout_upper'])
    prior_S = sum(gains) + sum(bads) + loading + route + 8 * bit
    if not prior_S / (1 - 2 * prior_S) == F(common['spatial_complete_error_upper']):
        raise ValueError("Required check at ancestor line 131: prior_S / (1 - 2 * prior_S) == F(common['spatial_complete_error_upper'])")
    S = prior_S + 8 * bit - bads[0] * bads[1]
    error = S / (1 - 2 * S)
    if not error == F(receipt['full_joint_complete_error_upper']):
        raise ValueError("Required check at ancestor line 133: error == F(receipt['full_joint_complete_error_upper'])")
    if not 2 * error == F(receipt['full_joint_halfturn_TV_upper']):
        raise ValueError("Required check at ancestor line 134: 2 * error == F(receipt['full_joint_halfturn_TV_upper'])")
    x = raw[0] + 512 * 8 * bit
    C = x / (1 - 2 * x) + gains[1] + loading + route + 16 * bit
    fa = read(I / 'FIELD_A_RESULT.json')
    delta = F(receipts[1]['certified_subblock_error_upper'])
    kappa = 4 * delta + 2 * delta * delta
    gamma = F(128, 2 ** 53 - 128)
    second_comparison = (512 * kappa + 512 * 12 * 4096 * gamma) / (1 - 512 * kappa) + F(1, 10 ** 12)
    control_error = C / (1 - 2 * C) + F(fa['comparison_512_error_upper']) + second_comparison
    lower = F(fa['dephasing_coordinate_halfturn_TV']) - 2 * control_error
    if not lower == F(receipt['active_dephasing_full_joint_TV_lower']):
        raise ValueError("Required check at ancestor line 141: lower == F(receipt['active_dephasing_full_joint_TV_lower'])")
    if not (recomputed['native'] < float(2 * error) and recomputed['dephasing'] > float(lower) > 0.03):
        raise ValueError("Required check at ancestor line 142: recomputed['native'] < float(2 * error) and recomputed['dephasing'] > float(lower) > 0.03")
    result = dict(producer_manifest_sha256=SEAL, verified_producer_files=count, exact_full_tensor_entries=12 ** 4, complete_joint_bins_per_case=65536, cases=4, independent_matrix_power_state_error=float(state_error), independent_joint_instrument_max_error=law_error, retained_old_chart_TV=str(oldtv), new_joint_TV_upper=float(2 * error), active_control_lower=float(lower), all_outcomes_and_costs_verified=True, producer_imports=False, review_group='Codex same-group', audit_kind='Independent complete reconstruction, adapted from prior same-group reviewer')
    return result

def strict_json(path):

    def pairs(rows):
        d = {}
        for k, v in rows:
            if k in d:
                raise ValueError('Duplicate JSON key: ' + k)
            d[k] = v
        return d

    def finite(v):
        x = float(v)
        if not math.isfinite(x):
            raise ValueError('Nonfinite JSON number')
        return x

    def forbidden(v):
        raise ValueError('Nonfinite JSON literal: ' + v)
    return json.loads(path.read_text(), object_pairs_hook=pairs, parse_float=finite, parse_constant=forbidden)

def compare(actual, saved, path='$'):
    if type(actual) is not type(saved):
        raise ValueError('Type mismatch ' + path)
    if isinstance(actual, dict):
        if actual.keys() != saved.keys():
            raise ValueError('Keys mismatch ' + path)
        for k in actual:
            compare(actual[k], saved[k], path + '.' + k)
    elif isinstance(actual, list):
        if len(actual) != len(saved):
            raise ValueError('Length mismatch ' + path)
        for i, (a, b) in enumerate(zip(actual, saved)):
            compare(a, b, path + '[' + str(i) + ']')
    elif isinstance(actual, float):
        if not math.isfinite(actual) or not math.isfinite(saved) or abs(actual - saved) > 1e-12 + 1e-10 * max(abs(actual), abs(saved)):
            raise ValueError('Numeric mismatch ' + path)
    elif actual != saved:
        raise ValueError('Value mismatch ' + path)

def require(value, message):
    if not value:
        raise ValueError(message)

def verify_seal(source, expected):
    manifest = source / 'MANIFEST.sha256'
    require(sha(manifest) == expected, 'manifest root ' + str(source))
    seen = set()
    for line in manifest.read_text().splitlines():
        digest, rel = line.split('  ', 1)
        require(rel not in seen and (not Path(rel).is_absolute()) and ('..' not in Path(rel).parts), 'manifest path')
        seen.add(rel)
        require(sha(source / rel) == digest, 'leaf hash ' + rel)
    return len(seen)

def common_join(common):
    count = verify_seal(common, '5cd7a135f811d1debfaeab274d5544937829d351604748f848c1b151ede3d521')
    r = read(common / 'RECEIPT.json')
    old = read(SOURCE / 'RECEIPT.json')
    require(r['all_gates_pass'] is True and len(r['gates']) == 15 and all((v is True for v in r['gates'].values())), 'complete common verdict')
    extraction = read(common / 'inputs/extraction/RECEIPT.json')
    eta = F(extraction['preparation_inconclusive_upper'])
    require(eta == F(8501, 2 ** 58), 'source extraction bound')

    def extend(x):
        original = x / (1 + 2 * x)
        return (original + eta) / (1 - 2 * (original + eta))

    def up(x):
        return F(-(-x.numerator * 10 ** 15 // x.denominator), 10 ** 15)
    e = up(extend(F(old['full_joint_complete_error_upper'])))
    c = up(extend(F(old['complete_control_error_upper'])))
    low = F(old['active_dephasing_full_joint_TV_lower']) - 2 * (c - F(old['complete_control_error_upper']))
    low = -up(-low)  # Published lower bound rounds outward, never inward.
    require(F(r['spatial_complete_error_upper']) == e and F(r['full_joint_halfturn_TV_upper']) == 2 * e, 'whole-source error propagation')
    require(F(r['active_dephasing_full_joint_TV_lower']) == low and F(r['complete_control_error_upper']) == c, 'whole-source control propagation')
    require(2 * e < F(11, 10000) and low > F(3, 100), 'common spatial scientific separation')
    a = np.load(SOURCE / 'PAIR_MAPS.npz')
    b = np.load(common / 'COMMON_MAPS.npz')
    for name in ['native', 'dephasing']:
        row = r['full_joint_spatial_records'][name]
        ancestor = old['laws'][name]
        for key in set(ancestor) & set(row):
            require(abs(np.asarray(row[key]) - np.asarray(ancestor[key])).max() < 1e-12, 'common full outcome ' + name + '/' + key)
    for name in ['native', 'dephasing']:
        require(np.max(abs(a[name+'_all65536_position_centres']-b[name+'_all65536_position_centres']))<1e-12,'all common joint bins '+name)
    require(np.max(abs(a['prepared_full144_density']-b['same_prepared_full144_density']))<1e-12,'common full prepared density')
    require(np.max(abs(a['selected_full144_density']-b['same_selected_full144_density']))<1e-14,'common full selected density')
    require(np.array_equal(a['global_Ry_address_permutation'],b['global_Ry_address_permutation']),'common global rotation')
    native=np.load(SOURCE/'inputs/common_MAPS.npz')
    pointer=np.load(SOURCE/'inputs/pointer_witness.npz')['pointer_effects_C']
    for name in ['native','dephasing']:
        local=[native['first_final_density' if name=='native' else 'dephasing_final_density'],native['second_final_density']]
        flags=[[float(np.trace(local_flag_state(native[f'loader{i}_clean_columns'],state,E)).real)for E in pointer]for i,state in enumerate(local)]
        require(np.max(abs(np.array(flags)-r['full_joint_spatial_records'][name]['all_region_flag_centres']))<1e-12,'complete common region flags '+name)
    return {'manifest_sha256': sha(common / 'MANIFEST.sha256'), 'files': count, 'joint_halfturn_upper': str(2 * e), 'active_control_lower': str(low), 'extraction_before_gain': str(eta)}

def main():
    global SOURCE
    parser = argparse.ArgumentParser(description='Read-only complete spatial replay and custody verification')
    parser.add_argument('--source', type=Path, default=SOURCE)
    parser.add_argument('--common', type=Path, default=HERE.parents[2] / 'shared/selected_source_common_v2')
    parser.add_argument('--reference', type=Path, default=HERE / 'RESULT.json')
    parser.add_argument('--record-to', type=Path)
    args = parser.parse_args()
    SOURCE = args.source.resolve()
    start = time.perf_counter()
    require(args.record_to is None or not args.record_to.exists(), 'exclusive first record')
    result = reconstruct()
    result['common_join'] = common_join(args.common.resolve())
    if args.record_to is None:
        compare(result, strict_json(args.reference))
    else:
        with args.record_to.open('x') as f:
            json.dump(result, f, indent=2, sort_keys=True, allow_nan=False)
            f.write('\n')
    print(json.dumps({'all_pass': True, 'seconds': time.perf_counter() - start, 'full_tensor_entries': result['exact_full_tensor_entries'], 'bins_per_case': 65536, 'common_joint_upper': result['common_join']['joint_halfturn_upper'], 'readonly': args.record_to is None}, sort_keys=True))
if __name__ == '__main__':
    main()
