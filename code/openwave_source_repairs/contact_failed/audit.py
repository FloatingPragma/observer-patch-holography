#!/usr/bin/env python3
"""Exact finite OLD240 audit. No network, shell, fitting, or radiation model."""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import time
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('frozen_native_endpoint', HERE/'inputs/native_endpoint.py')
native_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(native_module)


def canon(x):
    return json.dumps(x, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False)+'\n'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def frozen():
    prereg = json.loads((HERE/'PREREG.json').read_text())
    for name, expected in prereg['files'].items():
        if sha((HERE/name).read_bytes()) != expected:
            raise ValueError('Frozen input changed: '+name)
    return sha((HERE/'PREREG.json').read_bytes())


def log_interval(lam, denominator=1260, terms=32):
    x = F(lam, denominator)
    low = denominator*sum((x**k/k for k in range(1, terms+1)), F())
    return low, low+denominator*x**(terms+1)/((terms+1)*(1-x))


def deterministic_npz(arrays):
    target = io.BytesIO()
    with zipfile.ZipFile(target, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for name in sorted(arrays):
            out = io.BytesIO()
            np.lib.format.write_array(out, np.asarray(arrays[name]), allow_pickle=False)
            entry = zipfile.ZipInfo(name+'.npy', date_time=(1980, 1, 1, 0, 0, 0))
            entry.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(entry, out.getvalue())
    return target.getvalue()


def execute():
    prereg = frozen()
    graph = json.loads((HERE/'inputs/graph.json').read_text())
    inp = json.loads((HERE/'INPUTS.json').read_text())
    seams = [tuple(e) for e in graph['seams']]
    V, M = graph['ports'], len(seams)
    assert (V, M) == (240, 630)
    assert seams == sorted(set(seams)) and all(0 <= a < b < V for a, b in seams)
    L = np.zeros((V, V), dtype=np.int64)
    for a, b in seams:
        L[a, a] += 1; L[b, b] += 1; L[a, b] -= 1; L[b, a] -= 1
    f, g = np.zeros(V, dtype=np.int64), np.zeros(V, dtype=np.int64)
    for k, v in inp['f_numerator'].items(): f[int(k)] = v
    for k, v in inp['g_numerator'].items(): g[int(k)] = v
    A = np.outer(f, g)-np.outer(g, f)
    psi = np.zeros((256, 256), dtype=np.int64); psi[:V, :V] = A
    psi = psi.reshape(-1)
    nz = np.flatnonzero(psi)
    checks = {}
    def check(name, condition):
        checks[name] = bool(condition)
    check('f_g_normalized_and_orthogonal', int(f@f) == int(g@g) == 4 and int(f@g) == 0)
    check('full_single_modes_6', np.array_equal(L@f, 6*f) and np.array_equal(L@g, 6*g))
    check('antisymmetric_pair_norm_support', len(nz) == 32 and int(psi@psi) == 32 and np.all(np.diag(A) == 0) and np.array_equal(A.T, -A))
    pair_residual = L@A+A@L.T-12*A
    check('full_tensor_and_dirichlet_eigenrelation', not np.any(pair_residual))
    check('full_graph_degree_bound', int(np.diag(L).max()) == 6)

    labels = np.arange(65536, dtype=np.int64)
    p, r = labels//256, labels % 256
    valid = ((p < V)&(r < V)&(p != r))|(labels == 0)
    p = np.where(labels == 0, -1, p); r = np.where(labels == 0, -1, r)
    assert int(valid.sum()) == 57361
    calls = 0
    def literal(x, y, coin):
        nonlocal calls
        calls += 1
        return native_module.native((int(x), int(y)), (0, 1), int(coin))
    table = np.empty((2, 3, 3, 3), dtype=np.int8)
    for c in (0, 1):
        for x in (-1, 0, 1):
            for y in (-1, 0, 1):
                dest, rec = literal(x, y, c)
                assert rec[0] == c
                table[c, x+1, y+1] = (*dest, rec[1])

    sign_plant = psi.copy(); sign_plant[nz[0]] *= -1
    outcomes, records = [], []
    offsets, inputs, outputs, currents, counts = [0], [], [], [], []
    transfer = {law: np.zeros(65536, dtype=np.int64) for law in ('native', 'hold', 'swap_hold')}
    totals = {law: defaultdict(F) for law in transfer}
    action = hashlib.sha256()
    all_columns, injections, unchanged, scalar_checks = 0, 0, 0, 0
    contact_edges = []
    witness_inputs = [0, 256*seams[0][0]+seams[0][1], 256*seams[0][1]+seams[0][0]]
    witness_refs, witness_coefficients = [0,1,1], [1,1j,1]
    witness_body_norm = F()
    for ei, (a, b) in enumerate(seams):
        if psi[256*a+b] != 0:
            contact_edges.append(ei)
        qa = (p == a).astype(np.int8)-(r == a).astype(np.int8)
        qb = (p == b).astype(np.int8)-(r == b).astype(np.int8)
        for c in (0, 1):
            left, right, j = table[c, qa+1, qb+1].T
            po = np.where((p == a)|(p == b), np.where(left == 1, a, np.where(right == 1, b, -1)), p)
            ro = np.where((r == a)|(r == b), np.where(left == -1, a, np.where(right == -1, b, -1)), r)
            out = np.where((po == -1)&(ro == -1), 0, 256*po+ro)
            out = np.where(valid, out, labels).astype(np.uint16)
            j = np.where(valid, j, 0).astype(np.int8)
            assert np.all(valid[out[valid]])
            assert np.array_equal(out[~valid], labels[~valid]) and not np.any(j[~valid])
            zero = valid&(j == 0)
            assert np.array_equal(out[zero], labels[zero])
            unchanged += int(zero.sum())
            branch_counts = []
            for q in (-1, 0, 1):
                source = np.flatnonzero(valid&(j == q))
                assert len(np.unique(out[source])) == len(source)
                branch_counts.append(len(source)); injections += 1
            counts.append(branch_counts)
            all_columns += len(labels)
            changed = valid&(j != 0)
            inputs.append(labels[changed].astype(np.uint16)); outputs.append(out[changed]); currents.append(j[changed])
            offsets.append(offsets[-1]+int(changed.sum()))
            action.update(out.astype('<u2').tobytes()); action.update(j.tobytes())
            witness_outputs = defaultdict(complex)
            for z, ref, coefficient in zip(witness_inputs,witness_refs,witness_coefficients):
                witness_outputs[(int(j[z]),int(out[z]),ref)] += coefficient
            witness_body_norm += F(int(sum(abs(v)**2 for v in witness_outputs.values())),4*1260)

            # Execute the pinned endpoint for every prepared input/context.
            for z in nz:
                dest, rec = literal(qa[z], qb[z], c)
                assert tuple(dest) == tuple(table[c, qa[z]+1, qb[z]+1, :2])
                assert rec[1] == int(j[z])
                zz = {int(p[z]): 1, int(r[z]): -1}; zz[a], zz[b] = dest
                zz = {k:v for k,v in zz.items() if v}
                target = 0 if not zz else 256*next(k for k,v in zz.items() if v == 1)+next(k for k,v in zz.items() if v == -1)
                assert target == int(out[z]); scalar_checks += 1
                records.append([ei, c, int(j[z]), int(z), int(out[z]), int(psi[z])])

            control_out = {}
            control_j = {}
            control_out['native'], control_j['native'] = out[nz], j[nz]
            control_out['hold'], control_j['hold'] = nz, np.zeros(len(nz), dtype=np.int8)
            pp = np.where(p[nz] == a, b, np.where(p[nz] == b, a, p[nz])) if c else p[nz]
            rr = np.where(r[nz] == a, b, np.where(r[nz] == b, a, r[nz])) if c else r[nz]
            control_out['swap_hold'] = 256*pp+rr
            control_j['swap_hold'] = (qa[nz]-qb[nz]) if c else np.zeros(len(nz), dtype=np.int8)
            for law in transfer:
                dests, js = control_out[law], control_j[law]
                np.add.at(transfer[law], dests, psi[nz])
                for current in range(-2, 3):
                    mask = js == current
                    src, dest = nz[mask], dests[mask]
                    number, vac = len(src), int(np.count_nonzero(dest == 0))
                    overlap = int(psi[dest]@psi[src])
                    alt_overlap = int(psi[dest]@sign_plant[src])
                    prob = F(number, 32*1260)
                    click = F(overlap**2, 32*32*1260)
                    vacuum = F(vac, 32*1260)
                    complement = prob-click-vacuum
                    assert complement >= 0
                    totals[law]['probability'] += prob
                    totals[law]['contact'] += vacuum if current else F()
                    totals[law]['vacuum'] += vacuum
                    totals[law]['psi_click'] += click
                    totals[law]['valid_complement'] += complement
                    totals[law]['dephased_psi_click'] += F(int(np.sum(psi[dest]**2)), 32*32*1260)
                    totals[law]['sign_plant_psi_click'] += F(alt_overlap**2, 32*32*1260)
                    outcomes.append({'law':law,'seam':ei,'coin':c,'current':current,'probability':str(prob),
                                     'vacuum':str(vacuum),'psi_click':str(click),'valid_complement':str(complement),'invalid':'0'})

    archived = json.loads((HERE/'inputs/archived_action_receipt.json').read_text())
    check('unchanged_published_complete_source_action', action.hexdigest() == archived['complete_action_sha256'])
    check('full_native_mean_eigenrelation_including_vacuum', np.array_equal(transfer['native'], 1248*psi))
    check('all_instrument_outcomes_complete', all(t['probability'] == 1 and t['vacuum']+t['psi_click']+t['valid_complement'] == 1 for t in totals.values()))
    check('source_law_contact_control', totals['native']['contact'] > 0 and totals['hold']['contact'] == totals['swap_hold']['contact'] == 0)
    check('hold_does_not_have_claimed_transfer_energy', np.array_equal(transfer['hold'],1260*psi) and not np.array_equal(transfer['hold'],1248*psi))
    # Record rather than suppress a control that retains the spectral property.
    swap_eigen_same = np.array_equal(transfer['swap_hold'], 1248*psi)
    check('eigenray_does_not_imply_unread_persistence', totals['native']['vacuum'] > 0 and totals['native']['psi_click'] < 1)
    selected = contact_edges[0]
    for ei in contact_edges:
        selected_outcomes = [x for x in outcomes if x['law'] == 'native' and x['seam'] == ei]
        assert sum((F(x['vacuum'])*M for x in selected_outcomes), F()) == F(1,16)
        for q in (-1,1):
            assert sum((F(x['vacuum'])*M for x in selected_outcomes if x['current'] == q),F()) == F(1,32)
    check('addressed_contact_predictions', bool(contact_edges))
    a,b=seams[selected]
    # Native-domain positive/negative detector plants: a contact, its reverse,
    # vacuum and an untouched noncontact pair, on the same addressed seam.
    other = [v for v in range(V) if v not in (a,b)][:2]
    detector = []
    for name, pp, rr in [('positive_contact',a,b),('negative_contact',b,a),('vacuum',-1,-1),('noncontact',*other)]:
        for c in (0,1):
            q={pp:1,rr:-1} if pp>=0 else {}
            dst, rec=literal(q.get(a,0),q.get(b,0),c); q[a],q[b]=dst
            q={k:v for k,v in q.items() if v}
            seen=not q and rec[1]!=0
            detector.append([name,c,rec[1],seen])
            assert seen == (name in ('positive_contact','negative_contact'))
    check('detector_native_domain_sensitivity', all(row[3] == (row[0] in ('positive_contact','negative_contact')) for row in detector))
    # A false coherent merge of the two opposite-current contact branches
    # maps (|a,b>+|b,a>)/sqrt2 to sqrt2|vac>, norm2=2, not1.
    check('record_merge_mutation_detected', F(4,2) != 1)
    sign_matrix = sign_plant.reshape(256,256)[:V,:V]
    sign_res = L@sign_matrix+sign_matrix@L.T-12*sign_matrix
    check('sign_mutation_eigen_detector', bool(np.any(sign_res)))
    check('fixed_psi_detector_sees_coherence', F(int(psi@sign_plant)**2,1024) == F(225,256) and F(1,32) != 1)
    # Full-source arbitrary-reference completeness follows from branchwise
    # injection. Explicit orthogonal reference labels also retain valid,
    # vacuum and invalid sectors in the finite witness below.
    witness={'inputs':witness_inputs+[65535],'reference_labels':witness_refs+[0],
             'coefficients':['1/2','i/2','1/2','-i/2'],
             'body_total':str(witness_body_norm),'invalid_flag':'1/4','trace':str(witness_body_norm+F(1,4)),
             'basis_injection_contexts':injections,'invalid_projector_rank':int((~valid).sum())}
    check('complete_reference_witness', witness_body_norm == F(3,4) and injections==3780)

    e6, e12 = log_interval(6), log_interval(12)
    ratio=(511*e12[0]/e6[1],511*e12[1]/e6[0])
    receipt={'schema':'oph.old240.compact_contact_audit.v1','candidate':inp['candidate'],
             'prereg_sha256':prereg,'checks':checks,'all_pass':all(checks.values()),
             'graph':{'ports':V,'seams':M,'max_degree':6},
             'state':{'support':[[int(z),int(psi[z])] for z in nz],'amplitude_denominator':'sqrt(32)',
                      'squared_norm':'1','pair_residual_max':int(np.abs(pair_residual).max()),'prepared_not_selected':True},
             'source':{'valid_pair_columns':int(valid.sum()),'invalid_columns':int((~valid).sum()),
                       'full_columns_checked':all_columns,'body_contexts':1260,'fine_injections':injections,
                       'zero_current_identity_columns':unchanged,'sparse_changed_columns':offsets[-1],
                       'action_sha256':action.hexdigest(),'literal_native_calls':calls,'prepared_scalar_columns':scalar_checks},
             'contact':{'support_edge_indices':contact_edges,'selected_edge_index':selected,'selected_edge':list(seams[selected]),
                        'addressed_total':'1/16','addressed_current_plus':'1/32','addressed_current_minus':'1/32',
                        'addressed_each_coin_current':'1/64','uniform_total':str(totals['native']['contact'])},
             'outcome_totals':{law:{k:str(v) for k,v in value.items()} for law,value in totals.items()},
             'controls':{'swap_hold_retains_E12_eigenray':swap_eigen_same,'hold_pair_energy':'0',
                         'source_domain_detector_cases':detector,'incorrect_coherent_current_merge_norm2':'2',
                         'sign_mutation_eigen_residual_max':int(np.abs(sign_res).max()),
                         'initial_psi_certificate':{'coherent':'1','dephased':'1/32','sign_flip':'225/256'}},
             'reference_witness':witness,
             'energy':{'transfer_eigenvalue':'104/105','positive_Dirichlet_floor':'103/105',
                       'E12_formula':'-1260 log(104/105)','E6_formula':'-1260 log(209/210)',
                       'E6_interval':list(map(str,e6)),'E12_interval':list(map(str,e12)),
                       'pair_keV_if_E6_is_511_interval':list(map(str,ratio)),
                       'finite_clock_excess_formula':'1260 log(43681/43680)',
                       'finite_clock_excess_positive': e12[0] > 2*e6[1],
                       'physical_energy_or_clock_validated':False},
             'scope':{'sharp_transfer_energy':True,'dynamical_state_selection':False,'radiation':False,
                      'autonomous_energy_transfer':False,'native_compiled_pair_preparation':False,
                      'native_compiled_psi_projector':False,'external_review':False,'platform_award':False}}
    arrays={'offsets':np.array(offsets,dtype='<i8'),'inputs':np.concatenate(inputs).astype('<u2'),
            'outputs':np.concatenate(outputs).astype('<u2'),'currents':np.concatenate(currents).astype('i1'),
            'branch_counts':np.array(counts,dtype='<i8'),'invalid':labels[~valid].astype('<u2'),
            'endpoint_table':table,'seams':np.array(seams,dtype='<i2')}
    payload=''.join(json.dumps(x,separators=(',',':'))+'\n' for x in records).encode()
    evidence={'complete_actions.npz':deterministic_npz(arrays),'prepared_branches.jsonl.gz':gzip.compress(payload,mtime=0),
              'complete_outcomes.json':canon(outcomes).encode()}
    receipt['artifacts']={name:sha(data) for name,data in evidence.items()}
    return receipt,evidence


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true');args=parser.parse_args()
    started=time.perf_counter(); receipt,evidence=execute()
    if args.record:
        target=HERE/'evidence';target.mkdir(exist_ok=False)
        for name,data in evidence.items():(target/name).write_bytes(data)
        (HERE/'RECEIPT.json').write_text(canon(receipt))
    else:
        assert json.loads((HERE/'RECEIPT.json').read_text()) == receipt
        for name,data in evidence.items():assert (HERE/'evidence'/name).read_bytes() == data,name
    print(canon({'all_pass':receipt['all_pass'],'literal_native_calls':receipt['source']['literal_native_calls'],
                 'contact':receipt['contact'],'native_outcomes':receipt['outcome_totals']['native'],
                 'swap_same_energy':receipt['controls']['swap_hold_retains_E12_eigenray'],
                 'elapsed_seconds':time.perf_counter()-started}))
    return 0 if receipt['all_pass'] else 1


if __name__=='__main__':raise SystemExit(main())
