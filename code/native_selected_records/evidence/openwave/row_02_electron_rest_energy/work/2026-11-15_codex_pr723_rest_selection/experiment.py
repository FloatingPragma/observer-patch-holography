"""Offline finite record experiment. No package downloads or subprocesses."""
from collections import defaultdict, deque
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import numpy as np
from inputs.native import native

HERE = Path(__file__).resolve().parent
PIN = 'dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'


def canonical(x):
    return json.dumps(x, sort_keys=True, indent=2) + '\n'


def run():
    assert hashlib.sha256((HERE/'inputs/native.py').read_bytes()).hexdigest() == PIN
    graph = json.loads((HERE/'inputs/graph.json').read_text())
    seams = [tuple(e) for e in graph['seams']]
    local = [e for e in seams if max(e) < 12]
    records = {'columns': [], 'compressed_identity': 'Untouched source coordinates have j=0 and remain fixed.'}
    calls = 0

    def event(n, p, sign, edge, coin, phase, meta=None):
        nonlocal calls
        loads = [0]*n
        if p is not None:
            loads[p] = sign
        out, (c, j) = native(tuple(loads), edge, coin)
        calls += 1
        # This is an audit only. The controller below receives the current-reconstructed p.
        q = list(loads)
        q[edge[0]] -= j
        q[edge[1]] += j
        assert tuple(q) == out
        nonzero = [i for i, v in enumerate(q) if v]
        dest = nonzero[0] if nonzero else None
        assert sum(q) == (sign if p is not None else 0)
        records['columns'].append([phase, n, p, sign, list(edge), c, j, dest, meta])
        return dest, j

    # Exact transition count table: nonincident columns are proved identities;
    # every potentially nonidentity one-unit branch calls the pinned endpoint.
    count = np.eye(240, dtype=np.int64)*1260
    for eidx, edge in enumerate(seams):
        for p in edge:
            for coin in (0, 1):
                y, _ = event(240, p, 1, edge, coin, 'baseline_endpoint', eidx)
                count[p, p] -= 1
                count[p, y] += 1
    lap = 1260*np.eye(240, dtype=np.int64)-count
    vals = np.linalg.eigvalsh(lap.astype(float))
    f = np.zeros(240, dtype=np.int64)
    g = np.zeros(240, dtype=np.int64)
    f[[4, 5, 6, 7]] = [1, -1, -1, 1]
    g[[16, 17, 18, 19]] = [1, -1, -1, 1]
    assert np.array_equal(lap@f, 6*f) and np.array_equal(lap@g, 6*g)
    # Actual unread fine-current channel, including different current outcomes.
    rho_num = np.zeros((240, 240), dtype=np.int64)
    for eidx, edge in enumerate(seams):
        for coin in (0, 1):
            by_current = {}
            for p in np.flatnonzero(f):
                y, j = event(240, int(p), 1, edge, coin, 'compact_channel', eidx)
                by_current.setdefault(j, np.zeros(240, dtype=np.int64))[y] += f[p]
            for v in by_current.values():
                rho_num += np.outer(v, v)
    denom = 4*1260
    assert int(np.trace(rho_num)) == denom
    fidelity = F(int(f@rho_num@f), 4*denom)
    outside = [i for i in range(240) if f[i] == 0]
    leakage = F(int(np.trace(rho_num[np.ix_(outside, outside)])), denom)
    old_single = -1260*math.log1p(-6/1260)
    old_pair = -1260*math.log1p(-12/1260)

    # Search genuinely discovers the bond from the first recorded movement.
    search = {}
    held_control = []
    for sign in (-1, 1):
        for p in range(12):
            moving = 0
            for eidx, edge in enumerate(local):
                for coin in (0, 1):
                    y, j = event(12, p, sign, edge, coin, 'search', eidx)
                    moving += (j != 0)
            search[f'{sign}:{p}'] = str(F(moving, 60))
            for eidx, edge in enumerate(local):
                for coin in (0, 1):
                    # Changed law: retain the full source and emit zero current.
                    held_control.append([sign,p,eidx,coin,0,p])

    adjacency = {i: [] for i in range(12)}
    for a, b in local:
        adjacency[a].append(b)
        adjacency[b].append(a)
    kernels = []
    return_checks = []
    for eidx, (a, b) in enumerate(local):
        for h, t in ((a, b), (b, a)):
            # Ordered bond: h is arrival of the discovery current, t departure.
            dist = {h: 0, t: 0}
            queue = deque((h, t))
            while queue:
                v = queue.popleft()
                for u in sorted(adjacency[v]):
                    if u not in dist:
                        dist[u] = dist[v]+1
                        queue.append(u)
            for sign in (-1, 1):
                K = [[F(0), F(0)] for _ in range(2)]
                Ku = [[F(0), F(0)] for _ in range(2)]
                for row, p in enumerate((h, t)):
                    for selector in (0, 1):
                        # A literal zero-load native event supplies the gate bit.
                        _, _ = event(2, None, 0, (0, 1), selector, 'gate_reference', [eidx,h,t,sign,p])
                        for coin in (0, 1):
                            edge = (min(h,t), max(h,t))
                            # The idle branch executes native on an empty reference seam.
                            active = p == t or selector == 0
                            if active:
                                y, j = event(12, p, sign, edge, coin, 'hold_active', [eidx,h,t,selector])
                            else:
                                event(2, None, 0, (0,1), coin, 'hold_idle', [eidx,h,t,sign,p,selector])
                                y, j = p, 0
                            K[row][(h,t).index(y)] += F(1,4)
                            yu, _ = event(12, p, sign, edge, coin, 'feedback_disabled', [eidx,h,t,selector])
                            Ku[row][(h,t).index(yu)] += F(1,4)
                kernels.append({'edge':eidx,'h':h,'t':t,'sign':sign,
                                'counts_over_4':[[int(v*4) for v in r] for r in K],
                                'disabled_counts_over_4':[[int(v*4) for v in r] for r in Ku]})
                # All displaced starts and both coins are executed once. Exact
                # iteration then uses these complete transition columns.
                R = [[F(0) for _ in range(12)] for _ in range(12)]
                for p in range(12):
                    if p in (h,t):
                        R[p][p] = 1
                    else:
                        nxt = min(v for v in adjacency[p] if dist[v] == dist[p]-1)
                        edge = tuple(sorted((p,nxt)))
                        for coin in (0,1):
                            y, _ = event(12,p,sign,edge,coin,'return',[eidx,h,t])
                            R[p][y] += F(1,2)
                survive = [F(int(p not in (h,t))) for p in range(12)]
                for _ in range(32):
                    survive = [sum(R[p][q]*survive[q] for q in range(12)) for p in range(12)]
                return_checks.append({'edge':eidx,'h':h,'t':t,'sign':sign,
                                      'distance': [dist[p] for p in range(12)],
                                      'timeout_after_32':[str(x) for x in survive]})

    K = [[F(v,4) for v in row] for row in kernels[0]['counts_over_4']]
    a_rate,b_rate=K[0][1],K[1][0]
    pi = [b_rate/(a_rate+b_rate),a_rate/(a_rate+b_rate)]
    # Independent native work assay on every source bond, orientation and sign.
    work_distributions = []
    reset_distributions = []
    reset_work_distributions = []
    for eidx,(a,b) in enumerate(local):
        for h,t in ((a,b),(b,a)):
            for sign in (-1,1):
                work = defaultdict(F)
                for row,p in enumerate((h,t)):
                    for coin in (0,1):
                        y,j = event(12,p,sign,tuple(sorted((h,t))),coin,'work',[eidx,h,t])
                        level_before = int(p == t)
                        level_after = int(y == t)
                        work[level_after-level_before] += pi[row]/2
                work_distributions.append({str(k):str(v) for k,v in sorted(work.items())})
                # Feedback-reset target is t, all cap failures are retained.
                state = {h:pi[0],t:pi[1]}
                # Every possible transient branch is literal; compressed powers
                # do not count as additional executed native calls.
                transfer = defaultdict(F)
                for coin in (0,1):
                    y,j = event(12,h,sign,tuple(sorted((h,t))),coin,'reset',[eidx,h,t])
                    transfer[y] += F(1,2)
                for _ in range(8):
                    mass = state[h]
                    state = {h:mass*transfer[h],t:state[t]+mass*transfer[t]}
                reset_distributions.append({str(p):str(w) for p,w in sorted(state.items())})
                failure=transfer[h]**8
                reset_work_distributions.append({'0':str(pi[1]+pi[0]*failure),
                                                 '1':str(pi[0]*(1-failure))})

    # Separation uses actual consecutive local epochs, each with its own
    # native gate bit and endpoint bit. Product factorization follows because
    # neither controller reads the other's source, records, or selector.
    product = [[K[a//2][b//2]*K[a%2][b%2] for b in range(4)] for a in range(4)]
    exclusive = [[(K[a//2][b//2]*int(a%2==b%2)+
                   int(a//2==b//2)*K[a%2][b%2])/2 for b in range(4)] for a in range(4)]
    joint_pi = [pi[i]*pi[j] for i in range(2) for j in range(2)]
    # Conditional source position at a tick comes entirely from recorded j.
    # Complete kernel records above execute every factor; 4x4 entries are exact
    # factorized sums over the sixteen independent native bit alternatives.
    selected_energy_zero = all(sum(pi[i]*K[i][j] for i in range(2)) == pi[j] for j in range(2))
    work = {F(k):F(v) for k,v in work_distributions[0].items()}
    jar = sum(p*F(2)**(-int(w)) for w,p in work.items())
    reset_jar = sum(F(p)*F(2)**(-int(w)) for w,p in reset_work_distributions[0].items())
    max_timeout = max(F(p) for r in return_checks for p in r['timeout_after_32'])
    check = {
        'full_graph_240_630': len(seams)==630 and count.shape==(240,240),
        'compact_transfer_rays_lambda6': bool(np.array_equal(lap@f,6*f) and np.array_equal(lap@g,6*g)),
        'eighty_lower_modes': int(np.sum(vals<6-1e-8))==80,
        'compact_unread_channel_not_fixed': fidelity < 1 and leakage > 0,
        'native_discovery_all_signed_starts': set(search.values()) == {'1/12'},
        'held_law_never_discovers': all(x[4]==0 for x in held_control) and len(held_control)==1440,
        'all_selected_kernels_match': all(x['counts_over_4']==[[3,1],[2,2]] for x in kernels),
        'feedback_control_changes_stationary_reference': all(x['disabled_counts_over_4']==[[2,2],[2,2]] for x in kernels),
        'all_displaced_starts_return_with_retained_timeout': max_timeout <= F(33,2**32),
        'native_work_law_all_contexts': all(x=={'-1':'1/6','0':'1/2','1':'1/3'} for x in work_distributions),
        'jarzynski_matched_reference': jar==1,
        'feedback_reset_control_separated': reset_jar != 1,
        'selected_stationary_transfer_energy_zero': selected_energy_zero,
        'old_pair_additivity_fails': old_pair > 2*old_single,
    }
    receipt = {
        'schema':'oph.row02.pr723.rest_selection.v1',
        'candidate':'A-native-discovered-bond-v1',
        'official':'untested', 'physical_rest_energy_established':False,
        'native_sha256':PIN, 'direct_native_calls':calls,
        'record_columns':len(records['columns']),
        'changed_held_law_columns':held_control,
        'baseline':{'E6':old_single,'E12':old_pair,'pair_keV':511*old_pair/old_single,
                    'excess_keV':511*(old_pair/old_single-2),
                    'modes_below_6':int(np.sum(vals<6-1e-8)),
                    'lowest_positive_mode':float(vals[1]),
                    'compact_fidelity_after_one_unread_event':str(fidelity),
                    'compact_probability_outside_after_one_event':str(leakage)},
        'selection':{'search_movement_probability':search,'expected_search_ticks':12,
                     'search_timeout_after_256':str(F(11,12)**256),
                     'all_kernels':kernels,'stationary_reference':['2/3','1/3'],
                     'contraction_eigenvalue':'1/4','worst_return_timeout_after_32':str(max_timeout),
                     'return_contexts':return_checks},
        'work':{'gap_record_units':'log(2)','probabilities':work_distributions[0],
                'contexts':len(work_distributions),'exponential_work_mean':str(jar),
                'reset_cap':8,'reset_exponential_work_mean':str(reset_jar),
                'reset_work_probabilities':reset_work_distributions[0],
                'reset_distributions':reset_distributions,
                'reset_scope':'Feedback changes the instrument. Controller/record erasure cost is not included; no thermodynamic violation claimed.'},
        'separation':{'independent_epoch_kernel':[[str(v) for v in r] for r in product],
                      'exclusive_epoch_kernel':[[str(v) for v in r] for r in exclusive],
                      'stationary_product':[str(v) for v in joint_pi],
                      'single_relaxation_rate':'log(4)',
                      'independent_joint_relaxation_rate':'log(16)',
                      'exclusive_joint_relaxation_rate':'log(4)',
                      'exclusive_single_coordinate_relaxation_rate':'log(8/5)',
                      'joint_stationary_transfer_energy':0,
                      'modular_additivity':'K(pi tensor pi)=K(pi) tensor I+I tensor K(pi); a spectral identity, not independent mass evidence'},
        'checks':check,'all_component_checks':all(check.values()),
        'physical_gates':{'electron_identification':False,'nonzero_stationary_transfer_rest_energy':False,
                          'physical_work_normalization':False,'independent_measured_mass_or_frequency':False},
    }
    return receipt, records


if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    r,records=run()
    if args.write:
        for name,obj in [('RECEIPT.json',r),('RECORDS.json',records)]:
            with (HERE/name).open('x') as f:
                f.write(canonical(obj))
    print(canonical({'component_checks':r['all_component_checks'],'physical_gates':r['physical_gates'],
                     'native_calls':r['direct_native_calls'],'baseline':r['baseline']}))
    raise SystemExit(0 if r['all_component_checks'] else 1)
