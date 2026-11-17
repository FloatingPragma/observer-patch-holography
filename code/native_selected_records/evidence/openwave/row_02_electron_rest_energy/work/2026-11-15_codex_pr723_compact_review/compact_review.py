"""Independent finite verification only; no files, network, subprocesses or source imports."""
from fractions import Fraction as F
from itertools import product
import json
import math

EDGES = [(0,1),(0,5),(0,7),(0,10),(0,11),(1,5),(1,7),(1,8),(1,9),
         (2,3),(2,4),(2,6),(2,10),(2,11),(3,4),(3,6),(3,8),(3,9),
         (4,5),(4,9),(4,11),(5,9),(5,11),(6,7),(6,8),(6,10),
         (7,8),(7,10),(8,9),(10,11)]
ADJ = {p: sorted(b if a == p else a for a,b in EDGES if p in (a,b)) for p in range(12)}
EVALUATIONS = 0

def endpoint(q, edge, coin):
    """Independent floor/ceiling repair and the observed signed current."""
    global EVALUATIONS
    EVALUATIONS += 1
    a,b = edge
    total = q[a] + q[b]
    low, high = total // 2, total - total // 2
    y = list(q)
    y[a], y[b] = (low,high) if coin == 0 else (high,low)
    return y, q[a] - y[a]

def source_step(p, sign, edge, coin):
    q = [sign if k == p else 0 for k in range(12)]
    actual, current = endpoint(q, edge, coin)
    q[edge[0]] -= current
    q[edge[1]] += current
    assert q == actual
    destination, = [k for k,value in enumerate(q) if value]
    return destination, current

discovery = []
for sign,p in product((-1,1), range(12)):
    moved = sum(source_step(p,sign,e,c)[1] != 0 for e,c in product(EDGES,(0,1)))
    discovery.append(F(moved,60))
assert set(discovery) == {F(1,12)}

worst_return = F(0)
return_contexts = 0
for a,b in EDGES:
    distance = {a:0,b:0}
    queue = [a,b]
    for v in queue:
        for u in ADJ[v]:
            if u not in distance:
                distance[u] = distance[v] + 1
                queue.append(u)
    for h,t,sign in product((a,b),(a,b),(-1,1)):
        if h == t:
            continue
        destinations = []
        for p in range(12):
            if p in (h,t):
                destinations.append((p,p))
            else:
                nxt = min(u for u in ADJ[p] if distance[u] == distance[p] - 1)
                edge = tuple(sorted((p,nxt)))
                destinations.append(tuple(source_step(p,sign,edge,c)[0] for c in (0,1)))
        survival = [F(p not in (h,t)) for p in range(12)]
        for tick in range(32):
            survival = [(survival[x] + survival[y])/2 for x,y in destinations]
        worst_return = max(worst_return,max(survival))
        return_contexts += 1
assert return_contexts == 120 and worst_return == F(33,2**32)

results = {}
for gate_bits in (1,2):
    kernels, works = [], []
    for a,b in EDGES:
        for h,t,sign in product((a,b),(a,b),(-1,1)):
            if h == t:
                continue
            K = [[F(0),F(0)],[F(0),F(0)]]
            for p in (h,t):
                for bits in product((0,1),repeat=gate_bits+1):
                    observed = []
                    for c in bits[:-1]:
                        blank,j = endpoint([0,0],(0,1),c)
                        assert blank == [0,0] and j == 0
                        observed.append(c)
                    if p == t or all(c == 0 for c in observed):
                        y,j = source_step(p,sign,(a,b),bits[-1])
                    else:
                        blank,j = endpoint([0,0],(0,1),bits[-1])
                        assert blank == [0,0] and j == 0
                        y = p
                    K[int(p == t)][int(y == t)] += F(1,2**(gate_bits+1))
            forward,backward = K[0][1],K[1][0]
            pi = [backward/(forward+backward),forward/(forward+backward)]
            work = {-1:F(0),0:F(0),1:F(0)}
            for p,c in product((h,t),(0,1)):
                y,j = source_step(p,sign,(a,b),c)
                work[int(y == t)-int(p == t)] += pi[int(p == t)]/2
            kernels.append(K)
            works.append(work)
    assert all(k == kernels[0] for k in kernels) and all(w == works[0] for w in works)
    K,work = kernels[0],works[0]
    rate = F(1,2**gate_bits)
    assert pi == [1/(1+rate),rate/(1+rate)]
    assert sum(probability*rate**change for change,probability in work.items()) == 1
    independent = [[K[i//2][j//2]*K[i%2][j%2] for j in range(4)] for i in range(4)]
    exclusive = [[(K[i//2][j//2]*int(i%2 == j%2)+K[i%2][j%2]*int(i//2 == j//2))/2 for j in range(4)] for i in range(4)]
    assert independent != exclusive
    eigenvalue = K[0][0] + K[1][1] - 1
    assert all(sum(pi[i]*K[i][j] for i in range(2)) == pi[j] for j in range(2))
    results[str(rate)] = {'contexts':len(kernels),'kernel':K,'stationary':pi,
        'modular_gap':math.log(float(pi[0]/pi[1])),'work_in_gap_units':work,
        'independent_pair_kernel':independent,'exclusive_pair_kernel':exclusive,
        'relaxation_eigenvalue':eigenvalue,'independent_joint_eigenvalue':eigenvalue**2,
        'exclusive_joint_eigenvalue':eigenvalue,'stationary_transfer_energy':0,
        'native_events_per_local_epoch':gate_bits+1}

e6 = -1260*math.log1p(-6/1260)
e12 = -1260*math.log1p(-12/1260)
output = {'scope':'Verification of existing finite results, not discovery, full-record replay or physical mass.',
    'independent_endpoint_evaluations':EVALUATIONS,'discovery_probability':'1/12',
    'worst_return_failure_after_32':worst_return,'return_contexts':return_contexts,
    'scheduler_variants':results,'old_511keV_pair_excess_keV':511*(e12/e6-2),
    'physical_rest_energy_established':False,
    'supplied':'one-unit input, local graph, feedback support/gate policy, modular dictionary and local clock'}
print(json.dumps(output,default=str,indent=2,sort_keys=True))
