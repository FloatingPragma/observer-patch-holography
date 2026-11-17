"""Actual native unknown-total selection, extraction and neutral cleanup."""
import argparse
from collections import Counter, deque
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import random

from native_event import native
from observer import decode_census, transport

HERE = Path(__file__).resolve().parent
EDGES = [tuple(e) for e in json.loads((HERE/'graph.json').read_text())['seams']]
ADJ = {i: sorted(b if a == i else a for a,b in EDGES if i in (a,b)) for i in range(240)}
CAP = 64
COUNTS = Counter()


def encoded(x):
    return (json.dumps(x, sort_keys=True, separators=(',', ':'))+'\n').encode()


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def shortest(start, targets, allowed):
    todo = deque([start]); paths = {start: [start]}
    while todo:
        a = todo.popleft()
        if a in targets:
            return paths[a]
        for b in ADJ[a]:
            if b in allowed and b not in paths:
                paths[b] = paths[a]+[b]; todo.append(b)
    raise ValueError('Disconnected requested source route')


def spanning(nodes):
    root = min(nodes); seen = {root}; todo = [root]; edges = []
    for a in todo:
        for b in ADJ[a]:
            if b in nodes and b not in seen:
                seen.add(b); todo.append(b); edges.append((a,b))
    if seen != set(nodes):
        raise ValueError('Disconnected graph')
    return edges


CELLS = [set(range(12)), set(range(12,24))]
TREES = [spanning(c) for c in CELLS]


class Source:
    """The controller uses only event(); audit() is called outside decisions."""
    def __init__(self, initial):
        self.__state = tuple(initial)

    def event(self, edge, coin, kind, law='native'):
        if law == 'native':
            self.__state, (_,current) = native(self.__state, edge, coin)
            COUNTS[kind] += 1
        elif law == 'swap_hold':
            a,b = edge; current = coin*(self.__state[a]-self.__state[b])
            q = list(self.__state); q[a] -= current; q[b] += current
            self.__state = tuple(q); COUNTS['swap_hold'] += 1
        else:
            raise ValueError(law)
        return [*edge,coin,current]

    def audit(self):
        return list(self.__state)


class Controller:
    def __init__(self, source, seed, selection_law='native'):
        self.source = source; self.rng = random.Random(seed)
        self.selection_law = selection_law
        self.observed = [0]*240
        self.events = []; self.actions = []; self.movement_requests = 0

    def event(self, a, b, kind, law='native', forced_coin=None):
        coin = self.rng.randrange(2) if forced_coin is None else forced_coin
        record = self.source.event((a,b), coin, kind, law)
        self.events.append(dict(kind=kind, law=law, record=record))
        self.observed = transport(self.observed, record)
        return record

    def census(self, cell):
        tape = []
        for a,b in TREES[cell]:
            first = None
            for _ in range(CAP):
                r = self.event(a,b,'census'); tape.append(r)
                if first is None:
                    first = r[2]
                elif first != r[2]:
                    break
            else:
                return dict(cell=cell, tape=tape, decoded=None)
        decoded = decode_census(tape, TREES[cell], sorted(CELLS[cell]), (0,14)[cell])
        self.observed[12*cell:12*cell+12] = decoded['terminal']
        return dict(cell=cell, tape=tape, decoded=decoded)

    def exchange(self, a, b, kind, law='native', force_hold=False):
        before = [self.observed[a],self.observed[b]]
        if abs(before[0]-before[1]) != 1:
            self.actions.append(dict(kind=kind, request=[a,b], outcome='not_adjacent_levels', before=before))
            return False
        self.movement_requests += 1
        start = len(self.events)
        for _ in range(CAP):
            d = self.observed[a]-self.observed[b]
            r = self.event(a,b,kind,law, (d+1)//2 if force_hold else None)
            if r[3]:
                good = [self.observed[a],self.observed[b]] == before[::-1]
                self.actions.append(dict(kind=kind, request=[a,b], before=before,
                    span=[start,len(self.events)], outcome='moved' if good else 'wrong_move'))
                return good
        self.actions.append(dict(kind=kind, request=[a,b], before=before,
            span=[start,len(self.events)], outcome='timeout'))
        return False

    def select(self, cell):
        nodes = CELLS[cell]; drops = 0
        # Each native successful repair drops the integer potential by at least2.
        maximum_drops = 12*64*64//2
        for _ in range(maximum_drops+1):
            high = max(self.observed[v] for v in nodes)
            low = min(self.observed[v] for v in nodes)
            if high-low <= 1:
                return dict(outcome='selected', drops=drops,
                    state=[self.observed[v] for v in sorted(nodes)])
            start = min(v for v in nodes if self.observed[v] == high)
            route = shortest(start, {v for v in nodes if self.observed[v] <= high-2}, nodes)
            for a,b in zip(route[:-2],route[1:-1]):
                if self.observed[b] == high:
                    continue
                if not self.exchange(a,b,'selection_exchange',self.selection_law):
                    return dict(outcome='selection_move_failure', drops=drops)
            a,b = route[-2:]
            old = sum(self.observed[v]**2 for v in nodes)
            at = len(self.events)
            self.event(a,b,'selection_repair',self.selection_law)
            new = sum(self.observed[v]**2 for v in nodes)
            self.actions.append(dict(kind='strict_repair', request=[a,b], event=at,
                potential_before=old, potential_after=new))
            if new > old-2:
                return dict(outcome='nondecreasing_repair', drops=drops)
            drops += 1
        return dict(outcome='selection_bound_exceeded',drops=drops)

    def move(self, a, b, sign, kind, force_hold=False):
        if self.observed[a] != sign or self.observed[b] != 0:
            self.actions.append(dict(kind=kind, request=[a,b,sign],
                observed=[self.observed[a],self.observed[b]], outcome='occupied_or_wrong_sign'))
            return False
        return self.exchange(a,b,kind,force_hold=force_hold)

    def retain(self, target, sign, nodes):
        if self.observed[target] == sign:
            return True
        route = min((shortest(v,{target},nodes) for v in nodes if self.observed[v]==sign), key=lambda p:(len(p),p))
        for a,b in zip(route,route[1:]):
            if not self.move(a,b,sign,'retention_move'):
                return False
        return True

    def cleanup(self):
        nodes = set(range(240))-{0,14}; pairs = 0
        while any(self.observed[v] for v in nodes):
            positives = [v for v in nodes if self.observed[v] == 1]
            negatives = {v for v in nodes if self.observed[v] == -1}
            if not positives or not negatives:
                return dict(outcome='cleanup_unbalanced', pairs=pairs)
            route = min((shortest(v,negatives,nodes) for v in positives), key=lambda p:(len(p),p))
            if any(self.observed[v] != 0 for v in route[1:-1]):
                return dict(outcome='nonvacant_cleanup_path',pairs=pairs)
            for a,b in zip(route[:-2],route[1:-1]):
                if not self.move(a,b,1,'cleanup_move'):
                    return dict(outcome='cleanup_move_failure',pairs=pairs)
            a,b = route[-2:]; at = len(self.events)
            self.event(a,b,'cleanup_contact')
            if self.observed[a] or self.observed[b]:
                return dict(outcome='contact_failure',pairs=pairs)
            pairs += 1
            self.actions.append(dict(kind='annihilated_pair',request=[a,b],event=at))
        return dict(outcome='cleaned',pairs=pairs)

    def canonicalize(self, sign):
        if sign == -1:
            # Actual charge motion, never an algebraic sign flip.
            requests = [(0,1,-1,{14}), (14,0,1,{1}), (1,14,-1,{0})]
            for start,target,s,blocked in requests:
                route = shortest(start,{target},set(range(240))-blocked)
                for a,b in zip(route,route[1:]):
                    if not self.move(a,b,s,'orientation_move'):
                        return False
        return self.move(0,1,1,'frame_move') and self.move(1,0,1,'frame_move')


def initial_case(k):
    x = [0]*240
    for v,value in [(4,48),(5,k-48),(6,32),(7,-32),(8,16),(9,-16),
                    (16,-48),(17,48-k),(18,-32),(19,32),(20,-16),(21,16)]:
        x[v] = value
    return x


def execute(initial, seed, law):
    source = Source(initial); control = Controller(source,seed,law)
    censuses = [control.census(c) for c in (0,1)]
    selections = []
    outcome = 'census_timeout'; cleanup = None; orientation = None
    if all(r['decoded'] is not None for r in censuses):
        for cell in (0,1):
            selections.append(control.select(cell))
        outcome = 'selection_rejected'
        if all(r['outcome']=='selected' for r in selections):
            q = control.observed
            totals = [sum(q[v] for v in c) for c in CELLS]
            if totals[0] != -totals[1]:
                outcome = 'nonneutral_input'
            elif totals[0] == 0:
                outcome = 'no_unit'
            else:
                orientation = 1 if totals[0]>0 else -1
                outcome = 'retention_failure'
                if control.retain(0,orientation,CELLS[0]) and control.retain(14,-orientation,CELLS[1]):
                    cleanup = control.cleanup(); outcome = cleanup['outcome']
                    if outcome == 'cleaned':
                        outcome = 'handoff' if control.canonicalize(orientation) else 'orientation_failure'
    wanted = [0]*240; wanted[0]=1; wanted[14]=-1
    if outcome=='handoff':
        assert control.observed == wanted
    return dict(law=law, censuses=censuses, selections=selections,
        cleanup=cleanup, initial_orientation=orientation, outcome=outcome,
        movement_requests=control.movement_requests, events=control.events,
        actions=control.actions, observed_source=control.observed,
        source_audit_only=source.audit(), initial_audit_only=initial,
        pair_code=14 if outcome=='handoff' else None,
        positive_frame=[0,1] if outcome=='handoff' else None)


def challenges():
    initial=[0]*240;initial[0]=1;initial[14]=-1
    s=Source(initial);c=Controller(s,1);c.observed=list(initial)
    assert not c.move(0,1,1,'forced_native_hold',force_hold=True)
    held=dict(initial=initial,events=c.events,actions=c.actions,observed=c.observed,audit=s.audit())
    occupied=list(initial);occupied[1]=1;occupied[15]=-1
    s=Source(occupied);c=Controller(s,2);c.observed=list(occupied)
    assert not c.move(0,1,1,'occupied_target')
    assert not c.move(14,13,1,'wrong_sign')
    rejected=dict(initial=occupied,events=c.events,actions=c.actions,observed=c.observed,audit=s.audit())
    return dict(forced_native_timeout=held,invalid_route_requests=rejected)


def calculate():
    COUNTS.clear(); trials=[]
    for k in range(-12,13):
        for law in ('native','swap_hold'):
            r=execute(initial_case(k),72361000,law)
            r['harness_K_audit_only']=k;trials.append(r)
    challenge=challenges()
    native_rows=[r for r in trials if r['law']=='native']
    controls=[r for r in trials if r['law']=='swap_hold']
    graph_certificates={name:spanning(nodes) for name,nodes in {
        'cleanup_238':set(range(240))-{0,14},
        'orientation_avoid1':set(range(240))-{1},
        'orientation_avoid0':set(range(240))-{0},
        'orientation_avoid14':set(range(240))-{14}}.items()}
    census_ok=all(all(c['decoded'] is not None and c['decoded']['initial']==r['initial_audit_only'][12*c['cell']:12*c['cell']+12] for c in r['censuses']) for r in trials)
    gates=dict(all_record_only_censuses_exact=census_ok,
        native_all25_select=all(len(r['selections'])==2 and all(s['outcome']=='selected' for s in r['selections']) for r in native_rows),
        native24_actual_handoffs=sum(r['outcome']=='handoff' for r in native_rows)==24,
        native_zero_is_no_unit=[r['outcome'] for r in native_rows if r['harness_K_audit_only']==0]==['no_unit'],
        controls25_rejected=all(r['outcome']=='selection_rejected' for r in controls),
        all_terminal_states_agree=all(r['observed_source']==r['source_audit_only'] for r in trials),
        reverse12_physically_routed=sum(r['initial_orientation']==-1 and r['outcome']=='handoff' for r in native_rows)==12,
        all_controllers_below_public_move_bound=max(r['movement_requests'] for r in trials)<=544020,
        graphs_connected=all(len(v)==(237 if k=='cleanup_238' else 238) for k,v in graph_certificates.items()),
        native_timeout64_retained=len(challenge['forced_native_timeout']['events'])==64,
        occupied_and_wrong_sign_rejected=len(challenge['invalid_route_requests']['actions'])==2)
    preparation_failure=Fraction(22,2**63)+Fraction(544020,2**64)
    result=dict(candidate='A-old240-unknown-total-unit-extraction-v1',gates=gates,
        all_gates_pass=all(gates.values()),literal_native_calls_by_kind={k:v for k,v in COUNTS.items() if k!='swap_hold'},
        literal_native_calls=sum(v for k,v in COUNTS.items() if k!='swap_hold'),
        changed_law_calls=COUNTS['swap_hold'],preparation_inconclusive_upper=str(preparation_failure),
        preparation_inconclusive_upper_float=float(preparation_failure),
        complete_outcomes=dict(Counter(r['outcome'] for r in trials)),
        observed_move_request_max=max(r['movement_requests'] for r in native_rows),
        all_K_labels_hidden_from_controller=True,total_one_promise_used=False,
        controller_coin_seed=72361000,controller_seed_independent_of_K=True,
        known_initial_vacuum_datums=[0,14],reset_used=False,
        native_input_domain_max_abs=64,input_cell_total_domain=[-12,12],
        source_handoff_pair_code=14,positive_frame=[0,1],
        coherent_common_apparatus_freshly_executed=False,
        independent_topology_established=False,full_row_solved=False,other_group_G6=False)
    evidence=dict(trials=trials,challenges=challenge,graph_spanning_certificates=graph_certificates)
    return result,evidence


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true');args=parser.parse_args()
    prereg=json.loads((HERE/'PREREG.json').read_text())
    for n,h in prereg['sha256'].items():
        assert digest(HERE/n)==h,n
    result,evidence=calculate();raw=gzip.compress(encoded(evidence),mtime=0)
    result['records_sha256']=hashlib.sha256(raw).hexdigest();result['prereg_sha256']=digest(HERE/'PREREG.json')
    if args.record:
        for name,data in [('RECEIPT.json',encoded(result)),('RECORDS.json.gz',raw)]:
            with (HERE/name).open('xb') as f:f.write(data)
    else:
        assert json.loads((HERE/'RECEIPT.json').read_text())==result
        assert (HERE/'RECORDS.json.gz').read_bytes()==raw
    print(json.dumps(result,indent=2));assert result['all_gates_pass']


if __name__=='__main__':main()
