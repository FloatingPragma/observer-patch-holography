"""Persistent record-only controller with an opaque source and complete traces."""
from pathlib import Path
from collections import Counter, deque
import json
import random
from fractions import Fraction as F
from inputs.native import native

HERE=Path(__file__).resolve().parent
GRAPH=json.loads((HERE/'inputs/graph.json').read_text())
EDGES=[tuple(e) for e in GRAPH['seams'] if max(e)<12]
ADJ={i:[] for i in range(12)}
for a,b in EDGES:
    ADJ[a].append(b)
    ADJ[b].append(a)


class Source:
    def __init__(self,p,sign,held=False):
        self._loads=tuple(sign if i==p else 0 for i in range(14))
        self.held=held

    def step(self,edge,coin):
        if self.held:
            return [coin,0]
        self._loads,record=native(self._loads,edge,coin)
        return list(record)


class Observer:
    def __init__(self,p,sign):
        self.q=[sign if i==p else 0 for i in range(14)]
        self.sign=sign
        self.h=self.t=None
        self.mode='SEARCH'
        self.last_gate=None

    def position(self):
        occupied=[i for i,v in enumerate(self.q) if v]
        assert len(occupied)==1 and self.q[occupied[0]]==self.sign
        return occupied[0]

    def consume(self,edge,coin,current,tag):
        before=self.position()
        self.q[edge[0]]-=current
        self.q[edge[1]]+=current
        after=self.position()
        if tag=='search' and self.h is None and current!=0:
            self.t,self.h=before,after
            self.mode='HOLD'
        elif tag=='gate':
            self.last_gate=coin
        elif tag=='kick':
            self.mode='PERTURBED'
        elif tag=='return':
            self.mode='HOLD' if after in (self.h,self.t) else 'RETURN'
        return [after,self.h,self.t,self.mode,self.last_gate]

    def edge_toward(self,targets):
        dist={v:0 for v in targets}
        queue=deque(targets)
        while queue:
            v=queue.popleft()
            for u in ADJ[v]:
                if u not in dist:
                    dist[u]=dist[v]+1
                    queue.append(u)
        p=self.position()
        nxt=min(v for v in ADJ[p] if dist[v]==dist[p]-1)
        return tuple(sorted((p,nxt)))


class Trial:
    def __init__(self,p,sign,seed,held=False):
        self.source=Source(p,sign,held)
        self.observer=Observer(p,sign)
        self.rng=random.Random(seed)
        self.data={'preparation':[p,sign],'seed':seed,'held_control':held,
                   'events':[],'work':[],'failures':[]}

    def take(self,edge,tag):
        coin=self.rng.randrange(2)
        coin,current=self.source.step(edge,coin)
        checkpoint=self.observer.consume(edge,coin,current,tag)
        self.data['events'].append([list(edge),coin,current,tag,checkpoint])

    def discover(self,cap=256):
        for _ in range(cap):
            self.take(EDGES[self.rng.randrange(len(EDGES))],'search')
            if self.observer.h is not None:
                return True
        self.data['failures'].append('discovery_timeout')
        return False

    def epoch(self,feedback=True):
        o=self.observer
        self.take((12,13),'gate')
        active=(not feedback) or o.position()==o.t or o.last_gate==0
        self.take(tuple(sorted((o.h,o.t))) if active else (12,13),
                  'hold_active' if active else 'hold_idle')

    def displacement_and_return(self):
        o=self.observer
        target=min(p for p in range(12) if p not in (o.h,o.t))
        self.data['kick_target']=target
        for _ in range(32):
            if o.position()==target:
                break
            self.take(o.edge_toward([target]),'kick')
        if o.position()!=target:
            self.data['failures'].append('kick_timeout')
            return
        self.data['return_start']=len(self.data['events'])
        for _ in range(32):
            if o.position() in (o.h,o.t):
                o.mode='HOLD'
                self.data['return_end']=len(self.data['events'])
                return
            self.take(o.edge_toward([o.h,o.t]),'return')
        self.data['return_end']=len(self.data['events'])
        if o.position() not in (o.h,o.t):
            self.data['failures'].append('return_timeout')

    def assay(self,draws=32):
        o=self.observer
        for _ in range(draws):
            for _ in range(8):
                self.epoch()
            before=int(o.position()==o.t)
            start=len(self.data['events'])
            self.take(tuple(sorted((o.h,o.t))),'work')
            after=int(o.position()==o.t)
            self.data['work'].append([start,before,after,after-before])


def pair_run(exclusive):
    a=Trial(0,1,72302901+int(exclusive))
    b=Trial(7,-1,72302911+int(exclusive))
    found_a,found_b=a.discover(),b.discover()
    if not (found_a and found_b):
        return {'exclusive':exclusive,'a':a.data,'b':b.data,'joint':[],
                'selector_records':[],'failure':'discovery_timeout'}
    for _ in range(8):
        a.epoch()
        b.epoch()
    chooser=Source(0,0)
    rng=random.Random(72302921+int(exclusive))
    joint=[]
    selections=[]
    for k in range(2048):
        before=[int(a.observer.position()==a.observer.t),int(b.observer.position()==b.observer.t)]
        if exclusive:
            c,j=chooser.step((12,13),rng.randrange(2))
            selections.append([c,j])
            (a if c==0 else b).epoch()
        else:
            a.epoch()
            b.epoch()
        after=[int(a.observer.position()==a.observer.t),int(b.observer.position()==b.observer.t)]
        joint.append([before,after,len(a.data['events']),len(b.data['events'])])
    return {'exclusive':exclusive,'a':a.data,'b':b.data,'joint':joint,'selector_records':selections}


def run():
    trials=[]
    controls=[]
    for sign in (-1,1):
        for p in range(12):
            trial=Trial(p,sign,72302000+(sign+1)*100+p)
            if trial.discover():
                for _ in range(16):
                    trial.epoch()
                trial.displacement_and_return()
                if not trial.data['failures']:
                    trial.assay()
            trials.append(trial.data)
            held=Trial(p,sign,72302500+(sign+1)*100+p,held=True)
            held.discover(cap=16)
            controls.append(held.data)
    pairs=[pair_run(False),pair_run(True)]
    work=Counter(w[3] for trial in trials for w in trial['work'])
    total=sum(work.values())
    reference={-1:F(1,6),0:F(1,2),1:F(1,3)}
    deviation=max(abs(F(work[k],total)-p) for k,p in reference.items())
    # These finite-seed frequencies are diagnostics, separate from exact
    # population branch enumeration. All events and all failed runs stay stored.
    report={'schema':'oph.row02.trajectory.v1','candidate':'A-native-discovered-bond-v1',
            'trials':len(trials),'failed_trials':[i for i,t in enumerate(trials) if t['failures']],
            'discovery_bonds':[[next(e[4][1] for e in t['events'] if e[4][1] is not None),
                                next(e[4][2] for e in t['events'] if e[4][2] is not None)]
                               for t in trials if any(e[4][1] is not None for e in t['events'])],
            'native_event_records':sum(len(t['events']) for t in trials)+sum(len(p['a']['events'])+len(p['b']['events'])+len(p['selector_records']) for p in pairs),
            'held_control_event_records':sum(len(t['events']) for t in controls),
            'work_draws':total,'work_counts':{str(k):work[k] for k in (-1,0,1)},
            'max_work_frequency_deviation':str(deviation),
            'empirical_exponential_work_mean':str(sum(F(work[k],total)*F(2)**(-k) for k in (-1,0,1))),
            'finite_seed_scope':'declared pseudorandom trials; no independent physical measurement or fitted mass',
            'gates':{'all24_discover_and_return':all(not t['failures'] for t in trials),
                     'paired_discoveries_complete':all(len(p['joint'])==2048 for p in pairs),
                     'all24_held_controls_fail_discovery':all(t['failures']==['discovery_timeout'] for t in controls),
                     'all768_work_draws_retained':total==768,
                     'work_frequency_deviation_below_1_10':deviation<F(1,10)},
            'physical_rest_energy_established':False}
    return report,{'trials':trials,'held_controls':controls,'pairs':pairs}


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument('--write',action='store_true')
    args=p.parse_args()
    report,traces=run()
    for name,obj in [('TRAJECTORY_RECEIPT.json',report),('TRAJECTORIES.json',traces)]:
        encoded=json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n'
        if args.write:
            with (HERE/name).open('x') as f:f.write(encoded)
        elif (HERE/name).read_text()!=encoded:
            raise SystemExit('Mismatch: '+name)
    print(json.dumps(report,sort_keys=True,indent=2))
    raise SystemExit(0 if all(report['gates'].values()) else 1)
