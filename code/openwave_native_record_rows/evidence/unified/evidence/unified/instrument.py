"""One fixed native instrument, with a disclosed target-coupling preparation bit.

No event executes on import. This proposed candidate needs registration/pin;
its equivalence to old source observations is proved separately, not assumed.
"""
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path
import hashlib,importlib.util,itertools,sys

HERE=Path(__file__).resolve().parent
PATH=HERE/'inputs/native.py'
assert hashlib.sha256(PATH.read_bytes()).hexdigest()=='dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
spec=importlib.util.spec_from_file_location('unified_endpoint',PATH)
module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
native=module.native
FACES=((0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),
       (1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
       (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),
       (4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1))
INTERNAL=sorted({tuple(sorted((f[j],f[(j+1)%3]))) for f in FACES for j in range(3)})

@dataclass(frozen=True)
class Settings:
    period:int=3
    inverse_spacing:int=1
    inner_radius:int=2
    route_cap:int|None=128
    probe_cap:int|None=256

    def __post_init__(self):
        assert self.period>=3 and self.inverse_spacing>=1
        assert self.inner_radius>=self.period//2
        assert self.route_cap is None or self.route_cap>=1
        assert self.probe_cap is None or self.probe_cap>=1

@dataclass
class Memory:
    outer_position:int
    inner_position:int
    target_coupled:bool
    clock_origin:int
    clock_ticks:int=0
    macro_count:int=0
    native_calls:int=0

class Instrument:
    def __init__(self,settings=Settings()):
        self.settings=settings
        self.outer=self._graph(False)
        self.inner=self._graph(True)
        self.cumulative=[];total=0
        for _,_,axis in self.outer['edges']:
            total+=1 if axis==-1 else settings.inverse_spacing**2
            self.cumulative.append(total)
        self.proposal_weight=total

    def _graph(self,inner):
        s=self.settings
        spatial=range(-s.inner_radius,s.inner_radius+1) if inner else range(-(s.period//2),s.period-s.period//2)
        axes=[spatial]*3+([] if inner else [range(s.period)])
        cells=list(itertools.product(*axes));ci={x:i for i,x in enumerate(cells)};edges=[]
        for x in cells:
            c=ci[x];edges.extend((12*c+a,12*c+b,-1) for a,b in INTERNAL)
            for axis in range(len(x)):
                y=list(x);y[axis]+=1
                if inner:
                    if y[axis]>s.inner_radius:continue
                elif axis==3:y[axis]%=s.period
                elif y[axis]>=s.period-s.period//2:y[axis]=-(s.period//2)
                edges.extend((12*c+k,12*ci[tuple(y)]+k,axis) for k in range(12))
        adj=[[] for _ in range(12*len(cells))];pair={}
        for eid,(a,b,axis) in enumerate(edges):
            adj[a].append(eid);adj[b].append(eid);pair[frozenset((a,b))]=eid
        return dict(cells=cells,ci=ci,edges=edges,adj=adj,pair=pair)

    def prepare(self,outer_position,target_coupled):
        """Disclosed initial preparation, not a later hidden reset."""
        assert 0<=outer_position<12*len(self.outer['cells'])
        assert type(target_coupled) is bool
        return Memory(outer_position,12*self.inner['ci'][(0,0,0)],target_coupled,
                      self.outer['cells'][outer_position//12][3])

    @staticmethod
    def _event(g,position,eid,coin):
        a,b,_=g['edges'][eid]
        out,rec=native((int(position==a),int(position==b)),(0,1),coin)
        j=rec[1]
        if j:
            assert position==(a if j==1 else b)
            position=b if j==1 else a
        return position,j

    def step(self,memory,rng):
        """Complete macro-attempt; every inner event and hold remains returned."""
        s=self.settings
        eid=bisect_right(self.cumulative,rng.randrange(self.proposal_weight))
        coin=rng.randrange(2);a,b,axis=self.outer['edges'][eid]
        inner_events=[];prep_count=probe_count=0;outcome='skip';hold=False
        if axis==-1 and memory.outer_position in (a,b):
            target=self.outer['cells'][memory.outer_position//12][:3];failed=False
            for ax in range(3):
                while self.inner['cells'][memory.inner_position//12][ax]!=target[ax]:
                    z=list(self.inner['cells'][memory.inner_position//12]);z[ax]+=1 if target[ax]>z[ax] else -1
                    nxt=12*self.inner['ci'][tuple(z)]+memory.inner_position%12
                    peid=self.inner['pair'][frozenset((memory.inner_position,nxt))];attempt=0
                    while s.route_cap is None or attempt<s.route_cap:
                        cc=rng.randrange(2);memory.inner_position,j=self._event(self.inner,memory.inner_position,peid,cc)
                        inner_events.append([peid,cc,j]);prep_count+=1;attempt+=1
                        if j:break
                    else:failed=True
                    if failed:break
                if failed:break
            if failed:outcome='preparation_timeout'
            else:
                while True:
                    z=self.inner['cells'][memory.inner_position//12]
                    if z==(0,0,0):outcome='hit';break
                    if max(map(abs,z))==s.inner_radius:outcome='escape';break
                    if s.probe_cap is not None and probe_count>=s.probe_cap:outcome='probe_timeout';break
                    peid=rng.choice(self.inner['adj'][memory.inner_position]);cc=rng.randrange(2)
                    memory.inner_position,j=self._event(self.inner,memory.inner_position,peid,cc)
                    inner_events.append([peid,cc,j]);probe_count+=1
            hold=memory.target_coupled and outcome!='escape'
        if hold:
            out,rec=native((0,0),(0,1),coin)
            assert out==(0,0) and rec==(coin,0);j=0
        else:memory.outer_position,j=self._event(self.outer,memory.outer_position,eid,coin)
        memory.macro_count+=1;memory.native_calls+=len(inner_events)+1
        tick=None
        if self.outer['cells'][memory.outer_position//12][3]==(memory.clock_origin+memory.clock_ticks+1)%s.period:
            memory.clock_ticks+=1;tick=[memory.macro_count,memory.clock_ticks,memory.outer_position]
        return dict(proposal=eid,coin=coin,current=j,auxiliary_hold=hold,
                    target_coupled=memory.target_coupled,inner_events=inner_events,
                    preparation_calls=prep_count,probe_calls=probe_count,outcome=outcome,tick=tick)
