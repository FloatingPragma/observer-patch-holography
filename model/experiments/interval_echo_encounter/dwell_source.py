"""Native write, stationary dimer timebase and two-load restriction on the released boards."""
from collections import deque

import numpy as np
from scipy import sparse

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model.carrier import Carrier, federation


def require(condition, message):
    if not condition:
        raise ValueError(message)


def local_blocks(kappa):
    """Check physical columns on the full qutrit alphabet before restricting the federation."""
    require(candidate.native_law is law and candidate.KAPPA==8, "changed registered seam source")
    eye = np.eye(law.DIM)
    loads = np.eye(9)
    blank = eye[:,:9]
    def memory(sign):
        v = np.zeros(11)
        for coin in (0,1):
            v[1+law.LABELS.index((coin,sign))] = 1/np.sqrt(2)
        return v
    echo = law.echo()
    exchange,contrast = law.local_term_parts()
    native = exchange+kappa*contrast
    repaired = np.column_stack([np.kron(memory(s), loads[:,law.SRC.index((1,1))]) for s in (2,-2)])
    doubles = blank[:,[law.SRC.index((2,0)),law.SRC.index((0,2))]]
    require(np.linalg.norm(echo@doubles-repaired)<1e-13,"native write does not join the two loads")
    # The native pair reader first resolves (1,1). On a hit the two exterior
    # units are on the tag, so the neighbouring readout port is certainly empty.
    # Unwrite the tag, then echo (port 0, vacancy) and read the same pair
    # predicate there. It distinguishes the two old doubles with no sign probe.
    both = np.arange(law.DIM)%9==law.SRC.index((1,1))
    require(np.linalg.norm(echo@repaired-doubles)<1e-13,"tag unwrite loses the sign")
    rail = echo@blank[:,[law.SRC.index((n,0)) for n in range(3)]]
    require(np.linalg.norm(rail[both].T@rail[both]-np.diag([0,0,1]))<1e-13,
            "vacancy reader does not resolve a double")
    quiet = np.kron(memory(0)[:,None], loads)
    require(np.linalg.norm(native@quiet)<1e-12 and np.linalg.norm(echo@blank[:,0]-quiet[:,0])<1e-13,
            "vacuum-written boundary leaks")
    sign_term = 4*kappa*np.eye(2)+np.array([[.5,-.5],[-.5,.5]])
    always = native+repaired@sign_term@repaired.T
    for i,pair in enumerate(law.SRC):
        basis = np.column_stack([np.kron(memory(s),loads[:,i]) for s in (2,-2)])
        require(np.linalg.norm(native@basis-basis@(sign_term if pair!=(1,1) else 0*sign_term))<1e-12
                and np.linalg.norm(always@basis-basis@sign_term)<1e-12,"dirty contact leaks or ticks incorrectly")
    require(np.linalg.norm(native@blank-blank@law.seam_term(kappa))<1e-12,"blank contact changed")
    return native,always


class Station:
    def __init__(self,level,tag):
        self.host = host = federation(level)
        self.tag = tuple(tag)
        require(self.tag in host.seams and max(tag)<12,"tag must be an internal seam of the launch cell")
        neighbors = [[] for _ in range(host.cells)]
        for i,j in host.external:
            a,b = i//12,j//12
            neighbors[a].append(b)
            neighbors[b].append(a)
        distance = {0:0}
        queue = deque([0])
        while queue:
            a = queue.popleft()
            for b in neighbors[a]:
                if b not in distance:
                    distance[b]=distance[a]+1
                    queue.append(b)
        require(len(distance)==host.cells,"disconnected federation")
        self.clock_cell = min(c for c in distance if distance[c]==max(distance.values()))
        self.clock_distance = distance[self.clock_cell]
        room = set(range(12*self.clock_cell,12*self.clock_cell+12))
        options = [e for e in host.seams if set(e)<=room and all(host.degree[i]==5 for i in e)]
        self.dimers = [options[0],next(e for e in options if not set(e)&set(options[0]))]
        hands = set(self.dimers[0]+self.dimers[1])
        cuts = [e for e in host.seams if (len(set(e)&room)==1 or set(e)&hands) and e not in self.dimers]
        self.cuts = sorted(cuts)
        require(self.tag==(0,1) and (0,5) in host.seams and (0,5) not in self.cuts
                and host.degree[0]==host.degree[5]==5,"missing local vacancy readout contact")
        require(all(sum(i in e and e not in cuts for e in host.seams)==1 for i in hands),
                "dimer clock has unisolated contacts")
        outside = sorted(set(range(host.n))-room)
        ids = {p:i for i,p in enumerate(outside)}
        edges = [(ids[i],ids[j]) for i,j in host.seams if i in ids and j in ids and (i,j)!=self.tag]
        self.graph = Carrier("mobile ports",len(outside),edges)
        self.ends = tuple(ids[p] for p in tag)
        self.lost = host.degree[outside]-self.graph.degree
        self.d = float(sum(self.graph.degree[i] for i in self.ends)/2)
        # Bound the contrast saving on every complement word, without enumerating O(n²) words.
        self.saving = int(max(max(self.lost[a]+self.lost[b]+2 for a,b in edges),
                              sum(sorted(self.lost)[-2:])))
        self.tag_loss = int(sum(self.lost[i] for i in self.ends))
        require(self.tag_loss==2,"clock cut overlaps tag")

    def isolation(self,kappa,j):
        # Subtract 12*kappa and the common far-hand phase. The sharp word has
        # no tag energy. Off the word the sign parity contributes 4*kappa+j.
        seed = self.d-kappa*self.tag_loss
        if kappa>0:
            distance = (4-self.saving)*kappa+j-seed
        else:
            ceiling = (4-self.saving)*kappa+j+2*max(self.graph.degree)
            distance = seed-ceiling
        if distance<=0:
            return None
        ratio2 = (self.d/2)/distance**2
        w = 1/(1+ratio2)
        return dict(gap=float(distance),weight=float(w),survival=float(max(0,2*w-1)**2))

    def flight(self):
        sector = PairSector(self.graph)
        self.sector = sector
        self.p = sector.index(*self.ends)
        self.kinetic = sector.exchange()
        self.cost = -(self.lost[sector.aa]+self.lost[sector.bb]+2.*sector.adjacent)
        self.off_tag = np.ones(sector.P)
        self.off_tag[self.p]=0

    def generator(self,kappa,j,member):
        require(member in ("native","always_ticking") and j in (0,1),"unknown tag law or parity")
        clock = self.off_tag if member=="native" else 1.
        return self.kinetic+sparse.diags(kappa*self.cost+(4*kappa+j)*clock)
