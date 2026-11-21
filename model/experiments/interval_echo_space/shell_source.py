"""Interval Echo on signed clock collars and the remaining native load contacts."""
from math import pi

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.contrast_memory_spherical_v1 import law as native_contact
from model.executable_model import federation
from model.executable_model.carrier import CELL_SEAMS


def check(condition, message):
    if not condition:
        raise ValueError(message)


def local_source(exchange, contrast, amount=1):
    """Audit both coin orientations and both signs, including a load at either clock port."""
    law = candidate.native_law
    check(law is native_contact and candidate.KAPPA==8,"registered free interaction changed")
    h0,h1 = law.local_term_parts()
    h = exchange*h0+contrast*h1
    eye = np.eye(law.DIM)
    check(amount in (1,2),"one mobile load of amount one or two required")
    loads = [(0,0),(amount,0),(0,amount)]
    def ket(load, coin=None, sign=0):
        group = 0 if coin is None else 1+law.LABELS.index((coin,sign))
        return eye[group*law.N_SRC+law.SRC.index(load)]
    sign_term = 4*contrast*np.eye(2) + exchange*np.array([[.5,-.5],[-.5,.5]])
    residual = 0.
    for coin in (0,1):
        for load in loads:
            clock = np.column_stack([ket(load,coin,s) for s in (2,-2)])
            residual = max(residual,np.linalg.norm(h@clock-clock@sign_term))
    data = np.column_stack([ket(load) for load in loads])
    bare = np.zeros((3,3))
    bare[1:,1:] = contrast*amount**2*np.eye(2)+exchange*np.array([[.5,-.5],[-.5,.5]])
    residual = max(residual,np.linalg.norm(h@data-data@bare))
    # A native write from separated double fuel prepares either sign with the same outgoing loads and coin.
    for load,s in (((2,0),2),((0,2),-2)):
        written = law.echo()@ket(load)
        wanted = (ket((1,1),0,s)+ket((1,1),1,s))/np.sqrt(2)
        residual = max(residual,np.linalg.norm(written-wanted))
    check(residual < 1e-12,"native contact does not preserve declared source")
    return sign_term, bare[1:,1:]


def laplacian(n, edges):
    a,b = np.array(edges,dtype=int).T
    matrix = sparse.coo_matrix((np.ones(2*len(a)),(np.r_[a,b],np.r_[b,a])),shape=(n,n)).tocsr()
    return sparse.diags(np.asarray(matrix.sum(axis=1)).ravel())-matrix


class Board:
    def __init__(self, level, matching):
        self.level, self.h = level, 2.**(-level)
        self.car = federation(level)
        car = self.car
        pairs = [tuple(pair) for pair in matching]
        check(len(pairs)==6 and sorted(i for pair in pairs for i in pair)==list(range(12))
              and set(pairs)<=set(CELL_SEAMS),"a native perfect matching is required")
        self.side = np.tile([next(j for pair in pairs for j,p in enumerate(pair) if p==i)
                             for i in range(12)],car.cells)
        active = [e for e in car.seams if e[0]//12!=e[1]//12 or (e[0]%12,e[1]%12) not in pairs]
        internal = [e for e in active if e[0]//12==e[1]//12]
        self.active = active
        self.inside = laplacian(car.n,internal)/2
        self.outside = laplacian(car.n,car.external)/2
        self.gap = float(eigh(self.inside[:12,:12].toarray(),eigvals_only=True)[1])
        check(self.gap >= (3-np.sqrt(5))/2-1e-13, "internal gap lost")
        self.macros = np.arange(car.cells)//4**level
        self.bins = 2*np.repeat(self.macros,12)+self.side
        edges = [(a//12,b//12) for a,b in car.external]
        self.dual = laplacian(car.cells,edges)
        for step in range(level,0,-1):
            parent = federation(step-1)
            links = [tuple(sorted((a//4,b//4))) for a,b in edges if a//4!=b//4]
            adjacency = sorted((a//12,b//12) for a,b in parent.external)
            keys, counts = np.unique(links,axis=0,return_counts=True)
            check(keys.tolist()==[list(e) for e in adjacency] and np.all(counts==2),
                  "four-child quotient changes angular interfaces")
            # Each parent is a connected star of four edge-adjacent children.
            within = np.zeros(4*(parent.cells,),int)
            for a,b in edges:
                if a//4==b//4:
                    within[[a,b]]+=1
            check(np.all(np.sort(within.reshape(-1,4),axis=1)==[1,1,1,3]),"disconnected parent region")
            edges = adjacency

    def generator(self, member, amount=1):
        exchange,contrast = (.5,-8.) if member=="half_rate_contrast_reversed" else (1.,8.)
        check(member in ("native","half_rate_contrast_reversed","external_exchange_off","parked_signal"),
              "unknown law")
        _,blank = local_source(exchange,contrast,amount)
        # Reconstruct the native diagonal using actual active ends and the original host reference.
        kinetic = self.inside + self.outside
        active_degree = 2*kinetic.diagonal()
        diag = blank[0,0]*active_degree + contrast*amount**2*(6-self.car.degree)
        reference = diag-exchange*kinetic.diagonal()
        check(np.max(abs(reference-5*contrast*amount**2))<1e-11,
              "degree correction fails to be a common phase")
        if member=="external_exchange_off":
            kinetic = self.inside
        if member=="parked_signal":
            kinetic = 0*self.inside
        return (exchange*kinetic+sparse.diags(reference)).tocsr(), exchange, float(reference[0])

    def launch(self, name):
        check(name in ("smooth_0","smooth_1","polarized_0"),"unknown preparation")
        anchor = int(name[-1])
        adjacency = sparse.diags(self.dual.diagonal())-self.dual
        distance = shortest_path(adjacency,indices=np.flatnonzero(self.macros==anchor),
                                 unweighted=True,directed=False).min(axis=0)
        envelope = np.maximum(0.,1-self.h*distance/2)
        target = np.repeat(envelope,12)
        if name=="polarized_0":
            target[self.side==1]=0
        return target/np.linalg.norm(target)


def prepare(board, target, amount=1, root=None):
    """Execute the sharp-port launch with addressed native splitters and positive degree-reference phases."""
    car = board.car
    check(target.shape==(car.n,) and np.isfinite(target).all()
          and abs(np.linalg.norm(target)-1)<2e-13,"invalid requested source")
    edges = set(board.active)
    children = [[] for _ in range(car.n)]
    root = int(np.argmax(abs(target))) if root is None else root
    check(0 <= root < car.n, "unknown fuel contact")
    order = [root]
    visited = {root}
    for p in order:
        for q in car.neighbours[p]:
            if q not in visited and tuple(sorted((p,q))) in edges:
                visited.add(q)
                children[p].append(q)
                order.append(q)
    check(len(order)==car.n,"active source graph disconnected")
    mass = abs(target)**2
    for p in reversed(order):
        mass[p] += sum(mass[q] for q in children[p])
    psi = np.zeros(car.n,complex)
    psi[root]=1
    splits = phases = swaps = 0
    for p in order:
        for j,q in enumerate(children[p]):
            if mass[q]==0:
                continue
            keep = abs(target[p])**2+sum(mass[k] for k in children[p][j+1:])
            angle = 2*np.arctan2(np.sqrt(mass[q]),np.sqrt(keep))
            z = np.exp(-1j*(candidate.KAPPA*amount**2+.5)*angle)*psi[p]
            psi[p],psi[q] = z*np.cos(angle/2),1j*z*np.sin(angle/2)
            splits += 1
    for p in np.flatnonzero(target):
        angle = (np.angle(psi[p])-np.angle(target[p]))%(2*pi)
        q = p
        if car.degree[p]!=5:
            q = next(q for q in car.neighbours[p] if car.degree[q]==5
                     and tuple(sorted((p,q))) in edges)
            psi[[q,p]] = psi[[p,q]]
            swaps += 2
        psi[q] *= np.exp(-1j*angle)
        if q!=p:
            psi[[q,p]] = psi[[p,q]]
        phases += 1
    check(np.linalg.norm(psi-target)<1e-12,"native source failed to prepare requested wave")
    return psi,dict(root=root,splitters=splits,degree_pulses=phases,routing_swaps=swaps)
