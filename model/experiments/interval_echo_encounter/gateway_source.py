"""Two pulse phases at one triangular launch and sharp pair reads in the next cell."""
from collections import Counter
from itertools import combinations
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import presence
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import federation, seam_digest

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "gateway_protocol.json"
PIN = "30790962fccd5f90110219c943322744c2afbb158bf826820a144ef7aee0af8a"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def settings():
    data = PROTOCOL.read_bytes().replace(b"\r\n", b"\n")
    require(hashlib.sha256(data).hexdigest() == PIN, "registered gateway protocol changed")
    return json.loads(data)


def geometry(board):
    external = sorted((u,b) for u,b in board.external if u < 12)
    u, b = external[0]
    regular = sorted(set(board.neighbours[u]) - {g for g,_ in external} - {b})
    s, v = next((s,v) for s in regular for v in regular if v in board.neighbours[s])
    receiver = [edge for edge in board.seams if all(x//12 == b//12 for x in edge)]
    require(len(receiver) == 30 and s != v, "receiver or launch triangle missing")
    return s,u,v,b,receiver


def first_paths(board, initial, steps):
    """Integer walk counts on distinct unordered loads, without diagonal waiting steps."""
    paths = Counter({tuple(sorted(initial)): 1})
    first = {tuple(sorted(initial)): (0,1)}
    for n in range(1,steps+1):
        following = Counter()
        for (a,b), weight in paths.items():
            for move, fixed in ((a,b),(b,a)):
                for dest in board.neighbours[move]:
                    if dest != fixed:
                        following[tuple(sorted((dest,fixed)))] += weight
        for edge, count in following.items():
            first.setdefault(edge, (n,count))
        paths = following
    return first


def gate(kappa, duration):
    law = candidate.native_law
    h0,h1 = law.local_term_parts()
    full = expm(-1j*duration*(h0+kappa*h1))
    blank = full[:, :len(law.SRC)]
    require(np.linalg.norm(blank[len(law.SRC):]) < 1e-11, "collar leaked during isolated pulse")
    return blank[:len(law.SRC)]


def reference(bookings):
    law = candidate.native_law
    require(bool(bookings) and all(type(i) is int and 0<=i<128 for i in bookings)
            and bookings == sorted(set(bookings)), "invalid finite counter booking")
    block = gate(candidate.KAPPA,np.pi)
    indices = [law.SRC.index((2,0)),law.SRC.index((0,2))]
    probes = np.eye(law.DIM)[:,[law.SRC.index((0,0)),law.SRC.index((2,0))]]
    echo = law.echo()
    read = (echo @ probes)[np.arange(law.DIM)%len(law.SRC)==law.SRC.index((1,1))]
    require(np.allclose(read.T@read,np.diag([0,1]),atol=1e-13,rtol=0)
            and np.linalg.norm(echo@echo@probes-probes)<1e-13, "invalid nondestructive clock probe")
    state = np.eye(len(law.SRC),dtype=complex)[:,indices[0]]
    previous, word = 0, 0
    appointments = {0:0} if 0 in bookings else {}
    for exposure in range(1,bookings[-1]+1):
        state = block @ state
        weights = abs(state[indices])**2
        active = int(np.argmax(weights))
        require(abs(weights[active]-1)<2e-12, "clock pointer unresolved")
        if active != previous:
            word = (word+1)%128
            if word in bookings:
                appointments[word] = exposure
        previous = active
    require(list(appointments) == bookings, "missing clock booking")
    return appointments


def partition(size):
    # The binary sector permits at most one positive seam at a time.
    flags = np.vstack((np.eye(size,dtype=bool),np.zeros((1,size),bool)))
    records = presence.compile_records(flags,["nominal"]*size,True)
    for row,(status,words) in enumerate(records):
        expected = tuple(29 if j==row else 1 for j in range(size))
        require(words==expected and status==int(row<size), "incorrect public pair reading")


def source(sector, kappa, layout):
    s,u,v,_,_ = layout
    initial = next((a,b) for a,b in combinations(range(12),2) if b not in sector.car.neighbours[a])
    h = sector.hamiltonian(kappa)
    seed = np.zeros(sector.P)
    seed[sector.index(*initial)] = 1
    home, other = sector.index(s,u), sector.index(s,v)
    probability = float(abs(expm_multiply(-1j*np.pi*h,seed)[home])**2)
    require(0<probability<1, "encounter did not produce both first observations")
    partition(1)
    slots = [candidate.native_law.SRC.index(q) for q in ((1,0),(0,1))]
    launches = np.zeros((sector.P,2),complex)
    for j,time in enumerate((np.pi/2,3*np.pi/2)):
        pulse = gate(kappa,time)[:,slots[0]]
        z = pulse[slots]
        relative = z[1]/z[0]
        require(abs(relative-(1j if j==0 else -1j))<1e-12
                and abs(np.vdot(z,z)-1)<1e-12, "triangle pulse failed")
        launches[[home,other],j] = z
    return h, probability, launches


def high_projection(h, columns):
    """Projector onto the upper band via the binomial series for sign(H+5)."""
    x = ((h + 5*sparse.eye(h.shape[0],format="csr"))/16).tocsr()
    require(max(np.diff(x.indptr))<=13 and max(np.asarray(abs(x).sum(axis=1)).flat)<=2,
            "operator exceeds sparse arithmetic budget")
    term = columns.copy()
    result = term.copy()
    degree = 160
    for j in range(1,degree+1):
        term = (term-x@(x@term))*((2*j-1)/(2*j))
        result += term
    rho = 1-(5/16)**2
    error = (17/32)*rho**(degree+1)/(1-rho)+1e-8
    return (columns+x@result)/2, error


def band_budget(h, launches, receiver):
    # ||F Q||² is the largest eigenvalue of F Q F. The other factor is
    # ||Q psi||. Q is the complement of the isolated negative-energy band.
    sources = np.flatnonzero(np.any(launches!=0,axis=1))
    require(len(sources)==2, "band bound needs two emitting seams")
    sites = np.r_[receiver,sources]
    selected = np.zeros((h.shape[0],len(sites)))
    selected[sites,np.arange(len(sites))] = 1
    high, error = high_projection(h,selected)
    compressed = high[receiver,:len(receiver)]
    require(np.linalg.norm(compressed-compressed.T)<1e-11, "nonsymmetric band compression")
    cap = np.sqrt(np.maximum(0,np.r_[np.diag(compressed),np.linalg.eigvalsh(compressed)[-1]])
                  +len(receiver)*error)
    leaked = np.linalg.norm(high[:,len(receiver):]@launches[sources],axis=0)+2*error
    return cap[:,None]*leaked


def forecast(cfg):
    board = federation(cfg["level"])
    layout = geometry(board)
    sector = PairSector(board)
    receiver = [sector.index(*edge) for edge in layout[-1]]
    partition(len(receiver))
    bookings = reference(cfg["appointments"])
    arms = {}
    for kappa in (cfg["kappa"],0):
        h, admit, prepared = source(sector,kappa,layout)
        name = "native" if kappa else "exchange"
        arms[name] = (h,admit,prepared)
        if kappa:
            budget = band_budget(h,prepared,receiver)
            adjacent = sparse.diags(sector.adjacent.astype(float))
            separate = sparse.eye(sector.P,format="csr")-adjacent
            arms["blocked"] = ((adjacent@h@adjacent+separate@h@separate).tocsr(),admit,prepared)
    probabilities = {}
    for name in cfg["laws"]:
        h,admit,state = arms[name]
        previous = 0
        for word,tick in bookings.items():
            if tick>previous:
                state = expm_multiply(-1j*np.pi*(tick-previous)*h,state)
            previous = tick
            require(np.max(abs((abs(state)**2).sum(axis=0)-1))<2e-10, "lost free sector weight")
            local = abs(state[receiver,:])**2
            for j,menu in enumerate(cfg["menus"]):
                hits = local[:,j] if j<2 else local.mean(axis=1)
                p = np.r_[admit*hits, admit*(1-sum(hits)), 1-admit]
                require(min(p)>=0 and abs(sum(p)-1)<1e-13,"invalid joint detector law")
                probabilities[f"{name}_{menu}_{word}"] = p
    return probabilities,budget,dict(seams=layout[-1],source_triangle=layout[:3],
                                     gateway=layout[3],seam_digest=seam_digest(board))
