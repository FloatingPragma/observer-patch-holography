"""Recover the accessible material histories and apply the donor's spectral criterion."""

from hashlib import sha256
from pathlib import Path
import json

import numpy as np
import sympy as sp
from scipy.linalg import eigh

from .local import require
from .source import edge


DONOR = Path(__file__).resolve().parents[1]/"row_30_selective"
DONOR_HASHES = {
    "analyzers.json": "bca73380e3f4960b8b4aa1268023573845e9f40463c8545af35ed550e7f17a69",
    "results.json": "06f003e17f7ad04be6e76b4d48e9cbabc4afbd697711449d5cb7e33467d23447",
    "PROTOCOL.md": "6514822767043ce397a2f313bbfe9e2fcb93c4b1a309d0279298f6f47f3bdb9f",
    "../row_30_arrows/PROTOCOL.md": "5097dd73eb3da0f0aae98113162d14020c607af6b89a498e67b74adcdbf58ec9",
    "staggered.json": "8ea06d8a234a5a63c091ceb9d40a25f951c8209e09c04e8ff624f37cf7ed88bf",
}


class History:
    def __init__(self, board, heads):
        self.board = board
        self.heads = heads
        self.paths = [(0,)]
        parents = []
        for i, path in enumerate(self.paths):
            endpoint = path[-1]
            for q in board.neighbours[endpoint]:
                if heads[edge(q,endpoint)] == endpoint:
                    require(q not in path, "The donor is not an acyclic orientation")
                    parents.append(i)
                    self.paths.append((*path,q))
        self.index = {w:i for i,w in enumerate(self.paths)}
        self.sites = np.array([w[-1] for w in self.paths])
        dim = len(self.paths)
        child = np.arange(1,dim)
        a = np.zeros((dim,dim),dtype=int)
        a[child,parents] = a[parents,child] = 1
        self.lap = np.diag(a.sum(axis=1))-a
        self.twice_scaled = 2*np.diag(6-board.degree[self.sites])-np.diag(a.sum(axis=1))-a
        # Check the enabled moves and reconstruct the terminal position from
        # physical coin changes. Paths never enter the observation decoder.
        self.read_ports = []
        for j, path in enumerate(self.paths):
            changed = {edge(p,q) for p,q in zip(path,path[1:])}
            read = np.zeros(board.n,dtype=int)
            read[0] = 1
            for e in changed:
                start = heads[e]
                end = next(p for p in e if p!=start)
                read[start] -= 1
                read[end] += 1
            require(np.count_nonzero(read)==1 and read[path[-1]]==1, "Coin reports disagree on load position")
            self.read_ports.append(int(read.argmax()))
            p = path[-1]
            moves = set()
            for q in board.neighbours[p]:
                e = edge(p,q)
                if (heads[e]==p) != (e in changed):
                    new = path[:-1] if e in changed else (*path,q)
                    require(new and new[-1]==q, "A physical hop jumps outside the history tree")
                    moves.add(self.index[new])
            require(moves==set(np.flatnonzero(a[j])), "An allowed move is missing")


def groups(vals):
    return np.split(np.arange(len(vals)), np.where(np.diff(vals)>1e-8)[0]+1)


def donor_selection(tree):
    for name, expected in DONOR_HASHES.items():
        data = (DONOR/name).read_bytes().replace(b"\r\n",b"\n")
        require(sha256(data).hexdigest()==expected, f"Spectroscopy donor differs: {name}")
    archive = json.loads((DONOR/"results.json").read_text(encoding="utf-8"))
    la,z = eigh(tree.lap)
    en,v = eigh(4*tree.twice_scaled)  # historical spectrum at strength 8
    results = []
    for I in groups(en):
        best = []
        for J in groups(la):
            w = np.linalg.svd(v[:,I].T @ z[:,J],compute_uv=False)[0]**2
            if w>=.9:
                best.append(round(float(la[J[0]]),9))
        results.append((round(float(en[I[0]]),9),len(I),best))
    stored = archive["spectral_audit"]["path"]
    require(results==[(a["energy"],a["multiplicity"],a["compatible_k2"]) for a in stored],
            "Different result on the complete historical spectral menu")
    [(energy, multiplicity, [k2])] = [a for a in results if a[2]]
    require(k2==int(k2), "This exact compiler requires an integral donor wave number squared")
    level = sp.Matrix(tree.lap)-int(k2)*sp.eye(len(la))
    basis = sp.Matrix.hstack(*level.nullspace())
    projector = basis*(basis.T*basis).inv()*basis.T
    raw = []
    # The donor's column-order Gram-Schmidt, with no numerical choice on a null column.
    for column in projector.columnspace():
        w = column
        for old in raw:
            w -= old*(old.dot(w)/old.dot(old))
        if w.norm()!=0:
            raw.append(w)
    indices = np.flatnonzero(np.abs(la-k2)<1e-8)
    require(len(raw)==len(indices)==multiplicity,"Incomplete wavelength analyzer space")
    choices = []
    for index,w in zip(indices,raw,strict=True):
        i,j = [i for i,q in enumerate(w) if q]
        require(w[i]==-w[j] and tree.paths[i][:-1]==tree.paths[j][:-1],
                "The historical analyzer is not a signed two-child word")
        s = np.zeros(len(w),dtype=int)
        s[i],s[j] = 1,-1
        np.testing.assert_array_equal(tree.twice_scaled @ s,s)
        choices.append({"analyzer":int(index),"paths":[tree.paths[i],tree.paths[j]]})
    return {"energy_at_8":energy,"k_squared":k2,"distinct_energies_tested":len(results),
            "spatial_columns_tested":len(la),"choices":choices}


def donor_waits(tree, choice, k):
    spec = json.loads((DONOR/"analyzers.json").read_text(encoding="utf-8"))
    theta = np.asarray(spec["pulses"][choice["analyzer"]]["angles"])
    require(theta.shape==(spec["layers"],len(spec["ports"])), "Wrong donor pulse shape")
    masks = np.array([tree.sites==p for p in spec["ports"]],dtype=int)
    vals,vec = eigh(tree.twice_scaled)
    step = (vec*np.exp(-.5j*k*spec["step"]*vals)) @ vec.T
    psi = np.eye(len(vals))[0].astype(complex)
    for row in theta:
        # At the donor's strength, theta is the accumulated phase; the same
        # physical duration at another strength accumulates theta*k/8.
        psi = np.exp(-1j*(k/8)*(row @ masks))*(step @ psi)
    target = [tree.index[tuple(p)] for p in choice["paths"]]
    dark = abs(psi[target[0]]-psi[target[1]])**2/2
    dest = [p[-1] for p in choice["paths"]]
    at_37 = vec @ (np.exp(-.5j*k*37*vals)*(vec.T @ psi))
    return {"dark_weight":round(float(dark),12),
            "two_port_return_37":round(float(sum(abs(at_37[np.isin(tree.sites,dest)])**2)),12)}
