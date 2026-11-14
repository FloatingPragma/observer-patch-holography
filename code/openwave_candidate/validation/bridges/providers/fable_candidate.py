#!/usr/bin/env python3
"""candidate: reference implementation of the common candidate, one-load sector.

The candidate (identifier `ring-hm-mod6-local-coarse`, see ../CANDIDATE.md) on any graph (ports, seams):

  positions      ports 0..P-1 and transit seams P..P+S-1
  registers      one integer F[e,side] per (seam, endpoint); side 0 = first endpoint, 1 = second
  half-hops      port x -> seam e raises F[e,x] by one, seam e -> port x lowers it;
                 amplitude +1 at the first endpoint of the seam, -1 at the second
  generator      H_m = D + m Gamma, Gamma = +1 on ports, -1 on transit seams (m supplied)
  Gauss charges  G_x = [z = x] + sum_e F[e,x],  G_e = [z = e] - F[e,a] - F[e,b]   (constant on a sector)
  native event   seam e = (a, b) and a coin drawn uniformly; coin 0 moves a unit a -> b, coin 1 moves b -> a;
                 a transfer a -> b adds one to F[e,a] and removes one from F[e,b]

Three exact representations of the same model:

  Joint    positions x registers.  Regulator: `modulus=N` (registers in Z_N, no wall) or `cutoff=K` (|F| <= K,
           hard wall, the implementer's regulator).  Preparation: definite records, or the uniform fibre
           (zero-holonomy sector; modulus only).
  Gauged   one holonomy sector: positions only, half-hop phases exp(-2 pi i q[e,side]/N).  The joint model with
           definite records is the uniform average of these sectors over the holonomy group.
  Bare     the sector q = 0: positions only, the matrix H_m of the report.

Commit instruments (all act after the coherent interval):

  'MAX'   position and every register are measured, then the native jump is applied (implementer's instrument;
          supplied: the complete measurement, the transit no-op, the wall rejection)
  'C'     native C record law on the drawn seam: coin 1 asks "is the load at b", coin 0 "is it at a"
  'E'     native E record law on the drawn seam: resolves {a, b, elsewhere}
  'MIX'   C and E with equal weight

Protocol: rho' = (1 - p) U rho U* + p Commit(U rho U*), U = exp(-i t H).  Commit rate Gamma = p / t, rate per
seam Gamma / S.  The protocol is a first-order splitting of d rho/dt = -i[H, rho] + Gamma (Commit(rho) - rho).

Self-test: `python candidate.py`.  numpy only.  Nothing here writes outside the caller's directory.
"""
from __future__ import annotations

import itertools
import json
import math
from collections import deque
from pathlib import Path

import numpy as np

CANDIDATE_ID = "ring-hm-mod6-local-coarse"
VERSION = "F-1.0"


# ----------------------------------------------------------------------------------------------------------
# joint model
# ----------------------------------------------------------------------------------------------------------
class Joint:
    """Reachable Gauss sector of one load started on `start` with all registers zero."""

    def __init__(self, ports, seams, m=0.0, start=0, modulus=None, cutoff=None):
        if (modulus is None) == (cutoff is None):
            raise ValueError("give exactly one regulator: modulus=N or cutoff=K")
        self.P, self.seams, self.S = int(ports), tuple(tuple(e) for e in seams), len(seams)
        self.m, self.N, self.K = float(m), modulus, cutoff
        s0 = (start, (0,) * (2 * self.S))
        states, seen, queue = [s0], {s0}, deque([s0])
        while queue:
            st = queue.popleft()
            for tgt, _ in self.half_hops(st):
                if tgt not in seen:
                    seen.add(tgt); states.append(tgt); queue.append(tgt)
        self.states = states
        self.index = {s: i for i, s in enumerate(states)}
        self.n = len(states)
        self.pos = np.array([s[0] for s in states])

    # register arithmetic under the chosen regulator; returns None when the hard wall forbids the value
    def _shift(self, value, delta):
        v = value + delta
        if self.N is not None:
            return v % self.N
        return v if abs(v) <= self.K else None

    def half_hops(self, st):
        pos, f = st
        for e, ends in enumerate(self.seams):
            for side, port in enumerate(ends):
                if pos == port:
                    dest, delta = self.P + e, 1
                elif pos == self.P + e:
                    dest, delta = port, -1
                else:
                    continue
                v = self._shift(f[2 * e + side], delta)
                if v is None:
                    continue
                g = list(f); g[2 * e + side] = v
                yield (dest, tuple(g)), (1 if side == 0 else -1)

    def gauss(self, st):
        pos, f = st
        g = [int(pos == x) for x in range(self.P + self.S)]
        for e, (a, b) in enumerate(self.seams):
            g[a] += f[2 * e]; g[b] += f[2 * e + 1]; g[self.P + e] -= f[2 * e] + f[2 * e + 1]
        return tuple(x % self.N for x in g) if self.N is not None else tuple(g)

    def hamiltonian(self):
        h = np.zeros((self.n, self.n))
        for c, st in enumerate(self.states):
            h[c, c] = self.m * (1 if st[0] < self.P else -1)
            for tgt, sgn in self.half_hops(st):
                h[self.index[tgt], c] += sgn
        return h

    def move(self, st, e, direction):
        """Native transfer across seam e (direction +1: first -> second endpoint).  None when not applicable
        or rejected by the wall."""
        pos, f = st
        a, b = self.seams[e]
        if (direction == 1 and pos != a) or (direction == -1 and pos != b):
            return None
        va, vb = self._shift(f[2 * e], direction), self._shift(f[2 * e + 1], -direction)
        if va is None or vb is None:
            return None
        g = list(f); g[2 * e], g[2 * e + 1] = va, vb
        return (b if direction == 1 else a, tuple(g))

    def _moved(self, idx, e, direction):
        """index map of the native transfer; a rejected or inapplicable state stays where it is"""
        out = []
        for i in idx:
            tgt = self.move(self.states[i], e, direction)
            out.append(i if tgt is None else self.index[tgt])
        return np.array(out, dtype=int)

    def native_maps(self):
        """Kraus operators of the native record laws as (weight, source indices, target indices)."""
        idx = np.arange(self.n)
        maps = {"C": [], "E": []}
        w = 1.0 / (2 * self.S)
        for e, (a, b) in enumerate(self.seams):
            at_a, at_b = idx[self.pos == a], idx[self.pos == b]
            els = idx[(self.pos != a) & (self.pos != b)]
            ab, ba = self._moved(at_a, e, 1), self._moved(at_b, e, -1)
            maps["E"] += [(w, at_a, at_a), (w, at_b, ba), (w, els, els),      # coin 1
                          (w, at_a, ab), (w, at_b, at_b), (w, els, els)]      # coin 0
            not_b, not_a = idx[self.pos != b], idx[self.pos != a]
            maps["C"] += [(w, not_b, not_b), (w, at_b, ba), (w, not_a, not_a), (w, at_a, ab)]
        return maps

    def max_matrix(self):
        """Column-stochastic matrix of the MAX instrument on the diagonal (measure everything, native jump)."""
        t = np.zeros((self.n, self.n))
        w = 1.0 / (2 * self.S)
        for i, st in enumerate(self.states):
            for e, (a, b) in enumerate(self.seams):
                for coin in (0, 1):
                    direction = 1 if (st[0] == a and coin == 0) else (-1 if (st[0] == b and coin == 1) else 0)
                    tgt = self.move(st, e, direction) if direction else None
                    t[i if tgt is None else self.index[tgt], i] += w
        return t

    # preparations -------------------------------------------------------------------------------------
    def definite(self, state_index=0):
        rho = np.zeros((self.n, self.n), complex); rho[state_index, state_index] = 1
        return rho

    def uniform_fibre(self, position):
        """Zero-holonomy preparation at `position`: uniform superposition of every record state of the sector
        with the load there.  Normalisable only with a modulus."""
        if self.N is None:
            raise ValueError("the uniform fibre needs registers modulo N")
        v = (self.pos == position).astype(complex)
        v /= np.linalg.norm(v)
        return np.outer(v, v.conj())

    def intertwiner(self):
        """V: bare positions -> uniform fibres.  H_joint V = V H_bare holds with a modulus."""
        v = np.zeros((self.n, self.P + self.S))
        for z in range(self.P + self.S):
            col = (self.pos == z).astype(float)
            v[:, z] = col / np.linalg.norm(col)
        return v

    def position_probabilities(self, rho):
        d = np.real(np.diag(rho))
        return np.array([d[self.pos == z].sum() for z in range(self.P + self.S)])


# ----------------------------------------------------------------------------------------------------------
# one holonomy sector (positions only)
# ----------------------------------------------------------------------------------------------------------
class Gauged:
    """Positions only.  `phase[e][side]` is the angle theta of the half-hop port -> seam (amplitude
    sign * exp(-i theta)); for the sector q of registers modulo N, theta = 2 pi q[e,side] / N.
    phase=None is the bare model."""

    def __init__(self, ports, seams, m=0.0, phase=None):
        self.P, self.seams, self.S, self.m = int(ports), tuple(tuple(e) for e in seams), len(seams), float(m)
        self.n = self.P + self.S
        self.pos = np.arange(self.n)
        self.phase = np.zeros((self.S, 2)) if phase is None else np.asarray(phase, float).reshape(self.S, 2)

    def hamiltonian(self):
        h = np.zeros((self.n, self.n), complex)
        for e, ends in enumerate(self.seams):
            for side, port in enumerate(ends):
                amp = (1 if side == 0 else -1) * np.exp(-1j * self.phase[e, side])
                h[self.P + e, port] += amp
                h[port, self.P + e] += np.conj(amp)
        h += np.diag([self.m] * self.P + [-self.m] * self.S)
        return h

    def native_maps(self):
        maps = {"C": [], "E": []}
        idx = np.arange(self.n)
        w = 1.0 / (2 * self.S)
        for e, (a, b) in enumerate(self.seams):
            at_a, at_b = np.array([a]), np.array([b])
            els = idx[(idx != a) & (idx != b)]
            maps["E"] += [(w, at_a, at_a), (w, at_b, at_a), (w, els, els),
                          (w, at_a, at_b), (w, at_b, at_b), (w, els, els)]
            maps["C"] += [(w, idx[idx != b], idx[idx != b]), (w, at_b, at_a),
                          (w, idx[idx != a], idx[idx != a]), (w, at_a, at_b)]
        return maps

    def max_matrix(self):
        t = np.zeros((self.n, self.n))
        w = 1.0 / (2 * self.S)
        for z in range(self.n):
            for e, (a, b) in enumerate(self.seams):
                for coin in (0, 1):
                    tgt = b if (z == a and coin == 0) else (a if (z == b and coin == 1) else z)
                    t[tgt, z] += w
        return t

    def definite(self, position=0):
        rho = np.zeros((self.n, self.n), complex); rho[position, position] = 1
        return rho

    def position_probabilities(self, rho):
        return np.real(np.diag(rho)).copy()


# ----------------------------------------------------------------------------------------------------------
# protocol
# ----------------------------------------------------------------------------------------------------------
def unitary(h, t):
    w, v = np.linalg.eigh(h)
    return (v * np.exp(-1j * t * w)) @ v.conj().T


def commit(model, rho, instrument, cache=None):
    """One application of the commit channel."""
    cache = {} if cache is None else cache
    if instrument == "MAX":
        if "max" not in cache:
            cache["max"] = model.max_matrix()
        return np.diag(cache["max"] @ np.real(np.diag(rho))).astype(complex)
    if "maps" not in cache:
        cache["maps"] = model.native_maps()
    laws = {"C": ("C",), "E": ("E",), "MIX": ("C", "E")}[instrument]
    out = np.zeros_like(rho)
    for law in laws:
        for w, src, tgt in cache["maps"][law]:
            if len(src):
                out[np.ix_(tgt, tgt)] += (w / len(laws)) * rho[np.ix_(src, src)]
    return out


def run(model, rho, t, p, steps, instrument="MIX", observe=None):
    """Apply the protocol `steps` times.  Returns (final rho, list of observe(rho) after each step)."""
    u = unitary(model.hamiltonian(), t)
    cache, history = {}, []
    for _ in range(steps):
        ev = u @ rho @ u.conj().T
        rho = (1 - p) * ev + p * commit(model, ev, instrument, cache) if p else ev
        if observe is not None:
            history.append(observe(rho))
    return rho, history


def holonomy_sectors(ports, seams, modulus):
    """One register-conjugate label per holonomy class: q supported on the first half-seam of the seams outside
    a spanning forest.  Returns a list of phase arrays (S x 2)."""
    parent = list(range(ports))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    chords = []
    for e, (a, b) in enumerate(seams):
        ra, rb = find(a), find(b)
        if ra == rb:
            chords.append(e)
        else:
            parent[ra] = rb
    out = []
    for q in itertools.product(range(modulus), repeat=len(chords)):
        ph = np.zeros((len(seams), 2))
        for e, v in zip(chords, q):
            ph[e, 0] = 2 * math.pi * v / modulus
        out.append(ph)
    return out


# ----------------------------------------------------------------------------------------------------------
# graphs used by the self-test and by row experiments
# ----------------------------------------------------------------------------------------------------------
def cycle(k):
    return k, [(i, (i + 1) % k) for i in range(k)]


def dumbbell(k):
    """two k-cycles joined by a bridge between port 0 of each"""
    seams = [(i, (i + 1) % k) for i in range(k)] + [(k + i, k + (i + 1) % k) for i in range(k)] + [(0, k)]
    return 2 * k, seams


TRIANGLE = (3, ((0, 1), (1, 2), (2, 0)))          # the implementer's oriented triangle
DECLARED = dict(m=0.23, t=0.37, p=0.29, steps=8)  # the implementer's declared run
EXACT_TRIANGLE = (0.148516560376, 0.173545469750, 0.173545469750,
                  0.169718954539, 0.164954591047, 0.169718954539)   # K01: cutoff-free, instrument MAX


# ----------------------------------------------------------------------------------------------------------
def selftest(verbose=True):
    checks = []

    def ok(name, cond, detail=""):
        checks.append((name, bool(cond)))
        if verbose:
            print(("PASS " if cond else "FAIL ") + name + (" | " + detail if detail else ""))

    P, seams = TRIANGLE
    m, t, p, steps = (DECLARED[k] for k in ("m", "t", "p", "steps"))

    # 1. hard wall, instrument MAX: the implementer's numbers as stored by K01 (which replayed his code)
    ref = Path(__file__).resolve().parents[1] / "issues" / "K01_joint_instrument" / "results_cutoff.json"
    stored = {row["cutoff"]: row["position_probabilities"] for row in json.loads(ref.read_text())["hard_wall"]}
    for K in (1, 2):
        j = Joint(P, seams, m, 0, cutoff=K)
        rho, _ = run(j, j.definite(), t, p, steps, "MAX")
        err = float(np.max(np.abs(j.position_probabilities(rho) - np.array(stored[K]))))
        ok(f"hard wall K={K}, MAX: equals the stored implementer run", err < 1e-12, f"states {j.n}, err {err:.1e}")

    # 2. registers modulo 6, MAX: the cutoff-free probabilities
    j = Joint(P, seams, m, 0, modulus=6)
    rho, _ = run(j, j.definite(), t, p, steps, "MAX")
    err = float(np.max(np.abs(j.position_probabilities(rho) - np.array(EXACT_TRIANGLE))))
    ok("modulus 6, MAX: equals the cutoff-free probabilities", err < 1e-10, f"err {err:.1e}")
    ok("trace", abs(np.trace(rho).real - 1) < 1e-12)

    # 3. Gauss charges constant on the sector; Kraus completeness of the native laws
    g0 = j.gauss(j.states[0])
    ok("Gauss charges constant on the sector", all(j.gauss(s) == g0 for s in j.states))
    maps = j.native_maps()
    for law in ("C", "E"):
        cover = np.zeros(j.n)
        for w, src, tgt in maps[law]:
            cover[src] += w
            assert len(set(tgt.tolist())) == len(tgt)        # each Kraus operator is a partial isometry
        ok(f"native law {law}: sum K*K = 1", np.allclose(cover, 1))

    # 4. definite records = uniform average over holonomy sectors (dumbbell, two cycles, modulus 4)
    P2, seams2 = dumbbell(3)
    #    The native laws shift registers without reading them, so they conserve every holonomy sector.  MAX reads
    #    every register, which resets the holonomy distribution: it is not a sector-wise operation.
    for inst in ("MIX", "C", "E", "MAX"):
        j2 = Joint(P2, seams2, m, 0, modulus=4)
        rho, _ = run(j2, j2.definite(), t, p, 6, inst)
        direct = j2.position_probabilities(rho)
        sectors = holonomy_sectors(P2, seams2, 4)
        avg = np.zeros(P2 + len(seams2))
        for ph in sectors:
            g = Gauged(P2, seams2, m, ph)
            r, _ = run(g, g.definite(0), t, p, 6, inst)
            avg += g.position_probabilities(r) / len(sectors)
        err = float(np.max(np.abs(direct - avg)))
        if inst == "MAX":
            ok("MAX is not sector-wise: definite records differ from the sector average", err > 1e-5, f"diff {err:.1e}")
        else:
            ok(f"definite records = average over {len(sectors)} holonomy sectors ({inst})", err < 1e-12, f"err {err:.1e}")

    # 5. uniform fibre = bare model (zero-holonomy sector), with every instrument
    j2 = Joint(P2, seams2, m, 0, modulus=4)
    v = j2.intertwiner()
    bare = Gauged(P2, seams2, m)
    err = float(np.max(np.abs(j2.hamiltonian() @ v - v @ bare.hamiltonian().real)))
    ok("intertwiner H_joint V = V H_bare (modulus 4)", err < 1e-12, f"err {err:.1e}")
    for inst in ("MIX", "C", "E", "MAX"):
        rho, _ = run(j2, j2.uniform_fibre(0), t, p, 6, inst)
        r, _ = run(bare, bare.definite(0), t, p, 6, inst)
        err = float(np.max(np.abs(j2.position_probabilities(rho) - bare.position_probabilities(r))))
        if inst == "MAX":
            ok("MAX leaves the zero-holonomy sector (uniform fibre differs from the bare model)", err > 1e-4, f"diff {err:.1e}")
        else:
            ok(f"uniform fibre = bare model under the protocol ({inst})", err < 1e-12, f"err {err:.1e}")
    #    weight of the zero-holonomy sector after one commit: 1 under the native laws, 1/|sectors| under MAX
    vfib = j2.intertwiner()
    proj = vfib @ vfib.T
    u = unitary(j2.hamiltonian(), t)
    ev = u @ j2.uniform_fibre(0) @ u.conj().T
    w_native = float(np.trace(proj @ commit(j2, ev, "MIX")).real)
    w_max = float(np.trace(proj @ commit(j2, ev, "MAX")).real)
    ok("zero-holonomy weight after one commit: native 1", abs(w_native - 1) < 1e-12, f"{w_native:.12f}")
    ok("zero-holonomy weight after one commit: MAX 1/16", abs(w_max - 1 / 16) < 1e-12, f"{w_max:.6f}")

    # 6. definite and zero-holonomy preparations differ on a graph with cycles
    rho_d, _ = run(j2, j2.definite(), t, 0.0, 6, "MIX")
    rho_u, _ = run(j2, j2.uniform_fibre(0), t, 0.0, 6, "MIX")
    diff = float(np.max(np.abs(j2.position_probabilities(rho_d) - j2.position_probabilities(rho_u))))
    ok("definite and zero-holonomy preparations give different position laws", diff > 1e-3, f"max diff {diff:.4f}")

    # 7. H_m^2 = m^2 + D^2 and the port block of D^2 on the bare model is the graph Laplacian
    h = bare.hamiltonian().real
    d = h - np.diag(np.diag(h))
    lap = np.zeros((P2, P2))
    for a, b in seams2:
        lap[a, a] += 1; lap[b, b] += 1; lap[a, b] -= 1; lap[b, a] -= 1
    ok("H_m^2 = m^2 + D^2", np.allclose(h @ h, m * m * np.eye(len(h)) + d @ d))
    ok("port block of D^2 is the graph Laplacian", np.allclose((d @ d)[:P2, :P2], lap))

    passed = sum(c for _, c in checks)
    if verbose:
        print(f"{passed}/{len(checks)} checks passed")
    return passed == len(checks)


if __name__ == "__main__":
    raise SystemExit(0 if selftest() else 1)
