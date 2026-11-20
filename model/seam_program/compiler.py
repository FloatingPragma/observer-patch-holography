"""Compile slot amplitudes to local seam events with a stationary marker."""
from collections import deque
import math

import numpy as np

from .gate import require


class Compiler:
    def __init__(self, carrier, marker, shape):
        self.carrier = carrier
        self.marker = marker
        self.shape = tuple(shape)
        require(len(self.shape)==4 and all(type(n) is int and n>0 for n in self.shape),
                "four positive address sizes required")
        require(0<=marker<carrier.n, "marker outside the carrier")
        self.slots = [p for p in range(carrier.n) if p!=marker][:8*math.prod(shape)]
        require(len(self.slots)==8*math.prod(shape), "insufficient ports")
        self.seam = {pair:i for i,pair in enumerate(carrier.seams)}

    def exchange(self, a, b, quarter_turns):
        """Move one end up to its partner, apply the gate and restore the intermediates."""
        require(0<=a<self.carrier.n and 0<=b<self.carrier.n and self.marker not in (a,b)
                and quarter_turns in (1,4), "invalid seam instruction")
        if a==b:
            return []
        queue = deque([a])
        parent = {a:None}
        while b not in parent and queue:
            at = queue.popleft()
            for neighbour in self.carrier.neighbours[at]:
                if neighbour!=self.marker and neighbour not in parent:
                    parent[neighbour] = at
                    queue.append(neighbour)
        require(b in parent, "marker disconnects the requested ports")
        edges = []
        at = b
        while at!=a:
            prev = parent[at]
            edges.append(self.seam[tuple(sorted((at,prev)))])
            at = prev
        route = [(e,4) for e in reversed(edges[1:])]
        return route + [(edges[0],quarter_turns)] + route[::-1]

    def build(self, member):
        require(member in ("nominal","omit_mask_seven","shift_only"), "unknown gate program")
        instructions = []
        ports = np.asarray(self.slots).reshape((*self.shape,8))
        for address in np.ndindex(self.shape):
            block = ports[address]
            for mask in range(1,8):
                if member=="shift_only" or (mask==7 and member=="omit_mask_seven"):
                    continue
                for i in range(8):
                    if i<(i^mask):
                        instructions.extend(self.exchange(int(block[i]),int(block[i^mask]),1))
        # Every positive outgoing slot and its neighbour's negative slot form one transposition.
        for address in np.ndindex(self.shape):
            for axis in range(4):
                to = list(address)
                to[axis] = (to[axis]+1)%self.shape[axis]
                instructions.extend(self.exchange(int(ports[address][2*axis]),
                                                  int(ports[tuple(to)][2*axis+1]),4))
        return instructions

    def resources(self, instructions):
        count = np.zeros(self.carrier.n, dtype=int)
        depth = count.copy()
        for edge,q in instructions:
            a,b = self.carrier.seams[edge]
            require(self.marker not in (a,b), "marker used by an instruction")
            count[[a,b]] += 1
            depth[[a,b]] = 1+max(depth[a],depth[b])
        return dict(events=len(instructions),quarter_angle_units=sum(q for _,q in instructions),
                    depth=int(depth.max()),port_counts=count.tolist())

    def columns(self, instructions):
        """Exact invariant marked subspace in floating arithmetic; columns include all inputs."""
        state = np.zeros((self.carrier.n,len(self.slots)),complex)
        state[self.slots,np.arange(len(self.slots))] = 1
        for edge,q in instructions:
            a,b = self.carrier.seams[edge]
            pa,pb = state[[a,b]]
            z = np.exp(-1j*np.pi*q/4)
            difference = (1-z)*(pa-pb)/2
            state[a] = pa-difference
            state[b] = pb+difference
        return state
