"""Explicit path toggles and pointer projections, independent of the spectral shortcut."""

import numpy as np
import pytest
from scipy.linalg import eigh, expm

from ...row_08_annihilation import experiment
from . import measure as clock, routing


@pytest.mark.parametrize("pin", [0,2])
def test_control_target_memory_circuit(pin):
    h = experiment.wing(8,3,pin)
    e,v = eigh(h)
    rng = np.random.default_rng(310319)
    psi = rng.normal(size=8) + 1j*rng.normal(size=8)
    psi /= np.linalg.norm(psi)
    p1 = v[:,1:4] @ v[:,1:4].T
    flags = (clock.IDENTITY-p1,p1)
    axes = (np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]))
    for ref, matched, dephased in [(0,False,False),(2,False,False),(0,True,False),(0,False,True)]:
        reduced = routing.phase_bank(e,abs(v.T @ psi)**2,ref,matched,dephased)
        for i in (0,11,102,311):
            u = expm(-1j*routing.TIME[i]*(h-e[0]*clock.IDENTITY))
            m = u if matched else clock.IDENTITY
            for tag,q in enumerate(flags):
                # Real-space target and memory registers; swap axes on arm one, wait, undo.
                start = np.outer(q @ psi, v[:,ref])
                arm = np.stack([u @ start @ m.T,(u @ start.T @ m.T).T])/np.sqrt(2)
                for axis, pointer in enumerate(axes):
                    for j,s in enumerate(clock.SIGNS):
                        p = (np.eye(2)+s*pointer)/2
                        for end,terminal in enumerate(flags):
                            branches = np.stack([terminal @ a for a in arm])
                            if dephased:
                                weight = sum(np.linalg.norm(p[:,k,None,None]*branches[k])**2
                                             for k in (0,1))
                            else:
                                weight = np.linalg.norm(np.einsum("ij,jab->iab",p,branches))**2
                            np.testing.assert_allclose(weight,reduced[axis,i,tag,j,end],atol=1e-12)


def test_count_only_routing_decode():
    with np.load(clock.HERE/"routing_counts.npz",allow_pickle=False) as banks:
        routing.validate(routing.decode(banks))
