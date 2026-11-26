"""Compile the supplied normal-mode analyzer into adjacent native material holds."""

import numpy as np

from ...row_20_electric import experiment as prior
from ...row_20_electric.schedule import rotate_pairs


def analyzer(g):
    # The two real axes of each single mode are orthonormal. Match the physical
    # orientation with one complemented read bit if needed.
    frame = np.empty_like(g.r)
    frame[:,::2], frame[:,1::2] = g.modes.real, g.modes.imag
    vacuum_word = np.ones(g.n, dtype=int)
    if np.linalg.slogdet(frame)[0] < 0:
        frame[:,-1] *= -1
        vacuum_word[-1] = -1

    target = frame.copy()
    reader = np.eye(2*g.n)
    j, h = prior.couplings(8)
    holds = []
    # Physical row rotations multiply the running transfer from the left.
    # Reduce the vacuum frame to identity. Thus reader = frame.T.
    for col in range(len(frame)-1):
        for i in range(len(frame)-2,col-1,-1):
            x, y = target[i:i+2,col]
            if abs(x)+abs(y) < 1e-14:
                continue
            angle = np.arctan2(y,x)
            rate = 2*j if i%2 else -2*h
            dwell = (np.sign(rate)*angle % (2*np.pi))/abs(rate)
            rotate_pairs(target,np.array([i]),rate*dwell)
            rotate_pairs(reader,np.array([i]),rate*dwell)
            holds.append((i,dwell))
    np.testing.assert_allclose(target,np.eye(len(frame)),atol=3e-12)
    diagonal = reader @ g.c @ reader.T
    np.testing.assert_allclose(diagonal[::2,1::2],np.diag(vacuum_word),atol=3e-12)
    return reader, vacuum_word, holds
