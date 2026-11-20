"""A reversible travelling query for N consecutive binary loads and its public position bins."""
import numpy as np

from .. import reader
from .source import require


def query(window, width, scratch=0, pointer=0):
    """Collect one controlled bit per visit, write the tip pointer, and retrace the path."""
    for visit in range(width):
        scratch ^= window & (1 << visit)
    pointer ^= int(scratch == (1 << width)-1)
    for visit in reversed(range(width)):
        scratch ^= window & (1 << visit)
    return scratch,pointer


def position_bins(ring, control=None):
    """Compile the two-owner public read on every computational state of this sector."""
    if control is None:
        control = {}
    n = ring.number
    mask = (1 << n)-1
    truth = []
    for word in range(mask+1):
        restored, flag = query(word,n)
        require(restored==0,'query kept a private copy')
        public = reader.gate(flag,0,**control)
        require((public>>6)&3 == 1,'unavailable or disputed cluster position')
        truth.append(((public>>8)&7)==7)
    states = np.array(ring.states,dtype=np.int64)
    cyclic = states | (states << ring.size)
    bins = np.full(len(states),ring.size)
    for start in range(ring.size):
        accepted = np.take(truth,(cyclic>>start)&mask)
        require(np.all(bins[accepted]==ring.size),'multiple public cluster locations')
        bins[accepted] = start
    return bins
