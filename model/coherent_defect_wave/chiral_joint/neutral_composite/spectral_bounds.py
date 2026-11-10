"""The inherited rational projector bound with exact nearest-interval separation."""
from functools import lru_cache
import sympy as s
from ..periodic_spectrum.readout import approximation,clusters
from ..two_cell.readout import sqrt_upper


@lru_cache(maxsize=None)
def projector_errors(block):
    _,energies,residuals,delta=approximation(block)
    groups,intervals=clusters(block);n=len(intervals)
    assert len(energies)==n
    assert all(intervals[i][0]<=energies[i]<=intervals[i][1] for i in range(n))
    assert all(intervals[i][0]<=intervals[i+1][0] and intervals[i][1]<=intervals[i+1][1] for i in range(n-1))
    assert [i for indices,_,_ in groups for i in indices]==list(range(n))
    result={}
    for indices,lo,hi in groups:
        first,last=indices[0],indices[-1]
        assert indices==list(range(first,last+1))
        assert lo==intervals[first][0] and hi==intervals[last][1]
        if len(indices)==n:result[tuple(indices)]=delta;continue
        boundaries=[]
        if first:
            assert intervals[first-1][1]<lo
            boundaries.append(min(energies[i] for i in indices)-intervals[first-1][1])
        if last+1<n:
            assert hi<intervals[last+1][0]
            boundaries.append(intervals[last+1][0]-max(energies[i] for i in indices))
        separation=min(boundaries);assert separation>0
        sine=sqrt_upper(sum(residuals[i] for i in indices)/(1-delta))/separation
        result[tuple(indices)]=min(s.Integer(1),sine)+delta
    return result


def projector_error(block,indices):return projector_errors(tuple(block))[tuple(indices)]
