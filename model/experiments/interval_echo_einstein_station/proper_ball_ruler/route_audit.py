"""Replay selected longest regional legs with native addressed swaps and echo read/undo."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo.clock import Registers
from model.executable_model import federation
from ..equilibrium_refinement.measure import apparatus
from ..shape_stationarity.geometry import selected_shapes
from .measure import pulse_calibration


def run():
    pulse = pulse_calibration()
    experiments=[]
    for level in range(3):
        car=federation(level)
        regs=Registers(car,seed=351731)
        for name,cells in selected_shapes(car):
            ports,_,_,words=apparatus(car,cells)
            for index in (0,len(words)//2,len(words)-1):
                resident=words[index]
                active=np.setdiff1d(np.arange(car.n),resident)
                region=np.setdiff1d(ports,resident)
                station=int(region[0])
                graph=car.adjacency().tocsr()[active][:,active]
                distances,parents=shortest_path(graph,method='D',directed=False,unweighted=True,
                                                indices=np.searchsorted(active,station),
                                                return_predecessors=True)
                radius=distances[np.searchsorted(active,region)]
                target=int(region[radius.argmax()])
                node=np.searchsorted(active,target)
                path=[target]
                while path[-1]!=station:
                    node=parents[node]
                    path.append(int(active[node]))
                path.reverse()
                regs.loads.fill(0)
                regs.loads[list(resident)]=2
                regs.loads[station]=1
                hops=list(zip(path,path[1:]))
                pointers=[]
                for legs,dest in ((hops,target),([(b,a) for a,b in reversed(hops)],station)):
                    for a,b in legs:
                        regs.swap(a,b,{a//12,b//12})
                    auxiliary=next(p for p in car.neighbours[dest] if regs.loads[p]==0)
                    regs.loads[auxiliary]=1
                    pointers.append(regs.flag(dest,auxiliary,{dest//12,auxiliary//12}))
                    regs.loads[auxiliary]=0
                if pointers!=[True,True] or tuple(np.flatnonzero(regs.loads==2))!=resident:
                    raise ArithmeticError('Native tag or resident failed round-trip audit')
                experiments.append(dict(level=level,shape=name,resident=list(resident),routed_ports=path,
                                        terminal_tag_reads=list(map(bool,pointers)),
                                        round_trip_hand=float(2*len(hops)*pulse['decoded_single_step'])))
    return dict(replays=experiments)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--output',type=Path,required=True)
    target=cli.parse_args().output
    result=json.dumps(run(),indent=2)
    with target.open('x') as out:
        out.write(result+'\n')
