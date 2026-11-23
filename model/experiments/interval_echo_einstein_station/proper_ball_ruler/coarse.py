"""Clocked blank range growth on pinned spherical refinements, excluding simulated coordinates."""
import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path, dijkstra

from model.executable_model import cell, federation
from model.executable_model.refinement import load_refinement
from ..shape_stationarity.geometry import selected_shapes
from .measure import pulse_calibration, quartet_bound

N_TRIALS = 2**32


class RangeReceiver:
    def __init__(self,car):
        self.cells = car.cells
        self.external = car.external
        self.inside = shortest_path(cell().adjacency(),directed=False,unweighted=True).astype(int)
        self.gates = np.sort(car.external.ravel()).reshape(-1,3)
        ids = np.arange(3*self.cells).reshape(-1,3)
        i,j = np.array(list(combinations(range(3),2))).T
        a,b = ids[:,i].ravel(),ids[:,j].ravel()
        costs = self.inside[self.gates.ravel()[a]%12,self.gates.ravel()[b]%12]
        out = np.searchsorted(self.gates.ravel(),self.external)
        u,v = out.T
        self.routes = coo_matrix((np.r_[costs,costs,np.ones(len(out)*2)],
                                  (np.r_[a,b,u,v],np.r_[b,a,v,u])),
                                 shape=(3*self.cells,3*self.cells)).tocsr()

    def ranges(self,source):
        ids=3*(source//12)+np.arange(3)
        seeds=self.inside[source%12,self.gates[source//12]%12]
        routes=dijkstra(self.routes,indices=ids, directed=True)
        best=np.min(routes+seeds[:,None],axis=0).astype(np.int32).reshape(-1,3)
        go = np.full((self.cells,12),np.iinfo(np.int32).max,dtype=np.int32)
        for g in range(3):
            step = self.inside[self.gates[:,g]%12]
            go = np.minimum(go,best[:,g,None]+step)
        go[source//12] = np.minimum(go[source//12],self.inside[source%12])
        return go.ravel().astype(np.int32)


def run(directory):
    calibration = pulse_calibration()
    panels=[]
    for j in range(4):
        level=2*j+1
        get=lambda i: load_refinement(directory/f'FEDERATION_LEVEL{i}.npz',level=i)
        car,parent=get(level),get(j)
        reader=RangeReceiver(car)
        reference=reader.ranges(0)
        radius=3*2**j
        totals=[int(np.count_nonzero(reference<=r)) for r in (radius,2*radius)]
        inner,outer = totals
        degree = float(np.log2(outer/inner))
        sphere_counts=np.bincount(reference)
        shots=np.random.default_rng([351730,j]).multinomial(N_TRIALS,
                                np.array([inner,outer-inner,len(reference)-outer])/len(reference))
        shell=int(sphere_counts[radius])
        is_inside = reference<=radius
        tags = is_inside[car.external]
        external_contacts = int(np.count_nonzero(tags[:,0]^tags[:,1]))
        allseams=np.array(cell().seams)
        cells=np.flatnonzero(np.any(is_inside.reshape(-1,12),axis=1))
        bits=is_inside.reshape(-1,12)[cells]
        cut_faces=bits[:,allseams[:,0]]^bits[:,allseams[:,1]]
        cut_contacts=external_contacts+int(cut_faces.sum())
        division=4**(j+1)
        shape_radius=[]
        selected = selected_shapes(parent)
        landmark_groups=[[int(12*(i*division+division-1)) for i in sorted(cells_)]
                         for _,cells_ in selected]
        fixed=set(i for group in landmark_groups for i in group)
        distances={}
        for node in sorted(fixed):
            dist=reader.ranges(node)
            if level==1:
                native=federation(1)
                direct=shortest_path(native.adjacency(),directed=False,unweighted=True,indices=node)
                np.testing.assert_array_equal(dist,direct)
            distances[node] = {q:int(dist[q]) for q in fixed}
        for idx,(name,region) in enumerate(selected):
            per_cell=reference.reshape(-1,12)
            members=np.concatenate([np.arange(c*division,(c+1)*division) for c in sorted(region)])
            reached=int(per_cell[members].max())
            q=landmark_groups[idx]
            macro=np.array([[distances[x][y] for y in q] for x in q])
            bound,w=quartet_bound(macro,np.array(list(combinations(range(len(q)),4))))
            shape_radius.append(dict(shape=name,microcells=len(members),range_radius=reached,
                                     relative_radius=float(reached/reference.max()),
                                     fraction_of_host_ports=len(members)/car.cells,
                                     four_parent_cell_distortion_bound=bound,
                                     witness_ports=[q[i] for i in w],landmark_ranges=macro.tolist()))
        local=np.flatnonzero(is_inside)
        farthest = reference[local].copy()
        points=[0]
        rows=[reference[local].tolist()]
        for _ in range(11):
            site = int(local[np.argmax(farthest)])
            points.append(site)
            ranges=reader.ranges(site)[local]
            farthest=np.minimum(farthest,ranges)
            rows.append(ranges.tolist())
        d = np.array(rows,dtype=int)[:,np.searchsorted(local,points)]
        np.testing.assert_equal(d,d.T)
        eps,idx=quartet_bound(d,np.array(list(combinations(range(12),4))))
        panels.append(dict(parent=j,level=level,cells=car.cells,pinned_sha256=car.sha256,
                           home_max_range=int(reference.max()), microscopic_ball_radius=radius,
                           relative_ball_radius=float(radius/reference.max()),
                           ball_count=inner,doubled_ball_count=outer,doubling_dimension=degree,
                           sharp_sample_dimension=float(np.log2((shots[0]+shots[1])/shots[0])),
                           launch_classes=shots.tolist(),range_front=shell,cut_seams=cut_contacts,
                           ball_landmarks=points,ball_landmark_ranges=d.tolist(),
                           ball_distortion_lower_bound=eps,ball_quartet=[points[i] for i in idx],
                           shapes=shape_radius))
    return dict(seed=351730,shots_per_host=N_TRIALS,swap_hand=calibration['decoded_single_step'],
                panels=panels)


if __name__=='__main__':
    args=argparse.ArgumentParser(description=__doc__)
    args.add_argument('--directory',type=Path,required=True)
    args.add_argument('--output',type=Path,required=True)
    options=args.parse_args()
    data = json.dumps(run(options.directory),indent=2,allow_nan=False)
    with options.output.open('x') as out:
        out.write(data+'\n')
