"""Compare periodic and shortest native paths through the same 194 gateway pairs."""
import argparse
import json
from collections import Counter, deque
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.row_30_autonomous.channel import (
    column_checks, matrices, preparation,
)

HERE = Path(__file__).resolve().parent
SHOTS = 262_144


def shortest_internal(car, a, b):
    parent = {a: None}
    queue = deque([a])
    while queue:
        p = queue.popleft()
        if p == b:
            journey = []
            while p is not None:
                journey.append(p)
                p = parent[p]
            return journey[::-1]
        for other in sorted(car.neighbours[p]):
            if other//12 == a//12 and other not in parent:
                parent[other] = p
                queue.append(other)
    raise ArithmeticError("Disconnected cell gateways")


def portable(obj):
    if isinstance(obj, dict):
        return {str(k): portable(v) for k,v in obj.items()}
    if isinstance(obj, (tuple,list,np.ndarray)):
        return [portable(v) for v in obj]
    if isinstance(obj, (float,np.floating)):
        return round(float(obj), 10)
    if isinstance(obj, np.integer):
        return int(obj)
    return obj


def study(car, cells, route, seed):
    sites = np.concatenate(route)
    if len(set(sites)) != len(sites) or any(b not in car.neighbours[a] for a,b in zip(sites,sites[1:])):
        raise ArithmeticError("The native route visits a port twice or uses a non-neighbour")
    diag, off = matrices(car, sites)
    physical = diags((off,diag,off),(-1,0,1),format="csr")
    energies, waves = eigh_tridiagonal(diag,off)
    residual = np.linalg.norm(physical@waves - waves*energies, axis=0).max()
    inside = np.flatnonzero((energies > -5) & (energies < -3))
    if min(np.min(abs(energies-cut)) for cut in (-5,-3)) < 100*residual:
        raise ArithmeticError("Spectral count is ambiguous at an open-interval boundary")
    density = abs(waves)**2
    pathcenter = np.arange(len(sites))
    # Site centres strictly between the first and third quartiles in native path length.
    center = (4*pathcenter+2 >= len(sites)) & (4*pathcenter+2 < 3*len(sites))
    central_weight = density[center].sum(axis=0)
    end = np.isin(sites//12, [cells[0],cells[-1]])
    weights_at_ends = density[end].sum(axis=0)
    max_local = density.max(axis=0)
    if len(inside):
        local_projector = density[:,inside].sum(axis=1)
        spectral_weight = dict(central_trace=central_weight[inside].sum(),
                               central_trace_per_mode=central_weight[inside].mean(),
                               mean_central_local_density=local_projector[center].mean(),
                               maximum_at_port=local_projector.max(), port_of_maximum=int(sites[local_projector.argmax()]),
                               longest_central_fraction=float(central_weight[inside].max()),
                               interior_dominant_modes=int(np.sum(central_weight[inside]>.5)),
                               maximum_one_mode_one_port=max_local[inside].max(),
                               largest_end_cells_weight=weights_at_ends[inside].max())
    else:
        spectral_weight = dict(central_trace=0., central_trace_per_mode=None, mean_central_local_density=0.,
                               maximum_at_port=0.,port_of_maximum=None,longest_central_fraction=None,
                               interior_dominant_modes=0,maximum_one_mode_one_port=None,largest_end_cells_weight=None)
    def details(i):
        return dict(index=i, energy=energies[i], middle_half=central_weight[i], end_cells=weights_at_ends[i],
                    maximum_port_weight=max_local[i], maximum_port=sites[density[:,i].argmax()],
                    inverse_participation=np.dot(density[:,i],density[:,i]))
    below_center=np.flatnonzero(energies < -4)
    above_center=np.flatnonzero(energies >= -4)
    edge = dict(spectrum_floor=details(0),below_midpoint=details(below_center[-1]),
                at_or_above_midpoint=details(above_center[0]),
                below_test_interval=details(np.where(energies<=-5)[0][-1]),
                above_test_interval=details(np.where(energies>=-3)[0][0]))
    if len(inside):
        edge.update(test_low=details(inside[0]), test_high=details(inside[-1]),
                    largest_central=details(inside[np.argmax(central_weight[inside])]))
    k = len(cells)//2
    first = route[k][0]
    launch = np.zeros(len(sites))
    launch[np.where(sites==first)[0][0]]=1
    receiver = sites//12==cells[k+1]
    arrivals = []
    rng=np.random.default_rng(seed)
    for t in (.5,1.):
        psi=expm_multiply((-1j*t)*physical,launch)
        spectral = waves @ (np.exp(-1j*t*energies)*(waves.T @ launch))
        error=np.linalg.norm(psi-spectral)
        if error>1e-10:
            raise ArithmeticError("Direct evolution and spectral trajectory disagree")
        p=abs(psi)**2
        counts=rng.multinomial(SHOTS,p/p.sum())
        arrivals.append(dict(exposure=t, probability=p[receiver].sum(), received=counts[receiver].sum(),
                             quiet=counts[~receiver].sum(), per_port_counts=counts, solver_error=error))
    return dict(ports=sites, segment_lengths=[len(block) for block in route],
                degrees=car.degree[sites], stages=dict(total_ports=len(sites),middle_half_ports=sum(center),
                length_distribution=Counter(len(block) for block in route)), source_cell=cells[k],
                source_port=first, receiver_cell=cells[k+1], energies=energies, mode_central_weights=central_weight,
                mode_maximum_port_weights=max_local, inside_indices=inside, residual=residual, nearest_boundary_distance=min(
                    np.min(abs(energies+5)), np.min(abs(energies+3))), edge_states=edge,
                central_spectral_weight=spectral_weight, arrivals=arrivals)


def run():
    car,cells,periodic,_=preparation()
    original = periodic.reshape(-1,4).tolist()
    alternative = [shortest_internal(car,row[0],row[-1]) for row in original]
    return portable(dict(cells=cells, total_cells=car.cells, full_screen_seams=len(car.seams), shots=SHOTS,
                         local_checks=column_checks(), periodic=study(car,cells,original,970354),
                         shortest=study(car,cells,alternative,970355)))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    report=run()
    if args.output:
        with args.output.open("x") as handle:
            json.dump(report,handle,separators=(",",":"),allow_nan=False)
            handle.write("\n")
    elif report!=json.loads((HERE/"results.json").read_text()):
        raise SystemExit("Shortest-route sensitivity receipt differs from archive")
    for name in ("periodic","shortest"):
        row=report[name]
        print(name, row["stages"], len(row["inside_indices"]), row["central_spectral_weight"])
        print("edge", {k:v["energy"] for k,v in row["edge_states"].items()})
        print("arrival", [(a["exposure"],a["probability"]) for a in row["arrivals"]])
