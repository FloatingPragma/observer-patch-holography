"""Bounded-memory spectral refinement on the connected 56-cell gateway class."""

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components, shortest_path
from scipy.sparse.linalg import eigsh

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier, federation
from .classes import spatial_operators
from .experiment import HERE, SEED, dispersion, portable


def layout(car):
    graph = np.zeros((car.cells, car.cells))
    ports = defaultdict(list)
    for a, b in car.external:
        graph[a//12,b//12] = graph[b//12,a//12] = 1.
        for end in (a,b):
            ports[end//12].append(end%12)
    local = shortest_path(carrier.cell().adjacency(), directed=False, unweighted=True)
    groups = defaultdict(list)
    for cell, gateways in sorted(ports.items()):
        key = "".join(str(d) for d in sorted(int(local[a,b]) for i,a in enumerate(gateways) for b in gateways[i+1:]))
        groups[key].append(cell)
    counts = {}
    for key, cells in groups.items():
        _, labels = connected_components(graph[np.ix_(cells,cells)], directed=False)
        counts[key] = dict(cells=len(cells), component_sizes=sorted(Counter(map(int,labels)).values(), reverse=True))
    return graph, groups, counts


def lowest(matrix, count, seed, vectors=True):
    start = np.random.default_rng(seed).normal(size=matrix.shape[0])
    return eigsh(matrix, k=count, ncv=2*count+16, which="SA", tol=1e-12, v0=start,
                 return_eigenvectors=vectors)


def bracket(A, B, indices, energies, seed):
    boundary = abs(B).sum(axis=0).A1 > 0
    coupling2 = abs(B).sum(axis=0).max() * abs(B).sum(axis=1).max()
    limit = max(indices)
    # A single, conservative boundary penalty suffices for both tested cut energies.
    cuts = [(energies[k-1]+energies[k])/2 for k in indices]
    penalty = coupling2/min(abs(e) for e in cuts)
    damped = A-sparse.diags(penalty*boundary.astype(float))
    lower = np.sort(lowest(damped, limit+2, seed, vectors=False))
    return dict(boundary_configurations=int(sum(boundary)), lower_comparison_penalty=float(penalty),
                cuts=[dict(target_count=int(k), cut=float(e), ritz_count=int(np.sum(energies<e)),
                           penalized_count=int(np.sum(lower<e))) for k,e in zip(indices,cuts,strict=True)])


def analyze_radius(sector, graph, groups, h, distances, radius):
    inside = np.flatnonzero(distances <= radius)
    border = np.flatnonzero(distances == radius+1)
    A = h[inside][:,inside].tocsr()
    B = h[border][:,inside].tocsr()
    k = sector.car.cells
    energy, states = lowest(A,k+2,SEED+radius)
    # For distances > radius+1, a native single swap cannot touch the trial support.
    residuals = np.hypot(np.linalg.norm(A@states-states*energy,axis=0),
                         np.linalg.norm(B@states,axis=0))
    split = int(np.argmax(np.diff(energy[:k])))+1
    clusters = [slice(0,split), slice(split,k)]
    internal = sector.adjacent[inside] & (sector.aa[inside]//12 == sector.bb[inside]//12)
    rows = np.flatnonzero(internal)
    lift = sparse.csr_matrix((np.ones(len(rows))/np.sqrt(30),
                              (rows,sector.aa[inside[rows]]//12)), shape=(len(inside),k))
    # Compare all 80 orbitals with all 80 Ritz directions, without expanding to a full dense lift.
    overlaps = np.asarray(lift.T @ states[:,:k]).T
    weights = [float(np.sum(overlaps[band][:,groups["222"]]**2)/len(groups["222"])) for band in clusters]
    choice = int(np.argmax(weights))
    band, other = clusters[choice], clusters[1-choice]
    cells = groups["222"]
    complement = sorted(set(range(k))-set(cells))
    uniform_other = np.sum(overlaps[other][:,complement],axis=1)/np.sqrt(len(complement))
    other_rest = (uniform_other**2 @ energy[other]) / np.sum(uniform_other**2)
    panels = {}
    for metric, lap in spatial_operators(graph,cells).items():
        # dispersion needs states.T @ lift. All relevant overlaps were calculated on the full configuration basis.
        block = overlaps[band][:,cells]
        panel = dispersion(np.eye(len(cells)),np.linalg.eigh(lap),energy[band],block.T)
        centroid = panel["classes"][0]["centroid_above_bottom"] + min(energy[band])
        m = abs(centroid-other_rest)/2
        sign = np.sign(centroid-other_rest)
        first = panel["classes"][1]
        shift = sign*first["shift_from_rest"]
        speed2 = shift*(2*m+shift)/first["lambda_graph"]
        meaningful = m>0 and shift>0 and speed2>0
        panel.update(rest_half_gap=m, speed_squared=speed2, calibration_positive=bool(meaningful))
        for row in panel["classes"][1:]:
            row["shift_from_midgap"] = sign*row["shift_from_rest"]
            if meaningful:
                v = row["lambda_graph"]*speed2
                pred = v/(np.sqrt(m**2+v)+m)
                row["predicted_shift"] = pred
                row["relative_error"] = abs(pred-row["shift_from_midgap"])/pred
                row["tangent_separation"] = abs(v/(2*m)-pred)/pred
        panel["full_residual_over_first_shift"] = max(residuals[:k])/abs(shift)
        panels[metric] = panel
    result = dict(radius=radius, configurations=len(inside), ritz_bounds=[float(energy[0]),float(energy[k-1])],
                  next_gap=energy[k]-energy[k-1], internal_gap=energy[split]-energy[split-1],
                  split=split, upper_or_lower="upper" if choice else "lower",
                  mean_lower_upper_weights=weights, maximum_full_residual=max(residuals[:k]), panels=panels)
    if radius==7:
        result["eigenvalue_count_brackets"] = bracket(A,B,[split,k],energy,SEED+80)
    return result


def run():
    connectivity = {str(federation(level).cells): layout(federation(level))[2] for level in (0,1,2)}
    sector = PairSector(federation(1))
    graph, groups, _ = layout(sector.car)
    distances = shortest_path(sector.car.adjacency(), directed=False, unweighted=True)[sector.aa,sector.bb]
    h = sector.hamiltonian(candidate.KAPPA)
    results = [analyze_radius(sector,graph,groups,h,distances,r) for r in (6,7)]
    return portable(dict(candidate=candidate.NAME, full_binary_dimension=sector.P, active_seams=len(sector.car.seams),
                         carrier_sha256=carrier.seam_digest(sector.car), connectivity=connectivity, radii=results))


if __name__=="__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, help="write a separate receipt to a new path")
    path=cli.parse_args().output
    record=run()
    if path:
        with path.open("x") as f:
            json.dump(record,f,indent=2,allow_nan=False)
            f.write("\n")
    elif record != json.loads((HERE/"refinement.json").read_text()):
        raise SystemExit("free bound-pair refinement differs from the archive")
    for row in record["radii"]:
        print(row["radius"],row["maximum_full_residual"],
              [(k,v["narrow_at_every_momentum"]) for k,v in row["panels"].items()])
