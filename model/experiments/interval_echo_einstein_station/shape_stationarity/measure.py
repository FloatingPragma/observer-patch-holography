"""Equal-launch-volume regional shape changes and held stationary matter variations."""

import argparse
from collections import defaultdict
import json
from math import comb
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix

from model.executable_model import federation
from model.candidate_laws.interval_echo import candidate
from ..equilibrium_refinement.measure import apparatus, matched_energy_and_load, seam
from .geometry import geometry

SHOTS = 2**22
SEED = 616405


def partition(words, cut):
    inside_words, cut_words = {}, {}
    inside_ids, cut_ids, k = [], [], []
    for word in words:
        inside = tuple(v for v in word if v not in cut)
        exterior = tuple(v for v in word if v in cut)
        inside_ids.append(inside_words.setdefault(inside, len(inside_words)))
        cut_ids.append(cut_words.setdefault(exterior, len(cut_words)))
        k.append(len(exterior))
    return (np.array(inside_ids), np.array(cut_ids), np.array(k),
            list(inside_words), list(cut_words))


def native_sector(car, contacts, words):
    """Enumerate the actual two-load energy and certify native configuration connectivity."""
    block = seam()
    index = {w: i for i, w in enumerate(words)}
    energy = np.zeros(len(words))
    source, target = [], []
    for i, word in enumerate(words):
        energy[i] = 4*candidate.KAPPA*sum(6 - int(car.degree[p]) for p in word)
        for a, b in contacts:
            left, right = a in word, b in word
            state = 2*left+right
            energy[i] += block[state, state]
            if left != right and block[2*right+left, state]:
                source.append(i)
                target.append(index[tuple(sorted(set(word) ^ {a, b}))])
    graph = coo_matrix((np.ones(len(source)), (source, target)), shape=(len(words), len(words))).tocsr()
    if connected_components(graph, directed=False)[0] != 1:
        raise ValueError("Two-load native sharp-word graph is disconnected")
    if (graph != graph.T).nnz:
        raise ValueError("Native hopping graph breaks reciprocity")
    return energy, int(graph.nnz//2)


def read_counts(bulk_counts, edge_counts, interior_words, exterior_words, matched):
    b = np.asarray(bulk_counts, dtype=float)
    b /= b.sum()
    q = np.asarray(edge_counts, dtype=float)
    q /= q.sum()
    if min(b) <= 0 or min(q) <= 0:
        raise ArithmeticError("Meter support missing from acquisition")
    info = -np.log(b)
    i_k = np.array([2-len(w) for w in interior_words])
    e_k = np.array([len(w) for w in exterior_words])
    by_k, cut_law = [], []
    for j in range(3):
        block = q[e_k == j]
        cut_law.append(float(block.sum()))
        prob = block/block.sum()
        by_k.append(float(-prob @ np.log2(prob)))
    t = np.array(matched["sector_transfers"])
    means = np.array([info[i_k == j].mean() for j in range(3)])
    slope = float(t @ means + matched["compensation"]*(info[matched["high_group"]]-info[matched["low_group"]]))
    return dict(bulk_nats=float(b @ info), cut_bits=float(np.dot(cut_law, by_k)),
                strict_bulk_slope=slope, strict_cut_slope=float(t @ by_k))


def resident_tag_strength():
    law = candidate.native_law
    h, echo = law.seam_term(candidate.KAPPA), law.echo()
    pointer = np.arange(law.DIM) % len(law.SRC) == law.SRC.index((1, 1))
    for charge in (0, 1, 2):
        basis = np.eye(law.DIM)[:, law.SRC.index((charge, 1))]
        written = echo @ basis
        np.testing.assert_allclose(np.linalg.norm(written[pointer])**2, int(charge==1), atol=1e-14)
        np.testing.assert_allclose(echo @ written, basis, atol=1e-14)
    launch = law.SRC.index((1, 0))
    target = law.SRC.index((0, 1))
    column = h[:, launch].copy()
    amplitude = column[target]
    column[[launch, target]] = 0
    if np.any(column):
        raise ValueError("Tagged crossing leaks outside its native sharp seam")
    return float(4*abs(amplitude)**2)


def source(row, salt, support=None):
    car = federation(row["board_level"])
    ports, contacts, cut, words = apparatus(car, row["cells"])
    if sorted(cut) != sorted(row["cut_ports"]) or len(ports) != row["accessible_volume"]:
        raise ArithmeticError("Source geometry disagrees with independently observed addresses or volume")
    if support is not None:
        _, contacts, _, words = apparatus(car, support)
    inside_ids, cut_ids, k, interior_words, exterior_words = partition(words, cut)
    energies, transitions = native_sector(car, contacts, words)
    matched = matched_energy_and_load(energies, k, inside_ids)
    inside_mult = np.bincount(inside_ids)
    active_border = set(cut).intersection(set(v for word in words for v in word))
    for i, word in enumerate(interior_words):
        if inside_mult[i] != comb(len(active_border), 2-len(word)):
            raise ArithmeticError("Sharp completions disagree with observed cut size")
    edge_mult = np.bincount(cut_ids)
    population = read_counts(inside_mult.copy(), edge_mult.copy(), interior_words, exterior_words, matched)
    acquisitions = {}
    for label, number in (("training", 17), ("held", 59)):
        rng1, rng2 = [np.random.default_rng([SEED, salt, number, i]) for i in (0, 1)]
        b = rng1.multinomial(SHOTS, inside_mult/len(words))
        e = rng2.multinomial(SHOTS, edge_mult/len(words))
        acquisitions[label] = dict(read_counts(b, e, interior_words, exterior_words, matched),
                                   bulk_counts=b.tolist(), exterior_counts=e.tolist())
    return dict(**row, native_pair_states=len(words), native_pair_hopping_links=transitions,
                source_cells=sorted(row["cells"] if support is None else support),
                active_source_cut_ports=sorted(active_border),
                resident_tag_launch_volume=len(ports)-2,
                resident_tag_aperture=resident_tag_strength()*(len(cut)-k.mean()),
                interior_words=interior_words, exterior_words=exterior_words,
                matched_energy_and_cut_occupation=matched, population=population, acquisitions=acquisitions)


def change(cycle, path, group):
    first, second = [x["population"] if group == "population" else x["acquisitions"][group]
                     for x in (cycle, path)]
    if cycle["accessible_volume"] != path["accessible_volume"]:
        raise ValueError("Shape transfer changed the enrolled launch volume")
    return dict(bulk_nats=second["bulk_nats"]-first["bulk_nats"],
                cut_bits=second["cut_bits"]-first["cut_bits"],
                aperture=path["aperture"]-cycle["aperture"],
                blank_flight=path["blank_meter"]["sampled_extrapolated_area"]-cycle["blank_meter"]["sampled_extrapolated_area"],
                resident_tag=path["resident_tag_aperture"]-cycle["resident_tag_aperture"])


def score(observables, eta_aperture, eta_capacity, eta_flight, eta_tag):
    return dict(**observables, aperture_residual_nats=observables["bulk_nats"]+eta_aperture*observables["aperture"],
                cut_residual_nats=observables["bulk_nats"]+eta_capacity*observables["cut_bits"],
                blank_flight_residual_nats=observables["bulk_nats"]+eta_flight*observables["blank_flight"],
                resident_tag_residual_nats=observables["bulk_nats"]+eta_tag*observables["resident_tag"])


def run():
    surveyed = geometry()
    surveys_by_name = {(x["board_level"], x["shape"]): x for x in surveyed}
    regions = [source(row, i) for i, row in enumerate(surveyed)]
    pairs = defaultdict(dict)
    for region in regions:
        family = "path" if region["shape"].startswith("path") else "cycle"
        key = (region["board_level"], len(region["cells"]))
        if family in pairs[key]:
            raise ValueError("Ambiguous enrolled shapes")
        pairs[key][family] = region
    training = change(pairs[(0, 5)]["cycle"], pairs[(0, 5)]["path"], "training")
    eta_aperture = -training["bulk_nats"]/training["aperture"]
    eta_capacity = -training["bulk_nats"]/training["cut_bits"]
    eta_flight = -training["bulk_nats"]/training["blank_flight"]
    eta_tag = -training["bulk_nats"]/training["resident_tag"]
    coefficients = eta_aperture, eta_capacity, eta_flight, eta_tag
    cross, within, common, common_sources = [], [], [], []
    for (level, size), pair in sorted(pairs.items()):
        observations = {name: score(change(pair["cycle"], pair["path"], name), *coefficients)
                        for name in ("held", "population")}
        cross.append(dict(board_level=level, cells_per_shape=size, **observations))
        support = sorted(set(pair["cycle"]["cells"]) & set(pair["path"]["cells"]))
        shared = [source(surveys_by_name[level, row["shape"]], 300+level*20+size*2+i, support)
                  for i, row in enumerate((pair["cycle"], pair["path"]))]
        common_sources.extend(shared)
        observations = {name: score(change(*shared, name), *coefficients)
                        for name in ("held", "population")}
        common.append(dict(board_level=level, cells_per_shape=size, common_source_cells=support, **observations))
        for region in pair.values():
            observations = {}
            for name in ("held", "population"):
                meter = region["population"] if name == "population" else region["acquisitions"][name]
                observations[name] = score(dict(bulk_nats=meter["strict_bulk_slope"],
                                                cut_bits=meter["strict_cut_slope"], aperture=0., blank_flight=0., resident_tag=0.),
                                           *coefficients)
            within.append(dict(board_level=level, shape=region["shape"], **observations))
    return dict(shots_per_independent_meter=SHOTS, registered_exposure=.5, endpoint_mode="analytic stationary",
                calibration=training, fitted_nats_per_aperture=eta_aperture,
                fitted_nats_per_exterior_bit=eta_capacity, shape_changes=cross,
                fitted_nats_per_blank_flight_unit=eta_flight, fitted_nats_per_resident_tag_unit=eta_tag,
                within_shape_changes=within, same_source_shape_queries=common,
                regions=regions, same_source_regions=common_sources)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), separators=(",", ":"), allow_nan=False)+"\n")
