"""Measure uniform accessible volume and native crossing apertures on three screens."""
import argparse
import hashlib
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.contact import radar
from model.experiments.interval_echo_einstein_station import refinement
from model.experiments.interval_echo_einstein_station.regional_geometry.measure import hand_times
from model.experiments.interval_echo_null_balance.regional import verify

HERE = Path(__file__).resolve().parent
TICKS = [1, 2, 4]
FLIGHT_SHOTS = 2**32
MEMORY_SHOTS = 2**24
EPSILON = .25
SEED = 26112361


def propagate(h, launch):
    start = np.zeros((h.shape[0], len(launch)))
    start[launch, np.arange(len(launch))] = 1
    end = expm_multiply(-1j*h, start, start=0, stop=radar.DT*max(TICKS),
                        num=1+max(TICKS), traceA=-1j*h.diagonal().sum())
    return abs(end[TICKS])**2


def aperture(h, region, hand, dense=False, held=True):
    n = h.shape[0]
    exterior = np.setdiff1d(np.arange(n), region)
    crossing = h[exterior][:, region]
    leading_area = float(4*np.sum(abs(crossing.data)**2))
    flight = propagate(h, region)
    flux = flight[:, exterior].sum(axis=(1, 2))
    exterior_flux = flux.copy()
    dense_error = None
    if dense:
        exterior_flux = np.array([
            np.sum(abs(expm(-1j*t*radar.DT*h.toarray())[np.ix_(region, exterior)])**2)
            for t in TICKS
        ])
        dense_error = float(max(abs(flux-exterior_flux)))
        if dense_error > 1e-11:
            raise ArithmeticError("direct complementary preparation disagrees")
    held_rows = []
    if held:
        linked = np.asarray((abs(crossing)>1e-13).sum(axis=1)).ravel() > 0
        near = int(exterior[linked][0])
        # Choose a held address one outside-to-outside exchange from a crossing,
        # but with zero direct coupling to the region itself.
        neighbours = np.asarray(abs(h[exterior][:, [near]]).toarray()).ravel() > 1e-13
        far = int(exterior[~linked & neighbours][0])
        independent = propagate(h, [near, far])[:, region].sum(axis=1)
        predicted = flight[:, [near, far]].sum(axis=2)
        np.testing.assert_allclose(predicted, independent, atol=3e-13, rtol=0)
        for k, port in enumerate([near, far]):
            held_rows.append(dict(component_address=port, inward=independent[:, k].tolist(),
                                  reciprocal_prediction=predicted[:, k].tolist()))
    return dict(outward=flux/len(region), inward=exterior_flux/len(exterior),
                area=4*flux/hand**2, leading_area=leading_area, held=held_rows,
                full_inward_dense_max_error=dense_error)


def observe(outward, inward, rng):
    out_count = rng.binomial(FLIGHT_SHOTS, outward)
    in_count = rng.binomial(FLIGHT_SHOTS, inward)
    x, y = out_count/FLIGHT_SHOTS, in_count/FLIGHT_SHOTS
    fraction = y/(x+y)
    se = np.sqrt(y*y*x*(1-x)+x*x*y*(1-y))/(np.sqrt(FLIGHT_SHOTS)*(x+y)**2)
    return dict(out_counts=out_count.tolist(), in_counts=in_count.tolist(),
                volume_fraction=fraction.tolist(), standard_error=se.tolist())


def memory_records(ref, tangents, rng):
    ref_counts = rng.multinomial(MEMORY_SHOTS, ref.ravel()).reshape(3, 3)
    ref_read = ref_counts/MEMORY_SHOTS
    score = -np.log(ref_read)
    readings = []
    for tangent in tangents:
        lo, hi = (rng.multinomial(MEMORY_SHOTS, (ref+s*EPSILON*tangent).ravel()).reshape(3, 3)
                  for s in [-1, 1])
        slope = (hi-lo)/(2*EPSILON*MEMORY_SHOTS)
        variance = sum(np.sum(p*(score-np.sum(p*score))**2)/(4*EPSILON**2*MEMORY_SHOTS)
                       for p in [lo/MEMORY_SHOTS, hi/MEMORY_SHOTS])
        variance += (np.sum(slope**2/ref_read)-slope.sum()**2)/MEMORY_SHOTS
        readings.append(dict(minus_counts=lo.tolist(), plus_counts=hi.tolist(),
                             entropy_slope=float(np.sum(slope*score)), se=float(np.sqrt(variance))))
    return dict(reference_counts=ref_counts.tolist(), tangents=readings)


def balance(slopes, entropy):
    norm = float(slopes @ slopes)
    if norm < 1e-28:
        return dict(area_slopes=slopes.tolist(), coefficient=None, relative_residual=1.)
    coefficient = float(-(slopes @ entropy)/norm)
    residual = np.linalg.norm(entropy+coefficient*slopes)/np.linalg.norm(entropy)
    return dict(area_slopes=slopes.tolist(), coefficient=coefficient, relative_residual=float(residual))


def run():
    config = json.loads((HERE.parent/"refinement.json").read_text())
    source = dict(source_ports=config["source_ports"], hand_seam=config["hand_seam"],
                  source_word_energies=config["energies"])
    energies = config["energies"]
    ref = np.array(config["source_units"][1])/36
    rows = np.eye(3)[:2]-np.eye(3)[2]
    tangents = np.stack([np.outer(a,b)/36 for a,b in product(rows,repeat=2)])
    entropy = -np.einsum("kij,ij->k", tangents, np.log(ref))
    rng = np.random.default_rng(SEED)
    hand = hand_times()[TICKS]
    blank, occupied = [], []
    joint_tangent = None
    followup_rng = np.random.default_rng(SEED+1)
    for level in range(3):
        board = federation(level)
        region = np.arange(12*4**level)
        blank_read = aperture(candidate.generator(board).tocsr(), region, hand, dense=level==0)
        blank.append(dict(cells=board.cells, n_ports=board.n, n_region=len(region),
                          leading_area=blank_read["leading_area"], finite_areas=blank_read["area"].tolist(),
                          held=blank_read["held"], dense_inward_error=blank_read["full_inward_dense_max_error"],
                          counts=observe(blank_read["outward"], blank_read["inward"], rng),
                          areas_cubed_over_volume_squared=(np.r_[blank_read["leading_area"],
                                                            blank_read["area"]]**3/len(region)**2).tolist()))
        h, _, _ = refinement.station(level, source, tuple(config["protected_tag_seam"]))
        _, labels = connected_components(abs(h)>1e-13, directed=False)
        a, _, c = source["source_ports"]
        active = np.flatnonzero(labels==labels[a])
        if c not in active:
            raise ArithmeticError("matter endpoints are not in same mobile component")
        local = np.flatnonzero(active < 12)
        a, c = np.searchsorted(active, [a, c])
        base = h[active][:, active]
        records = []
        for i,j in product(range(3), repeat=2):
            word_h = base.tolil()
            word_h[a,a] -= energies[i]
            word_h[c,c] -= energies[j]
            result = aperture(word_h.tocsr(), local, hand, dense=level==0, held=(i==j==0))
            records.append(result)
        areas = np.array([np.r_[r["leading_area"], r["area"]] for r in records]).reshape(3,3,4)
        area_slopes = np.einsum("kij,ijt->tk", tangents, areas)
        if level==0:
            normalized = area_slopes[1:]/np.linalg.norm(area_slopes[1:],axis=1)[:,None]
            _, singular_values, right = np.linalg.svd(normalized, full_matrices=True)
            if singular_values[-1] < 1e-5:
                raise ArithmeticError("three appointment gradients are numerically unresolved")
            joint_tangent = right[-1]*np.sign(right[-1] @ entropy)
        interference = np.einsum("k,kij->ij",joint_tangent,tangents)
        joint_record = memory_records(ref, [interference], followup_rng)
        trials = {}
        out = np.array([r["outward"] for r in records]).reshape(3,3,3)
        inc = np.array([r["inward"] for r in records]).reshape(3,3,3)
        crosscheck = None
        if level==0:
            alternate_areas = 4*(len(active)-len(local))*inc/hand**2
            independent = np.einsum("kij,ijt->tk", tangents, alternate_areas)
            crosscheck = float(np.max(abs(independent-area_slopes[1:])))
            if crosscheck >= 1e-13:
                raise ArithmeticError("independent dense aperture gradients disagree")
        for label,p in [("reference",ref)]+[(f"{k}_{s}",ref+s*EPSILON*d)
                                            for k,d in enumerate(tangents) for s in [-1,1]]:
            outward, inward = [np.einsum("ij,ijt->t", p, v) for v in [out,inc]]
            np.testing.assert_allclose(inward/(inward+outward), len(local)/len(active), atol=1e-15)
            trials[label] = observe(outward,inward,rng)
        blind = np.einsum("kij,ijt->tk", tangents, np.broadcast_to(areas[0,0],areas.shape))
        occupied.append(dict(cells=board.cells, accessible_volume=len(active), home_volume=len(local),
                             home_port_addresses=active[local].tolist(), conditional_word_areas=areas.tolist(),
                             exact_volume_fraction=len(local)/len(active), counts=trials,
                             balance=[balance(s,entropy) for s in area_slopes],
                             dense_inward_max_slope_discrepancy=crosscheck,
                             transferred_joint_tangent=dict(coefficients=joint_tangent.tolist(),
                                 finite_area_slopes=(area_slopes[1:] @ joint_tangent).tolist(),
                                 entropy_slope=float(joint_tangent @ entropy), memory=joint_record),
                             blind_maximum_area_slope=float(np.max(abs(blind))),
                             matter_memory=memory_records(ref,tangents,rng)))
    baseline = np.array(blank[0]["areas_cubed_over_volume_squared"])
    for record in blank:
        record["ratio_to_20_cell_calibration"] = (record["areas_cubed_over_volume_squared"]/baseline).tolist()
    return dict(followup_sha256=hashlib.sha256((HERE/"FOLLOWUP.md").read_bytes()).hexdigest(),
                protocol_sha256=hashlib.sha256((HERE/"DESIGN.md").read_bytes()).hexdigest(), seed=SEED,
                aperture_order=["short_hand_limit","appointment_1","appointment_2","appointment_4"],
                material_hand=hand.tolist(), flight_shots=FLIGHT_SHOTS, memory_shots=MEMORY_SHOTS,
                entropy_slopes=entropy.tolist(), memory_tangents=tangents.tolist(), blank=blank,
                occupied=occupied)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE/"results.json"
    if args.write:
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    else:
        verify(json.loads(path.read_text()),result)
    for b,o in zip(result["blank"],result["occupied"],strict=True):
        print(b["cells"], b["n_region"], b["leading_area"], b["ratio_to_20_cell_calibration"],
              [x["relative_residual"] for x in o["balance"]])


if __name__ == "__main__":
    main()
