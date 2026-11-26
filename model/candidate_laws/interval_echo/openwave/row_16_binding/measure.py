"""Energy/position banks for six ordinary loads and two separated face stocks."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy import linalg

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.executable_model import carrier
from model.experiments.count_intervals import clopper_pearson
from .certificates import compactness, face_band, partitions
from .apparatus import ARMS, NATIVE_LINE, STARTS, STEP, Shape, calibrate, sample, time_port

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "records.json.gz"
SEED = 161433


def boards():
    result = []
    geometry = compactness()
    for name in ("cell", 0, 1, 2):
        graph = carrier.carrier(name)
        neighbourhoods = list(map(set, graph.neighbours))
        for a, b in graph.seams:
            common = neighbourhoods[a] & neighbourhoods[b]
            if any(neighbourhoods[x] & common for x in common):
                raise ValueError("Found a four-clique")
        external = [site for edge in graph.external for site in edge]
        if len(external) != len(set(external)):
            raise ValueError("Intercell edges must form a matching")
        same_cell = [(a % 12, b % 12) for a, b in graph.seams if a//12 == b//12]
        if len(same_cell) != 30*graph.cells or set(same_cell) != set(carrier.cell().seams):
            raise ValueError("The owners must share the registered cell")
        cap = set([0, *carrier.cell().neighbours[0]])
        internal = sum(a in cap and b in cap for a, b in graph.seams)
        exits = sum((a in cap) != (b in cap) for a, b in graph.seams)
        upper = 8 * (36-2*internal) + exits / 2
        result.append(dict(carrier=name, ports=graph.n, maximum_degree=int(max(graph.degree)),
                           cap_cut=exits, cap_edges=internal, cap_energy_upper=upper,
                           proper_split_lower=geometry["disconnected_lower"], two_trimers_lower=192,
                           free_connected_probability_lower=1-(upper-geometry["potential_lower"])
                           /(geometry["disconnected_lower"]-geometry["potential_lower"])))
    return result


def finite_bank(rng, shape, amplitude):
    counts = sample(rng, abs(amplitude)**2)
    return dict(counts=counts.tolist(), dark=int(counts[-1]), **shape.read(counts[:-1]))


def run():
    rng = np.random.default_rng(SEED)
    sectors = {n: BinarySector(carrier.cell(), n) for n in range(1, 7)}
    six = sectors[6]
    shape = Shape(six)
    trimer = sectors[3]
    face = trimer.index[carrier.cell().faces()[0]]
    record = dict(starts=STARTS, seed=SEED, probability_grid_bits=40, timing_arms=ARMS,
                  timing_step=STEP, native_line=NATIVE_LINE, boards=boards(), compactness=compactness(),
                  word_order="lexicographic combinations of twelve ports; last bin dark",
                  separate_order="lexicographic product of the two three-load registers; last bin dark",
                  axes=shape.axes, source_words={k: six.states[v] for k, v in shape.sources.items()},
                  signed_and_timing_calibration=calibrate(rng), sharp_standards={}, laws={})
    for name, idx in (("axis_max", np.argmax(shape.cell_q[:, 0].round(10))),
                      ("axis_min", np.argmin(shape.cell_q[:, 0].round(10))),
                      *shape.sources.items()):
        pulse = np.zeros(len(six.states))
        pulse[idx] = 1
        record["sharp_standards"][name] = dict(word=six.states[idx], **finite_bank(rng, shape, pulse))
    for name, contrast in (("native", 8), ("exchange_only", 0), ("reversed", -8)):
        values, vectors = linalg.eigh(six.hamiltonian(contrast, False))
        three, basis = linalg.eigh(trimer.hamiltonian(contrast, False))
        floors = {n: float(linalg.eigvalsh(sectors[n].hamiltonian(contrast, False),
                                            subset_by_index=[0, 0])[0])
                  for n in (1, 2, 4, 5)}
        floors.update({3: float(three[0]), 6: float(values[0])})
        fragments = {"+".join(map(str, p)): sum(floors[n] for n in p)
                     for p in partitions(6) if len(p) > 1}
        result = dict(floors=floors, fragment_floors=fragments, compact={}, separated={})
        for source, index in shape.sources.items():
            weights = vectors[index]
            runs = {}
            for time in (0, 1, 8):
                phase = np.exp(-1j*time*values)
                runs[f"free_{time}"] = finite_bank(rng, shape, vectors @ (phase*weights))
                if time != 1:
                    references = [NATIVE_LINE, 0.] if name == "exchange_only" else [NATIVE_LINE]
                    for ref in references:
                        amplitude = vectors @ (phase * time_port(values, ref) * weights)
                        runs[f"filter_{ref}_{time}"] = finite_bank(rng, shape, amplitude)
            result["compact"][source] = runs
        references = [NATIVE_LINE, 0.] if name == "exchange_only" else [NATIVE_LINE]
        for ref in references:
            joint = three[:, None] + three[None, :]
            weights = basis[face, :, None] * basis[face, None, :]
            amplitude = basis @ (time_port(joint, ref) * weights) @ basis.T
            counts = sample(rng, abs(amplitude)**2)
            result["separated"][str(ref)] = dict(counts=counts.tolist(), dark=int(counts[-1]))
        if name == "native":
            # Beyond the cap configurations the positive potential is >=144.
            # The cap subspace has constant compression 133 I.
            band = values < 144
            v = vectors[:, band]
            joint_norm = linalg.svdvals(v[shape.pairs])[0]**2
            result["cap_band"] = dict(certificate=face_band(six, shape), dimension=int(band.sum()), minimum=float(min(values[band])),
                                      maximum=float(max(values[band])),
                                      next=float(min(values[~band])), two_face_operator_norm=float(joint_norm))
            high_lo, high_hi = 160, 6*6*8+30
            edge = np.array([high_lo, high_hi])
            result["high_energy_aperture_probability_bound"] = float(
                1/(ARMS*min(np.sin((edge-NATIVE_LINE)*STEP/2)))**2)
        result["incoherent_bright"] = int(rng.binomial(STARTS, 1/ARMS))
        record["laws"][name] = result
    native = record["laws"]["native"]
    cap = native["compact"]["cap"][f"filter_{NATIVE_LINE}_0"]
    pair = native["compact"]["two_faces"][f"filter_{NATIVE_LINE}_0"]
    separated = STARTS - native["separated"][str(NATIVE_LINE)]["dark"]
    observations = {"cap_bright": (cap["hits"], STARTS), "separated_bright": (separated, STARTS),
                    "bright_cap_two_faces": (cap["two_faces"], cap["hits"]),
                    "two_face_input_bright": (pair["hits"], STARTS),
                    "bright_two_face_input_two_faces": (pair["two_faces"], pair["hits"])}
    record["sampler_intervals_99_family"] = {
        k: list(map(float, clopper_pearson(*obs, tail=0.001)))
        for k, obs in observations.items()
    }
    return record


def encode(record):
    payload = json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return gzip.compress(payload, compresslevel=9, mtime=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ARCHIVE)
    args = parser.parse_args()
    args.output.write_bytes(encode(run()))
