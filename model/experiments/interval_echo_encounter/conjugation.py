"""Excess and vacancy records above one common uniform load background."""
import argparse
import json
from functools import lru_cache

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector, pair_index
from model.executable_model.carrier import federation
from .transport import (
    HERE, PROTOCOL as RELEASE_PROTOCOL, census, probability_vector, read_protocol, release, write_receipt,
)

PROTOCOL = HERE / "conjugation_protocol.json"
RECEIPT = HERE / "conjugation.json"

@lru_cache(maxsize=1)
def singleton_profile():
    """Read |written contrast| = 1 through the native echo on a blank reference."""
    labels = np.repeat([abs(d) == 1 for _,d in law.LABELS], law.N_SRC)
    write = law._write_pattern()  # integer sqrt(2) times the native isometry
    selected = np.array([int(abs(a-b) == 1) for a,b in law.SRC])
    if not np.array_equal(labels[:,None]*write, write*selected):
        raise ArithmeticError("The material label projector fails to intertwine the native write")
    echo = law.echo()
    embedded = np.eye(law.DIM, law.N_SRC)
    selector = np.r_[np.zeros(law.N_SRC), labels]
    if not np.allclose(echo @ (selector[:,None] * (echo @ embedded)),
                       embedded*selected, rtol=0, atol=2e-14):
        raise ArithmeticError("The inverse echo leaks outside the blank reference code")
    return tuple(int(selected[law.SRC.index((n,0))]) for n in range(3))


def reference_profile(reference):
    if reference == "singleton":
        return singleton_profile()
    if reference in (0,1):
        return tuple((n-reference)**2 for n in range(3))
    raise ValueError("Unknown reference service")


@lru_cache(maxsize=None)
def free_seam(strength):
    """Check the blank-collar reduction and all transitions omitted by the sector models."""
    seam = law.seam_term(strength)
    exchanged = [law.SRC.index((b, a)) for a, b in law.SRC]
    columns = np.arange(law.N_SRC)
    supported = np.eye(law.N_SRC, dtype=bool)
    supported[exchanged, columns] = True
    if np.any(seam[~supported]):
        raise ArithmeticError("Free seam changes the load multiset")
    kinetic, contrast = law.local_term_parts()
    embedded = np.eye(law.DIM, law.N_SRC)
    if not np.allclose((kinetic + strength * contrast) @ embedded,
                       embedded @ seam, atol=3e-13, rtol=0):
        raise ArithmeticError("Native generator differs on the blank collar")
    return seam


def defect_generator(car, sign, strength, reference, background=1):
    """Assemble native seam columns relative to the energy of the uniform medium."""
    amount = background + sign
    if (background, amount) not in ((1,0), (1,2), (0,2)):
        raise ValueError("Unregistered one-position defect")
    profile = reference_profile(reference)
    seam = free_seam(strength)
    medium = law.SRC.index((background, background))
    left = law.SRC.index((amount, background))
    right = law.SRC.index((background, amount))
    a, b = np.array(car.seams).T
    diag = np.zeros(car.n)
    np.add.at(diag, a, seam[left,left] - seam[medium,medium])
    np.add.at(diag, b, seam[right,right] - seam[medium,medium])
    diag += strength*(6-car.degree)*(profile[amount]-profile[background])
    exchange = sparse.coo_matrix((np.r_[np.full(len(a), seam[right,left]),
                                       np.full(len(a), seam[left,right])],
                                  (np.r_[b,a], np.r_[a,b])), shape=(car.n,car.n))
    return (exchange + sparse.diags(diag)).tocsr()


def local_pulse(strength, sign):
    u = expm(-0.5j*np.pi*free_seam(strength))
    inputs = [law.SRC.index((1+sign,1)), law.SRC.index((1,1+sign))]
    spectator = law.SRC.index((1,1))
    if not np.allclose(u[:,spectator], np.eye(9)[:,spectator], atol=2e-14, rtol=0):
        raise ValueError("A spectator changed during the pulse")
    pulse = u[np.ix_(inputs, inputs)]
    if not np.allclose(u[:,inputs], np.eye(9)[:,inputs] @ pulse, atol=2e-14, rtol=0):
        raise ValueError("Pulse leaves the one-defect description")
    return pulse


def run_arm(car, ports, home, exposures, strength, reference, sign):
    a,b,c = ports
    pulse = local_pulse(strength, sign)
    h = defect_generator(car, sign, strength, reference)
    if reference == "singleton":
        potential = strength*(2*car.degree-6)
    else:
        potential = strength*(6 + 2*sign*(1-reference)*(6-car.degree))
    difference = h - (car.laplacian()/2 + sparse.diags(potential, dtype=float))
    if difference.count_nonzero():
        raise ArithmeticError("Qutrit columns differ from the defect Hamiltonian")
    state = np.zeros(car.n, complex)
    state[a] = 1
    for bond in ([a,b], [b,c]):
        state[bond] = pulse @ state[bond]
    generator = -1j*h
    tangent = pulse @ (generator @ state)[[a,b]]
    initial = pulse @ state[[a,b]]
    slope = 2*(initial[0].conj()*tangent[0]).real
    inside = slice(12*home,12*(home+1))
    readings = {}
    for t in exposures:
        evolved = expm_multiply(t*generator, state) if t else state.copy()
        evolved[[a,b]] = pulse @ evolved[[a,b]]
        weights = abs(evolved)**2
        quiet = weights[:inside.start].sum() + weights[inside.stop:].sum()
        record = np.r_[weights[inside], quiet]
        if not abs(record.sum() - 1) <= 2e-12:
            raise ArithmeticError("Defect read lost probability")
        readings[str(t)] = record.round(12).tolist()
    return {"slope": round(float(slope),12), "readings": readings}


def campaign():
    spec = read_protocol(PROTOCOL)
    boards = {}
    for level in spec["levels"]:
        car = federation(level)
        home = spec["cell"]
        local = range(12*home,12*(home+1))
        a = next(p for p in local if car.degree[p] == 6)
        b = min(p for p in car.neighbours[a] if p in local and car.degree[p] == 5)
        c = min(p for p in car.neighbours[b] if p in local and p != a)
        arms = []
        if list(singleton_profile()) != spec["singleton_reference"]["profile"]:
            raise ValueError("The registered reference profile differs from the native projector")
        for setting in [spec["native"], *spec["controls"], {"kappa": 8, "reference": "singleton"}]:
            records = [run_arm(car, (a,b,c), home, spec["exposures"], setting["kappa"],
                               setting["reference"], sign) for sign in spec["defects"]]
            arms.append({**setting, "excess": records[0], "vacancy": records[1]})
        boards[str(level)] = {"ports": [a,b,c], "arms": arms}
    return boards


def verify(result, archive):
    if result.keys() != archive.keys():
        raise ValueError("Carrier levels differ")
    for level, board in result.items():
        other = archive[level]
        if board["ports"] != other["ports"]:
            raise ValueError("Preparation ports differ")
        for arm, stored in zip(board["arms"], other["arms"], strict=True):
            if (arm["kappa"],arm["reference"]) != (stored["kappa"],stored["reference"]):
                raise ValueError("Control law differs")
            for sign in ("excess", "vacancy"):
                readings = arm[sign]["readings"]
                if (arm[sign]["slope"] != stored[sign]["slope"] or
                        readings.keys() != stored[sign]["readings"].keys()):
                    raise ValueError("Tangent or exposure menu differs")
                for t, bins in readings.items():
                    target = stored[sign]["readings"][t]
                    if not np.allclose(bins, probability_vector(target, 13), rtol=0, atol=2e-11):
                        raise ValueError(f"Defect read differs at {level}, {sign}, {t}")


def binary_generator(sector, strength, reference):
    """Assemble by physical seam, relative to vacuum and the common 12κ shift."""
    seam = free_seam(strength)
    empty, left, right, both = [law.SRC.index(s) for s in ((0,0), (1,0), (0,1), (1,1))]
    missing = 6 - sector.car.degree
    profile = reference_profile(reference)
    diagonal = strength * float(profile[1] - profile[0]) * (missing[sector.aa] + missing[sector.bb])
    diagonal -= 12 * strength
    ports = np.arange(sector.n)
    left_positions, right_positions = [], []
    for a, b in sector.car.seams:
        other = ports[(ports != a) & (ports != b)]
        ia = pair_index(sector.n, np.minimum(a, other), np.maximum(a, other))
        ib = pair_index(sector.n, np.minimum(b, other), np.maximum(b, other))
        diagonal[ia] += seam[left,left] - seam[empty,empty]
        diagonal[ib] += seam[right,right] - seam[empty,empty]
        diagonal[sector.index(a,b)] += seam[both,both] - seam[empty,empty]
        left_positions.append(ia)
        right_positions.append(ib)
    ia, ib = np.concatenate(left_positions), np.concatenate(right_positions)
    amplitudes = np.r_[np.full(len(ia), seam[right,left]), np.full(len(ib), seam[left,right])]
    moves = sparse.coo_matrix((amplitudes, (np.r_[ib,ia], np.r_[ia,ib])), shape=(sector.P,sector.P))
    return (moves + sparse.diags(diagonal)).tocsr()


def pair_regression(reference):
    spec = read_protocol(PROTOCOL)["pair_regression"]
    launch = read_protocol(RELEASE_PROTOCOL)["release"]
    sector = PairSector(federation(launch["level"]))
    seams, source = release(sector, launch["cell"])
    generator = -1j * binary_generator(sector, spec["kappa"], reference)
    state = np.zeros(sector.P, complex)
    state[sector.index(*source)] = 1
    readings = {}
    previous = 0
    for t in spec["exposures"]:
        state = expm_multiply((t-previous)*generator, state)
        readings[str(t)] = census(sector, launch["cell"], seams, state)
        previous = t
    return readings


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    comparisons = parser.add_mutually_exclusive_group()
    comparisons.add_argument("--pair-regression", action="store_true")
    comparisons.add_argument("--singleton-pair", action="store_true")
    args = parser.parse_args()
    if args.pair_regression or args.singleton_pair:
        selected = "singleton" if args.singleton_pair else read_protocol(PROTOCOL)["pair_regression"]["reference"]
        readings = pair_regression(selected)
        path = HERE / "centered_pair.json"
        if args.write and not args.singleton_pair:
            write_receipt(path, readings)
        else:
            stored = (json.loads((HERE/"transport.json").read_text())["arms"]["8"]
                      if args.singleton_pair else json.loads(path.read_text()))
            if (readings.keys() != stored.keys() or
                any(np.abs(v - probability_vector(stored[t], 34)).sum() > 1e-9
                    for t, v in readings.items())):
                raise ValueError("Reference-service pair census differs")
        native = json.loads((HERE / "transport.json").read_text())["arms"]["8"]
        for t, bins in readings.items():
            hit = sum(bins[:30])
            print(t, "home", round(hit, 9), "change from native", round(hit - sum(native[t][:30]), 9))
        raise SystemExit(0)
    result = campaign()
    if args.write:
        write_receipt(RECEIPT, result)
    else:
        verify(result, json.loads(RECEIPT.read_text()))
    for level, board in result.items():
        port = board["ports"][0] % 12
        print(level, [(arm["excess"]["readings"]["0.0625"][port],
                       arm["vacancy"]["readings"]["0.0625"][port]) for arm in board["arms"]])
