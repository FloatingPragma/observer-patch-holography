"""Excess and vacancy records above one common uniform load background."""
import argparse
import json

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model.carrier import federation
from .transport import HERE, PROTOCOL as RELEASE_PROTOCOL, census, release, write_receipt

PROTOCOL = HERE / "conjugation_protocol.json"
RECEIPT = HERE / "conjugation.json"


def defect_generator(car, sign, strength, reference):
    """Assemble native seam columns relative to the energy of the uniform medium."""
    if sign not in (-1, 1) or reference not in (0, 1):
        raise ValueError("Unregistered defect or reference")
    seam = law.seam_term(strength)
    medium = law.SRC.index((1, 1))
    left = law.SRC.index((1 + sign, 1))
    right = law.SRC.index((1, 1 + sign))
    a, b = np.array(car.seams).T
    diag = np.zeros(car.n)
    np.add.at(diag, a, seam[left,left] - seam[medium,medium])
    np.add.at(diag, b, seam[right,right] - seam[medium,medium])
    diag += strength*(6-car.degree)*((1+sign-reference)**2 - (1-reference)**2)
    exchange = sparse.coo_matrix((np.r_[np.full(len(a), seam[right,left]),
                                       np.full(len(a), seam[left,right])],
                                  (np.r_[b,a], np.r_[a,b])), shape=(car.n,car.n))
    return (exchange + sparse.diags(diag)).tocsr()


def local_pulse(strength, sign):
    u = expm(-0.5j*np.pi*law.seam_term(strength))
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
        if abs(record.sum() - 1) > 2e-12:
            raise ArithmeticError("Defect read lost probability")
        readings[str(t)] = record.round(12).tolist()
    return {"slope": round(float(slope),12), "readings": readings}


def campaign():
    spec = json.loads(PROTOCOL.read_text())
    boards = {}
    for level in spec["levels"]:
        car = federation(level)
        home = spec["cell"]
        local = range(12*home,12*(home+1))
        a = next(p for p in local if car.degree[p] == 6)
        b = min(p for p in car.neighbours[a] if p in local and car.degree[p] == 5)
        c = min(p for p in car.neighbours[b] if p in local and p != a)
        arms = []
        for setting in [spec["native"], *spec["controls"]]:
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
                    if np.shape(target) != (13,) or not np.allclose(bins, target, rtol=0, atol=2e-11):
                        raise ValueError(f"Defect read differs at {level}, {sign}, {t}")


def pair_regression():
    spec = json.loads(PROTOCOL.read_text())["pair_regression"]
    launch = json.loads(RELEASE_PROTOCOL.read_text())["release"]
    sector = PairSector(federation(launch["level"]))
    seams, source = release(sector, launch["cell"])
    missing = 6 - sector.car.degree
    potential = -2.0*spec["kappa"]*spec["reference"]*(missing[sector.aa]+missing[sector.bb])
    generator = -1j*(sector.hamiltonian(spec["kappa"]) + sparse.diags(potential))
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
    parser.add_argument("--pair-regression", action="store_true")
    args = parser.parse_args()
    if args.pair_regression:
        readings = pair_regression()
        path = HERE / "centered_pair.json"
        if args.write:
            write_receipt(path, readings)
        else:
            stored = json.loads(path.read_text())
            if (readings.keys() != stored.keys() or
                any(np.shape(stored[t]) != (34,) or not np.allclose(v, stored[t], atol=1e-9, rtol=0)
                    for t,v in readings.items())):
                raise ValueError("Centered-reference pair census differs")
        print({t: sum(bins[:30]) for t,bins in readings.items()})
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
