"""Read free dispersion on the exact confinement controls and an end-connector intervention."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from ...row_22_confinement.mobile_string import commission, equal_degree, material
from .. import experiment
from ..confirmation import acquire as confirmation, spatial


HERE = Path(__file__).resolve().parent
N = 16


def free_diagonal(level, sites):
    cal = material.calibration_energies(8)
    contact, gap, private = cal
    board = material.federation(level)
    degree = np.array(board.degree)[sites].reshape(-1,2)
    left = material.bits(N)
    right = 1-left
    z = 2*left-1
    energy = contact*(1-right[:,:-1]*left[:,1:]).sum(axis=1)
    energy += private*(left @ (6-degree[:,0]) + right @ (6-degree[:,1]))
    expected = experiment.K*((z[:,:-1]*z[:,1:]).sum(axis=1) + z[:,0]-z[:,-1])
    defect = float(np.ptp(energy-expected))
    if defect>1e-10 or abs(gap/2-experiment.H)>1e-12 or np.any(degree[:,0]!=degree[:,1]):
        raise ArithmeticError("Physical free diagonal or dimer flips differ")
    return dict(orientation_words=len(z), residual=defect, dimer_h=gap/2,
                contact_cost=contact, end_fields=[experiment.K,-experiment.K])


def edge_control(route, level):
    owners = set(v//12 for v in route[0])
    dest = next(i for i in range(3,N) if owners.isdisjoint((route[i-1][1]//12, route[i][0]//12)))
    starts = []
    for reverse in (False, True):
        z = np.ones(N+2,dtype=int)
        z[-1] = -1
        if reverse:
            z[1] = -1
        starts.append((1+z[1:]*z[:-1])/2)
    shots = 2**22
    radius = 2*np.sqrt(np.log(2*96/.01)/(2*shots))
    rng = np.random.default_rng(223041+level)
    results = []
    # With the first connector idle, dimers 1..15 are an independent chain.
    # Their left virtual sign equals +1 regardless of the isolated dimer 0.
    cut_start = np.ones(N)
    cut_start[-1] = 0
    for t in (0,2,4,8):
        kernel = experiment.walls(N,8,t)
        native = [.5+kernel[dest]@(s-.5) for s in starts]
        cut = .5+experiment.walls(N-1,8,t)[dest-1]@(cut_start-.5)
        bank = rng.binomial(shots,np.clip([*native,cut,cut],0,1)).tolist()
        results.append(dict(t=t, counts=bank, native=(bank[0]-bank[1])/shots, cut=(bank[2]-bank[3])/shots))
    static = {}
    for length in material.LENGTHS:
        first = (N-length)//2
        first_connector = first-1
        if first_connector<=0:
            raise ArithmeticError("Static interior reaches altered connector")
        static[length] = [first_connector, first+length-1]
    return dict(receiver=dest, shots=shots, joint_99_radius=float(radius), rows=results,
                participating_static_connector_ranges=static)


def run(directory):
    directory.mkdir(exist_ok=True,parents=True)
    # Serialize the temporary CountStream capture if the package is called from another thread.
    meters = equal_degree.decode()
    study = []
    with np.load(equal_degree.RECEIVERS,allow_pickle=False) as archive:
        for level in range(3):
            sites = equal_degree.ports(level)
            sock = commission.material_stock(level,N,sites)
            route = sock["dimers"]
            if (not np.array_equal(archive[f"ports_{level}"],sites) or
                    sock["ordered_address_sha256"]!=str(archive[f"stock_sha_{level}"])):
                raise ValueError("Different row 22 material")
            diagonal = free_diagonal(level,sites)
            folder = directory/str(level)
            folder.mkdir()

            def keep(stream, name, counts):
                words = np.asarray(counts,dtype="<u8",order="C")
                np.save(folder/(name+".npy"),words,allow_pickle=False)
                stream.hashes[name] = sha256(memoryview(words).cast("B")).hexdigest()

            with confirmation._ACQUISITION, patch.object(experiment.CountStream,"keep",keep):
                obs = experiment.spectroscopy(N,np.random.default_rng(223031+level))
            physical = dict(route=route,observer=obs)
            held = spatial.score(physical,folder)
            calibrations = meters[level,8][0]
            h, j = calibrations[1]/2, calibrations[0]/4
            meter = dict(h=h,j=j,m=2*abs(j-h),c_squared=4*j*h)
            disp = None
            if "even_trained_k" in held:
                disp = confirmation.even_hold_dispersion(physical,held,meter)
                rows = disp["modes"]
                a,b = rows[0],rows[-1]
                slope = (b["omega"]-a["omega"])/(np.cos(b["k"])-np.cos(a["k"]))
                for r in rows[1:-1]:
                    r["refit_hopping_residual"] = r["omega"]-(a["omega"]+slope*(np.cos(r["k"])-np.cos(a["k"])))
            study.append(dict(level=level,route=route,stock_sha256=sock["ordered_address_sha256"],
                              native_diagonal=diagonal,held=held,dispersion=disp,count_arrays=obs["count_sha256"],
                              pair_sum_closure=obs["withheld_pair_sum_max_error"],
                              changed_contrast=obs["same_spatial_source_changed_contrast"],
                              end_control=edge_control(route,level)))
    return dict(protocol_sha256=sha256((HERE/"PROTOCOL.md").read_bytes()).hexdigest(),
                static_archive_sha256=sha256(equal_degree.RECEIVERS.read_bytes()).hexdigest(),
                boards=study)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--acquire",type=Path)
    cli.add_argument("--save-counts",type=Path)
    args = cli.parse_args()
    if args.save_counts and args.save_counts.exists():
        cli.error("The output count directory must be unused")
    with TemporaryDirectory() as temp:
        result = run(args.save_counts or Path(temp))
    if args.acquire:
        with args.acquire.open("x") as dest:
            json.dump(result,dest,indent=1,allow_nan=False)
            dest.write("\n")
    else:
        expected = json.loads((HERE/"results.json").read_text())
        if not confirmation.same_readings(json.loads(json.dumps(result)),expected):
            raise ArithmeticError("Free-release receipt differs")
    for b in result["boards"]:
        print(b["level"],b["held"]["status"],b["dispersion"]["gates"] if b["dispersion"] else None,
              len(b["count_arrays"]),"stored array hashes")


if __name__=="__main__":
    main()
