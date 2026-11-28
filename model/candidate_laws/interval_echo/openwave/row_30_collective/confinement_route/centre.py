"""Check the confinement work and free propagation with dimer 7 held by an idle write."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
from scipy.linalg import svd

from ...row_22_confinement.mobile_string import measure as static
from . import measure as study


CENTRE = 7


def native_columns():
    law = study.experiment.candidate.native_law
    memory = law.echo()[:,law.SRC.index((1,1))].reshape(-1,9)[:,law.SRC.index((1,1))]
    inputs = np.kron(memory[:,None],np.eye(9)[:,[law.SRC.index(w) for w in ((0,1),(1,0))]])
    errors = [float(abs(h@inputs).max()) for h in law.local_term_parts()]
    if max(errors)>1e-13:
        raise ValueError("The centre memory does not suppress both generator parts")
    return errors


def altered_matrix(level,length,pin):
    site = study.equal_degree.ports(level)
    board = study.material.federation(level)
    cal = study.material.calibration_energies(8)
    h = study.equal_degree.pinned_matrix(board,site,length,pin,cal)
    local = CENTRE - ((16-length)//2 + 1)
    if 0<=local<length-2:
        words = np.arange(len(h))
        h[words,words^(1<<local)] += cal[1]/2
    return h


def endpoint_arms():
    rng = np.random.default_rng(223075)
    bits = []
    for flip in (False,True):
        z = np.ones(18)
        z[-1] = -1
        if flip:
            z[6] = -1  # physical dimer 5, after the virtual left end
        bits.append((1+z[:-1]*z[1:])/2)
    rows = []
    for t in (0,2,4,8):
        native = study.experiment.walls(16,8,t)[10]
        with patch.object(study.experiment,"dual",split):
            cut = study.experiment.walls(16,8,t)[10]
        p = [.5+w@(a-.5) for w in (native,cut) for a in bits]
        counts = rng.binomial(2**22,np.clip(p,0,1)).tolist()
        rows.append(dict(t=t,counts=counts,native_delta=(counts[0]-counts[1])/2**22,
                         idle_delta=(counts[2]-counts[3])/2**22,
                         exact_idle_difference=float(p[2]-p[3])))
    return rows


def split(n,j):
    matrix = 2*j*np.eye(n+1)-np.eye(n+1,k=-1)
    matrix[CENTRE+1,CENTRE] = 0
    u,e,vh = svd(matrix)
    return u,e,vh.T


def run():
    failures = native_columns()
    sources = {}
    for length in study.material.LENGTHS:
        lower,upper = (altered_matrix(0,length,pin) for pin in (0,1))
        error = float(max(abs(upper-lower[::-1,::-1]-32*np.eye(len(upper))).ravel()))
        if error>1e-10:
            raise ArithmeticError("Centre-idle work conjugacy failed")
        _, panels = study.material.pinned_laws(lower,upper)
        sources[length] = panels
    original = study.equal_degree.decode()
    boards = []
    with np.load(study.equal_degree.RECEIVERS,allow_pickle=False) as bank:
        for level in range(3):
            rng = np.random.default_rng(223071+level)
            clock = np.polyfit(static.TICK,np.unwrap(static.port_phase(bank[f"clock_{level}_8"])),1)[0]
            rows = []
            for length,panels in sources.items():
                for name,(signal,_) in zip(("sharp_intact","sharp_screened"),panels,strict=True):
                    counts = static.sample_ports(rng,signal)
                    decoded = static.ramsey(counts,clock)
                    rows.append(dict(length=length,source=name,counts=counts.tolist(),read=decoded,
                                     original_visibility=original[level,8][1][length][name]["min_visibility"]))
            boards.append(rows)
    route = study.commission.material_stock(0,16,study.equal_degree.ports(0))["dimers"]
    with TemporaryDirectory() as temp:
        folder = Path(temp)
        def keep(stream,name,counts):
            words = np.asarray(counts,dtype="<u8",order="C")
            np.save(folder/(name+".npy"),words,allow_pickle=False)
            stream.hashes[name] = sha256(memoryview(words).cast("B")).hexdigest()
        with study.confirmation._ACQUISITION, patch.object(study.experiment.CountStream,"keep",keep), patch.object(study.experiment,"dual",split):
            obs = study.experiment.spectroscopy(16,np.random.default_rng(223074))
        held = study.spatial.score(dict(route=route,observer=obs),folder)
    return dict(protocol_sha256=sha256((study.HERE/"CENTRE.md").read_bytes()).hexdigest(),
                source_columns=failures, boards=boards, held=held, spectra_sha=obs["count_sha256"],
                cross_centre=endpoint_arms())


if __name__=="__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--acquire",type=Path)
    args=cli.parse_args()
    result=run()
    if args.acquire:
        with args.acquire.open("x") as dest:
            json.dump(result,dest,indent=1,allow_nan=False)
            dest.write("\n")
    elif not study.confirmation.same_readings(json.loads(json.dumps(result)),json.loads((study.HERE/"centre.json").read_text())):
        raise ArithmeticError("Centre-idle evidence differs")
    print(result["held"]["status"],result["cross_centre"][-1])
