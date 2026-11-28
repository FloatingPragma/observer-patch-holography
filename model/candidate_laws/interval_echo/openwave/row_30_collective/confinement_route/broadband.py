"""Predicted and freely refitted held waveforms on the 16-dimer passive material."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.optimize import least_squares

from ..impulse import measure as impulse
from . import measure


def run(save=None):
    times = .06*np.arange(181)
    stock = measure.commission.material_stock(0,16,measure.equal_degree.ports(0))
    cal = measure.equal_degree.decode()[0,8][0]
    j,h = cal[0]/4,cal[1]/2
    with measure.confirmation._ACQUISITION, patch.object(impulse,"TIMES",times):
        raw,hand = impulse.acquire(16)
    if save:
        with save.open("xb") as out:
            np.savez_compressed(out,signs=raw,hand=hand)
    rows=[]
    # The independently calibrated square-gradient operator includes the physical end field.
    lap = np.diag(np.full(17,2.))-np.eye(17,k=1)-np.eye(17,k=-1)
    lap[0,0] -= h/j
    lam,basis = np.linalg.eigh(lap)
    for tag, counts, clock in zip(("native","half"),raw,hand,strict=True):
        cols,observed,overlap,gap = impulse.read(counts,clock)
        field=cols[:,::2]
        kg=impulse.forecasts(16,j,h,observed)[0]
        def hopping(parameters, holds):
            low,high=parameters
            e=low+(high-low)*lam/4
            return (np.cos(np.outer(holds,e))*basis[impulse.SOURCE[0]//2])@basis.T
        def residual(parameters):
            return (hopping(parameters,observed[:61])-field[:61]).ravel()
        starts=((15,17),(17,15),(15,15),(10,20),(20,10),(1,45))
        fits=[least_squares(residual,s,bounds=(0,np.pi/.06)) for s in starts]
        fit=min(fits,key=lambda ret:ret.cost)
        predicted=hopping(fit.x,observed)
        misses=[np.linalg.norm(field-v,axis=1) for v in (kg,predicted)]
        rows.append(dict(material=tag,free_hopping_edges=fit.x.tolist(),min_overlap=overlap,
                         min_projector_gap=gap,space_time=settings(counts),
                         clock_hash=sha256(clock.tobytes()).hexdigest(),
                         training_max=[float(max(v[:61])) for v in misses],
                         held_max=[float(max(v[61:])) for v in misses],
                         held_rms=[float(np.sqrt(np.mean(v[61:]**2))) for v in misses],
                         held_vectors=120, late_errors=[v[61:].tolist() for v in misses]))
    return dict(protocol=sha256((measure.HERE/"IMPULSE.md").read_bytes()).hexdigest(),
                stock=stock["ordered_address_sha256"],calibrated_j=j,calibrated_h=h,rows=rows)


def settings(counts):
    return dict(shape=list(counts.shape),trials_per_setting=impulse.SHOTS,
                array_sha256=sha256(memoryview(counts).cast("B")).hexdigest())


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--save-counts",type=Path)
    cli.add_argument("--acquire",type=Path)
    args=cli.parse_args()
    result=run(args.save_counts)
    if args.acquire:
        with args.acquire.open("x") as out:
            json.dump(result,out,indent=1,allow_nan=False)
            out.write("\n")
    elif not measure.confirmation.same_readings(json.loads(json.dumps(result)),json.loads((measure.HERE/"broadband.json").read_text())):
        raise ArithmeticError("Broadband receipt differs")
    print([(r["material"],r["held_max"]) for r in result["rows"]])
