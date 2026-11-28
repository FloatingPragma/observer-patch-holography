"""The surviving 15-dimer branch and exact invisible connector perturbations."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import sympy as sp

from . import measure


def clamped_nullspace():
    rows = set()
    for length in measure.material.LENGTHS:
        first = (measure.N-length)//2
        inside = measure.material.bits(length-2)
        coefficients = []
        for pin in (0,1):
            ends = np.full(len(inside),pin)
            zeros = np.zeros(len(inside),int)
            upper = np.column_stack((zeros,ends,inside,ends,zeros))
            locations = np.arange(first-1,first+length+1)
            left = np.where(locations%2,upper,1-upper)
            active = 1-(1-left[:,:-1])*left[:,1:]
            full = np.zeros((len(inside),measure.N-1),int)
            full[:,first-1:first+length-1+1] = active
            coefficients.extend(full)
        for v in np.array(coefficients)-coefficients[0]:
            rows.add(tuple(int(x) for x in v))
    matrix = sp.Matrix(sorted(rows))
    null = [[int(x) for x in vec] for vec in matrix.nullspace()]
    independent = [i for i in range(measure.N-1) if not any(matrix[:,i])]
    return dict(distinct_rows=matrix.rows, rank=int(matrix.rank()), nullspace=null,
                individually_invisible=independent)


def run():
    original = measure.equal_degree.ports(0)
    stock = measure.commission.material_stock(0,measure.N,original)
    route = stock["dimers"][1:]
    with TemporaryDirectory() as temp:
        folder = Path(temp)
        def keep(stream, name, counts):
            words = np.asarray(counts,dtype="<u8",order="C")
            np.save(folder/(name+".npy"),words,allow_pickle=False)
            stream.hashes[name]=sha256(memoryview(words).cast("B")).hexdigest()
        with measure.confirmation._ACQUISITION, patch.object(measure.experiment.CountStream,"keep",keep):
            obs = measure.experiment.spectroscopy(len(route),np.random.default_rng(223061))
        physical = dict(route=route,observer=obs)
        readings = measure.spatial.score(physical,folder)
        cal = measure.equal_degree.decode()[0,8][0]
        h,j = cal[1]/2,cal[0]/4
        meter = dict(h=h,j=j,m=2*abs(j-h),c_squared=4*j*h)
        dispersion = None
        if "even_trained_k" in readings:
            dispersion = measure.confirmation.even_hold_dispersion(physical,readings,meter)
    return dict(protocol_sha256=sha256((measure.HERE/"SUBCHAIN.md").read_bytes()).hexdigest(),
                original_stock_sha=stock["ordered_address_sha256"], route=route,
                clamped=clamped_nullspace(),held=readings, dispersion=dispersion, hashes=obs["count_sha256"])


if __name__=="__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--acquire",type=Path)
    args = cli.parse_args()
    result = run()
    if args.acquire:
        with args.acquire.open("x") as dest:
            json.dump(result,dest,indent=1,allow_nan=False)
            dest.write("\n")
    elif not measure.confirmation.same_readings(json.loads(json.dumps(result)),json.loads((measure.HERE/"shortened.json").read_text())):
        raise ArithmeticError("Shortened material evidence differs")
    print(result["held"]["status"],result["clamped"])
