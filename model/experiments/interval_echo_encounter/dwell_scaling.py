"""Fixed memory charges and the rate of the clock written by a native pair source."""
from fractions import Fraction
from hashlib import sha256
import json

import numpy as np
from scipy import sparse
from scipy.optimize import brentq
from scipy.sparse.linalg import cg

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.experiments.seam_program.gate import receipt
from .dwell import HERE, plan
from .dwell_source import Station, require


def memory_charges():
    """Integer algebra, with twice the written-range and memory projectors."""
    write = law._write_pattern()
    permutation, squared = law._label_ops()
    sign = np.kron(permutation,np.eye(9,dtype=int))
    square = np.kron(squared,np.eye(9,dtype=int))
    written = write@write.T
    twice_complement = 2*np.eye(90,dtype=int)-written
    charges = {}
    # The two generator parts act inside the memory labels before the echo.
    for m in (0,1,2):
        select = np.diag(np.diag(squared)==m*m).astype(int)
        for parity in (1,-1):
            charge = select@(np.eye(10,dtype=int)+parity*permutation)
            if not charge.any():
                continue
            require(np.array_equal(charge@charge,2*charge),"invalid sign projector")
            full=np.kron(charge,np.eye(9,dtype=int))
            for G in (np.eye(90,dtype=int)-sign,square):
                mixed=write.T@G@twice_complement
                dirty=twice_complement@G@twice_complement
                require(not mixed.any() and np.array_equal(dirty@full,full@dirty),
                        "a local free term changes the memory charge")
            charges[f"{m}_{parity}"]=int(np.trace(charge)//2)
    require(1+sum(charges.values())==11,"missing seam sector")
    return charges


def coefficient(station):
    row=station.kinetic.getrow(station.p)
    exits=row.indices[row.indices!=station.p]
    require(np.all(row.data[row.indices!=station.p]==-.5),"unexpected pair move")
    gaps=(4+station.tag_loss+station.cost[exits]).astype(int)
    require(np.all(gaps>0),"degenerate contrast minimum")
    return sum(Fraction(1,4*int(g)**2) for g in gaps)


def secular(matrix,p,gap):
    """Solve the scalar Feshbach equation below the complementary band."""
    keep = np.arange(matrix.shape[0])!=p
    seed = matrix[p,p]
    coupling = matrix[keep,p].toarray().ravel()
    complement = matrix[keep][:,keep]
    eye=sparse.eye(complement.shape[0],format="csr")
    def tail(energy):
        x,info=cg(complement-energy*eye,coupling,rtol=5e-14,atol=0,maxiter=200)
        require(info==0,"complement resolvent did not converge")
        return x
    energy=brentq(lambda z: z-seed+np.dot(coupling,tail(z)),seed-np.linalg.norm(coupling),seed,
                  xtol=2e-14,rtol=1e-15)
    x=tail(energy)
    v=np.empty(matrix.shape[0])
    v[p],v[keep]=1.,-x
    residual=np.linalg.norm(matrix@v-energy*v)/np.linalg.norm(v)
    require(residual<1e-9 and abs(energy-seed)<gap/2,"secular root not certified in isolated interval")
    return energy,x


def run(spec, base):
    b=Station(spec["level"],base["tag"])
    b.flight()
    C=coefficient(b)
    M=max(b.d,2*max(b.graph.degree)+1-b.d+np.sqrt(b.d/2))
    delta_min=4+b.tag_loss-b.saving
    out=dict(local_charge_ranks=memory_charges(),coefficient=str(C),levels={})
    for kappa in spec["couplings"]:
        branches=[secular(b.generator(kappa,j,"native"),b.p,b.isolation(kappa,j)["gap"]) for j in (0,1)]
        tail_overlap=np.dot(branches[0][1],branches[1][1])
        rate=tail_overlap/(1+tail_overlap)
        require(abs(rate-(branches[1][0]-branches[0][0]))<2e-9,"overlap identity fails")
        rho=M/(kappa*delta_min)
        require(0<rho<1,"coupling outside the finite enclosure")
        lo,hi=[float(C)/(kappa**2*(1+sign*rho)**2+float(C)) for sign in (1,-1)]
        require(lo<=rate<=hi,"observed rate outside finite strong-contrast bounds")
        out["levels"][str(kappa)]=dict(rate=float(rate),scaled_rate=float(kappa**2*rate),
            relative_to_limit=float(kappa**2*rate/float(C)-1),bounds=[float(lo),float(hi)])
    return out


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record",action="store_true")
    args=parser.parse_args()
    files=("dwell_scaling.py","dwell_scaling_protocol.json","dwell_sources.json")
    pins=json.loads((HERE/"dwell_scaling_sources.json").read_text(encoding="utf8"))
    require(set(pins)==set(files) and all(pins[n]==sha256((HERE/n).read_bytes().replace(b"\r\n",b"\n")).hexdigest()
                                          for n in files),"changed scaling source")
    spec=json.loads((HERE/"dwell_scaling_protocol.json").read_text(encoding="utf8"))
    out=run(spec,plan())
    receipt(HERE/"dwell_scaling_result.json",out,args.record)
    print("Strong-contrast rate coefficient:",out["coefficient"])


if __name__=="__main__":
    main()
