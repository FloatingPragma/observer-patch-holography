"""Full two-parity head tomography from charge-hidden radar return timestamps."""
import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from .. import certificate as c
from .. import station
from . import clock_history as clock
from . import measure

# Numerical pivoting proposed the unweighted rows; only the rational inverse proves rank.
# Here an index is 11 * route index + (first completed return tick - 1).
BINS = (6,38,129,2,10,97,124,43,4,131,93,54,9,123,90,35,7,41,98,125,42,37,
        94,86,40,39,130,126,85,84,83)
AXES = ("IX","IY","IZ","XI","XX","XY","XZ","YI","YX","YY","YZ","ZI","ZX","ZY","ZZ")


def affine_rows():
    """Return all 132 pooled timestamp rows and tracial offsets with one denominator."""
    labels = 2*station.iso.basis(False)
    re,im = labels.real.astype(int), labels.imag.astype(int)
    assert np.array_equal(labels,re+1j*im)
    # Row ordering is route-major, ticks 1..11. Tick 1 cannot contain two successive hits.
    tab=np.zeros((132,31),dtype=object)
    offsets=np.zeros(132,dtype=object)
    for _,parity,i,t,effect in clock.native_events():
        row=11*i+t-1
        moments=(re*effect[0]+im*effect[1]).sum(axis=(1,2))
        tr=int(np.trace(effect[0]))
        # Charge average 1/2, Pauli coefficient 1/2, parity-coordinate coefficient
        # +/-1/4 and tracial reference weight 1/8. The common denominator is 16*F.
        tab[row,15*parity:15*(parity+1)] += 4*moments
        tab[row,-1] += (2 if parity else -2)*tr
        offsets[row] += tr
    return tab,offsets,16*clock.F


def run():
    c.iv.dps=30
    tab,offsets,denom=affine_rows()
    num=tab[list(BINS),:]
    a=sp.Matrix(num)/denom
    eps,_=clock.timing_error()
    # The Pauli/2 operators have trace norm 2, as does the last full two-block
    # operator (-I/4,+I/4). Each unit-trace tracial origin has offset error <= eps.
    bounds=clock.square_certificate(a,c.ri(62*sp.Rational(eps)))
    assert bounds["inverse_operator_norm_upper"] < 2e7
    selected=[[*measure.ROUTES[i//11],i%11+1] for i in BINS]
    return dict(normalized_head_algebra="M4 direct sum M4",candidate_pooled_pre_cap_bins=132,
                rank=31,free_real_state_parameters=31,source_reference="I8/8",
                kept_clock_bins=selected,
                coordinate_order=dict(even_then_odd_Pauli_over_two=list(AXES),
                    last_coordinate="odd probability minus 1/2",
                    last_direction_even_odd=["-I4/4","+I4/4"]),
                affine_clock=dict(denominator=str(denom),
                    numerators=[[str(int(v)) for v in row] for row in num],
                    offset_numerators=[str(int(offsets[i])) for i in BINS],
                    entry_error_upper=float(np.nextafter(float(2*eps),np.inf)),
                    offset_error_upper=float(np.nextafter(float(eps),np.inf))),
                certificate=bounds)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    s=json.dumps(run(),indent=2,allow_nan=False)+"\n"
    if parser.parse_args().write:
        Path(__file__).with_suffix(".json").write_text(s)
    else:
        print(s,end="")
