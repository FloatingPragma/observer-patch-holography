"""Native-frequency transmission on screened and unscreened geometric routes."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

from ..row_30_autonomous.channel import apply, matrices
from ..row_30_route.receiver import PULSES, SHOTS
from .experiment import HERE, portable
from model.executable_model import federation


def reference(car, cells):
    ordinary = [v for v in range(len(car.degree)) if v//12 not in cells and car.degree[v]==5]
    for centre in ordinary:
        ns = sorted(v for v in car.neighbours[centre] if car.degree[v]==5)
        if len(ns) >= 2:
            sites = [ns[0], centre, ns[1]]
            d, h = matrices(car, sites)
            ref = np.array([1.,2.,1.])/np.sqrt(6)
            np.testing.assert_allclose(apply(d,h,ref),-4*ref,atol=1e-14)
            return sites
    raise ArithmeticError('No disjoint stationary reference')


def run():
    census = json.loads((HERE/'junction.json').read_text())
    rng = np.random.default_rng(320932)
    rows = []
    for level in (0,1,2):
        car = federation(level)
        for tie,variants in census[str(car.cells)].items():
            for selector in ('longest', 'longest_whole_trimer'):
                route = variants[selector]
                sites=np.array(route['sites'])
                cells=route['cells']
                starts=np.r_[0,np.flatnonzero(np.diff(sites//12))+1]
                middle=len(cells)//2
                ref=reference(car,set(cells))
                d,off=matrices(car,sites)
                energies, v=eigh_tridiagonal(d,off)
                for freq in (-10,-1,-6):
                    # Phase and amplitude of the accepted read history. All other
                    # histories are one rejected outcome of the very same trial.
                    delta=energies-freq
                    f=np.cos(delta/16)**PULSES*np.exp(-1j*PULSES*delta/16)
                    psi=v @ (v[starts[middle]] * f)
                    percell=np.add.reduceat(abs(psi)**2,starts)
                    p=percell.sum()
                    assert p <= 1 and abs(p-np.sum(abs(v[starts[middle]]*f)**2))<1e-11
                    counts=rng.multinomial(SHOTS,np.r_[percell,1-p])
                    rows.append(dict(board=car.cells,tie=tie,selector=selector,frequency=freq,
                        reference=ref,source=cells[middle],receiver=cells[middle+3],
                        accepted=counts[:-1].sum(),rejected=counts[-1],
                        all_accepted_cells=counts[:-1],receiver_clicks=counts[middle+3],
                        receiver_probability=percell[middle+3]))
    return portable(dict(pulses=PULSES,shots=SHOTS,rows=rows))


if __name__ == '__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--output',type=Path)
    args=cli.parse_args()
    r=run()
    if args.output:
        with args.output.open('x') as dest:
            json.dump(r,dest,separators=(',',':'),allow_nan=False)
            dest.write('\n')
    elif r != json.loads((HERE/'filter.json').read_text()):
        raise SystemExit('Frequency-resolved route receipt differs from archive')
    for row in r['rows']:
        if row['tie']=='low':
            print(row['board'],row['selector'],row['frequency'],row['accepted'],row['receiver_clicks'])
