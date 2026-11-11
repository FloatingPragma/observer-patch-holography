"""One stream of quenched SU(3) configurations (Wilson action). Saves configurations and the sweep series.

    python -m model.quenched_hadron_binding.production.generate NS NT BETA STREAM NTHERM NSEP NCFG OUTDIR [SEEDBASE]

Even stream numbers start cold, odd ones hot. One compound sweep = 1 heat bath + 4 overrelaxation sweeps.
Files: OUTDIR/cfg_SS_K.npy (links, complex128, shape (V, 4, 3, 3)), OUTDIR/series_SS.npz (plaquette and Polyakov
loop after every compound sweep). The gauge kernel is model/colour_deconfinement/su3_numba.py. Needs numba.
"""
import json
import os
import sys
import time

import numpy as np

from ...colour_deconfinement import su3_numba as S

N_OR = 4


def main():
    if not S.AVAILABLE:
        raise SystemExit('the production needs numba')
    a = sys.argv
    ns, nt, beta, stream = int(a[1]), int(a[2]), float(a[3]), int(a[4])
    ntherm, nsep, ncfg, out = int(a[5]), int(a[6]), int(a[7]), a[8]
    seedbase = int(a[9]) if len(a) > 9 else 7000
    os.makedirs(out, exist_ok=True)
    dims = (ns, ns, ns, nt)
    up, dn = S.neighbours(*dims)
    V = ns ** 3 * nt
    seed = seedbase + 97 * ns + 13 * nt + stream
    S.seed_rng(seed)
    start = 'cold' if stream % 2 == 0 else 'hot'
    U = S.cold_start(V) if start == 'cold' else S.hot_start(V)
    nsw = ntherm + nsep * ncfg
    ps = np.zeros(nsw); pt = np.zeros(nsw); L = np.zeros(nsw, np.complex128)
    t0 = time.time()
    k = 0
    for n in range(nsw):
        S.compound_sweeps(U, up, dn, beta, 1, N_OR)
        if (n + 1) % 50 == 0:
            S.reunitarize(U)
        ps[n], pt[n] = S.plaquette(U, up)
        L[n] = S.polyakov(U, up, nt)
        if n + 1 > ntherm and (n + 1 - ntherm) % nsep == 0:
            S.reunitarize(U)
            tmp = os.path.join(out, f'cfg_{stream:02d}_{k}.tmp.npy')
            np.save(tmp, U)
            os.replace(tmp, os.path.join(out, f'cfg_{stream:02d}_{k}.npy'))
            k += 1
            meta = dict(ns=ns, nt=nt, beta=beta, stream=stream, seed=seed, start=start, ntherm=ntherm, nsep=nsep,
                        ncfg=ncfg, done=k, n_or=N_OR, seconds=time.time() - t0)
            tmp = os.path.join(out, f'series_{stream:02d}.tmp.npz')
            np.savez(tmp, ps=ps[:n + 1], pt=pt[:n + 1], L=L[:n + 1], meta=json.dumps(meta))
            os.replace(tmp, os.path.join(out, f'series_{stream:02d}.npz'))
    print('done', out, stream, time.time() - t0, flush=True)


if __name__ == '__main__':
    main()
