"""Read all four blank endpoint occupations on the historical local law, before any free release."""
import argparse
import json
from math import comb

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from .record_action import addresses, changes, signed_table
from model.executable_model import cell


SIZE = 256
HISTORICAL = (0., .5, 1., 2., 4., 8., 16., 64., 128., -8.)
ADDITIONAL = (7., 9., 10.)
WORDS = ((0,0), (0,1), (1,0), (1,1))


def read(h, coherent=False):
    law = candidate.native_law
    columns = [law.SRC.index(word) for word in WORDS]
    inputs = np.eye(law.DIM, dtype=complex)[:,columns]
    if coherent:
        left, right = inputs[:,1], inputs[:,2]
        inputs = np.column_stack((inputs, (left+right)/np.sqrt(2), (left+1j*right)/np.sqrt(2)))
    routes = np.zeros((SIZE, law.DIM, inputs.shape[1]), dtype=complex)
    routes[0] = inputs
    pulse = expm(-2j*np.pi*h/SIZE)
    span = 1
    while span < SIZE:
        routes[span:2*span] = pulse @ routes[:span]
        span *= 2
        if span < SIZE:
            pulse = pulse @ pulse
    # Trace over every one of the 99 final seam configurations on each output.
    probabilities = (abs(np.fft.ifft(routes, axis=0))**2).sum(axis=1).T
    np.testing.assert_allclose(probabilities.sum(axis=1), 1., atol=2e-13, rtol=0)
    return probabilities


def kernel(energy):
    # Complete digital clock response, including fractional and negative energies.
    times = np.arange(SIZE)*2*np.pi/SIZE
    return abs(np.fft.ifft(np.exp(-1j*energy*times)))**2


def panel(arm, exchange, contrast):
    if arm == "pair-only":
        # The rival costs +16 precisely on 11 in the occupied binary material.
        # Its action outside the blank binary preparation is not used.
        occupied = np.eye(len(exchange))[:,candidate.native_law.SRC.index((1,1))]
        h = exchange + candidate.KAPPA*contrast + 16*np.outer(occupied, occupied)
        return read(h)
    return read(exchange + float(arm)*contrast)


def coherent_check(exchange, contrast, cost):
    """Complete tomography of each clock-bin effect on the blank one-load block."""
    native = exchange + candidate.KAPPA*contrast
    law = candidate.native_law
    left, right = (law.SRC.index(w) for w in WORDS[1:3])
    selected = [left, right]
    results = []
    nominal = read(native, coherent=True)
    for name, phase in (("native", 1.), ("reversed_exchange", -1.), ("imaginary_exchange", 1j)):
        h = native.astype(complex)
        h[left, right] *= phase
        h[right, left] = h[left, right].conjugate()
        p = read(h, coherent=True)
        mean = (p[1]+p[2])/2
        off = p[4]-mean + 1j*(mean-p[5])
        recovered = np.array([[p[1], off], [off.conj(), p[2]]]) @ np.arange(SIZE)
        np.testing.assert_allclose(recovered, h[np.ix_(selected,selected)], atol=3e-12, rtol=0)
        forecast = cost*np.eye(2) + np.array([[.5,-.5],[-.5,.5]])
        agrees = np.linalg.norm(recovered-forecast, ord=2) < 1e-10
        sharp_TV = float(np.max(np.sum(abs(p[:4]-nominal[:4]), axis=1))/2)
        if sharp_TV > 1e-10 or agrees != (name == "native"):
            raise ArithmeticError("Coherent contact tomography missed an invisible sharp-input rival")
        results.append(dict(arm=name, sharp_TV=sharp_TV,
                            recovered_off_diagonal=[float(recovered[0,1].real), float(recovered[0,1].imag)],
                            coherent_TV=(np.sum(abs(p[4:]-nominal[4:]), axis=1)/2).tolist()))
    return results


def echo_reference():
    """Two ordinary free exposures with a stationary reference in the same load sector."""
    law = candidate.native_law
    seeds = np.eye(law.DIM)[:, [law.SRC.index(word) for word in WORDS]]
    reference = law.echo() @ seeds
    exchange, contrast = law.local_term_parts()
    np.testing.assert_allclose(exchange@reference, 0, atol=3e-15, rtol=0)
    np.testing.assert_allclose(contrast@reference, 0, atol=3e-15, rtol=0)
    np.testing.assert_allclose(seeds.T@reference, 0, atol=3e-15, rtol=0)
    inputs = (seeds+reference)/np.sqrt(2.)
    exposures = (1/3, np.sqrt(2)/3)
    data = {}
    for k in (*HISTORICAL, *ADDITIONAL, 264., "pair-only"):
        strength = candidate.KAPPA if k == "pair-only" else k
        h = exchange+strength*contrast
        if k == "pair-only":
            h = h+16*np.outer(seeds[:,3],seeds[:,3])
        rows = []
        for t in exposures:
            final = expm(-1j*t*h)@inputs
            # Each read has two interfering output rays and the full orthogonal complement.
            reads = []
            for phase in (1, 1j):
                ports = []
                for sign in (1,-1):
                    detector = (seeds+sign*phase*reference)/np.sqrt(2.)
                    ports.append(abs(np.sum(detector.conj()*final, axis=0))**2)
                other = final.copy()
                for sign in (1,-1):
                    detector = (seeds+sign*phase*reference)/np.sqrt(2.)
                    other -= detector*np.sum(detector.conj()*final, axis=0)
                probabilities = np.array([*ports, (abs(other)**2).sum(axis=0)])
                np.testing.assert_allclose(probabilities.sum(axis=0), 1., atol=2e-13, rtol=0)
                reads.append(probabilities)
            amplitude = reads[0][0]-reads[0][1] - 1j*(reads[1][0]-reads[1][1])
            one = np.exp(-1j*t*strength)*(1+np.exp(-1j*t))/2
            expected = [1.,one,one,np.exp(-16j*t) if k=="pair-only" else 1.]
            np.testing.assert_allclose(amplitude, expected, atol=2e-13, rtol=0)
            rows.append(np.array(reads).tolist())
        data[str(k)] = rows
    return dict(exposures=list(exposures), words=WORDS, reference="native echo of each sharp input", outcomes=data)


def compare():
    board = cell()
    edge = board.seams[0]
    probe = np.zeros(board.n, dtype=int)
    probe[edge[0]] = 1
    cost = sum(bool(d) for d in changes(addresses(board, edge), edge, probe, signed_table()))
    expected = np.stack([kernel(0), (kernel(cost)+kernel(cost+1))/2,
                         (kernel(cost)+kernel(cost+1))/2, kernel(0)])
    exchange, contrast = candidate.native_law.local_term_parts()
    results = []
    for arm in (*HISTORICAL, *ADDITIONAL, "pair-only"):
        raw = panel(arm, exchange, contrast)
        k = candidate.KAPPA if arm == "pair-only" else arm
        equal = [0, 16 if arm == "pair-only" else 0]
        lines = [[equal[0]], [k,k+1], [k,k+1], [equal[1]]]
        reconstruction = np.array([np.mean([kernel(e) for e in row], axis=0) for row in lines])
        defect = float(np.max(np.sum(abs(raw-reconstruction), axis=1))/2)
        if not defect < 1e-10:
            raise ArithmeticError("Stored line representation loses native detector outcomes")
        distances = (np.sum(abs(raw-expected), axis=1)/2).tolist()
        # Sixteen independent 01 preparations: all in the two predicted bins
        # with both seen. One 11 preparation must also read zero. Sum mutually
        # exclusive counts to avoid cancellation on near-zero controls.
        x, y = raw[1, cost:cost+2]
        acceptance = raw[3,0]*sum(comb(16,m)*x**m*y**(16-m) for m in range(1,16))
        results.append(dict(arm=arm, equal_weight_lines=lines, four_word_TV=distances,
                            seventeen_read_acceptance=float(acceptance),
                            max_reconstruction_TV=defect, meets_source=max(distances) < 1e-10))
    matching = [r["arm"] for r in results if r["meets_source"]]
    if matching != [candidate.KAPPA]:
        raise ArithmeticError("The local record forecast did not single out the historical setting in its sweep")
    alias = panel(cost+SIZE, exchange, contrast)
    alias_tv = float(np.max(np.sum(abs(alias-expected), axis=1))/2)
    if not alias_tv < 1e-10:
        raise ArithmeticError("The out-of-range period shift should be invisible to this digital clock")
    return dict(bin_count=SIZE, words=WORDS, predicted_cost=cost, alternatives=results,
                outside_simple_graph=dict(strength=cost+SIZE, max_TV=alias_tv),
                coherent_probes=coherent_check(exchange, contrast, cost))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    modes = p.add_mutually_exclusive_group()
    modes.add_argument("--echo-reference", action="store_true", help="ordinary free exposures, both quadratures and quiet")
    modes.add_argument("--raw", choices=(*map(str,HISTORICAL), *map(str,ADDITIONAL), "pair-only"), help="print all 4 × 256 probabilities")
    args = p.parse_args()
    if args.echo_reference:
        result = echo_reference()
    elif args.raw is not None:
        result = panel(args.raw, *candidate.native_law.local_term_parts()).tolist()
    else:
        result = compare()
    print(json.dumps(result, allow_nan=False))
