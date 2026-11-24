"""Can an independently addressed terminal cell receive load at the old gap centre?"""

import argparse
import json

import numpy as np
from scipy.linalg import eigh_tridiagonal

from ..row_30_autonomous.channel import apply, matrices, preparation
from .experiment import HERE, shortest_internal

PULSES = 8192
SHOTS = 2**22
SOURCE_CELLS = (48, 97, 145)


def run():
    car, cells, original, _ = preparation()
    four = original.reshape(-1, 4).tolist()
    shortest = [shortest_internal(car, sites[0], sites[-1]) for sites in four]
    reference = [2, 11, 0]
    stationary = np.array([1.0, 2.0, 1.0]) / np.sqrt(6)
    ref_d, ref_off = matrices(car, reference)
    np.testing.assert_allclose(
        apply(ref_d, ref_off, stationary), -4 * stationary, atol=2e-14
    )
    rng = np.random.default_rng(301198)
    records = []
    for route, blocks in (("four", four), ("shortest", shortest)):
        sites = np.concatenate(blocks)
        assert not set(reference).intersection(sites)
        starts = np.cumsum([0] + [len(b) for b in blocks[:-1]])
        for degree in (True, False):
            d, off = matrices(car, sites, degree)
            e, waves = eigh_tridiagonal(d, off)
            assert min(e) >= -16.00000001 and max(e) <= 8.00000001
            for freq in (-4, -6):
                delta = e - freq
                # Exact success operator of forward holds, including their phase.
                filt = np.cos(delta / 16) ** PULSES * np.exp(-1j * PULSES * delta / 16)
                for idx in SOURCE_CELLS:
                    spectrum = waves[starts[idx]] * filt
                    for delay in (0, 16):
                        psi = waves @ (spectrum * np.exp(-1j * delay * e))
                        p = np.add.reduceat(abs(psi) ** 2, starts)
                        ok = np.vdot(spectrum, spectrum).real
                        assert abs(ok - sum(p)) < 2e-12
                        # Draw the complete address distribution, with every
                        # rejected filter history grouped into its own outcome.
                        counts = rng.multinomial(SHOTS, np.r_[p, 1 - ok])
                        records.append(
                            dict(
                                route=route,
                                degree=degree,
                                frequency=freq,
                                source_index=idx,
                                source_port=int(sites[starts[idx]]),
                                extra_hold=delay,
                                accepted_by_cell=counts[:-1].tolist(),
                                rejected=int(counts[-1]),
                            )
                        )
    return dict(
        cells=list(map(int, cells)), shots=SHOTS, pulses=PULSES, records=records
    )


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record", action="store_true")
    opts = cli.parse_args()
    result = run()
    path = HERE / "receiver.json"
    if opts.record:
        with path.open("x") as dest:
            json.dump(result, dest, separators=(",", ":"))
            dest.write("\n")
    elif result != json.loads(path.read_text()):
        raise SystemExit("The complete frequency-gated receiver counts differ.")
    else:
        print(
            "Verified receiver, quiet and filter-reject records for both 194-cell paths."
        )


if __name__ == "__main__":
    main()
