"""Held predictions at three radii, including a rescaled local charge vertex."""

import argparse
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigvalsh
from scipy.sparse import diags
from model.executable_model import EMITTER_SEAM, federation
from ..free_bath import native
from . import readback, virial

OUTPUT = Path(__file__).with_name("confirmation.json.gz")


def run():
    rng = np.random.Generator(np.random.PCG64(23042917).jumped(4))
    training = json.loads(gzip.decompress(virial.OUTPUT.read_bytes()))
    blank = training["calibration"]["8"]
    gain = blank["gains"]
    rows = {}
    for level in (0, 1):
        board = federation(level)
        home = EMITTER_SEAM[level][0]
        for shell in (1, 2, 3):
            probe = training["panels"][f"{board.cells},8,{shell}"]["probe"]
            signs = {}
            for s in (-1, 1):
                for t in (-1, 1):
                    h, free = native.vacancy_operator(
                        board, 8, {home: s + 1, probe: t + 1}
                    )
                    n = len(free)
                    enc = virial.bounds(h)
                    band = (int(np.floor(enc[0])), int(np.ceil(enc[1])))
                    lines = eigvalsh(h.toarray(), check_finite=False)
                    bank = virial.trace_bank(lines, rng, 256 * virial.SHOTS)
                    baseline = virial.predicted_generator(
                        board, home, probe, 8, blank, gain, s, t
                    )
                    field = gain[str(t)] * np.isin(free, board.neighbours[probe])
                    pred = {}
                    for eta in (0.9, 1.0, 1.1):
                        mat = baseline + diags((eta - 1) * field)
                        pred[str(eta)] = eigvalsh(mat.toarray(), check_finite=False)
                    readings = {}
                    for beta in virial.BETAS:
                        trace, error = virial.thermal(bank, beta, n, enc, band=band)
                        dh = 15.1 * blank["line_error"]
                        readings[str(beta)] = dict(
                            trace=trace,
                            error=error,
                            forecasts={
                                eta: float(np.exp(-beta * spec).sum())
                                for eta, spec in pred.items()
                            },
                            forecast_error=float(
                                n
                                * beta
                                * dh
                                * np.exp(
                                    -beta * (enc[0] - max(abs(field)) * 0.1 - 2 * dh)
                                )
                            ),
                        )
                    signs[f"{s},{t}"] = dict(bank=bank, band=band, readings=readings)
            summaries = {}
            for beta in virial.BETAS:
                c = {key: row["readings"][str(beta)] for key, row in signs.items()}
                value = (
                    -readback.signed({key: r["trace"] for key, r in c.items()}) / beta
                )
                error = sum(b["error"] for b in c.values()) / (4 * beta)
                tol = sum(b["forecast_error"] for b in c.values()) / (4 * beta)
                forecasts = {
                    eta: -readback.signed(
                        {key: r["forecasts"][eta] for key, r in c.items()}
                    )
                    / beta
                    for eta in ("0.9", "1.0", "1.1")
                }
                summaries[str(beta)] = dict(
                    value=value,
                    error=error,
                    forecasts=forecasts,
                    forecast_error=tol,
                    agrees={
                        eta: bool(abs(target - value) <= error + tol)
                        for eta, target in forecasts.items()
                    },
                )
            rows[f"{board.cells},{shell}"] = dict(
                probe=probe, signs=signs, readings=summaries
            )
    return dict(rows=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.dumps(
        run(), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    if args.output:
        args.output.write_bytes(gzip.compress(data, mtime=0))
    elif gzip.decompress(OUTPUT.read_bytes()) != data:
        raise ValueError("Confirmation differs from deposited readout")
    else:
        print("The held contact-strength experiment replayed.")


if __name__ == "__main__":
    main()
