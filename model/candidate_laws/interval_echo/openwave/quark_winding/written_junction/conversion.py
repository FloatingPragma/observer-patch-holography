"""Charge-changing core conversion and native sign/coin reads on a three-wall preparation."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import currents
from model.executable_model import carrier
from ..belt_walls import belt
from . import exchange


def run():
    board = carrier.federation(0)
    tips = [(a, b) for a, b in board.external if a // 12 == 0]
    layouts = []
    for _, b in tips:
        order = next(o for o in belt.belts().values() if b % 12 in o)
        layouts.append([12 * (b // 12) + p for p in order])
    site = 1 << layouts[0].index(tips[0][1])
    minus = next(
        w
        for w in sorted(exchange.walls(3))
        if not w & site and (w | site) in exchange.walls(4)
    )
    spectators = [min(exchange.walls(n)) for n in (3, 4)]
    inputs = [("capture", 1, 0), ("reverse", 0, 1)]
    x, c = law.local_term_parts()
    echo = law.echo()
    panels = []
    for form, a, b in inputs:
        for reflected in (False, True):
            pair = (b, a) if reflected else (a, b)
            source = np.eye(99)[:, law.SRC.index(pair)]
            for kappa in (8, 0):
                h = x + kappa * c
                for service, psi0 in [
                    ("blank", source),
                    ("write", echo @ source),
                    ("erase", echo @ echo @ source),
                ]:
                    for t in (0.0, 0.25, 1.0, np.pi):
                        psi = expm(-1j * t * h) @ psi0
                        probabilities = abs(psi.reshape(11, 9)) ** 2
                        output = {}
                        for m, i in zip(*np.where(probabilities > 1e-20)):
                            loads = law.SRC[i]
                            if reflected:
                                loads = loads[::-1]
                            if m:
                                coin, sign = law.LABELS[m - 1]
                                if reflected:
                                    coin, sign = 1 - coin, -sign
                                tag = f"{coin},{sign}"
                            else:
                                tag = "blank"
                            key = f"{loads[0]},{loads[1]}|{tag}"
                            output[key] = float(probabilities[m, i])
                        assert abs(sum(output.values()) - 1) < 1e-11
                        conv = sum(
                            p
                            for key, p in output.items()
                            if key.startswith(f"{b},{a}|")
                        )
                        panels.append(
                            dict(
                                input=form,
                                reflected=reflected,
                                kappa=kappa,
                                service=service,
                                time=t,
                                converted=conv,
                                records=output,
                            )
                        )
    return currents._portable(
        dict(
            central_owner=0,
            bridges=tips,
            layouts=layouts,
            negative_word=minus,
            positive_word=minus | site,
            spectator_words=spectators,
            panels=panels,
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    encoded = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as f:
        f.write(encoded)
