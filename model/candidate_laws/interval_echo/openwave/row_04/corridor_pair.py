"""A two-load corridor sealed by native zero-contrast records in the scored cell."""
from functools import lru_cache
from itertools import combinations
import json
from math import isfinite
from numbers import Real

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation


PATH = (156, 157, 164)
BASIS = tuple(combinations(PATH, 2))
LIVE = (BASIS[0], BASIS[2])
LAW = candidate.native_law


@lru_cache(maxsize=1)
def geometry():
    """Installed neutral edges and the thirty ordered M10.3 receiver edges."""
    board = federation(0)
    ports = set(PATH)
    home = tuple((a, b) for a, b in board.seams if a//12 == b//12 == 13)
    live = tuple(e for e in board.seams if set(e) <= ports)
    neutral = tuple(e for e in board.seams if len(set(e) & ports) == 1)
    if (len(home) != 30 or live != LIVE or len(neutral) != 12
            or tuple(board.degree[list(PATH)]) != (5, 6, 5)):
        raise ValueError("The scored home does not contain the declared corridor")
    return neutral, home


def neutral_column(loads):
    """Native record at a cut seam, with its two endpoint loads in seam order."""
    state = np.zeros(LAW.DIM)
    for coin in (0, 1):
        row = LAW.N_SRC * (1 + LAW.LABELS.index((coin, 0))) + LAW.SRC.index(loads)
        state[row] = 1/np.sqrt(2)
    return state


def _close(actual, expected, name):
    if not np.allclose(actual, expected, atol=5e-14, rtol=0):
        raise ArithmeticError(name)


@lru_cache(maxsize=1)
def certify_native_collar():
    """Prove annihilation from the native write support and check the physical queries."""
    # U=[[0,W*],[W,I-WW*]]. A d=0 record has overlap only with writes
    # whose entire support has d=0. U sends it to blank plus d=0 records.
    # Both parts of G annihilate this space, hence U G U* kills it for all κ.
    zero = np.repeat([d == 0 for _, d in LAW.LABELS], LAW.N_SRC)
    pattern = LAW._write_pattern()
    touched = pattern[zero].any(axis=0)
    if (pattern[~zero][:, touched].any()
            or set(np.flatnonzero(touched)) != {i for i, (a, b) in enumerate(LAW.SRC) if a == b}):
        raise ArithmeticError("Neutral write columns mix different absolute contrasts")
    r, d2 = LAW._label_ops()
    kept = [i for i, (_, d) in enumerate(LAW.LABELS) if d == 0]
    if not np.array_equal(r[:, kept], np.eye(LAW.N_LAB, dtype=int)[:, kept]) or d2[:, kept].any():
        raise ArithmeticError("The record generator moves a neutral mark")
    echo = LAW.echo()
    neutral_rows = np.r_[np.zeros(LAW.N_SRC, dtype=bool), zero]
    if neutral_rows.sum() != 18:
        raise ArithmeticError("Expected both neutral coin states for all nine load words")
    for part in LAW.local_term_parts():
        _close(part[:, neutral_rows], 0, "A written neutral collar exchanges or phases a load")

    blank = np.arange(LAW.DIM) < LAW.N_SRC
    pair = np.arange(LAW.DIM) % LAW.N_SRC == LAW.SRC.index((1, 1))
    vacuum = np.eye(LAW.DIM)[:, LAW.SRC.index((0, 0))]
    _close(echo @ vacuum, neutral_column((0, 0)), "Empty native write differs from the neutral source")
    # Copying the ordinary M10.3 material pair flag between echo and undo
    # has the same two selective operations on blank and neutral collars.
    for loads in ((0, 0), (0, 1), (1, 0), (1, 1)):
        marked = neutral_column(loads)
        sharp = np.eye(LAW.DIM)[:, LAW.SRC.index(loads)]
        for state in (sharp, marked):
            written = echo @ state
            hit = loads == (1, 1)
            _close(echo.T @ (pair*written), hit*state, "Terminal pair click changed the binary state")
            _close(echo.T @ (~pair*written), (not hit)*state, "Terminal quiet changed the binary state")
        if sum(loads) <= 1:
            # On a cut with a known vacant exterior endpoint this measures
            # whether the interior is empty while restoring the sealed collar.
            empty = loads == (0, 0)
            written = echo @ marked
            _close(echo.T @ (blank*written), empty*marked, "Incorrect vacancy branch")
            _close(echo.T @ (~blank*written), (not empty)*marked, "Incorrect occupied branch")


def seal_path(path, neutral):
    """Write each cut with both ends vacant, using the original two units."""
    certify_native_collar()
    a, b, c = path
    live = {tuple(sorted(edge)) for edge in ((a,b), (b,c))}
    loads = {a, b}
    written, swaps = [], []
    route = (c, b, a, b, c)
    for i, hole in enumerate(route):
        if set(path)-loads != {hole}:
            raise ArithmeticError("The sealing head did not reach a vacancy")
        for edge in neutral:
            if hole in edge and edge not in written:
                if loads.intersection(edge):
                    raise ArithmeticError("A boundary echo would write a nonzero contrast")
                written.append(edge)
        if i+1 < len(route):
            following = route[i+1]
            edge = tuple(sorted((hole, following)))
            if edge not in live or following not in loads:
                raise ArithmeticError("The sealing swap misses a live contact")
            # On 10/01 a π pulse is SWAP times exp(-iπκ); only one sharp
            # configuration is present at each preparation contact.
            swaps.append(edge)
            loads.symmetric_difference_update(edge)
    if set(written) != set(neutral) or loads != {a, b}:
        raise ArithmeticError("Incomplete cut or changed pair launch")
    return {"neutral_edges": tuple(written), "live_swaps": tuple(swaps), "final_pair": tuple(sorted((a, b)))}


def seal_from_sharp():
    return seal_path(PATH, geometry()[0])


@lru_cache(maxsize=1)
def corridor_parts():
    """Native exchange and contrast reduced on X=ab, Y=ac, Z=bc, including references."""
    board = federation(0)
    neutral, _ = geometry()
    certify_native_collar()
    blank = tuple(edge for edge in board.seams if edge not in neutral)
    exchange = np.zeros((3, 3))
    contrast = np.zeros((3, 3))
    for column, occupied in enumerate(BASIS):
        for edge in blank:
            word = tuple(int(v in occupied) for v in edge)
            if sum(word) != 1:
                continue
            after = tuple(sorted(set(occupied).symmetric_difference(edge)))
            if after not in BASIS:
                raise ArithmeticError("An unsealed blank seam carries a load outside the path")
            exchange[column, column] += .5
            exchange[BASIS.index(after), column] -= .5
            contrast[column, column] += 1
        contrast[column, column] += sum(6-board.degree[v] for v in occupied)
    exchange.setflags(write=False)
    contrast.setflags(write=False)
    return exchange, contrast


def corridor_generator(kappa):
    if not isinstance(kappa, Real) or isinstance(kappa, (bool, np.bool_)) or not isfinite(kappa):
        raise ValueError("Expected a finite real agreement strength")
    exchange, contrast = corridor_parts()
    return exchange + kappa*contrast


def receiver_bins(populations):
    """Thirty enrolled home seams, in the original order, followed by quiet."""
    _, home = geometry()
    weights = np.asarray(populations)
    if (weights.shape != (3,) or weights.dtype.kind not in "fiu"
            or not np.isfinite(weights).all() or np.any(weights < 0) or np.any(weights > 1+2e-12)
            or weights.sum() > 1+2e-12):
        raise ValueError("Expected three subnormalized material weights")
    bins = np.zeros(31)
    for word, weight in zip(BASIS, weights, strict=True):
        bins[home.index(word) if word in home else -1] += weight
    return bins


def exact_flight():
    """Both parts of the generator preserve the antisymmetric detected pair."""
    setup = seal_from_sharp()
    exchange, contrast = corridor_parts()
    dark = np.array([1, 0, -1])
    if not (np.array_equal(2*exchange @ dark, dark) and np.array_equal(contrast @ dark, 2*dark)):
        raise ArithmeticError("The full two-load law does not preserve the dark pair")
    return {**setup, "twice_exchange": (2*exchange).astype(int).tolist(),
            "contrast": contrast.astype(int).tolist(), "dark_terminal": receiver_bins(dark**2/2).tolist()}


if __name__ == "__main__":
    print(json.dumps(exact_flight(), indent=2, allow_nan=False))
