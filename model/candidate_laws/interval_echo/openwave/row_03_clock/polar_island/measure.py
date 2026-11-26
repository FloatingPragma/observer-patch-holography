"""Free-release joint dimer reads on the polar ring. Times use native exchange units."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from model.experiments.material_inertia.pair_spectrum.neutral_guide.polar_ring import (
    fringes,
    measure as ring,
)
from model.candidate_laws.interval_echo.openwave.row_22_confinement.mobile_string import (
    material,
)

SEED = 3031927
SHOTS = 2**38
STEP = 1 / 128
WINDOWS = (32768, 524288)
SIDES = (-1 / 16, 0, 1 / 16)
SOURCES = ("one_upper", "two_alternative_sites", "all_lower")
CONTROLS = ("free", "ground_dephased", "selected_island")
BRANCHES = ("ground_island", "ground_other", "within_island", "remaining_excited")


def reference(panel):
    """Decode the previous instrument's independent, sharp-word return hands."""
    pairs = zip(
        panel["sharp_lower_return"], panel["clocks"][0]["plus_counts"], strict=True
    )
    signals = [
        fringes.amplitude(h) * fringes.amplitude(r).conjugate() for h, r in pairs
    ]
    return fringes.unwrap(ring.HOLDS, signals)


def kernel(delta, length):
    """Finite two-draw delay characteristic, continuous at all aliases."""
    reduced = np.remainder(delta * STEP / (2 * np.pi) + 0.5, 1) - 0.5
    return np.square(np.sinc(length * reduced) / np.sinc(reduced))


def effects(n):
    bits = material.bits(n)
    words = np.arange(1 << n)
    tag = sparse.diags(bits.mean(axis=1), format="csc")
    neighbor = sparse.csc_matrix(tag.shape)
    tagged = sparse.csc_matrix(tag.shape)
    for j in range(n):
        # bits[:, j] is bit j. Flipping the right neighbor leaves the tag unchanged.
        flipped = words ^ (1 << ((j + 1) % n))
        neighbor += sparse.csc_matrix(
            (np.full(len(words), 1 / n), (words, flipped)), shape=tag.shape
        )
        tagged += sparse.csc_matrix((bits[:, j] / n, (words, flipped)), shape=tag.shape)
    return tag, neighbor, tagged


def source_words(n):
    starts = np.zeros((1 << n, 3))
    starts[1, 0] = 1
    starts[[1, 2], 1] = 1 / np.sqrt(2)
    starts[0, 2] = 1
    return starts


def block_measures(n, kappa, rest):
    _, blocks = ring.spectrum(n, kappa)
    starts = source_words(n)
    operators = effects(n)
    # Indices are delay window, then DC or signed frequency, then transition
    # class, source and observed moment (tag, sign, tagged sign).
    total = np.zeros((len(WINDOWS), 1 + len(SIDES), len(BRANCHES), 3, 3))
    ground_diag = np.zeros_like(total[:, :, 0])
    band_weight = np.zeros(3)
    levels = []
    for k, (trans, e, v, _, band) in enumerate(blocks):
        c = v.conj().T @ (trans.conj().T @ starts)
        band_weight += abs(c[band]) ** 2
        levels.append(e[band])
        diff = e[:, None] - e
        # The real input, real effects and time-reversal-paired momentum
        # sectors cancel all sine terms in the total expectation at |r-s|.
        filters = np.array(
            [
                [kernel(diff, w), *[kernel(diff - (rest + off), w) for off in SIDES]]
                for w in WINDOWS
            ]
        ).reshape(-1, e.size**2)
        for oi, op in enumerate(operators):
            reduced = (trans.conj().T @ op @ trans).toarray()
            observable = v.conj().T @ (reduced @ v)
            for s in range(3):
                product = (c[:, s].conj()[:, None] * observable * c[:, s]).real
                all_terms = (filters @ product.ravel()).reshape(len(WINDOWS), -1)
                shaped = filters.reshape(len(WINDOWS), -1, e.size, e.size)
                within = shaped[:, :, band, band] * product[band, band]
                if k == 0:
                    ground_band = (
                        shaped[:, :, 0, band] * product[0, band]
                        + shaped[:, :, band, 0] * product[band, 0]
                    )
                    diag = shaped[:, :, 0, 0] * product[0, 0]
                    grounds = (
                        shaped[:, :, 0, :] @ product[0]
                        + shaped[:, :, :, 0] @ product[:, 0]
                        - diag
                    )
                    ground_diag[:, :, s, oi] = diag
                else:
                    grounds = ground_band = 0
                total[:, :, 0, s, oi] += ground_band
                total[:, :, 1, s, oi] += grounds - ground_band
                total[:, :, 2, s, oi] += within
                total[:, :, 3, s, oi] += all_terms - grounds - within
    plain = total.sum(axis=2)
    dephased = ground_diag + total[:, :, 2:].sum(axis=2)
    selected = total[:, :, 2] / band_weight[None, None, :, None]
    return np.stack([plain, dephased, selected]), total, band_weight, np.array(levels)


def joint(identity, moments):
    tag, sign, tagged = moments
    return (
        np.array(
            [
                tag + tagged,
                tag - tagged,
                identity - tag + sign - tagged,
                identity - tag - sign + tagged,
            ]
        )
        / 2
    )


def bank_counts(mean, signed, freq, length, rng):
    ordinary = joint(1, mean)
    contrast = joint(kernel(freq, length), signed)
    p = np.stack([(ordinary + contrast) / 2, (ordinary - contrast) / 2], axis=1)
    if p.min() < -1e-12 or abs(p.sum() - 1) > 1e-12:
        raise ValueError(f"Invalid retained distribution {p}")
    p = np.maximum(p.ravel(), 0)
    return rng.multinomial(SHOTS, p / p.sum()).reshape(4, 2)


def tagged_contrast(counts):
    """Original-denominator signed P(tag) times the neighboring sign."""
    a = np.asarray(counts)
    return float((a[0, 0] - a[0, 1] - a[1, 0] + a[1, 1]) / SHOTS)


def acquire(panel):
    rest = reference(panel)
    n, kappa = panel["n"], panel["kappa"]
    data, parts, weight, levels = block_measures(n, kappa, rest)
    rng = np.random.default_rng([SEED, n, int(kappa)])
    banks = []
    for ci, control in enumerate(CONTROLS):
        for s, source in enumerate(SOURCES):
            for wi, w in enumerate(WINDOWS):
                for ti, off in enumerate(SIDES, 1):
                    counts = bank_counts(
                        data[ci, wi, 0, s], data[ci, wi, ti, s], rest + off, w, rng
                    )
                    banks.append(
                        dict(
                            control=control,
                            source=source,
                            length=w,
                            offset=off,
                            counts=counts.tolist(),
                            tagged_contrast=tagged_contrast(counts),
                        )
                    )
    return dict(
        n=n,
        kappa=kappa,
        rest_hand=rest,
        band_width=float(np.ptp(levels)),
        band_weights=dict(zip(SOURCES, weight.tolist(), strict=True)),
        center_tagged_components={
            s: dict(zip(BRANCHES, parts[-1, 2, :, i, 2].tolist(), strict=True))
            for i, s in enumerate(SOURCES)
        },
        banks=banks,
    )


def run():
    prior = json.loads(Path(fringes.__file__).with_name("fringes.json").read_text())
    panels = [p for p in prior["panels"] if p["cells"] == 80]
    assert {(p["n"], p["kappa"]) for p in panels} == {
        (n, g) for n in [10, 12] for g in [4, 8]
    }
    return dict(
        cells=80,
        seed=SEED,
        shots_per_bank=SHOTS,
        step=STEP,
        outcomes=["upper,+", "upper,-", "lower,+", "lower,-"],
        reference_signs=[1, -1],
        panels=[acquire(p) for p in panels],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    output = parser.parse_args().output
    output.write_text(json.dumps(run(), indent=1, allow_nan=False) + "\n")
