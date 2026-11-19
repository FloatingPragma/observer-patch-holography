"""Replay both observer instruments and verify the retained counts."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ...polynomial_sv import collar
from ..carrier import require
from ..cage.instruments import draw
from .law import collision, ramsey

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def protocol():
    hashes = json.loads((HERE / 'sources.json').read_text(encoding='utf-8'))
    for name, expected in hashes.items():
        # Sources are UTF-8 text. Git may check out CRLF on Windows.
        data = (ROOT / name).read_bytes().replace(b'\r\n', b'\n')
        require(hashlib.sha256(data).hexdigest() == expected, f'source differs: {name}')
    return json.loads((HERE / 'protocol.json').read_text(encoding='utf-8'))


def experiment(spec):
    joint = collision(spec)
    full = np.array([ramsey(spec['kappa'], g) for g in spec['couplings_in_exchange_units']])
    return {'collision': joint, 'ramsey': np.stack([1 - full, full], axis=-1)}


def observe(counts, shots):
    """Read only setting-indexed integer counts of write, capture and Ramsey pointers."""
    joint = counts['collision'] / shots
    full = joint[..., 3] + joint[..., 7]
    rates = np.array([[[collar.interpret(settings)[1] for settings in signs]
                       for signs in fillings] for fillings in counts['ramsey'][..., 1] / shots])
    reference = rates[:, 0, 1] - rates[:, 0, 0]
    require(np.all(reference > 0), 'collar reference has no positive period')
    read_g = (rates[:, 1, 0] - rates[:, 0, 0]) / reference
    contrast = full - full[0]
    return {'coupling_from_collar': read_g.tolist(),
            'largest_scalar_pointer_change': np.max(abs(contrast), axis=(1, 2, 3)).tolist(),
            'largest_sharp_minus_undriven_change': np.max(abs(contrast[:, 1] - contrast[:, 2]), axis=(1, 2)).tolist()}


def assess(spec, laws, counts):
    gates = spec['gates']
    shots = spec['shots_per_setting']
    require(set(counts) == set(laws), 'wrong instrument names')
    for key, law in laws.items():
        recorded = counts[key]
        tol = gates['maximum_absolute_probability_roundoff']
        require(np.isfinite(law).all() and np.min(law) >= -tol and np.max(abs(law.sum(axis=-1) - 1)) <= tol,
                f'incomplete {key} record law')
        require(recorded.shape == law.shape and np.issubdtype(recorded.dtype, np.integer)
                and np.all((recorded >= 0) & (recorded <= shots))
                and np.all(recorded.sum(axis=-1) == shots), f'invalid {key} counts')
        require(np.all(recorded[law == 0] == 0), f'forbidden {key} outcome')
        require(np.max(abs(recorded / shots - law)) < gates['maximum_probability_sampling_error'],
                f'{key} counts lie outside the frozen sampling band')
    readings = observe(counts, shots)
    couplings = spec['couplings_in_exchange_units']
    require(np.max(abs(np.array(readings['coupling_from_collar']) - couplings)) <
            gates['maximum_error_in_ramsey_coupling'], 'collar missed the magnitude or sign')
    for key, gate in [('largest_scalar_pointer_change', 'minimum_joint_pointer_contrast_against_additive'),
                      ('largest_sharp_minus_undriven_change', 'minimum_sharp_minus_undriven_interaction_contrast')]:
        require(min(readings[key][1:]) > gates[gate], f'collision failed {gate}')
    return readings


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--record', type=Path, help='draw into an unused npz file and score every outcome')
    args = cli.parse_args()
    spec = protocol()
    laws = experiment(spec)
    if args.record:
        rng = np.random.Generator(np.random.PCG64(spec['seed']))
        counts = {key: draw(rng, law, spec['shots_per_setting']) for key, law in laws.items()}
        with args.record.open('xb') as f:
            np.savez_compressed(f, **counts)
    else:
        receipt = HERE / 'counts.npz'
        expected = (HERE / 'counts.sha256').read_text(encoding='ascii').strip()
        require(digest(receipt) == expected, 'count receipt differs')
        with np.load(receipt, allow_pickle=False) as saved:
            counts = {name: saved[name] for name in saved.files}
    print(json.dumps(assess(spec, laws, counts), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
