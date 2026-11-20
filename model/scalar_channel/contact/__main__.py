"""Replay both observer instruments and verify the retained counts."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from tools.strict_json import load as load_json

from ...contrast_memory_spherical_v1 import collar
from ..carrier import require
from ..cage.instruments import PREPARATIONS, draw
from .law import collision, ramsey

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE_FILES = {
    'model/__init__.py', 'model/requirements.txt', 'tools/strict_json.py',
    *('model/contrast_memory_spherical_v1/' + f for f in ('__init__.py', 'carrier.py', 'collar.py', 'common.py', 'law.py')),
    *('model/scalar_channel/' + f for f in ('__init__.py', 'carrier.py', 'compiler.py', 'instruments.py',
                                            'cage/__init__.py', 'cage/construction.py', 'cage/instruments.py',
                                            'contact/__init__.py', 'contact/__main__.py', 'contact/law.py',
                                            'contact/protocol.json')),
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def protocol():
    check_sources(HERE / 'sources.json', SOURCE_FILES)
    return load_json(HERE / 'protocol.json')


def check_sources(manifest, expected_files):
    hashes = load_json(manifest)
    require(type(hashes) is dict and set(hashes) == expected_files, 'incomplete or unexpected source list')
    for name, expected in hashes.items():
        # Sources are UTF-8 text. Git may check out CRLF on Windows.
        data = (ROOT / name).read_bytes().replace(b'\r\n', b'\n')
        require(hashlib.sha256(data).hexdigest() == expected, f'source differs: {name}')


def experiment(spec):
    joint = collision(spec)
    full = np.array([ramsey(spec['kappa'], g) for g in spec['couplings_in_exchange_units']])
    return {'collision': joint, 'ramsey': np.stack([1 - full, full], axis=-1)}


def draw_counts(spec, laws):
    rng = np.random.Generator(np.random.PCG64(spec['seed']))
    return {name: draw(rng, law, spec['shots_per_setting']) for name, law in laws.items()}


def read_counts(receipt, spec, laws):
    expected = receipt.with_suffix('.sha256').read_text(encoding='ascii').strip()
    require(digest(receipt) == expected, 'count receipt differs')
    with np.load(receipt, allow_pickle=False) as saved:
        require(sorted(saved.files) == sorted(laws), 'missing, extra or repeated count array')
        counts = {name: saved[name] for name in saved.files}
    for name, generated in draw_counts(spec, laws).items():
        require(np.array_equal(counts[name], generated), f'{name} differs from the seeded acquisition')
    return counts


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
    coherent, sharp, undriven = PREPARATIONS
    require(coherent[2] == undriven[2] and undriven[1] == 0 and sharp[2] != undriven[2],
            'preparation timing differs from the declared comparison')
    matched = np.max(abs(contrast[:, 0] - contrast[:, 2]), axis=(1, 2))
    # Each difference of differences uses four binary full-fuel frequencies.
    # Union over every setting before taking a maximum over waits and analyzers.
    radius = 4 * math.sqrt(math.log(2 * full.size / 1e-8) / (2 * shots))
    return {'coupling_from_collar': read_g.tolist(),
            'largest_scalar_pointer_change': np.max(abs(contrast), axis=(1, 2, 3)).tolist(),
            'largest_sharp_minus_undriven_change': np.max(abs(contrast[:, 1] - contrast[:, 2]), axis=(1, 2)).tolist(),
            'largest_equal_duration_change': matched.tolist(),
            'simultaneous_difference_error_at_1e_minus_8': radius}


def assess(spec, laws, counts):
    gates = spec['gates']
    shots = spec['shots_per_setting']
    require(set(counts) == set(laws), 'wrong instrument names')
    # Simultaneous Bernstein bounds also resolve the rare missed writes.
    log_tail = math.log(2 * sum(law.size for law in laws.values()) / 1e-8)
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
        variance = shots * np.maximum(0, law * (1-law))
        limit = np.sqrt(2 * variance * log_tail) + 2 * log_tail / 3
        require(np.all(abs(recorded - shots * law) <= limit), f'{key} counts fail the variance bound')
    readings = observe(counts, shots)
    couplings = spec['couplings_in_exchange_units']
    require(np.max(abs(np.array(readings['coupling_from_collar']) - couplings)) <
            gates['maximum_error_in_ramsey_coupling'], 'collar missed the magnitude or sign')
    for key, gate in [('largest_scalar_pointer_change', 'minimum_joint_pointer_contrast_against_additive'),
                      ('largest_sharp_minus_undriven_change', 'minimum_sharp_minus_undriven_interaction_contrast')]:
        require(min(readings[key][1:]) > gates[gate], f'collision failed {gate}')
    require(min(readings['largest_equal_duration_change'][1:]) >
            readings['simultaneous_difference_error_at_1e_minus_8'], 'equal-duration effect unresolved')
    return readings


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--record', type=Path, help='draw into an unused npz file and score every outcome')
    args = cli.parse_args()
    spec = protocol()
    laws = experiment(spec)
    if args.record:
        counts = draw_counts(spec, laws)
        with args.record.open('xb') as f:
            np.savez_compressed(f, **counts)
    else:
        counts = read_counts(HERE / 'counts.npz', spec, laws)
    print(json.dumps(assess(spec, laws, counts), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
