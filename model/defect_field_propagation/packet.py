"""Replay complete finite native-seam propagation and source readout."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from tools import strict_json
from ..defect_field_bridge.dynamics_packet import compare_arrays, compare_json
from ..defect_field_benchmark import quasiparticles
from ..common_source import calibrated_matter
from ..common_source.chiral_candidate import gaussian_source
from ..defect_protection import packet as protection_packet
from . import propagation

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def bindings():
    paths = list(HERE.glob('*.py'))+[HERE/'README.md',
        Path(quasiparticles.__file__), Path(calibrated_matter.__file__),
        Path(gaussian_source.__file__),
        Path(gaussian_source.car_model.__file__),
        ROOT/'model/common_source/calibrated_matter_parameters.json']
    values = protection_packet.source_bindings()
    values.update({p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in sorted(paths)})
    return dict(sorted(values.items()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    pins = bindings()
    report, arrays = propagation.build()
    report['source_sha256'] = pins
    if bindings() != pins:
        raise RuntimeError('source changed during propagation replay')
    if args.write:
        (HERE/'results.json').write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
        np.savez_compressed(HERE/'results.npz', **arrays)
    elif args.check:
        expected = strict_json.load(HERE/'results.json')
        compare_json(expected, report)
        with np.load(HERE/'results.npz', allow_pickle=False) as values:
            compare_arrays({key: values[key] for key in values.files}, arrays)
    print('Verified native-seam calibrated propagation, complete modes and source current readout.')


if __name__ == '__main__':
    main()
