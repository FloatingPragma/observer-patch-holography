"""Replay fixed-output native owners through their unchanged flag-only CLIs.

The receipt name is an exact allowlisted path, checked against the producer's
actual OUTPUT before invocation. A caller cannot select another file or module.
Every producer executes its complete existing --check operation.
"""
from __future__ import annotations

import argparse
import importlib
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
OWNERS = {
    'native_fair_local_preparation': 'model.native_fair_local_preparation.packet',
    'native_fair_local_readout': 'model.native_fair_local_readout.packet',
    'native_fair_local_clock': 'model.native_fair_local_clock.packet',
    'native_local_charge_history': 'model.native_local_charge_history.packet',
    'native_energy_identification': 'model.native_energy_identification.packet',
    'native_gravity_identification': 'model.native_gravity_identification.certificate',
    'native_noncentral_reference': 'model.native_noncentral_reference.packet',
    'native_energy_measurement': 'model.native_energy_measurement.packet',
    'native_local_energy_clock': 'model.native_local_energy_clock.packet',
    'native_regional_modular_clock': 'model.native_regional_modular_clock.packet',
    'native_fair_ground_readout': 'model.native_fair_ground_readout.packet',
    'native_coherent_occupation_control': 'model.native_coherent_occupation_control.packet',
    'native_source_chsh': 'model.native_source_chsh.packet',
    'native_dilute_ground_limit': 'model.native_dilute_ground_limit.packet',
}


def command_for(owner: str, receipt: str) -> list[str]:
    if owner not in OWNERS:
        raise ValueError('unregistered fixed-output native owner')
    expected = f'model/{owner}/results.json'
    if type(receipt) is not str or receipt != expected:
        raise ValueError('require the exact canonical receipt of this native owner')
    path = ROOT / expected
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.resolve() != path:
        raise ValueError('fixed-output receipt may not be redirected by a symlink')
    module = importlib.import_module(OWNERS[owner])
    actual = getattr(module, 'OUTPUT', None)
    if not isinstance(actual, Path) or actual.resolve() != path:
        raise ValueError('producer OUTPUT differs from its registered receipt')
    return [sys.executable, '-m', OWNERS[owner], '--check']


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('owner', choices=tuple(OWNERS))
    parser.add_argument('--check', required=True)
    args = parser.parse_args()
    subprocess.run(command_for(args.owner, args.check), cwd=ROOT, check=True)


if __name__ == '__main__':
    main()
