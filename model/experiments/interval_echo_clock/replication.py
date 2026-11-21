"""The prespecified one-look replication on fresh uniform-journey data."""
import argparse
import hashlib
import json
from pathlib import Path

from . import experiment

SEED = 2026112300
PER_BLOCK = 262144


def run():
    pattern_index = 2
    pattern = experiment.PATTERNS[pattern_index]
    if pattern != (4, 4):
        raise ValueError("the frozen uniform fixture has changed")
    probe = experiment.native_probe()
    handoffs, boards = experiment.material_ensemble(pattern, SEED, probe)
    return {
        "schema": "oph.interval_echo.uniform_clock_replication.v1",
        "protocol_sha256": hashlib.sha256((experiment.HERE / "REPLICATION.md").read_bytes()).hexdigest(),
        "seed": SEED, "prepared_services": list(pattern),
        "preparations_per_analyzer_and_rule": experiment.BLOCKS * PER_BLOCK,
        "boards": boards,
        "readouts": {rule: experiment.cohort(pattern_index, ri, handoffs, probe, seed=SEED, block_size=PER_BLOCK)
                     for ri, rule in enumerate(experiment.RULES)},
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='write to a fresh file')
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open('x') as f:
            f.write(json.dumps(result, separators=(',', ':'), allow_nan=False) + '\n')
    elif result != json.loads((experiment.HERE / 'replication.json').read_text()):
        raise SystemExit('uniform replication differs from registered receipt')
    for name, arm in result['readouts'].items():
        print(name, arm['verdict'], arm['readings']['dimer_minus_registered'])
