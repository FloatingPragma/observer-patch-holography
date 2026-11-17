"""Verify complete copied scientific packets, with optional offline replay."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent

BASE_RUNS = [
    ('row_01_charge_quantization/work/2026-11-15_codex_pr723_vacuum_reservoir', 'headline_strict.py', []),
    ('row_02_electron_rest_energy/work/2026-11-15_codex_native_work_vacuum', 'headline.py', []),
    ('row_02_electron_rest_energy/work/2026-11-15_codex_native_work_vacuum', 'check.py', []),
    ('row_02_electron_rest_energy/work/2026-11-15_codex_native_work_vacuum', 'detector_check.py', []),
    ('row_02_electron_rest_energy/work/2026-11-16_codex_native_word_continuum', 'check.py', []),
    ('row_07_double_cover/work/2026-11-15_codex_subdivided_rotation', 'headline.py', []),
    ('row_07_double_cover/work/2026-11-15_codex_subdivided_rotation', 'independent_check.py', []),
    ('row_07_double_cover/work/2026-11-15_codex_subdivided_rotation', 'independent_check.py', ['FIXED_REFERENCE_RECEIPT.json']),
    ('row_07_double_cover/work/2026-11-15_codex_subdivided_rotation', 'independent_error.py', []),
    ('row_07_double_cover/work/2026-11-15_codex_subdivided_rotation', 'detector_check.py', []),
    ('row_01_charge_quantization/work/2026-11-15_codex_common_selected_source', 'headline_v2.py', []),
    ('row_01_charge_quantization/work/2026-11-16_codex_bounded_cosine_reader', 'check.py', []),
    ('row_01_charge_quantization/reviews/2026-11-16_codex_bounded_reader_review', 'audit.py', []),
    ('row_01_charge_quantization/work/2026-11-16_codex_common_refinement_family', 'headline.py', []),
    *[('row_01_charge_quantization/work/2026-11-16_codex_common_refinement_family',
       'independent.py', ['--k', str(k)]) for k in range(7)],
    ('row_01_charge_quantization/work/2026-11-16_codex_common_refinement_family', 'selector_check.py', []),
    ('row_01_charge_quantization/work/2026-11-16_codex_pair_stride_readout', 'experiment.py', []),
    ('row_01_charge_quantization/work/2026-11-16_codex_pair_stride_readout', 'headline.py', []),
    ('row_01_charge_quantization/reviews/2026-11-16_codex_pair_stride_audit', 'audit.py', []),
]

# A separately registered source/reference and probability-cost amendment is
# added here only after its producer packet and adversarial review are sealed.
AMENDMENT_RUNS = [
    ('row_01_charge_quantization/work/2026-11-16_codex_common_reference_menu', 'headline.py', []),
    ('row_01_charge_quantization/reviews/2026-11-16_codex_reference_menu_adversarial',
     'audit.py', ['../../work/2026-11-16_codex_common_reference_menu', '6', 'stationary']),
    ('row_07_double_cover/reviews/2026-11-16_codex_common_stationary_audit', 'spectral_check.py', []),
]
CURRENT_RUNS = BASE_RUNS + AMENDMENT_RUNS


def contract_module():
    import importlib.util
    # The archived descendant is immutable and reviewed separately. Its schema
    # contracts bind the historical 28-job output, not successor experiments.
    folder = HERE/'evidence/openwave/row_01_charge_quantization/work/2026-11-16_codex_verification_contract_v3/descendants'
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location('strict_selected_source_contract', folder/'verify_source.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, folder


def verify():
    module, _ = contract_module()
    checked = module.custody(HERE.resolve())
    roots = json.loads((HERE/'ROOT_INPUTS.json').read_text())
    root_paths = set()
    for entry in roots['files']:
        path = (HERE/entry['path']).resolve()
        module.require(path.is_relative_to(HERE.resolve()) and entry['path'] not in root_paths,
                       'invalid root-input dependency path')
        root_paths.add(entry['path'])
        module.require(hashlib.sha256(path.read_bytes()).hexdigest() == entry['sha256'],
                       'changed frozen root input')
        frozen = json.loads((HERE/entry['frozen_by']).read_text())
        module.require(frozen['source_inputs'][Path(entry['path']).name] == entry['sha256'],
                       'root input differs from scientific freeze')
    module.require(len(root_paths) == 3, 'incomplete frozen root inputs')
    return {**checked, 'frozen_root_inputs': len(root_paths),
            'custody_passed': True, 'official_awards_added': 0}



def replay_current(runs=None):
    """Strict historical schema comparison with unchanged evidence metadata."""
    module, folder = contract_module()
    all_contracts = json.loads((folder/'SOURCE_JOB_CONTRACTS.json').read_text())['jobs']
    indexed = {(c['directory'], c['script'], tuple(c['arguments'])): c for c in all_contracts}
    selected = CURRENT_RUNS if runs is None else runs
    contracts = [indexed[(directory, script, tuple(arguments))] for directory, script, arguments in selected]
    module.require(len(contracts) == len(selected), 'Missing replay contract')
    before = module.snapshot(HERE)
    checked = verify()
    start = time.perf_counter()
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs = list(pool.map(lambda contract: module.run_one(HERE, contract), contracts))
    unchanged = before == module.snapshot(HERE)
    return {'utc': datetime.now(timezone.utc).isoformat(),
            'harness_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'contracts_sha256': hashlib.sha256((folder/'SOURCE_JOB_CONTRACTS.json').read_bytes()).hexdigest(),
            'scope': 'Historical source checks under exact saved-output contracts. Every current verdict is checked; runtime-only exclusions are prebound. Successor experiments have separate contracts.',
            'jobs': jobs, 'seconds': time.perf_counter()-start,
            'saved_evidence_unchanged': unchanged,
            'passed': unchanged and all(job['passed'] for job in jobs),
            'custody_before': checked, 'custody_after': verify(),
            'external_G6': False, 'official_awards_added': 0}



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', action='store_true', help='Alias for strict --replay-current')
    parser.add_argument('--replay-current', action='store_true')
    parser.add_argument('--replay-amendments', action='store_true')
    args = parser.parse_args()
    result = verify()
    if args.replay or args.replay_current:
        result['current_replay'] = replay_current()
    if args.replay_amendments:
        if not AMENDMENT_RUNS:
            raise SystemExit('No sealed amendment replay is configured')
        result['amendment_replay'] = replay_current(AMENDMENT_RUNS)
    print(json.dumps(result, indent=2))
    if (args.replay or args.replay_current) and not result['current_replay']['passed']:
        raise SystemExit(1)
    if args.replay_amendments and not result['amendment_replay']['passed']:
        raise SystemExit(1)
