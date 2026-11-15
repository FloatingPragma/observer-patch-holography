"""Replay the unchanged independent audit in a temporary original-path layout."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import tempfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((HERE/'ADDENDUM_MANIFEST.json').read_text())
    for name, expected in manifest['files_sha256'].items():
        if sha(HERE/name) != expected:
            raise ValueError('changed addendum input: '+name)
    provenance = HERE/'provenance'
    origins = json.loads((provenance/'ORIGINS.json').read_text())['files']
    for original, record in origins.items():
        if sha(provenance/record['path']) != record['sha256']:
            raise ValueError('changed independent audit input: '+original)
    # No native or scientific function is replaced. Only the audit's explicit
    # filesystem roots are rebound to a temporary copy of its original layout.
    with tempfile.TemporaryDirectory(prefix='oph-law-independent-') as directory:
        root = Path(directory).resolve()
        rer = root/'reverse-engineering-reality'
        code = rer/'code'; code.mkdir(parents=True)
        shutil.copytree(HERE.parent/'openwave_row1_native', code/'openwave_row1_native',
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        protocol = code/'openwave_row1_protocol'; protocol.mkdir()
        for name in ('law_controls.py', 'LAW_CONTROLS.json'):
            shutil.copyfile(HERE/name, protocol/name)
        for original, record in origins.items():
            destination = root/original
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(provenance/record['path'], destination)
        original = root/'openwave-tuning/next/CHAIN2/law_gate_audit/audit.py'
        spec = importlib.util.spec_from_file_location('independent_law_auditor', original)
        auditor = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(auditor)
        auditor.main()
        actual = json.loads((original.parent/'AUDIT_RESULT.json').read_text())
        expected = json.loads((provenance/'AUDIT_RESULT.json').read_text())
        if actual != expected:
            raise ValueError('portable independent replay differs from frozen audit')
        args.output.write_text(json.dumps(actual, indent=2)+'\n')
    print('PASS: unchanged independent audit reproduced exactly')


if __name__ == '__main__':
    main()
