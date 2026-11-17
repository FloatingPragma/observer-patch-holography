"""Export the signed-source particle-record claim ledger from archived evidence."""
import argparse
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'code/openwave_particle_records'
OUTPUT = ROOT / 'claims/openwave_particle_records.json'


def build(workspace=None):
    spec = importlib.util.spec_from_file_location('particle_records_verify', PACKAGE / 'verify.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    verification = module.verify(workspace)
    verification.pop('original_leaves_verified')
    assessment = json.loads((PACKAGE / 'ASSESSMENT.json').read_text())
    inventory = json.loads((PACKAGE / 'EVIDENCE.json').read_text())
    return dict(assessment, evidence=inventory, verification=verification,
                physical_premises_closed=[], official_OpenWave_awards_added=0,
                same_group_review_is_external_acceptance=False,
                scientific_rerun_performed_by_export=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--workspace', type=Path)
    args = parser.parse_args()
    content = json.dumps(build(args.workspace), indent=2, ensure_ascii=False) + '\n'
    if args.check:
        if not OUTPUT.is_file() or OUTPUT.read_text() != content:
            raise SystemExit('particle-record claim ledger differs; run this exporter without --check')
        print('particle-record claim ledger matches archived evidence')
    else:
        OUTPUT.write_text(content)
        print('wrote ' + str(OUTPUT.relative_to(ROOT)))
