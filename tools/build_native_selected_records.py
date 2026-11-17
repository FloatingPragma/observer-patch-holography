"""Project finite selection claims and optionally archive exact sealed packets.

--check requires no external checkout. --workspace and --add-packet explicitly
copy only manifest-listed files; an existing archived byte is never replaced.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT/'code/native_selected_records'
OUTPUT = ROOT/'claims/native_selected_records.json'
IDS = {'OPH-NATIVE-DYNAMIC-UNIT-SELECTION', 'OPH-NATIVE-LOCALIZED-LOWEST-MODE',
       'OPH-NATIVE-TRANSPORTED-CLOCK-LIMIT', 'OPH-NATIVE-PLANAR-ROTATION-HOLONOMY',
       'OPH-NATIVE-SELECTED-STATE-JOIN'}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def archive(workspace, name):
    relative = Path(name)
    if relative.is_absolute() or '..' in relative.parts or relative.parts[0] != 'openwave':
        raise ValueError('Use a manifest-bearing openwave/ relative path')
    source = workspace/relative
    target = PACKAGE/'evidence'/relative
    manifest = source/'MANIFEST.sha256'
    entries = manifest.read_text().splitlines()
    for line in entries:
        expected, leaf = line.split('  ', 1)
        original = (source/leaf).resolve()
        if not original.is_relative_to(source.resolve()) or digest(original) != expected:
            raise ValueError('Unbound source leaf: '+leaf)
        copied = target/leaf
        copied.parent.mkdir(parents=True, exist_ok=True)
        if copied.exists():
            if copied.read_bytes() != original.read_bytes():
                raise ValueError('Refusing to replace archived bytes: '+str(copied))
        else:
            shutil.copyfile(original,copied)
    copied_manifest = target/'MANIFEST.sha256'
    if copied_manifest.exists() and copied_manifest.read_bytes() != manifest.read_bytes():
        raise ValueError('Refusing changed source manifest')
    copied_manifest.write_bytes(manifest.read_bytes())
    return {'source':str(relative),'snapshot':'evidence/'+str(relative),
            'manifest_sha256':digest(manifest),'files':len(entries)}


def build():
    spec = importlib.util.spec_from_file_location('selected_record_custody',PACKAGE/'verify.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    checked = module.verify()
    registry = json.loads((ROOT/'claims/claim_registry.yaml').read_text())
    claims = [row for row in registry['claims'] if row['claim_id'] in IDS]
    assert {row['claim_id'] for row in claims} == IDS
    return {'schema':1,'claims':claims,
            'source_inventory':'code/native_selected_records/EVIDENCE.json',
            'custody':checked,'export_performs_scientific_run':False,
            'official_awards_added':0}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--workspace',type=Path)
    parser.add_argument('--add-packet',action='append',default=[])
    args = parser.parse_args()
    if args.add_packet:
        if args.check or args.workspace is None:
            parser.error('Archiving requires --workspace and cannot use --check')
        inventory_path = PACKAGE/'EVIDENCE.json'
        inventory = json.loads(inventory_path.read_text())
        providers = {row['source']:row for row in inventory['providers']}
        for relative in args.add_packet:
            row = archive(args.workspace.resolve(),relative)
            if relative in providers and providers[relative] != row:
                raise ValueError('Source provider changed: '+relative)
            providers[relative] = row
        inventory['providers'] = list(providers.values())
        inventory_path.write_text(json.dumps(inventory,indent=2)+'\n')
    content = json.dumps(build(),indent=2,ensure_ascii=False)+'\n'
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text() != content:
            raise SystemExit('Native selection claim projection differs')
        print('Native selection claims and archived custody agree')
    else:
        OUTPUT.write_text(content)
        print('Wrote claims/native_selected_records.json')
