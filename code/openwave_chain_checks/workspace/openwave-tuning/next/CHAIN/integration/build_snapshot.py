"""Maintainer construction of an immutable original-layout source bundle."""
from pathlib import Path
import hashlib,json,shutil
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
SNAPSHOT=HERE/'frozen_workspace'
LOGICAL='openwave-push/issues/N3D_cubic_routing/n3dlib.py'
ARCHIVED='openwave-push/issues/N3D_cubic_routing/consumed/n3dlib.py'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    admitted=json.loads((HERE/'INPUT_PINS.json').read_text())['files']
    discovered=json.loads((HERE/'IMPORT_CLOSURE.json').read_text())['files']
    files=dict(discovered)
    for rel,digest in admitted.items():
        if rel in files and files[rel]!=digest:raise ValueError('discovery disagrees with frozen input: '+rel)
        files[rel]=digest
    names=['integrate.py','strict_facade.py','gates.py','seal.py','run.py','snapshot_guard.py','build_snapshot.py',
           'inspect_imports.py','freeze_inputs.py','PROOF.md','README.md','INPUT_PINS.json','IMPORT_CLOSURE.json',
           'SOURCE_CLOCK_MAP.json','EXTERNAL_SOURCE_MAP.json','ROW1_SCOPE.json']
    for name in names:files[str((HERE/name).relative_to(ROOT))]=sha(HERE/name)
    metadata={'schema':'oph.chain.historical_source_snapshot.v1','original_live_root':str(ROOT),
        'scope':'The declared audited source versions are executed in an isolated original-relative-path workspace. This result does not validate the concurrently edited live tree.',
        'original_reviewed_input_count':len(admitted),'additional_discovered_local_files':len(set(discovered)-set(admitted)),
        'substitution':{'logical_path':LOGICAL,'archived_origin':ARCHIVED,'executed_sha256':admitted[LOGICAL],
            'live_sha256_observed_at_packaging':sha(ROOT/LOGICAL),'content_rewritten':False},
        'source_layout_preserved':True,'symlink_fallback':False,'mutable_root_fallback':False,
        'independent_audit_pins_changed':False}
    (HERE/'HISTORICAL_SOURCE.json').write_text(json.dumps(metadata,indent=2)+'\n')
    files[str((HERE/'HISTORICAL_SOURCE.json').relative_to(ROOT))]=sha(HERE/'HISTORICAL_SOURCE.json')
    if SNAPSHOT.exists():shutil.rmtree(SNAPSHOT)
    SNAPSHOT.mkdir()
    total=0
    for rel,digest in sorted(files.items()):
        source=ROOT/(ARCHIVED if rel==LOGICAL else rel)
        if sha(source)!=digest:raise ValueError('source cannot satisfy frozen snapshot: '+rel)
        destination=SNAPSHOT/rel;destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,destination)
        if destination.is_symlink() or sha(destination)!=digest:raise ValueError('snapshot copy failed: '+rel)
        total+=destination.stat().st_size
    manifest={'schema':'oph.chain.historical_workspace_manifest.v1','files':dict(sorted(files.items())),
        'file_count':len(files),'total_bytes':total,'original_reviewed_input_count':len(admitted),
        'additional_discovered_local_files':len(set(discovered)-set(admitted)),
        'scope':'File count describes source/import custody, not physical claims or scientific criteria.'}
    (HERE/'SNAPSHOT_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    inside=SNAPSHOT/HERE.relative_to(ROOT)
    shutil.copyfile(HERE/'SNAPSHOT_MANIFEST.json',inside/'SNAPSHOT_MANIFEST.json')
    print(json.dumps({'snapshot_files':len(files),'bytes':total,'reviewed_inputs':len(admitted),
        'extra_local_dependency_files':len(set(discovered)-set(admitted))},indent=2))


if __name__=='__main__':main()
