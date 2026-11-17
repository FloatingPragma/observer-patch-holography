#!/usr/bin/env python3
"""Copy the frozen dependency closure. Never modifies any provider file."""
from pathlib import Path
import hashlib,json,shutil
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
DEST=HERE/'publication'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    if DEST.exists():raise RuntimeError('refuse to overwrite publication artifact')
    pins=json.loads((HERE/'PREREG.json').read_text())['files']
    for path,h in pins.items():assert sha(ROOT/path)==h,path
    packet=[p for p in HERE.rglob('*') if p.is_file() and p.name not in ('MANIFEST.sha256','EXPORT_RECEIPT.json') and '__pycache__' not in p.parts]
    DEST.mkdir()
    inputs={ROOT/path for path in pins}|set(packet)
    for source in sorted(inputs):
        target=DEST/source.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        assert sha(source)==sha(target)
    files=sorted(p for p in DEST.rglob('*') if p.is_file())
    (DEST/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.relative_to(DEST)}\n' for p in files))
    report={'schema':'oph.static_source_export.v1','source_inputs':len(pins),'copied_files':len(files),
        'copied_bytes':sum(p.stat().st_size for p in files),'manifest_sha256':sha(DEST/'MANIFEST.sha256'),
        'live_workspace_imports_required':False,'original_provider_files_modified':False}
    (HERE/'EXPORT_RECEIPT.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps(report,sort_keys=True))
if __name__=='__main__':main()
