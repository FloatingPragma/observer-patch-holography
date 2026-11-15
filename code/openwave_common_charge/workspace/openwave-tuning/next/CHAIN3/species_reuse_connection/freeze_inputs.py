"""Copy the consumed runtime closure and traced data inputs without edits."""
from pathlib import Path
import hashlib
import json
import shutil
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
reads=set()


def hook(event,args):
    if event!='open' or not args or not isinstance(args[0],(str,bytes)):
        return
    p=Path(args[0]).resolve()
    mode=args[1] if len(args)>1 else None
    if isinstance(mode,str)and any(c in mode for c in 'wax+'):
        return
    if p.is_file()and p.is_relative_to(ROOT)and not p.is_relative_to(HERE):
        if p.suffix=='.pyc':
            import importlib.util
            try:p=Path(importlib.util.source_from_cache(str(p)))
            except ValueError:pass
        if p.is_file():reads.add(p)


sys.addaudithook(hook)
import source
from model.common_source.chiral_candidate import kinetic
data=json.loads((HERE/'RUNTIME_INPUTS.json').read_text())['files_sha256']
reads.update(ROOT/p for p in data)
manifest=json.loads((source.FIELD/'MANIFEST.json').read_text())
reads.add((source.FIELD/'MANIFEST.json').resolve())
reads.update((source.FIELD/p).resolve()for p in manifest['files_sha256'])
reads.update((source.FIELD/'inputs'/p).resolve()for p in manifest['producer_files'])
for rel in ('openwave-push/blockers/X2b_mass_spectrum/PREMISES_P1_P9.md',
            'openwave-push/blockers/X2b_mass_spectrum/VERDICT.md'):
    reads.add(ROOT/rel)
pins={}
for path in sorted(reads):
    relative=str(path.relative_to(ROOT))
    value=hashlib.sha256(path.read_bytes()).hexdigest()
    if relative in data:assert data[relative]==value,relative
    target=HERE/'inputs'/relative
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(path,target)
    assert hashlib.sha256(target.read_bytes()).hexdigest()==value
    pins[relative]=value
aliases={}
for relative,value in manifest['producer_files'].items():
    source_file=source.FIELD/'inputs'/relative
    assert hashlib.sha256(source_file.read_bytes()).hexdigest()==value
    target=HERE/'inputs'/relative
    if target.exists():assert hashlib.sha256(target.read_bytes()).hexdigest()==value,relative
    else:
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source_file,target)
    pins[relative]=value
    aliases[relative]=str(source_file.relative_to(ROOT))
(HERE/'INPUTS.json').write_text(json.dumps({'schema':'oph.species_connection.inputs.v1',
     'files_sha256':pins,'traced_runtime_and_data_files':len(pins),
     'same_byte_top_level_provider_aliases':aliases},indent=2)+'\n')
print(json.dumps({'frozen_files':len(pins),'bytes':sum((HERE/'inputs'/p).stat().st_size for p in pins)}))
