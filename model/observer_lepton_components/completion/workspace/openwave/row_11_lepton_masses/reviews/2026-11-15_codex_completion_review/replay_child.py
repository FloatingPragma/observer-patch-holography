#!/usr/bin/env python3
"""Isolated child: original byte replay with an audit of repository imports."""
from pathlib import Path
import argparse
import contextlib
import hashlib
import io
import json
import platform
import runpy
import sys
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    assert sys.flags.isolated and not sys.flags.optimize and sys.dont_write_bytecode
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--packet',type=Path,required=True)
    ap.add_argument('--extra',type=Path,required=True)
    args=ap.parse_args()
    root=args.root.resolve();packet=args.packet.resolve()
    sources=json.loads((packet/'SOURCES.json').read_text())
    extra=json.loads(args.extra.read_text())
    # The original manifest is not rewritten when a source has drifted.
    assert all(sha(root/p)==value for p,value in sources.items())
    for line in (packet/'PREREG.txt').read_text().splitlines():
        if '  ' in line:
            value,path=line.split('  ',1)
            assert sha(packet/path)==value, path
    output=io.StringIO()
    start=time.perf_counter()
    sys.argv=[str(packet/'headline.py'),'--check',str(packet/'RECEIPT.json')]
    with contextlib.redirect_stdout(output):
        runpy.run_path(str(packet/'headline.py'),run_name='__main__')
    imported=[]
    external={}
    for name,module in sorted(sys.modules.items()):
        path=getattr(module,'__file__',None)
        if not path:
            continue
        path=Path(path).resolve()
        if path.suffix not in ('.py','.so','.dylib'):
            continue
        try:
            relative=str(path.relative_to(root))
        except ValueError:
            if 'site-packages' in str(path):
                external[name]={'file':str(path),'sha256':sha(path)}
            continue
        digest=sha(path)
        assert (sources|extra).get(relative)==digest, f'unbound repository import: {name}: {relative}'
        imported.append({'module':name,'path':relative,'sha256':digest,
                         'binding':'original SOURCES' if relative in sources else 'review supplementary custody'})
    result={'schema':'row11.isolated_replay.v1','all_pass':True,
            'packet':packet.name,'original_receipt_sha256':sha(packet/'RECEIPT.json'),
            'original_headline_sha256':sha(packet/'headline.py'),
            'original_source_count':len(sources),'repository_imports':imported,
            'third_party_imports':external,
            'python':sys.version,'platform':platform.platform(),
            'isolated':bool(sys.flags.isolated),'optimization':sys.flags.optimize,
            'bytecode_writes_disabled':bool(sys.dont_write_bytecode),
            'headline_stdout':json.loads(output.getvalue()),
            'seconds':time.perf_counter()-start}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
