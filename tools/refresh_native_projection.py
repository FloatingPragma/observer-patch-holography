#!/usr/bin/env python3
"""Refresh official-model navigation while preserving every scientific ledger row."""
from __future__ import annotations
import argparse
from copy import deepcopy
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
LEDGER='code/particles/runs/status/postdiction_ledger.json'
MARKDOWN='docs/POSTDICTION_LEDGER.md'
METADATA='executable_completion_program'


def require(condition,message):
    if not condition:raise ValueError(message)


def scientific_projection(ledger):
    return {k:v for k,v in ledger.items() if k!=METADATA}


def validate_projection(before,after):
    require(json.dumps(scientific_projection(before),sort_keys=True,separators=(',',':'),allow_nan=False)==
            json.dumps(scientific_projection(after),sort_keys=True,separators=(',',':'),allow_nan=False),
            'scientific ledger values changed during navigation refresh')
    require(after.get(METADATA,{}).get('program',{}).get('official_model')=='CCF-LT-2',
            'wrong official model projection')


def project(ledger,completion):
    result=deepcopy(ledger);result[METADATA]=deepcopy(completion)
    validate_projection(ledger,result)
    return result


def canonical_renderer(root):
    spec=importlib.util.spec_from_file_location('oph_postdiction_renderer',Path(root)/'code/particles/scripts/build_postdiction_ledger.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module._render_md


def refresh(root=ROOT,write=False):
    from tools.build_ccf_lt2_claims import generate
    from tools.build_completion_program import ledger_view
    root=Path(root).resolve()
    require(root==ROOT,'projection uses the actual local scientific producers')
    expected=generate()
    require(json.loads((root/'claims/ccf_lt2.json').read_text())==expected,'official finite claim register is stale')
    raw=(root/LEDGER).read_bytes();before=json.loads(raw)
    result=project(before,ledger_view(root))
    text=json.dumps(result,sort_keys=True,indent=2)+'\n'
    markdown=canonical_renderer(root)(result)
    require((root/LEDGER).read_bytes()==raw,'ledger changed concurrently')
    if write:
        (root/LEDGER).write_text(text);(root/MARKDOWN).write_text(markdown)
    else:
        require((root/LEDGER).read_text()==text,'official ledger metadata is stale')
        require((root/MARKDOWN).read_text()==markdown,'ledger Markdown is stale')
    validate_projection(before,json.loads((root/LEDGER).read_text()))
    return {'official_model':'CCF-LT-2','all_scientific_values_preserved':True,
            'scientific_owners_reexecuted':False,'physical_comparisons_added':0,
            'scientific_projection_sha256':sha256(json.dumps(scientific_projection(before),sort_keys=True,separators=(',',':')).encode()).hexdigest(),
            'written':write}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args();print(json.dumps(refresh(write=args.write),sort_keys=True,indent=2))


if __name__=='__main__':main()
