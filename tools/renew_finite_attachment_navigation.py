#!/usr/bin/env python3
"""Regenerate current navigation from explicit scientific domains and claims."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'tools'))
from tools import build_model_theorem_map as theorem_map
from tools import build_completion_program as programme


def require(condition,message):
    if not condition:raise ValueError(message)


def build(root=ROOT):
    root=Path(root).resolve()
    premises=json.loads((root/theorem_map.REGISTER).read_text())
    claims=json.loads((root/theorem_map.CLAIMS).read_text())
    mapping=premises[theorem_map.SECTION]
    errors=theorem_map.validate(mapping,premises,claims,root)
    require(not errors,'invalid theorem navigation: '+repr(errors))
    counts=theorem_map.inventory(mapping,premises,claims,root)
    program,hypothesis=programme.load_program(root)
    return {
        str(theorem_map.OUTPUT):theorem_map.render(mapping,counts),
        programme.OUTPUT:programme.render(program,hypothesis),
    }


def execute(root=ROOT,write=False):
    root=Path(root).resolve();outputs=build(root)
    changed=[p for p,value in outputs.items() if not (root/p).exists() or (root/p).read_text()!=value]
    if write:
        for path,value in outputs.items():(root/path).write_text(value)
    return {'official_model':'CCF-LT-2','changed_paths':changed,'written':write,
            'gate_statuses_changed':False,'physical_premises_discharged':[]}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,default=ROOT)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true');mode.add_argument('--preview',action='store_true')
    args=parser.parse_args();result=execute(args.root,args.write)
    if args.check:require(not result['changed_paths'],'stale official navigation')
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':main()
