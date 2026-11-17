#!/usr/bin/env python3
"""Print one local acceptance contract; no external issue is edited."""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import build_completion_program_current as builder
import build_completion_work_queue as work_queue


def issue_body(stage,program,*,include_work_queue=True,root=builder.ROOT):
    body=builder.render_issue(stage,program)
    if include_work_queue:
        section=work_queue.render_section(work_queue.load_queue(root,program),stage['id'])
        if section:body+='\n'+section
    return body


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('gate');parser.add_argument('--without-work-queue',action='store_true');args=parser.parse_args()
    program,_=builder.load_program();stages={s['id']:s for s in program['stages']}
    if args.gate not in stages:raise SystemExit('unknown acceptance condition: '+args.gate)
    print(issue_body(stages[args.gate],program,include_work_queue=not args.without_work_queue),end='')


if __name__=='__main__':main()
