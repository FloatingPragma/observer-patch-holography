#!/usr/bin/env python3
"""Validate the bounded same-source acceptance tasks without promoting physics."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
QUEUE='claims/completion_work_queue.json'
PROGRAM='claims/completion_program.json'
BEGIN='<!-- completion-work-queue:start -->'
END='<!-- completion-work-queue:end -->'


def require(condition,message):
    if not condition:raise ValueError(message)


def validate_queue(queue,root=ROOT,program=None):
    root=Path(root).resolve()
    if program is None:program=json.loads((root/PROGRAM).read_text())
    require(queue.get('schema')=='oph.completion_work_queue.v2' and queue.get('official_model')=='CCF-LT-2','wrong queue model')
    require(queue.get('gate_acceptance_changed') is False and queue.get('physical_comparison_count_added')==0,'planning cannot accept physical gates')
    stages={s['id']:s for s in program['stages']};tasks=queue.get('tasks',[])
    require(bool(tasks),'missing tasks')
    ids=[t['id'] for t in tasks];require(len(set(ids))==len(ids),'duplicate task')
    completed=set()
    for task in tasks:
        require(task['gate'] in stages,'unknown acceptance gate')
        stage=stages[task['gate']]
        require(task['status']=='required','unsupported completion status')
        require(set(task['depends_on'])<=completed,'cyclic or unknown task dependency')
        require(task['acceptance']==stage['acceptance'] and task['scope_boundary']==stage['domain'],'acceptance or domain substitution')
        require(task['source_paths']==stage['evidence'],'source substitution')
        for name in task['source_paths']:
            target=(root/name).resolve()
            require(target.is_relative_to(root) and target.is_file(),'missing or escaping source')
        completed.add(task['id'])
    return queue


def load_queue(root=ROOT,program=None):
    queue=json.loads((Path(root)/QUEUE).read_text());validate_queue(queue,root,program);return queue


def render_section(queue,gate=None):
    tasks=[t for t in queue['tasks'] if gate is None or t['gate']==gate]
    if not tasks:return ''
    lines=[BEGIN,'','## Same-source calculations','']
    for task in tasks:
        lines += ['### '+task['title'],'',task['scope_boundary'],'']+['- '+s for s in task['acceptance']]+['']
    return '\n'.join(lines+[END,''])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');parser.parse_args()
    queue=load_queue();print('PASS_CCF_LT2_TASK_DOMAINS: '+str(len(queue['tasks']))+' required calculations')


if __name__=='__main__':main()
