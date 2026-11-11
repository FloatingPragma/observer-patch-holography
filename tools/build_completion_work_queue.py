#!/usr/bin/env python3
"""Validate bounded research tasks and render their completion-guide section.

Validation checks ownership, prerequisites and evidence custody. It does not
execute the listed commands, verify a scientific proof or accept a CP gate.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path, PurePosixPath
import re

if __package__:
    from tools import strict_json
else:
    import strict_json

ROOT = Path(__file__).resolve().parents[1]
QUEUE = 'claims/completion_work_queue.json'
PROGRAM = 'claims/completion_program.json'
PREMISES = 'tracking/premise_register.json'
REPO = 'https://github.com/FloatingPragma/observer-patch-holography'
BEGIN = '<!-- completion-work-queue:start -->'
END = '<!-- completion-work-queue:end -->'
TASK_ID = re.compile(r'WQ-[0-9]{2,}\Z')
OUTCOMES = {'positive', 'negative', 'inconclusive'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def nonempty(value, label):
    require(isinstance(value, str) and bool(value.strip()), f'missing {label}')


def strings(value, label, *, required=False):
    require(isinstance(value, list) and (value or not required), f'invalid {label}')
    for item in value:
        nonempty(item, label)
    require(len(value) == len(set(value)), f'duplicate {label}')


def local_path(root, name, *, existing=True):
    nonempty(name, 'repository path')
    path = PurePosixPath(name)
    require(not path.is_absolute() and '..' not in path.parts
            and str(path) == name and not any(c in name for c in '*?[]{}\\'),
            f'nonliteral or escaping repository path: {name}')
    target = (root / name).resolve()
    require(target.is_relative_to(root), f'escaping repository path: {name}')
    require(target.is_file() if existing else not target.is_dir(),
            f'expected {"existing" if existing else "planned"} file: {name}')
    return target


def validate_queue(queue, root=ROOT, program=None):
    root = Path(root).resolve()
    if program is None:
        program = strict_json.load(local_path(root, PROGRAM))
    premises = strict_json.load(local_path(root, PREMISES))
    require(isinstance(queue, dict) and queue.get('schema') == 'oph.completion_work_queue.v1',
            'wrong completion work-queue schema')
    nonempty(queue.get('queue_id'), 'queue identifier')
    nonempty(queue.get('purpose'), 'queue purpose')
    require(queue.get('gate_acceptance_changed') is False
            and type(queue.get('physical_comparison_count_added')) is int
            and queue['physical_comparison_count_added'] == 0,
            'task planning cannot accept gates or add physical comparisons')
    for key in ('execution_policy', 'premise_policy'):
        strings(queue.get(key, []), key)
    stages = {stage['id']: stage for stage in program['stages']}
    premise_ids = {row['id'] for row in premises['rows']}
    tasks = queue.get('tasks')
    require(isinstance(tasks, list) and bool(tasks), 'missing bounded tasks')
    ids = [task.get('id') for task in tasks if isinstance(task, dict)]
    require(len(ids) == len(tasks) and all(isinstance(i, str) and TASK_ID.fullmatch(i) for i in ids)
            and len(ids) == len(set(ids)), 'missing, invalid or duplicate task identifier')
    by_id = dict(zip(ids, tasks, strict=True))
    allowed = {'id', 'title', 'gate', 'issue', 'depends_on', 'status', 'blocked_by',
               'source_paths', 'reused_theorems', 'premise_ids', 'exact_prompt',
               'outputs', 'validation', 'closure', 'scope_boundary', 'completion_evidence',
               'completion_outcome', 'completion_bindings', 'accepted_dependency_outcomes'}
    for task in tasks:
        ident = task['id']
        require(set(task) <= allowed, f'{ident}: unknown task fields {sorted(set(task) - allowed)}')
        for key in ('title', 'exact_prompt', 'scope_boundary'):
            nonempty(task.get(key), f'{ident} {key}')
        gate = task.get('gate')
        require(gate in stages and type(task.get('issue')) is int
                and task['issue'] == stages[gate]['issue'], f'{ident}: wrong CP gate/issue owner')
        require(task.get('status') in {'ready', 'blocked', 'completed'}, f'{ident}: invalid task status')
        for key in ('depends_on', 'blocked_by', 'source_paths', 'premise_ids'):
            strings(task.get(key), f'{ident} {key}', required=key == 'source_paths')
        require(set(task['depends_on']) <= set(ids) and ident not in task['depends_on'],
                f'{ident}: unknown or self dependency')
        admitted = task.get('accepted_dependency_outcomes', {})
        require(isinstance(admitted, dict) and set(admitted) <= set(task['depends_on']),
                f'{ident}: outcome admission requires a direct dependency')
        for dependency, outcomes in admitted.items():
            strings(outcomes, f'{ident} accepted outcomes for {dependency}', required=True)
            require(set(outcomes) <= OUTCOMES, f'{ident}: unknown dependency outcome')
        require(set(task['premise_ids']) <= premise_ids, f'{ident}: unknown canonical premise')
        for name in task['source_paths']:
            local_path(root, name)
        theorems = task.get('reused_theorems')
        require(isinstance(theorems, list), f'{ident}: missing theorem-reuse list')
        for theorem in theorems:
            require(isinstance(theorem, dict) and set(theorem) == {'path', 'anchor', 'role'},
                    f'{ident}: untyped theorem reuse')
            source = local_path(root, theorem['path'])
            nonempty(theorem['role'], f'{ident} theorem role')
            require(isinstance(theorem['anchor'], str), f'{ident}: invalid theorem anchor')
            if theorem['anchor']:
                require(theorem['anchor'] in source.read_text(encoding='utf-8'),
                        f'{ident}: theorem anchor absent: {theorem["anchor"]}')
        outputs = task.get('outputs')
        require(isinstance(outputs, list) and bool(outputs), f'{ident}: missing bounded outputs')
        output_paths = []
        for output in outputs:
            require(isinstance(output, dict) and set(output) == {'path', 'description'},
                    f'{ident}: untyped bounded output')
            local_path(root, output['path'], existing=task['status'] == 'completed')
            nonempty(output['description'], f'{ident} output description')
            output_paths.append(output['path'])
        require(len(output_paths) == len(set(output_paths)), f'{ident}: duplicate output path')
        commands = task.get('validation')
        require(isinstance(commands, list) and bool(commands), f'{ident}: missing validation')
        for command in commands:
            require(isinstance(command, dict) and set(command) == {'command', 'expected'},
                    f'{ident}: untyped validation command')
            nonempty(command['command'], f'{ident} validation command')
            nonempty(command['expected'], f'{ident} validation expectation')
        closure = task.get('closure')
        require(isinstance(closure, dict) and set(closure) == {'pass', 'fail', 'inconclusive'},
                f'{ident}: closure requires pass, fail and inconclusive criteria')
        for outcome, criteria in closure.items():
            strings(criteria, f'{ident} {outcome} closure', required=True)
        evidence = task.get('completion_evidence', [])
        require(isinstance(evidence, list), f'{ident}: untyped completion evidence')
        require(bool(evidence) == (task['status'] == 'completed'),
                f'{ident}: completion status/evidence mismatch')
        outcome = task.get('completion_outcome')
        require((isinstance(outcome, str) and outcome in OUTCOMES)
                if task['status'] == 'completed' else outcome is None,
                f'{ident}: completed task requires a typed scientific outcome')
        bindings = task.get('completion_bindings', {})
        require(isinstance(bindings, dict) and bool(bindings) == (task['status'] == 'completed'),
                f'{ident}: completion status/source bindings mismatch')
        required_bindings = set(output_paths)
        for row in evidence:
            require(isinstance(row, dict) and set(row) == {'path', 'sha256', 'independent_check', 'scope'},
                    f'{ident}: untyped completion evidence')
            source = local_path(root, row['path'])
            checker = local_path(root, row['independent_check'])
            require(not source.samefile(checker),
                    f'{ident}: an evidence file cannot be its own independent checker')
            nonempty(row['scope'], f'{ident} completion scope')
            require(hashlib.sha256(source.read_bytes()).hexdigest() == row['sha256'],
                    f'{ident}: completion evidence changed: {row["path"]}')
            required_bindings.update((row['path'], row['independent_check']))
            require(bindings.get(row['path']) == row['sha256'],
                    f'{ident}: completion manifests disagree: {row["path"]}')
        if task['status'] == 'completed':
            require(required_bindings <= set(bindings),
                    f'{ident}: completion bindings omit an output, result or checker')
        for name, digest in bindings.items():
            source = local_path(root, name)
            require(hashlib.sha256(source.read_bytes()).hexdigest() == digest,
                    f'{ident}: completed source changed: {name}')
    visiting, visited = set(), set()

    def visit(ident):
        require(ident not in visiting, f'work-queue dependency cycle at {ident}')
        if ident in visited:
            return
        visiting.add(ident)
        for dependency in by_id[ident]['depends_on']:
            visit(dependency)
        visiting.remove(ident)
        visited.add(ident)

    for task in tasks:
        visit(task['id'])
        incomplete = [d for d in task['depends_on'] if by_id[d]['status'] != 'completed'
                      or by_id[d].get('completion_outcome') not in
                      task.get('accepted_dependency_outcomes', {}).get(d, ['positive'])]
        blocked = bool(incomplete or task['blocked_by'])
        require(blocked == (task['status'] == 'blocked'),
                f'{task["id"]}: readiness disagrees with prerequisites or explicit blockers')
    order = queue.get('priority_order', ids)
    strings(order, 'priority_order', required=True)
    require(set(order) == set(ids), 'priority_order must contain every task exactly once')
    first = queue.get('first_batch', [])
    strings(first, 'first_batch')
    require(set(first) <= set(ids) and all(by_id[i]['status'] == 'ready' for i in first),
            'first_batch requires ready task identifiers')
    return queue


def load_queue(root=ROOT, program=None):
    root = Path(root).resolve()
    return validate_queue(strict_json.load(local_path(root, QUEUE)), root, program)


def code_block(value, language=''):
    """Retain copyable prompts even when they contain fenced examples."""
    runs = [len(m.group()) for m in re.finditer(r'`+', value)]
    fence = '`' * max(3, max(runs, default=0) + 1)
    return f'{fence}{language}\n{value}\n{fence}'


def render_section(queue, gate=None):
    """Render a validated queue, optionally restricted to one CP owner."""
    by_id = {task['id']: task for task in queue['tasks']}
    order = queue.get('priority_order', list(by_id))
    tasks = [by_id[i] for i in order if gate is None or by_id[i]['gate'] == gate]
    if not tasks:
        return ''
    def source_link(path):
        return f'[{path}]({REPO}/blob/main/{path})'
    def task_link(ident):
        target = f'#work-{ident.lower()}'
        if gate is not None and by_id[ident]['gate'] != gate:
            target = f'{REPO}/blob/main/docs/COMPLETION_PROGRAM.md{target}'
        return f'[{ident}]({target})'
    def cell(value):
        return value.replace('|', '\\|').replace('\n', ' ')
    lines = [BEGIN, '## Bounded completion tasks', '', queue['purpose'], '',
             f'Source: {source_link(QUEUE)}. Task completion records a bounded result; '
             'it does not accept a physical completion gate, discharge an unproved '
             'premise or add an independent physical comparison. Passing software tests '
             'alone does not establish a positive scientific outcome.', '']
    if gate is None:
        for key, title in [('execution_policy', 'Execution policy'), ('premise_policy', 'Premise policy')]:
            if queue.get(key):
                lines.extend([f'### {title}', ''])
                for paragraph in queue[key]:
                    lines.extend([paragraph, ''])
    else:
        lines += [f'The [execution and premise policies]({REPO}/blob/main/docs/'
                  'COMPLETION_PROGRAM.md#bounded-completion-tasks) apply to every task.', '']
    first = [i for i in queue.get('first_batch', []) if gate is None or by_id[i]['gate'] == gate]
    if first:
        lines.extend(['First batch: ' + ', '.join(task_link(i) for i in first) + '.', ''])
    lines += ['| Task | Bounded result | Owner | Status | Prerequisites |',
              '| --- | --- | --- | --- | --- |']
    for task in tasks:
        dependencies = []
        for ident in task['depends_on']:
            accepted = task.get('accepted_dependency_outcomes', {}).get(ident, ['positive'])
            label = task_link(ident)
            if accepted != ['positive']:
                label += ' (' + '/'.join(accepted) + ')'
            dependencies.append(label)
        deps = ', '.join(dependencies) or 'none'
        if task['blocked_by']:
            deps += '; ' + '; '.join(cell(s) for s in task['blocked_by'])
        status = task['status']
        if status == 'completed':
            status += ': ' + task['completion_outcome']
        lines.append(f'| {task_link(task["id"])} | {cell(task["title"])} | '
                     f'{task["gate"]} / [#{task["issue"]}]({REPO}/issues/{task["issue"]}) | '
                     f'{status} | {deps} |')
    lines.append('')
    for task in tasks:
        status = task['status']
        if status == 'completed':
            status += ': ' + task['completion_outcome']
        lines += [f'<a id="work-{task["id"].lower()}"></a>',
                  f'### {task["id"]}: {task["title"]}', '',
                  f'Owner: {task["gate"]} / [#{task["issue"]}]({REPO}/issues/{task["issue"]}). '
                  f'Status: `{status}`.', '', task['scope_boundary'], '',
                  '**Exact prompt**', '', code_block(task['exact_prompt'], 'text'), '',
                  '**Inputs and theorem reuse**', '']
        lines += [f'- {source_link(path)}' for path in task['source_paths']]
        for theorem in task['reused_theorems']:
            anchor = f' Anchor: `{theorem["anchor"]}`.' if theorem['anchor'] else ' Whole proof.'
            lines.append(f'- {source_link(theorem["path"])}:{anchor} {theorem["role"]}')
        lines += ['', 'Canonical premise references: ' +
                  (', '.join(task['premise_ids']) or 'none assigned') + '.', '', '**Bounded outputs**', '']
        lines += [f'- `{row["path"]}`: {row["description"]}' for row in task['outputs']]
        lines += ['', '**Validation**', '']
        for row in task['validation']:
            lines += [code_block(row['command'], 'sh'), '', row['expected'], '']
        lines += ['**Task closure**', '']
        for outcome in ('pass', 'fail', 'inconclusive'):
            lines += [f'- {outcome.capitalize()}: {criterion}' for criterion in task['closure'][outcome]]
        if task.get('completion_evidence'):
            lines += ['', '**Bounded completion evidence**', '']
            lines += [f'- {source_link(row["path"])}: {row["scope"]} '
                      f'Independent check: {source_link(row["independent_check"])}.'
                      for row in task['completion_evidence']]
        lines.append('')
    return '\n'.join([*lines, END, ''])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='validate without rendering')
    parser.add_argument('--gate', help='render only tasks owned by one CP gate')
    args = parser.parse_args()
    queue = load_queue()
    if args.gate:
        require(args.gate in {t['gate'] for t in queue['tasks']}, 'gate has no queued tasks')
    if args.check:
        print(f'Completion work queue verified: {len(queue["tasks"])} bounded tasks; no gate accepted')
    else:
        print(render_section(queue, args.gate), end='')


if __name__ == '__main__':
    main()
