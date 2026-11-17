#!/usr/bin/env python3
"""Validate and render the official CCF-LT-2 acceptance programme."""
from __future__ import annotations
import argparse
import hashlib
import ast
import importlib.util
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PROGRAM='claims/completion_program.json'
HYPOTHESIS='claims/executable_model_hypothesis.json'
OUTPUT='docs/COMPLETION_PROGRAM.md'
JOINT_DEPENDENCY_PATH='model/particle_joint_reconstruction/dependencies.json'


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def local_file(root: Path, name: str) -> Path:
    require(isinstance(name, str) and bool(name), "missing local artifact")
    path = (root / name).resolve()
    require(path.is_relative_to(root.resolve()) and path.is_file(),
            f"missing or escaping artifact: {name}")
    return path


def load_program(root=ROOT,*,allow_unlinked=False):
    root=Path(root).resolve()
    program=json.loads(local_file(root,PROGRAM).read_text())
    hypothesis=json.loads(local_file(root,HYPOTHESIS).read_text())
    graph=json.loads(local_file(root,JOINT_DEPENDENCY_PATH).read_text())
    require(program.get('schema')=='oph.completion_program.v2','wrong program schema')
    require(program.get('program_id')==program.get('official_model')==graph.get('primary_program')=='CCF-LT-2','wrong official model')
    require(program.get('specification')=='model/SPEC.md' and program.get('dependency_map')==JOINT_DEPENDENCY_PATH,'wrong official source contract')
    require(program.get('finite_evidence')=='claims/ccf_lt2.json','wrong finite evidence')
    for path in ('model/SPEC.md',program['finite_evidence']):local_file(root,path)
    require(program.get('hypothesis')==HYPOTHESIS and hypothesis.get('official_model')=='CCF-LT-2','wrong hypothesis parent')
    require(hypothesis.get('claim_class')=='working_hypothesis' and hypothesis.get('status')=='PROPOSED_NOT_ESTABLISHED','physical hypothesis promotion')
    for flag in ['global_positive_mass_gap_claimed','modular_flow_automatically_physical_time','einstein_premises_realized']:
        require(hypothesis.get(flag) is False,'unsupported physical promotion: '+flag)
    require(program.get('decision')=='NOT_DECIDED' and program.get('whole_program_falsified') is False and program.get('physical_comparison_count_added')==0,'unsupported programme verdict')
    require(program.get('boundary_status')=='SPECIFIED_FINITE_LAW','finite law specification changed')
    require(not any(k in program for k in ['predecessor_program','legacy_issue_policy','superseded_issues','adopted_class']),'history or different law in official programme')
    stages=program.get('stages',[])
    require([s['id'] for s in stages]==[n['id'] for n in graph['nodes']],'acceptance inventory differs')
    completed=set()
    for index, (stage,node) in enumerate(zip(stages,graph['nodes'],strict=True)):
        require(stage['requires']==node['depends_on'] and set(stage['requires'])<=completed,'cyclic or substituted dependency')
        expected = 'finite_constructed' if index < 3 else 'conditional_criterion' if index == 8 else 'physical_proof_required'
        require(stage['status']==node['status']==expected,'status promotion or changed finite domain')
        require(stage['domain']==node['domain'],'source domain changed')
        require(stage['acceptance']==node['acceptance'],'acceptance condition omitted')
        require(stage['evidence']==[e['path'] for e in node['evidence']],'evidence substitution')
        for path in stage['evidence']:local_file(root,path)
        completed.add(stage['id'])
    return program,hypothesis


def issue_link(number,label=None):
    return f'[{label or str(number)}](https://github.com/FloatingPragma/observer-patch-holography/issues/{number})' if number else 'Local acceptance condition'


def readiness(stage,program):
    return 'Checked in the stated finite domain' if stage['status']=='finite_constructed' else 'Physical source attachment required' if stage['status']=='physical_proof_required' else 'Conditional implication'


def render(program,hypothesis):
    lines=['# Common-source physical acceptance','',
      'The official executable model is [CCF-LT-2](../model/SPEC.md). Its',
      '[finite evidence](../claims/ccf_lt2.json) establishes the complete native',
      'operator and named diagnostic domains. Physical correspondence uses the',
      'same Hamiltonian, preparation, retained reference, clock and instrument.',
      '',program['acceptance_policy'],'',
      'The [machine dependency map](../'+JOINT_DEPENDENCY_PATH+') fixes every',
      'acceptance condition. Specification of a construction is separate from',
      'proof that its records reproduce nature. No physical comparison or',
      'Einstein premise is awarded by this programme.','']
    for stage in program['stages']:
      lines += ['## '+stage['title'],'',readiness(stage,program)+'.','',stage['domain'],'']
      lines += ['- '+item for item in stage['acceptance']]
      lines += ['','Evidence: '+', '.join('['+Path(path).name+'](../'+path+')' for path in stage['evidence'])+'.','']
    lines += ['The [Einstein bridge](../model/ccf_lt2/EINSTEIN_BRIDGE.md) reuses the',
      'established matched tensor and Newton implications. Its regional data',
      'and common physical refinement must be supplied by this source.','']
    return '\n'.join(lines)


def render_issue(stage,program):
    return '\n'.join(['## Acceptance condition','',stage['title'],'',stage['domain'],'',
        '## Dependencies','',', '.join(stage['requires']) or 'None.','',
        '## Acceptance','',*['- '+x for x in stage['acceptance']],'',
        'Finite checks retain their stated domains; physical acceptance requires same-source evidence.',''])


def ledger_view(root=ROOT):
    program,hypothesis=load_program(root)
    parents={p:sha256(local_file(Path(root),p).read_bytes()).hexdigest() for p in [PROGRAM,HYPOTHESIS,JOINT_DEPENDENCY_PATH,'model/SPEC.md']}
    return {'verification':'Source contract and dependency structure only; not scientific gate acceptance or a physical comparison','parents':parents,'hypothesis':hypothesis,'program':program,'retained_foundations':[]}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');parser.add_argument('--allow-unlinked',action='store_true');args=parser.parse_args()
    program,hypothesis=load_program();content=render(program,hypothesis)
    if args.check:require((ROOT/OUTPUT).read_text()==content,'completion programme view is stale')
    else:(ROOT/OUTPUT).write_text(content)
    print('PASS_CCF_LT2_COMPLETION_DOMAINS')


if __name__=='__main__':main()


def integration_source_paths(root: Path, contract: dict) -> set[str]:
    """Repository Python closure, including function imports and eager parents."""
    document = local_file(root, contract.get('domain_and_proof'))
    paths = {document.relative_to(root).as_posix(), 'model/requirements.txt'}
    paths.update(p.relative_to(root).as_posix() for p in document.parent.glob('*.py'))
    paths.update(contract.get('independent_checks', []))
    receipts = contract.get('source_receipts', [])
    require(isinstance(receipts, list) and bool(receipts),
            'integration source custody requires owning source receipts')
    for name in receipts:
        receipt = json.loads(local_file(root, name).read_text())
        inherited = receipt.get('source_bindings', receipt.get('source_sha256'))
        require(isinstance(inherited, dict) and bool(inherited),
                f'integration source receipt has no source manifest: {name}')
        paths.add(name)
        paths.update(inherited)
    pending = list(paths)
    visited = set()

    def add_module(name):
        if not name.startswith(('model.', 'tools.')):
            return
        stem = name.replace('.', '/')
        candidates = [stem+'.py', stem+'/__init__.py']
        parts = stem.split('/')
        candidates.extend('/'.join(parts[:end])+'/__init__.py'
                          for end in range(1, len(parts)))
        for candidate in candidates:
            if (root/candidate).is_file() and candidate not in paths:
                paths.add(candidate)
                pending.append(candidate)

    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        module = name.removesuffix('.py').replace('/', '.')
        package = (module.removesuffix('.__init__') if module.endswith('.__init__')
                   else module.rpartition('.')[0])
        add_module(module)
        for node in ast.walk(ast.parse(local_file(root, name).read_text())):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add_module(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add_module(target)
                for alias in node.names:
                    add_module(target+'.'+alias.name)
    return paths



def validate_integration_sources(root: Path, contract: dict) -> None:
    paths = integration_source_paths(root, contract)
    bindings = contract.get('source_bindings')
    require(isinstance(bindings, dict) and paths <= set(bindings),
            'integration source custody omits a transitive dependency')
    for name, expected in bindings.items():
        require(hashlib.sha256(local_file(root, name).read_bytes()).hexdigest() == expected,
                f'integration source custody changed: {name}')
    for name in contract['source_receipts']:
        receipt = json.loads(local_file(root, name).read_text())
        inherited = receipt.get('source_bindings', receipt.get('source_sha256'))
        require(all(bindings.get(path) == expected for path, expected in inherited.items()),
                f'integration source custody differs from owning receipt: {name}')
    require(contract.get('integrated_gate_accepted') is False
            and contract.get('source_law_changed') is False,
            'integration simplification cannot accept a gate or change its source law')
