"""Validate the current finite-attachment review without promoting its scope.

This checks custody, complete canonical ownership and explicit model maps.
It is not a proof checker for the mathematical arguments it records.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
PATH = Path('model/completion/finite_attachment_map.json')
CLASSES = {'realized_on_claimed_domain','required_but_conditional',
           'alternative_branch','redundant','irrelevant_to_target','needs_correction'}
INTERFACES = {'state','algebra','operation','clock','refinement'}
REALIZED = {f'FA-{n:02d}' for n in range(1,8)} | {f'FA-{n:02d}' for n in range(13,40)}


def canonical_contract(row):
    return {key:row[key] for key in ('id','name','type','disposition','statement')}


def digest_object(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def predecessor_review(program, issue):
    """Derive dependency status from the owning completion registry."""
    stages = {stage['id']: stage for stage in program['stages']}
    stage = next(stage for stage in stages.values() if stage['issue'] == int(issue))
    required = list(stage['requires'])
    return {'required_predecessors': required,
            'unaccepted_predecessors': [gate for gate in required
                                       if stages[gate]['status'] != 'ACCEPTED']}


def validate(packet, root=ROOT):
    root=Path(root).resolve()
    register=strict_json.load(root/'tracking/premise_register.json')
    program=strict_json.load(root/'claims/completion_program.json')
    owners={row['id']:row['issue'] for row in program['stages']}
    require(packet.get('schema')=='oph.finite_attachment_map.v2','wrong attachment schema')
    require(packet.get('gate_acceptance_owner') == 'claims/completion_program.json' and
            packet.get('physical_comparison_count_added')==0 and
            packet.get('canonical_premises_discharged')==[], 'scope promotion forbidden')
    require(packet.get('selected_issues')==[1138,1139,1140,1141,1142,1143,1147,1137], 'selected issue ownership changed')
    bindings=packet.get('source_bindings')
    require(isinstance(bindings,dict) and bindings,'missing source bindings')
    for name,expected in bindings.items():
        p=(root/name).resolve()
        require(p.is_relative_to(root) and p.is_file(),f'missing bound source: {name}')
        require(hashlib.sha256(p.read_bytes()).hexdigest()==expected,f'stale source hash: {name}')
    reviewed=packet.get('reviewed_sources',{})
    for name,source in reviewed.items():
        require(source['path'] in bindings,f'unbound reviewed source: {name}')
        text=(root/source['path']).read_text()
        require(source.get('anchors') and all(a in text for a in source['anchors']),
                f'missing reviewed theorem anchor: {name}')
    rows=packet.get('canonical_ownership',[])
    require([r.get('premise') for r in rows]==[r['id'] for r in register['rows']],
            'canonical ownership must cover every current premise exactly once')
    for row,canonical in zip(rows,register['rows'],strict=True):
        require(row.get('contract_sha256')==digest_object(canonical_contract(canonical)),
                f'canonical premise contract changed: {canonical["id"]}')
        require(row.get('preserved_discharge_status')==canonical['discharge']['status'],
                f'canonical discharge was silently changed: {canonical["id"]}')
        require(row.get('classification') in CLASSES and row.get('reason'),
                f'unclassified canonical premise: {canonical["id"]}')
        gates=row.get('owners')
        require(isinstance(gates,list) and gates and len(gates)==len(set(gates)) and
                all(g in owners for g in gates),f'missing canonical owner: {canonical["id"]}')
        require(row.get('issues')==[owners[g] for g in gates],f'wrong owner issue: {canonical["id"]}')
    existing=register['model_supply']['contract_coverage']
    snapshot=packet.get('existing_navigation_coverage',{})
    require(snapshot.get('contract_coverage_sha256')==digest_object(existing),
            'existing claim/Lean navigation coverage changed; review and refresh it')
    require(snapshot.get('claim_count')==len(existing['claim_ids']) and
            snapshot.get('lean_module_count')==len(existing['lean_files']), 'navigation inventory mismatch')
    applications=packet.get('attachments',[])
    identifiers=[a.get('id') for a in applications]
    require(len(identifiers)==len(set(identifiers)) and set(identifiers)=={f'FA-{n:02d}' for n in range(1,40)},
            'missing or duplicate finite attachment')
    for a in applications:
        require(a.get('classification') in CLASSES,'unclassified finite attachment')
        require(set(a.get('interface',{}))==INTERFACES and all(a['interface'].values()),
                f'incomplete state/algebra/operation/clock/refinement map: {a["id"]}')
        require(a.get('owners') and all(g in owners for g in a['owners']),f'missing attachment owner: {a["id"]}')
        require(a.get('remaining') and a.get('witnesses') and a.get('countercontrols'),
                f'missing witness or boundary: {a["id"]}')
        require(set(a.get('reviewed_sources',[])) <= set(reviewed),f'unknown theorem source: {a["id"]}')
        for p in a['witnesses']+a['countercontrols']:
            require((root/p).is_file(),f'missing attachment witness: {p}')
    actual={a['id'] for a in applications if a['classification']=='realized_on_claimed_domain'}
    require(actual==REALIZED and packet.get('finite_obligations_realized')==len(REALIZED),
            'finite realization count or domain changed')
    review=packet.get('issue_closure_review',{})
    require(set(review)=={'1138','1139','1140','1141','1142','1143','1147','1137'},'missing issue closure review')
    for n,row in review.items():
        require(row.get('close_issue') is False, f'map cannot close its own issue: {n}')
        expected = predecessor_review(program, n)
        require(all(row.get(key) == value for key, value in expected.items()),
                f'stale prerequisite status or ownership: {n}')
        if n == '1137':
            require(row.get('required_predecessors') == ['CP-00','CP-01'] and
                    row.get('finite_map_review_complete') is True and row.get('closure_policy'),
                    'CP13 closure requires independently accepted predecessors and review')
        else:
            require(row.get('unfulfilled_acceptance'),
                    f'unsupported construction closure: {n}')
    accepted = next(s for s in program['stages'] if s['id']=='CP-13')['status']=='ACCEPTED'
    return {'canonical_premises':len(rows),'finite_obligations_realized':len(REALIZED),
            'canonical_premises_discharged':0,'classifications':dict(Counter(r['classification'] for r in rows)),
            'full_gate_accepted':accepted, 'gate_acceptance_owner':packet['gate_acceptance_owner']}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true',help='validate retained review and source custody')
    parser.parse_args()
    print(json.dumps(validate(strict_json.load(ROOT/PATH)),indent=2))


if __name__=='__main__':
    main()
