"""One-way paper custody: presentation binds accepted science, never conversely."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
PAPER='extra/yang_mills_gap_clay_problem.tex'
OWNERS=('model/particle_yang_mills_gap_contract/results.json',
        'model/particle_joint_reconstruction/candidate_results.json',
        'model/particle_observer_outputs/results.json',
        'model/particle_joint_reconstruction/results.json')
KEYS=('source_sha256','input_sha256','source_bindings')


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def binding_maps(value):
    if isinstance(value,dict):
        for key,child in value.items():
            if key in KEYS and isinstance(child,dict):yield child
            yield from binding_maps(child)
    elif isinstance(value,list):
        for child in value:yield from binding_maps(child)


def check(root=ROOT):
    owners=[];bindings={}
    for name in OWNERS:
        path=root/name;value=json.loads(path.read_text());count=0
        for mapping in binding_maps(value):
            assert PAPER not in mapping,('paper-to-science backedge',name)
            for relative,expected in mapping.items():
                if relative in bindings:assert bindings[relative]==expected,('conflicting source binding',relative)
                bindings[relative]=expected;count+=1
        owners.append({'receipt':name,'receipt_sha256':sha(path),'authenticated_binding_edges':count})
    # Authenticate the exact union once per invocation. Every repeated edge
    # above must agree; there is no cross-run cache or stat-only substitute.
    for relative,expected in bindings.items():
        source=(root/relative).resolve()
        assert source.is_relative_to(root.resolve())
        assert sha(source)==expected,('scientific source changed',relative)
    # The mathematical criterion remains a hard scientific dependency.
    assert 'model/particle_yang_mills_gap_contract/PROOF.md' in bindings
    for path in ('Lean/ObserverPatchHolography/RepairGapChain.lean',
                 'Lean/Screen/KogutSusskindFiberRateComparison.lean'):
        assert path in bindings
    return {'schema':'oph.one_way_paper_dependency.v1',
        'verdict':'PASS_CURRENT_PAPER_BOUND_TO_ACCEPTED_SCIENCE_WITHOUT_REVERSE_EDGE',
        'accepted_scientific_receipts':owners,
        'presentation':{'path':PAPER,'sha256':sha(root/PAPER)},
        'scientific_source_pin_count':len(bindings),
        'paper_edit_requires_this_publication_check':True,
        'paper_edit_alone_requires_expensive_science_replay':False,
        'operative_source_change_invalidates_custody':True,
        'mathematical_proofs_Lean_and_data_remain_bound':True,
        'scientific_calculation_or_physical_claim_added':False}


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    a=p.parse_args();result=check();encoded=json.dumps(result,sort_keys=True,indent=2)+'\n'
    output=HERE/'paper_results.json'
    if a.write:output.write_text(encoded)
    else:assert output.read_text()==encoded,'paper publication custody is stale'
    print(result['verdict'])


if __name__=='__main__':main()
