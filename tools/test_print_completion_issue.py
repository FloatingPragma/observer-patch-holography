"""Printing acceptance contracts does not publish or award their conclusions."""
from tools import build_completion_program as program
from tools import print_completion_issue as printer


def test_printed_contract_has_exact_domain_and_acceptance():
    p,_=program.load_program()
    for stage in p['stages']:
        body=printer.issue_body(stage,p,include_work_queue=False)
        assert stage['domain'] in body
        assert all(item in body for item in stage['acceptance'])
        assert body==program.render_issue(stage,p)


def test_scoped_tasks_are_attached_only_to_owner():
    p,_=program.load_program()
    stage=next(s for s in p['stages'] if s['id']=='geometric_gravity')
    body=printer.issue_body(stage,p)
    assert 'Same-source calculations' in body
    assert 'OPH.EinsteinBranch.Quantitative.small_ball_tensor_bound' in body
