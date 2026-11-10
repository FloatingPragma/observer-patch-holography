"""Reject status changes, lost evidence and unlisted narrative or scope edits."""
import copy

import pytest

from .field_navigation import digest, normalize, project


def fixture():
    baseline = {'register': {'claims': [{'statement': 'retained prose', 'status': 'conditional',
                 'evidence': ['original.py'], 'assumptions': ['old']}],
                 'bindings': [{'assumptions': ['old'], 'scope': 'finite only'}]}}
    prose = {'file': 'register', 'json_path': ['claims', 0, 'statement'],
             'before_value_sha256': digest('retained prose'), 'after': 'Reviewed concise prose',
             'after_value_sha256': digest('Reviewed concise prose')}
    approved = {'prose_replacements': [prose], 'ordered_appends': [
        {'file': 'register', 'json_path': ['claims', 0, 'evidence'], 'items': ['new.py']},
        {'file': 'register', 'json_path': ['claims', 0, 'assumptions'], 'items': ['declared']}],
        'assumption_unions': [{'file': 'register', 'json_path': ['bindings', 0, 'assumptions'],
                              'before': ['old'], 'after': ['declared', 'old']}]}
    comparison = copy.deepcopy(baseline)
    comparison['register']['claims'][0]['statement'] = '<REVIEWED_PROSE_FIELD>'
    approved['normalized_predecessor_semantics_sha256'] = {'register': digest(comparison['register'])}
    return baseline, approved


def test_exact_reviewed_projection_preserves_inputs():
    before, approved = fixture()
    untouched = copy.deepcopy(before)
    current = project(before, approved)
    normalized = normalize(current, approved)
    assert before == untouched
    assert normalized['register']['claims'][0]['evidence'] == ['original.py']
    assert normalized['register']['claims'][0]['status'] == 'conditional'


@pytest.mark.parametrize('mutation', ['status', 'old_evidence', 'new_evidence', 'duplicate',
                                     'prose', 'assumptions', 'binding', 'scope', 'extra_claim'])
def test_unreviewed_change_rejected(mutation):
    before, approved = fixture()
    current = project(before, approved)
    claim = current['register']['claims'][0]
    if mutation == 'status':
        claim['status'] = 'proved'
    elif mutation == 'old_evidence':
        claim['evidence'].remove('original.py')
    elif mutation == 'new_evidence':
        claim['evidence'].append('unreviewed.py')
    elif mutation == 'duplicate':
        claim['evidence'].insert(0, 'new.py')
    elif mutation == 'prose':
        claim['statement'] += ' A physical mass is derived.'
    elif mutation == 'assumptions':
        claim['assumptions'].remove('old')
    elif mutation == 'binding':
        current['register']['bindings'][0]['assumptions'].remove('declared')
    elif mutation == 'scope':
        current['register']['bindings'][0]['scope'] = 'continuum'
    elif mutation == 'extra_claim':
        current['register']['claims'].append(copy.deepcopy(claim))
    with pytest.raises(ValueError):
        normalize(current, approved)
