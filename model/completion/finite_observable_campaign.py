"""Check source custody and the domains of three finite observable theorems.

The mathematical checks live in the owning packages. This validator binds
their conclusions to the finite premise map and five completion issues.
"""
import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
PATH = Path('model/completion/finite_observable_campaign.json')
CAMPAIGN = 'OPH-FINITE-OBSERVABLE-CLOSURE-2026-11-05'
COMPONENTS = {
    'charged_composite_search': ('FA-150', ['CP-03', 'CP-07']),
    'quadratic_conservation': ('FA-151', ['CP-04', 'CP-06']),
    'prepared_regional_entropy': ('FA-152', ['CP-09']),
}
ISSUES = [1139, 1140, 1142, 1143, 1145]
OWNERS = dict(zip(['CP-03', 'CP-04', 'CP-06', 'CP-07', 'CP-09'], ISSUES))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate(packet, root=ROOT):
    root = Path(root).resolve()
    require(type(packet) is dict and set(packet) == {
        'schema', 'campaign_id', 'selected_issues', 'scope', 'whole_program_frozen',
        'integrated_gate_accepted', 'canonical_premises_discharged',
        'physical_comparison_count_added', 'specification_readout_allowed',
        'component_results', 'source_bindings'}, 'wrong campaign fields')
    require(packet['schema'] == 'oph.finite_observable_campaign.v1'
            and packet['campaign_id'] == CAMPAIGN, 'wrong campaign identity')
    require(packet['selected_issues'] == ISSUES, 'wrong five-issue ownership')
    require(packet['scope'] == 'FINITE_COMPONENT_FAMILY'
            and packet['whole_program_frozen'] is False
            and packet['integrated_gate_accepted'] is False
            and packet['canonical_premises_discharged'] == []
            and type(packet['physical_comparison_count_added']) is int
            and packet['physical_comparison_count_added'] == 0
            and packet['specification_readout_allowed'] is True,
            'unsupported physical promotion or specification restriction')
    rows = packet['component_results']
    require(type(rows) is list and all(type(r) is dict for r in rows)
            and [r.get('id') for r in rows] == list(COMPONENTS),
            'missing or duplicated component')
    bindings = packet['source_bindings']
    require(type(bindings) is dict and bindings, 'missing source custody')

    def check_file(name, expected):
        require(type(name) is str and type(expected) is str, 'invalid source binding')
        path = (root/name).resolve()
        require(path.is_relative_to(root) and path.is_file()
                and not name.startswith('docs/'), 'missing or external source')
        require(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
                'stale source: '+name)

    for name, digest in bindings.items():
        check_file(name, digest)
    finite = strict_json.load(root/'model/completion/finite_attachment_map.json')
    program = strict_json.load(root/'claims/completion_program.json')

    def unique(items, identifier):
        matches = [item for item in items if item.get('id') == identifier]
        require(len(matches) == 1, 'missing or duplicate ledger row: '+identifier)
        return matches[0]

    required_sources = set()
    receipts = {f'model/{name}/results.json' for name in COMPONENTS}
    for row in rows:
        require(set(row) == {'id', 'attachment', 'issues', 'result', 'interface', 'remaining'},
                'unexpected component fields')
        name = row['id']
        attachment_id, owners = COMPONENTS[name]
        source = unique(finite['attachments'], attachment_id)
        require(source['owners'] == owners and source['reviewed_sources'] == [name]
                and source['classification'] == 'realized_on_claimed_domain'
                and source['physical_identification'] is False,
                'finite attachment ownership or scope changed')
        require(row['attachment'] == attachment_id
                and row['issues'] == [OWNERS[owner] for owner in owners], 'wrong component owner')
        for target, key in [('result', 'conclusion'), ('interface', 'interface'),
                            ('remaining', 'remaining')]:
            require(row[target] == source[key], 'component domain differs from finite theorem map')
        evidence = set(source['witnesses'] + source['countercontrols'])
        required_sources.update(evidence)
        require(evidence <= set(bindings) and all(
            bindings[p] == finite['source_bindings'].get(p) for p in evidence),
            'missing or conflicting reviewed evidence')
        receipt_path = f'model/{name}/results.json'
        require(receipt_path in evidence, 'missing component receipt')
        receipt = strict_json.load(root/receipt_path)
        dependencies = receipt.get('source_bindings', receipt.get('source_sha256'))
        require(type(dependencies) is dict and dependencies, 'missing receipt dependencies')
        for dependency, digest in dependencies.items():
            check_file(dependency, digest)
        for owner in owners:
            stage = unique(program['stages'], owner)
            require(stage['issue'] == OWNERS[owner], 'stage issue changed')
            progress = unique([r for r in stage.get('progress_evidence', [])
                               if r.get('campaign_id') == CAMPAIGN], name)
            require(progress['path'] == receipt_path
                    and progress['sha256'] == bindings[receipt_path]
                    and progress['independent_check'] in source['countercontrols']
                    and progress['scope'] == row['result']+' '+row['remaining'][0],
                    'unbound progress evidence')
    require(set(bindings) == required_sources, 'campaign evidence coverage changed')
    for stage in program['stages']:
        for accepted in stage.get('accepted_evidence', []):
            require(accepted.get('campaign_id') != CAMPAIGN
                    and accepted.get('path') not in receipts,
                    'finite theorem used as integrated gate acceptance')
        for progress in stage.get('progress_evidence', []):
            if progress.get('campaign_id') == CAMPAIGN:
                require(progress.get('id') in COMPONENTS
                        and stage['id'] in COMPONENTS[progress['id']][1],
                        'progress assigned outside its theorem domain')
    registered = unique(program['component_campaigns'], CAMPAIGN)
    require(registered['supports'] == list(OWNERS)
            and registered['scope'] == 'FINITE_COMPONENT_FAMILY'
            and registered['whole_program_frozen'] is False, 'registered domain changed')
    contract = registered['contract']
    require(contract['path'] == str(PATH)
            and contract['independent_check'] == 'model/completion/test_finite_observable_campaign.py',
            'wrong registered contract')
    check_file(str(PATH), contract['sha256'])
    return {'selected_issues': ISSUES, 'finite_components': len(rows),
            'canonical_premises_discharged': 0, 'integrated_gate_accepted': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.parse_args()
    print(json.dumps(validate(strict_json.load(ROOT/PATH)), indent=2))


if __name__ == '__main__':
    main()
