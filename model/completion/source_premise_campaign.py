"""Check the bounded five-issue campaign and its explicit non-composition.

This verifies custody and scope; mathematical proofs remain in the source
packages. Separate finite algebras and references are never joined by name.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json


ROOT = Path(__file__).resolve().parents[2]
PATH = Path('model/completion/source_premise_campaign.json')
CAMPAIGN = 'OPH-SOURCE-PREMISE-CLOSURE-2026-11-04'
OWNERS = {
    'adaptive_sweep_reference': ('FA-110', ['CP-02'], [1138]),
    'one_particle_search': ('FA-111', ['CP-03', 'CP-07'], [1139, 1143]),
    'bipartite_spin_limit': ('FA-112', ['CP-06'], [1142]),
    'finite_modular_inclusion': ('FA-113', ['CP-09'], [1145]),
    'finite_area_response': ('FA-114', ['CP-09'], [1145]),
    'moving_contact_scattering': ('FA-115', ['CP-06', 'CP-07'], [1142, 1143]),
}
ISSUES = [1138, 1139, 1142, 1143, 1145]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate(packet, root=ROOT):
    root = Path(root).resolve()
    require(packet.get('schema') == 'oph.source_premise_campaign.v1'
            and packet.get('campaign_id') == CAMPAIGN, 'wrong source-premise contract')
    require(packet.get('scope') == 'FINITE_COMPONENT_FAMILY'
            and packet.get('whole_program_frozen') is False,
            'finite campaign cannot freeze the whole program')
    require(packet.get('selected_issues') == ISSUES, 'five-issue ownership changed')
    require(packet.get('physical_comparison_count_added') == 0
            and type(packet.get('physical_comparison_count_added')) is int
            and packet.get('canonical_premises_discharged') == []
            and packet.get('integrated_gate_accepted') is False
            and packet.get('common_interacting_QFT_and_Einstein_model') is False,
            'component evidence cannot accept physical completion')
    require(packet.get('specification_readout_allowed') is True,
            'the declared simulator and observer specifications are admitted inputs')
    rows = packet.get('component_results', [])
    require([r.get('id') for r in rows] == list(OWNERS), 'missing or duplicate component')
    finite = strict_json.load(root/'model/completion/finite_attachment_map.json')
    attachments = {a['id']: a for a in finite['attachments']}
    program = strict_json.load(root/'claims/completion_program.json')
    stages = {stage['id']: stage for stage in program['stages']}
    for row in rows:
        identifier = row['id']
        attachment, owners, issues = OWNERS[identifier]
        require(row.get('attachment') == attachment and row.get('issues') == issues,
                f'component owner changed: {identifier}')
        source = attachments[attachment]
        require(source['owners'] == owners
                and source['classification'] == 'realized_on_claimed_domain'
                and source['physical_identification'] is False,
                f'finite map scope changed: {identifier}')
        require(row.get('result') == source['conclusion']
                and row.get('remaining') == source['remaining'],
                f'component result or boundary differs from its reviewed attachment: {identifier}')
        require(row.get('interface') == source['interface'],
                f'component state/algebra/operation/clock/refinement changed: {identifier}')
        evidence = row.get('evidence', {})
        require(isinstance(evidence, dict) and bool(evidence), f'missing source custody: {identifier}')
        for name, digest in evidence.items():
            path = (root/name).resolve()
            require(path.is_relative_to(root) and path.is_file(), f'missing source: {name}')
            require(hashlib.sha256(path.read_bytes()).hexdigest() == digest,
                    f'stale component source: {name}')
        receipt = f'model/{identifier}/results.json'
        require(receipt in evidence, f'unbound receipt: {identifier}')
        for owner in owners:
            progress = [p for p in stages[owner].get('progress_evidence', [])
                        if p.get('campaign_id') == CAMPAIGN and p.get('id') == identifier]
            require(len(progress) == 1 and progress[0]['path'] == receipt
                    and progress[0]['sha256'] == evidence[receipt],
                    f'completion progress does not bind this component: {identifier}/{owner}')
            require(not any(p.get('campaign_id') == CAMPAIGN
                            for p in stages[owner].get('accepted_evidence', [])),
                    f'component campaign entered accepted evidence: {owner}')
    boundary = packet.get('composition_boundary', {})
    require(boundary.get('automatic_cross_domain_join') is False
            and boundary.get('common_state_reference_constructed') is False
            and boundary.get('common_clock_and_refinement_constructed') is False
            and boundary.get('distinct_domains') == list(OWNERS)
            and bool(boundary.get('reason')),
            'distinct source states cannot be silently composed')
    return {'selected_issues': ISSUES, 'reviewed_components': len(rows),
            'canonical_premises_discharged': 0, 'integrated_gate_accepted': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.parse_args()
    print(json.dumps(validate(strict_json.load(ROOT/PATH)), indent=2))


if __name__ == '__main__':
    main()
