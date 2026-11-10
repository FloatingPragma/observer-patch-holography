"""Typed field/current-record evidence projection from the sealed field checkpoint.

Only four enumerated prose fields, primary evidence, four architecture premises
and their existing binding unions change. Physical statuses and every unlisted
semantic value are conserved. The predecessor validators remain immutable.
"""
from __future__ import annotations
import argparse
import copy
import json
import re
from pathlib import Path
from tools import strict_json
from . import field_navigation as retained

ROOT=Path(__file__).resolve().parents[2]
AUTHORITY='model/completion/field_current_records_navigation_authority.json'
AUTHORITY_SHA256='f9b296a2fa4419c64ee46becb95a3a284eb498d88cd3669a7bba0a8a3a58bb24'
REGISTRY=retained.REGISTRY
TRACKING=retained.TRACKING
DICTIONARY=retained.DICTIONARY
JSON_FILES=retained.JSON_FILES
canonical=retained.canonical
digest=retained.digest
bytes_digest=retained.bytes_digest
require=retained.require
at=retained.at


def release_transition(approved):
    transition = approved['release_identity_transition']
    require(transition == {'file': REGISTRY, 'json_path': ['release_id'],
                           'before': 'r2063', 'after': 'r2064'},
            'unreviewed release identity transition')
    require(approved['release_id_transition_included'] is True,
            'release identity transition is not explicitly admitted')
    return transition


def normalize(documents, approved):
    transition = release_transition(approved)
    value = copy.deepcopy(documents)
    require(value[REGISTRY]['release_id'] == transition['after'], 'current release identity differs')
    value[REGISTRY]['release_id'] = transition['before']
    return retained.normalize(value, approved)


def project(baseline, approved):
    transition = release_transition(approved)
    require(baseline[REGISTRY]['release_id'] == transition['before'], 'retained release identity differs')
    value = retained.project(baseline, approved)
    value[REGISTRY]['release_id'] = transition['after']
    normalize(value, approved)
    return value


def authority(root=ROOT):
    root=Path(root)
    raw=(root/AUTHORITY).read_bytes()
    require(bytes_digest(raw)==AUTHORITY_SHA256,'field/current-record navigation authority changed')
    approved=strict_json.loads(raw.decode())
    require(approved['immutable_checkpoint']=='a0367dc26b07ad959c915054dbfcb3f9828af295',
            'retained field checkpoint changed')
    for key in ('predecessor_custody_authority','predecessor_navigation_authority'):
        row=approved[key]
        require(bytes_digest((root/row['path']).read_bytes())==row['sha256'],
                'retained field authority changed')
    for name,row in approved['baseline_files'].items():
        require(bytes_digest((root/row['retained_path']).read_bytes())==row['sha256'],
                'retained field ledger bytes changed: '+name)
    for name,expected in approved['new_evidence_sha256'].items():
        path=(root/name).resolve()
        require(path.is_relative_to(root.resolve()) and path.is_file()
                and bytes_digest(path.read_bytes())==expected,
                'reviewed field/current-record evidence changed: '+name)
    release_transition(approved)
    require(set(approved['actual_release_identity_sha256']) == {
        'paper/release_info.tex', 'paper/paper_release_manifest.json'},
        'release identity input scope changed')
    for name, expected in approved['actual_release_identity_sha256'].items():
        require(bytes_digest((root/name).read_bytes()) == expected, 'actual release identity source changed: '+name)
    require(strict_json.load(root/'paper/paper_release_manifest.json')['release_id'] == 'r2064',
            'publication manifest release differs')
    require(re.search(r'\\newcommand\{\\OPHPaperReleaseID\}\{r2064\}',
                      (root/'paper/release_info.tex').read_text()) is not None,
            'paper release macro differs')
    return approved


def validate(root=ROOT,approved=None):
    root=Path(root);approved=authority(root) if approved is None else approved
    documents={name:strict_json.load(root/name) for name in JSON_FILES}
    normalize(documents,approved)
    raw=(root/DICTIONARY).read_bytes();suffix=approved['dictionary_append'].encode()
    require(bool(suffix) and raw.endswith(suffix),'reviewed architecture definitions changed')
    require(bytes_digest(raw[:-len(suffix)])==approved['baseline_files'][DICTIONARY]['sha256'],
            'retained architecture dictionary changed')
    return {'schema':'oph.field_current_records_navigation_validation.v1',
            'retained_checkpoint':approved['immutable_checkpoint'],
            'prose_fields_replaced':len(approved['prose_replacements']),
            'added_evidence_files':len(approved['new_evidence_sha256']),
            'added_assumptions':list(approved['new_assumption_definitions']),
            'claim_count':len(documents[REGISTRY]['claims']),
            'physical_statuses_and_nonlisted_semantics_exactly_preserved':True,
            'scientific_producer_replay_claimed_by_this_projection':False,
            'release_identity_transition':approved['release_identity_transition'],
            'current_files':{name:bytes_digest((root/name).read_bytes())
                             for name in (*JSON_FILES,DICTIONARY)}}


def apply(root=ROOT):
    root=Path(root);approved=authority(root)
    raw={name:(root/name).read_bytes() for name in (*JSON_FILES,DICTIONARY)}
    baseline={name:strict_json.load(root/approved['baseline_files'][name]['retained_path']) for name in JSON_FILES}
    candidate={name:strict_json.loads(raw[name].decode()) for name in JSON_FILES}
    transition=release_transition(approved)
    if candidate[REGISTRY].get('release_id') in (transition['before'], transition['after']):
        candidate[REGISTRY]['release_id']=transition['before']
    if candidate==baseline and bytes_digest(raw[DICTIONARY])==approved['baseline_files'][DICTIONARY]['sha256']:
        result=project(baseline,approved)
        outputs={name:(json.dumps(result[name],indent=2,ensure_ascii=False)+'\n').encode() for name in JSON_FILES}
        outputs[DICTIONARY]=raw[DICTIONARY]+approved['dictionary_append'].encode()
        for name in raw:require((root/name).read_bytes()==raw[name],'ledger changed during projection')
        for name,value in outputs.items():(root/name).write_bytes(value)
    return validate(root,approved)


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--apply',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();print(json.dumps(apply() if args.apply else validate(),indent=2,sort_keys=True))


if __name__=='__main__':main()
