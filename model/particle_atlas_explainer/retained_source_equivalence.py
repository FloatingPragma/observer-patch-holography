"""Authenticate retained primary evidence without rewriting its source IDs.

Only six enumerated live receipt files may differ from the historical input
view, and only by their reviewed custody fields. Every other input byte and
every original contract guard is retained. Current navigation has its own
digest; it is not a replacement physical source, state, reference or clock ID.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path

from tools import strict_json
from tools.renew_canonical_live_receipts import canonical, carry_metadata, metadata_view
from tools.replay_model_receipts import bindings

ROOT = Path(__file__).resolve().parents[2]
AUTHORITY = 'model/particle_atlas_explainer/retained_source_equivalence_authority.json'
AUTHORITY_SHA256 = 'a5cad7a312b7588507bff6472fe5bf54b685a11e3e8cf831cc2ec4499cdbd527'
REPLAY = 'model/particle_atlas_explainer/retained_source_replay.json'
OUTPUT = 'model/particle_atlas_explainer/retained_source_navigation.json'
REPLAY_SHA256 = 'f6b8dcc383f2deaaf5375e6c9015c2466a6b21ad4c919db835ac47b54eb3f533'
RETAINED_NAMES = {
    'model/particle_primary_regional_fermion/results.json',
    'model/particle_primary_composite_spin/results.json',
    'model/particle_primary_reduced_sector/results.json',
    'model/particle_primary_loop_fusion/results.json',
    'model/standard_model_source_map/candidates/correspondence.json',
    'model/particle_candidate_confirmation/results.json',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_authority(root=ROOT):
    root = Path(root)
    if sha(root / AUTHORITY) != AUTHORITY_SHA256:
        raise ValueError('retained source equivalence authority changed')
    value = strict_json.load(root / AUTHORITY)
    if set(value['retained_receipts']) != RETAINED_NAMES:
        raise ValueError('retained proof inventory changed')
    return value


def validate_transition(name, spec, *, root=ROOT):
    root = Path(root)
    archive = root / spec['retained_path']
    if sha(archive) != spec['retained_sha256']:
        raise ValueError('retained source receipt changed: ' + name)
    old, live = strict_json.load(archive), strict_json.load(root / name)
    pointers = spec['reviewed_nested_receipt_pointers']
    admitted = carry_metadata(name, old, live, pointers=pointers, root=root)
    if canonical(admitted) != canonical(live):
        raise ValueError('retained/live source science differs: ' + name)
    if canonical(metadata_view(name, old, pointers, root)) != canonical(metadata_view(name, live, pointers, root)):
        raise ValueError('retained/live source scope differs: ' + name)
    for _, relative, expected, path in bindings(root, name, live):
        if sha(path) != expected:
            raise ValueError('live source custody is stale: ' + name + ':' + relative)
    return {'retained_sha256': spec['retained_sha256'], 'current_sha256': sha(root / name),
            'retained_path': spec['retained_path'], 'non_provenance_payload_equal': True}


def historical_input_paths(*, root=ROOT, approved=None):
    root = Path(root).resolve()
    approved = load_authority(root) if approved is None else approved
    result = {}
    for name, expected in approved['retained_input_sha256'].items():
        transition = approved['allowed_live_receipt_transitions'].get(name)
        path = (root / (transition['retained_path'] if transition else name)).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError('retained input missing or escaping: ' + name)
        if transition and transition['retained_sha256'] != expected:
            raise ValueError('retained input differs from transition authority: ' + name)
        if sha(path) != expected:
            raise ValueError('unreviewed historical source change: ' + name)
        result[name] = path
    return result


def validate_inputs(*, root=ROOT, approved=None):
    approved = load_authority(root) if approved is None else approved
    historical_input_paths(root=root, approved=approved)
    transitions = {name: validate_transition(name, spec, root=root)
                   for name, spec in approved['allowed_live_receipt_transitions'].items()}
    return transitions


def physical_specification(member):
    """Remove only the explicitly identified content-address fields."""
    result = copy.deepcopy(member)
    del result['source']['implementation_sha256']
    del result['state']['source_id']
    del result['reference']['source_id']
    del result['reference']['state_id']
    del result['clock']['source_id']
    del result['ids']
    del result['member_id']
    return result


def navigation(*, root=ROOT):
    root = Path(root)
    approved = load_authority(root)
    transitions = validate_inputs(root=root, approved=approved)
    if sha(root / REPLAY) != REPLAY_SHA256:
        raise ValueError('recorded historical execution changed')
    replay = strict_json.load(root / REPLAY)
    if replay['authority_sha256'] != AUTHORITY_SHA256 or not replay['all_retained_receipts_exactly_replayed']:
        raise ValueError('historical replay does not bind this retained source view')
    if set(replay['receipts']) != set(approved['retained_receipts']):
        raise ValueError('historical replay omitted a retained owner')
    for name, spec in approved['retained_receipts'].items():
        if (sha(root / name) != spec['sha256'] or sha(root / spec['path']) != spec['sha256']
                or replay['receipts'][name]['retained_sha256'] != spec['sha256']):
            raise ValueError('retained scientific evidence changed: ' + name)
    value = {'schema': 'oph.retained_primary_navigation.v1',
             'authority_sha256': AUTHORITY_SHA256, 'historical_replay_sha256': sha(root / REPLAY),
             'projection_source_sha256': sha(Path(__file__)),
             'retained_physical_ids': replay['members'], 'live_custody_transitions': transitions,
             'physical_ids_rewritten': False, 'physical_specification_changed': False,
             'identity_intertwiner_domain': 'same complete finite basis, Hamiltonian, prepared state, reference and clock formulas',
             'live_producer_receipts_claimed_for_retained_primary_evidence': False}
    value['navigation_digest'] = hashlib.sha256(canonical(value).encode()).hexdigest()
    return value


def authenticate_pins(pins, *, root=ROOT):
    """Resolve only authenticated historical exceptions; return current file pins."""
    root = Path(root).resolve()
    approved = load_authority(root)
    result = {}
    needed = False
    for name, expected in pins.items():
        path = (root / name).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError('missing or escaping atlas source: ' + name)
        actual = sha(path)
        if actual != expected:
            spec = approved['allowed_live_receipt_transitions'].get(name)
            if spec is None or expected != spec['retained_sha256']:
                raise ValueError('unreviewed stale atlas input: ' + name)
            needed = True
        result[name] = actual
    if needed:
        navigation(root=root)
        result.update({AUTHORITY: AUTHORITY_SHA256, REPLAY: sha(root / REPLAY),
                       'model/particle_atlas_explainer/retained_source_equivalence.py': sha(Path(__file__))})
        for spec in approved['allowed_live_receipt_transitions'].values():
            result[spec['retained_path']] = spec['retained_sha256']
    return result


def read_receipt(path, *, root=ROOT):
    root = Path(root).resolve()
    path = Path(path).resolve()
    name = str(path.relative_to(root))
    approved = load_authority(root)
    if name not in approved['retained_receipts']:
        raise ValueError('receipt is outside retained primary admission: ' + name)
    navigation(root=root)
    value = strict_json.load(path)
    authenticate_pins(value['source_sha256'], root=root)
    return value


def validate_receipt_value(name, value, *, root=ROOT):
    """Use the exact result of the unchanged historical owner and validator."""
    approved = load_authority(root)
    if name not in approved['retained_receipts']:
        raise ValueError('unknown retained proof owner')
    navigation(root=root)
    expected = strict_json.load(Path(root) / approved['retained_receipts'][name]['path'])
    if value != expected:
        raise ValueError('retained proof science, identity or gate changed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--renew', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = navigation()
    if args.renew:
        (ROOT / OUTPUT).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    elif strict_json.load(ROOT / OUTPUT) != value:
        raise ValueError('retained source navigation projection is stale')
    print('Verified retained physical identities and distinct current navigation digest ' + value['navigation_digest'])


if __name__ == '__main__':
    main()
