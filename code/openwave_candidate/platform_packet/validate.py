"""Validate preparation completeness; never award physical criteria."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FROZEN_SHA256 = 'd6e36d90bc253f5af641760de2e6efc50e181f16e17c8f5334c31d0d7928085f'
COMMIT = '5d686cb049b82cc55f8822210a254d5833be23a7'
MODELS_SHA256 = 'a69895e4af4f6bc180748c4e9a8ae11d745c75462650aea44c3a34bd687a139d'
CANDIDATE = 'ring-hm-mod6-local-coarse-v1.4-proposed'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate(packet, frozen_bytes=None):
    """Reject altered literals, invalid IDs and any status promotion."""
    raw = frozen_bytes if frozen_bytes is not None else (ROOT / 'criteria_frozen.json').read_bytes()
    require(hashlib.sha256(raw).hexdigest() == FROZEN_SHA256, 'frozen criteria digest mismatch')
    frozen = json.loads(raw)
    require(frozen['openwave_commit'] == COMMIT and frozen['models_sha256'] == MODELS_SHA256,
            'upstream provenance mismatch')
    require(frozen['upstream_url'] == f'https://github.com/openwave-labs/openwave/blob/{COMMIT}/MODELS.md',
            'upstream URL mismatch')
    require(packet['schema'] == 'oph.openwave.platform_packet.v1', 'unknown packet schema')
    require(packet['platform_commit'] == COMMIT, 'platform commit mismatch')
    require(packet['candidate_id'] == CANDIDATE, 'candidate mismatch')
    require(packet['official_status_change'] is False, 'official promotion prohibited')
    rows = packet['rows']
    require(type(rows) is list and len(rows) == 31, 'exactly 31 rows required')
    for expected, row, literal in zip(range(1, 32), rows, frozen['criteria']):
        require(type(row['id']) is int and row['id'] == expected, 'row IDs must be exact ordered integers')
        require(row['candidate_id'] == CANDIDATE, 'row candidate mismatch')
        require(row['criterion'] == literal['criterion'] and row['literal_gate'] == literal['simplest_test'],
                'criterion literal mismatch')
        require(type(row['required_before_run']) is list and bool(row['required_before_run'])
                and all(type(item) is str and bool(item.strip()) for item in row['required_before_run']),
                'missing execution prerequisites')
        require(row['status'] == 'blocked' and row['result'] is None
                and row['official_status'] == 'untested', 'row promotion prohibited')
    return {'preparation_valid': True, 'rows': 31, 'physical_gates_ready': 0,
            'official_status_change': False}


if __name__ == '__main__':
    print(json.dumps(validate(json.loads((ROOT / 'packet.json').read_text()))))
