"""Fixed-policy closing observer; v2 remains an immutable dependency."""
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent.parent/'row1_v2'))
import closing_reader as v2

ORDER = tuple(range(7, 30))+tuple(range(7))


def check_order(records):
    if not isinstance(records, list) or any(not isinstance(r, dict) for r in records):
        raise ValueError('closing record list required')
    if tuple(r.get('seam') for r in records) != ORDER:
        raise ValueError('records do not follow the registered closing order')


def observe(records, region):
    check_order(records)
    result = v2.observe_closing(records, region)
    result['source_channel_premise'] = 'Every source movement since sweep start is one of the thirty recorded events.'
    result['channel_completeness_authenticated'] = False
    result['measurement_time'] = 'End of the registered sweep, conditional on the source-channel premise.'
    return result


def admit(records, region, pointer):
    try:
        check_order(records)
    except ValueError as error:
        return {'admitted': False, 'reason': str(error)}
    result = v2.admit_closing(records, region, pointer, retries=2)
    if result['admitted']:
        result['channel_completeness_authenticated'] = False
    return result
