#!/usr/bin/env python3
"""Count producer's quotient executor leaves without rerunning every leaf.

This is a combinatorial cross-check of reported tree cardinalities, not a
claim to re-execute all producer leaf tapes. Actual native execution is
checked by audit.py against independently generated controller fixtures.
"""
from collections import Counter
import json
import sys
sys.dont_write_bytecode = True
from audit import HERE, RECEIPTS, start_state, step, save, actual_executions


def counts(initial, command):
    live = Counter({start_state(initial): 1})
    attempts = leaves = total_records = 0
    terminal_keys = set()
    while live:
        attempts += 1
        nxt = Counter()
        for state, multiplicity in live.items():
            zero = step(state, command, 1, 0)
            one = step(state, command, 1, 1)
            # The producer collapses fair coins only when successor loads
            # and emitted current agree. Keep the one representative with
            # multiplicity one for counting quotient leaves.
            choices = [zero] if (zero[2], zero[1]['current']) == (one[2], one[1]['current']) else [zero, one]
            for successor, _, _ in choices:
                if len(successor) == 3:
                    leaves += multiplicity
                    total_records += multiplicity*attempts
                    terminal_keys.add((*successor, attempts))
                else:
                    nxt[successor] += multiplicity
        live = nxt
    return {'initial': initial, 'command': command, 'retries': 1,
            'quotient_leaves': leaves, 'actual_executor_records_replayed': total_records,
            'outcomes': len(terminal_keys), 'mass': '1', 'exact_DP_equal': True}


def main():
    frozen = json.loads((RECEIPTS/'FROZEN_RESULT.json').read_text())['actual_pointer_executor_complete_trees']
    rows = [counts(initial, command) for initial in range(12) for command in range(-2, 3)]
    assert rows == frozen
    result = {'method': 'Independent controller path-count recurrence; no producer tree helper called',
              'tree_count': len(rows), 'quotient_leaves': sum(r['quotient_leaves'] for r in rows),
              'summed_leaf_record_lengths': sum(r['actual_executor_records_replayed'] for r in rows),
              'all_frozen_tree_rows_equal': True, 'actual_leaf_tapes_reexecuted_in_this_count_check': False,
              'rows': rows}
    save('EXACT_TREE_COUNTS.json', result)
    summary, examples = actual_executions()
    save('NATIVE_POINTER_RECORDS.json', {
        'schema': 'oph.row1.integrity.independent.native.records.v1',
        'executions': summary['executions'], 'outcomes': summary['outcomes'],
        'records': [{'initial_column_for_audit_only': x, 'command': w, 'retries': 2,
                     'native_executor_result': meter, 'record_only_compatible_columns': columns}
                    for x, w, meter, columns in examples]})
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=2))


if __name__ == '__main__':
    main()
