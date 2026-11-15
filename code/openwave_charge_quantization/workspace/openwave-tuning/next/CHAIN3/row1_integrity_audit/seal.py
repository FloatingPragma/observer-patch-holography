#!/usr/bin/env python3
"""Verify cross-artifact consistency and seal this independent audit."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
sys.dont_write_bytecode = True
from audit import HERE, PRODUCER, EXPECTED_MANIFEST, pins, digest, save


def main():
    result = json.loads((HERE/'RESULT.json').read_text())
    assert result['verdict'] == 'VALIDATED_WITH_RECORD_ONLY_FACADE'
    assert result['independent_pointer_laws']['every_effect_sum_equals_identity']
    assert result['complete_joint_laws']['independently_regenerated_joint_maps'] == 768
    assert result['observer_only_facade']['pointer_snapshot_defect']['facade_verdict']['compatible_pointer_initials'] == [0, 9]
    measured, _ = pins()
    assert measured == json.loads((HERE/'PINNED_INPUTS.json').read_text())
    native = json.loads((HERE/'NATIVE_POINTER_RECORDS.json').read_text())
    checks = result['actual_native_pointer_executions']
    assert native['executions'] == checks['executions'] == 276
    assert native['outcomes'] == checks['outcomes']
    for observed, expected in zip(native['records'], checks['rows']):
        assert digest(observed['native_executor_result']) == expected['meter_sha256']
        assert observed['record_only_compatible_columns'] == expected['tape_only_columns']
    trees = json.loads((HERE/'EXACT_TREE_COUNTS.json').read_text())
    assert trees['tree_count'] == 60 and trees['quotient_leaves'] == 2970 and trees['summed_leaf_record_lengths'] == 65910

    # Isolated standard-library execution, not the producer's expensive run.
    with tempfile.TemporaryDirectory(prefix='portable-check-', dir=HERE) as temporary:
        isolated = Path(temporary)
        manifest = json.loads((PRODUCER/'MANIFEST.json').read_text())
        for name in manifest['files_sha256']:
            target = isolated/name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(PRODUCER/name, target)
        shutil.copyfile(HERE/'observer_only.py', isolated/'observer_only.py')
        code = '''import json,sys
from pathlib import Path
root=Path(sys.argv[1]).resolve()
sys.path.insert(0,str(root))
import integrity as p
from observer_only import admit_record_only
out=[]
for total in (1,2):
    G=(total,)+(0,)*11
    claim,audit=p.execute(G,[(0,0,'C')],(0,),11,
        lambda w:p.meter.favorable_coins(11,(w,),2),[0]*30,2)
    del claim['meter']['final_loads']
    result=admit_record_only(claim,G,(0,),2,producer=p)
    assert result['admitted'] and result['readout']==audit['source_L1'][0]-G[0]
    out.append({'total':total,'admission':result})
for name in ('integrity','native_kernel','native_cycle_meter','closing_sweep','strict_tape_validator','observer_only'):
    assert Path(sys.modules[name].__file__).resolve().is_relative_to(root)
print(json.dumps(out))
'''
        process = subprocess.run([sys.executable, '-I', '-B', '-S', '-c', code, str(isolated)],
                                 cwd=isolated, capture_output=True, text=True, check=True)
        isolated_results = json.loads(process.stdout)
    save('PORTABILITY_CHECK.json', {'isolated_standard_library_execution': True,
        'python_flags': ['-I', '-B', '-S'], 'workspace_imports': False,
        'producer_manifest_sha256': EXPECTED_MANIFEST, 'executions': isolated_results})
    after, _ = pins()
    assert after == measured
    names = ['README.md', 'PROOF.md', 'audit.py', 'observer_only.py', 'branch_tree_counts.py', 'seal.py',
             'PINNED_INPUTS.json', 'RESULT.json', 'SUMMARY.json', 'EXACT_TREE_COUNTS.json',
             'NATIVE_POINTER_RECORDS.json', 'PORTABILITY_CHECK.json', 'run.log', 'tree_counts.log']
    manifest = {'schema': 'oph.row1.integrity.independent.audit.manifest.v1',
        'producer_manifest_sha256': EXPECTED_MANIFEST, 'verdict': result['verdict'],
        'files_sha256': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
        'reproduce': ['python3 audit.py', 'python3 branch_tree_counts.py', 'python3 seal.py'],
        'runtime': 'Python 3 standard library', 'source_and_receipts_unchanged': True,
        'reserved_or_load_four_executed': False, 'official_rows_changed': []}
    save('MANIFEST.json', manifest)
    print(json.dumps({'verdict': result['verdict'], 'manifest_sha256': hashlib.sha256((HERE/'MANIFEST.json').read_bytes()).hexdigest(),
                      'sealed_files': len(names), 'isolated_executions': len(isolated_results)}, indent=2))


if __name__ == '__main__':
    main()
