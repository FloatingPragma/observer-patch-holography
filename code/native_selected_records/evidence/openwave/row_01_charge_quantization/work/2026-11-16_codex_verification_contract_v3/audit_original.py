"""Process-level contract probes; all altered code/data live in temporary copies."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INTEGRATION = ROOT/'row_01_charge_quantization/reviews/2026-11-16_codex_resubmission_integration'
SOURCE = INTEGRATION/'source_export_v2/source'
SCRIPTS = INTEGRATION/'outgoing_v2/openwave/xperiments/m10_oph/research/scripts'
ENV = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
ENV.pop('PYTHONOPTIMIZE', None)


def invoke(command, cwd, label, env=ENV):
    begin = time.monotonic()
    try:
        p = subprocess.run(command, cwd=cwd, env=env, capture_output=True, text=True, timeout=60)
        result = dict(label=label, command=command, exit_code=p.returncode,
                      stdout=p.stdout, stderr=p.stderr, seconds=time.monotonic()-begin)
    except subprocess.TimeoutExpired as error:
        result = dict(label=label, command=command, timeout=True, seconds=time.monotonic()-begin,
                      stdout=str(error.stdout), stderr=str(error.stderr))
    (HERE/'original_processes'/f'{label}.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('stdout','stderr')}), flush=True)
    return result


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise RuntimeError(f'Expected unique mutation site: {old!r}')
    return text.replace(old, new)


def main():
    (HERE/'original_processes').mkdir(exist_ok=False)
    results = []
    # Change real scientific rejection conditions; no passed label is substituted.
    mutants = {
        'charge': ('assert entries == {"native": 120, "swap_hold": 0}', 'assert entries == {"native": 121, "swap_hold": 0}'),
        'energy': ('assert np.max(abs(R @ R - np.eye(33))) == 0', 'assert np.max(abs(R @ R - np.eye(33))) > 1'),
        'spin': ('assert visibility < 0 if law == "native" and turns == 1 else visibility > 0', 'assert visibility > 0 if law == "native" and turns == 1 else visibility < 0'),
        'common': ('assert impurity < .002 and pair_impurity < .01', 'assert impurity < -1 and pair_impurity < -1'),
        'pair': ('assert energy[1]<OLD_PAIR', 'assert energy[1]<1000'),
    }
    with tempfile.TemporaryDirectory(prefix='verification-contract-original-') as temp:
        folder = Path(temp)
        for name, (old, new) in mutants.items():
            script = folder/f'm10_4_{name}.py'
            script.write_text(replace_once((SCRIPTS/script.name).read_text(), old, new))
            for optimized in (False, True):
                command = [sys.executable, '-B']+(['-O'] if optimized else [])+[str(script)]
                results.append(invoke(command, folder, name+('_optimized' if optimized else '_normal')))
        packet = SOURCE/'evidence/openwave/row_01_charge_quantization/work/2026-11-15_codex_pr723_vacuum_reservoir'
        copy = folder/'reservoir';shutil.copytree(packet, copy)
        script = copy/'experiment.py'
        script.write_text(replace_once(script.read_text(),
            "'conditional_transport_prediction': abs(observed-float(expected)) <= threshold}",
            "'conditional_transport_prediction': abs(observed-float(expected)) < -1}"))
        results.append(invoke([sys.executable, '-B', str(script), '--output', str(folder/'failed_record')], copy, 'reservoir_actual_failed_gate'))
        # A scientific value is changed without changing any code or custody hash.
        # The direct independent reader must reject its value comparison.
        receipt = copy/'RECEIPT.json';data=json.loads(receipt.read_text())
        data['heldout_nonzero_observed'] += 1;receipt.write_text(json.dumps(data))
        before=(receipt.read_bytes(), receipt.stat().st_mtime_ns)
        results.append(invoke([sys.executable, '-B', str(copy/'verify.py')], copy, 'reservoir_changed_saved_observation'))
        if before != (receipt.read_bytes(), receipt.stat().st_mtime_ns):
            raise RuntimeError('Verification modified the saved mutated expectation')
        # Harness-only parser probe, explicitly not scientific evidence.
        for name, output in [('earlier_failure', [{'passed':False},{'passed':True}]),
                             ('nested_gate_failure',[{'gates':{'source_selection':False}}])]:
            probe=folder/(name+'.py');probe.write_text('import json\n'+''.join('print(json.dumps('+repr(x)+'))\n' for x in output))
            driver=folder/'driver.py'
            driver.write_text('import runpy,json\nn=runpy.run_path('+repr(str(SOURCE/'verify.py'))+')\nf=n["replay_current"]\nf.__globals__["verify"]=lambda:{}\nr=f([('+repr(str(folder))+','+repr(probe.name)+',[])])\nprint(json.dumps(r))\nraise SystemExit(0 if r["passed"] else 1)\n')
            results.append(invoke([sys.executable,'-B',str(driver)],folder,'harness_'+name))
    summary={'processes':[{k:v for k,v in r.items() if k not in ('stdout','stderr')} for r in results],
             'integrity_or_gate_mutations_only':True,'scientific_emergence_credit':False,
             'sealed_source_modified':False}
    (HERE/'ORIGINAL_AUDIT.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__ == '__main__':
    main()
