"""Inspect retained process evidence and prove mechanical calculation preservation."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
INTEGRATION=HERE.parents[2]/'row_01_charge_quantization/reviews/2026-11-16_codex_resubmission_integration'


class Calculation(ast.NodeTransformer):
    def visit_FunctionDef(self,node):
        if node.name in ('_require','_verify_or_record'):return None
        return self.generic_visit(node)

    def visit_If(self,node):
        if any(isinstance(x,ast.Call) and isinstance(x.func,ast.Name) and x.func.id=='_verify_or_record' for x in ast.walk(node)):
            return None
        return self.generic_visit(node)

    def visit_Assign(self,node):
        if len(node.targets)==1 and isinstance(node.targets[0],ast.Name) and node.targets[0].id=='_verification_payload':
            return ast.Expr(value=ast.Call(func=ast.Name(id='print',ctx=ast.Load()),args=[node.value],keywords=[]))
        return self.generic_visit(node)

    def visit_Expr(self,node):
        if isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Name) and node.value.func.id=='_require':
            args=node.value.args;message=args[1] if len(args)>1 else None
            if isinstance(message,ast.Constant) and str(message.value).startswith('Scientific check at original line '):message=None
            return ast.Assert(test=args[0],msg=message)
        return self.generic_visit(node)


def calculation_equivalence():
    rows=[]
    for path in sorted((HERE/'descendants/scripts').glob('*.py')):
        original=INTEGRATION/'outgoing_v2/openwave/xperiments/m10_oph/research/scripts'/path.name
        expected=ast.dump(ast.parse(original.read_text()),include_attributes=False)
        actual=ast.dump(Calculation().visit(ast.parse(path.read_text())),include_attributes=False)
        rows.append({'script':path.name,'calculation_AST_identical_after_explicit_contract_and_assert_normalization':actual==expected,
                     'original_sha256':hashlib.sha256(original.read_bytes()).hexdigest(),
                     'descendant_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    return rows


def main():
    rows=calculation_equivalence()
    if not all(r['calculation_AST_identical_after_explicit_contract_and_assert_normalization'] for r in rows):
        raise AssertionError('A scientific calculation changed')
    verdicts=[]
    for folder in ['validation_payloads','validation_producers','validation_harness','validation_supplement','validation_pair_receipt']:
        summary=json.loads((HERE/folder/'SUMMARY.json').read_text())
        if not summary['passed']:raise AssertionError(folder+' failed')
        for result in summary['results']:
            record=next((p for p in (HERE/folder).glob('*.json') if p.stem==result['label']),None)
            if record is None:raise AssertionError('Missing real process record '+str(result))
            full=json.loads(record.read_text());exit_code=full['exit_code']
            if result.get('expected_success',True):
                if exit_code!=0:raise AssertionError('Unexpected process failure')
            else:
                if exit_code==0 or full.get('timeout'):raise AssertionError('Mutant did not reject')
                if result['label']=='failed_gate_record_nonzero':
                    receipt=json.loads(full['stdout'])
                    if receipt['passed'] or receipt['checks']['conditional_transport_prediction']:
                        raise AssertionError('Wrong failed scientific gate')
                elif result['label']=='real_job_changed_saved_baseline':
                    receipt=json.loads(full['stdout'])
                    if receipt['passed'] or not receipt['saved_evidence_unchanged'] or receipt['jobs'][0]['error']!='Saved scientific output differs':
                        raise AssertionError('Wrong saved-baseline rejection')
                elif not any(word in full['stderr'] for word in ('AssertionError','FileExistsError','FileNotFoundError')):
                    raise AssertionError('Non-scientific/unexpected process rejection: '+result['label'])
            verdicts.append({'group':folder,'label':result['label'],'exit_code':exit_code,'passed':True,
                             'receipt_sha256':hashlib.sha256(record.read_bytes()).hexdigest()})
    wrappers=json.loads((HERE/'WRAPPER_AUDIT.json').read_text())
    for row in wrappers['results']:
        record=HERE/'wrapper_processes'/(row['label']+'.json');full=json.loads(record.read_text())
        if full.get('exit_code')!=1 or 'AssertionError' not in full.get('stderr','') or full.get('timeout'):
            raise AssertionError('Wrong wrapper mutation rejection')
    source=json.loads((HERE/'validation_harness/SOURCE_REPLAY.json').read_text())
    if not source['passed'] or not source['saved_evidence_unchanged'] or len(source['jobs'])!=28:
        raise AssertionError('Historical full replay failed')
    original=json.loads((HERE/'original_processes/reservoir_actual_failed_gate.json').read_text())
    receipt=json.loads(original['stdout'])
    if original['exit_code']!=0 or receipt['passed'] or receipt['checks']['conditional_transport_prediction']:
        raise AssertionError('Original direct producer defect not demonstrated')
    result={'passed':True,'same_group':True,'external_G6':False,'scientific_emergence_credit':False,
            'payload_calculation_equivalence':rows,'process_verdicts':verdicts,
            'distinct_historical_entry_predicate_mutations':len(wrappers['results']),
            'full_replay_jobs':len(source['jobs']),'max_job_seconds':max(j['seconds'] for j in source['jobs']),
            'saved_historical_source_unchanged':True}
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--verify',type=Path,default=HERE/'FINAL_VERDICT.json')
    mode.add_argument('--record-to',type=Path)
    args=parser.parse_args()
    if args.record_to:
        with args.record_to.open('x') as stream:stream.write(json.dumps(result,indent=2)+'\n')
    elif result!=json.loads(args.verify.read_text()):
        raise AssertionError('Saved audit verdict differs; verification does not replace it')
    print(json.dumps({k:v for k,v in result.items() if k not in ('payload_calculation_equivalence','process_verdicts')}))


if __name__=='__main__':main()
