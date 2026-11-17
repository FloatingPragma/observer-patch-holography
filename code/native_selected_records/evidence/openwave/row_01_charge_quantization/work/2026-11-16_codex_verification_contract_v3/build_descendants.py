"""Mechanical descendants: preserve calculations; enforce checks and record contract."""
import ast
import hashlib
import json
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent
INTEGRATION = HERE.parents[2]/'row_01_charge_quantization/reviews/2026-11-16_codex_resubmission_integration'
ORIGINAL = INTEGRATION/'outgoing_v2/openwave/xperiments/m10_oph/research'

PREFIX = '''
# V3 verification contract: these checks are unconditional, also under python -O.
def _require(condition, message="Scientific check failed"):
    if not condition:
        raise AssertionError(message)

'''
SUFFIX = '''

def _verify_or_record(payload):
    """Default: compare saved JSON. Explicit first record: exclusive file creation."""
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--verify', type=Path)
    mode.add_argument('--record-to', type=Path)
    args = parser.parse_args()
    actual = json.loads(payload)
    if args.record_to is not None:
        # Never overwrite. Scientific checks have already completed successfully.
        with args.record_to.open('x') as stream:
            stream.write(payload+'\\n')
    else:
        saved = args.verify or Path(__file__).resolve().parent.parent/'data'/(__NAME__+'_summary.json')
        # Compare the complete typed JSON result. No tolerance or ignored science fields.
        expected = json.loads(saved.read_text())
        _require(json.dumps(actual, sort_keys=True, allow_nan=False) ==
                 json.dumps(expected, sort_keys=True, allow_nan=False),
                 'Saved evidence differs; verification never replaces it: '+str(saved))
    print(payload)

if __name__ == '__main__':
    _verify_or_record(_verification_payload)
'''


def unconditional(text):
    """Replace each assert with the same always-evaluated predicate and message."""
    tree = ast.parse(text)
    lines = text.splitlines(keepends=True)
    offsets = [0]
    for line in lines:offsets.append(offsets[-1]+len(line))
    changes=[]
    for node in ast.walk(tree):
        if isinstance(node, ast.Assert):
            test = ast.get_source_segment(text, node.test)
            message = ast.get_source_segment(text,node.msg) if node.msg else repr('Scientific check at original line '+str(node.lineno))
            changes.append((offsets[node.lineno-1]+node.col_offset,
                            offsets[node.end_lineno-1]+node.end_col_offset,
                            f'_require({test}, {message})'))
    for start,end,new in sorted(changes, reverse=True):text=text[:start]+new+text[end:]
    # Keep original docstring and future imports before new definitions.
    tree=ast.parse(text);at=0
    for node in tree.body:
        if isinstance(node,ast.Expr) and isinstance(node.value,ast.Constant) and isinstance(node.value.value,str):
            at=node.end_lineno
        elif isinstance(node,ast.ImportFrom) and node.module=='__future__':at=node.end_lineno
        else:break
    lines=text.splitlines(keepends=True)
    return ''.join(lines[:at])+PREFIX+''.join(lines[at:]),len(changes)


def build():
    target=HERE/'descendants';target.mkdir(exist_ok=False)
    (target/'scripts').mkdir();shutil.copytree(ORIGINAL/'data',target/'data')
    records=[]
    for original in sorted((ORIGINAL/'scripts').glob('m10_4_*.py')):
        text,count=unconditional(original.read_text())
        tree=ast.parse(text)
        prints=[n for n in ast.walk(tree) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call)
                and isinstance(n.value.func,ast.Name) and n.value.func.id=='print']
        if len(prints)!=1:raise RuntimeError('Expected exactly one final JSON print: '+original.name)
        node=prints[0];old=ast.get_source_segment(text,node)
        text=text.replace(old,'_verification_payload = '+ast.get_source_segment(text,node.value.args[0]),1)
        text+=SUFFIX.replace('__NAME__',repr(original.stem))
        compile(text,original.name,'exec')
        path=target/'scripts'/original.name;path.write_text(text)
        records.append({'name':original.name,'original_sha256':hashlib.sha256(original.read_bytes()).hexdigest(),
                        'descendant_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'unconditional_checks':count,
                        'calculation_changed':False,'saved_data_changed':False})
    (HERE/'DESCENDANT_BUILD.json').write_text(json.dumps(records,indent=2)+'\n')


if __name__=='__main__':build()
