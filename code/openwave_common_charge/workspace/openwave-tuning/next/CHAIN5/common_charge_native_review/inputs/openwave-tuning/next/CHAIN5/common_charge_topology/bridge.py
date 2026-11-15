"""Pinned g=1 common source and exact Gauss-fiber geometry."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,importlib.util,json,sys
sys.dont_write_bytecode=True
import numpy as np
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
COMMON=ROOT/'openwave-tuning/next/CHAIN4/dirac_kg_connection'
ROW1=ROOT/'openwave-tuning/row_completion/ROW_01'
EXPECTED='398482e1a7e05fcba687cabe6056b9f870c8e364fad47abd4f8244db4a2187f8'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    assert sha(COMMON/'MANIFEST.json')==EXPECTED
    result={str((COMMON/'MANIFEST.json').relative_to(ROOT)):EXPECTED}
    for section in ('files','source_pins'):
        for name,digest in read(COMMON/'MANIFEST.json')[section].items():
            assert sha(ROOT/name)==digest,name;result[name]=digest
    for name in ('CRITERION.json','PROOF.md','CONTRACT.json','MANIFEST.json','upstream/MODELS.md','upstream/m5_1_winding.py'):
        p=ROW1/name;result[str(p.relative_to(ROOT))]=sha(p)
    criterion=read(ROW1/'CRITERION.json')
    assert criterion['criterion_sha256']=='3af50af0955c64a0aa2019fa84af3cf0dd36e6c63a5abc76eb2c60c9b413cb6f'
    assert hashlib.sha256(criterion['criterion'].encode()).hexdigest()==criterion['criterion_sha256']
    return result
verify()
sys.path.insert(0,str(COMMON))
import connection as common
native=common.native;words=common.words;g=common.g;source=g.source;matter=common.matter
execute=common.execute;events=common.events

def save(name,obj):
    (HERE/name).write_text(json.dumps(obj,indent=2,default=str)+'\n');print(name,'PASS',flush=True)

def exact_geometry():
    graph,B,A,periods,winding,labels=source.field.cell_complex.build(2)
    assert B*A==s.zeros(96,208)
    d=g.data();Q=s.Matrix(d['Q']);K=A*A.T
    # B K=0 proves B Y=0 exactly for Y=K Q Lambda^-1. The positive
    # polar normalization does not alter its range; no rounded Q is used
    # as an exact physical transverse frame.
    assert B*K==s.zeros(96,288)
    harmonics=s.Matrix.hstack(*(s.Matrix(x) for x in d['harmonics']))
    assert B*harmonics==s.zeros(96,3) and A.T*harmonics==s.zeros(208,3)
    particulars=s.Matrix.hstack(*(s.Matrix(row['particular']) for row in d['sectors']))
    charges=s.zeros(96,8)
    for cell in range(8):charges[0,cell]+=1;charges[12*cell,cell]-=1
    assert B*particulars==charges
    return {'graph':graph,'B':B,'A':A,'K':K,'Q':Q,'harmonics':harmonics,
            'particulars':particulars,'charges':charges,'periods':periods,'winding':winding,'source_data':d}
