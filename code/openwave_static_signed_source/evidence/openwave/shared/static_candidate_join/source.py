"""One signed source. Protocol context chooses addressed seams, not a row law."""
from pathlib import Path
import importlib.util,json,sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
INTERFACE=HERE.parent/'compact_static_adjudication'
spec=importlib.util.spec_from_file_location('static_signed_coordinate_source',INTERFACE/'encoded_source.py')
S=importlib.util.module_from_spec(spec);sys.modules[spec.name]=S;spec.loader.exec_module(S)
GRAPH=json.loads((HERE.parent/'localized_rest_join/evidence/graph.json').read_text())
SEAMS=tuple(map(tuple,GRAPH['seams']))
CELL=tuple(e for e in SEAMS if e[1]<12)

def context(name):
    if name=='cell_compare':return CELL
    if name=='federation_energy':return SEAMS
    raise ValueError('unregistered observer program')

def endpoint(label,seam,coin):
    return S.apply_pair_column(None,seam,coin,label)

def represented(label):
    if label==0:return {}
    p,r=divmod(int(label),256)
    if p<240 and r<240 and p!=r:return {p:1,r:-1}
    raise ValueError('invalid source label')

def basis(q):
    q={k:v for k,v in q.items() if v}
    if not q:return 0
    p,=[k for k,v in q.items() if v==1]
    r,=[k for k,v in q.items() if v==-1]
    assert len(q)==2
    return 256*p+r
