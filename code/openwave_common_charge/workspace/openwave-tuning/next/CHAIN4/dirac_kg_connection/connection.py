"""Free routed-Dirac factor of the pinned native oscillator-field source."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from collections import deque
import hashlib,importlib.util,json,math,sys
sys.dont_write_bytecode=True
import numpy as np
from scipy.linalg import expm

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
SPECIES=ROOT/'openwave-tuning/next/CHAIN3/species_reuse_connection'
FIELD=ROOT/'openwave-tuning/next/CHAIN3/native_field_connection'
PINNED={SPECIES/'MANIFEST.json':'c36680ec7e0021158a66b643a15c7b74549a0d4ba460fd40722c2ec83bf9c9a2',
 FIELD/'MANIFEST.json':'719036fe46769e68202731f765045dc3d75ebc8fe6e901e85319ce435776779f',
 ROOT/'openwave-tuning/common_dynamics/spin/MANIFEST.json':'65075e0bc14f883c5cdb877171650b08ee57c8d1bdffce281e68937bb07eb401',
 ROOT/'openwave-tuning/row_completion/ROW_29/MANIFEST.json':'993744d14c2da8d66a93dc76422f81799074febad9914ae33bb60bf455cf00b6'}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    result={}
    for p,digest in PINNED.items():
        assert sha(p)==digest,str(p);result[str(p.relative_to(ROOT))]=digest
        d=read(p)
        groups=[(p.parent,d.get('files_sha256',{})),(p.parent/'inputs',d.get('producer_files',{})),
                (ROOT,d.get('source_sha256',{})),(ROOT,d.get('files',{})),(ROOT,d.get('source_pins',{}))]
        for base,files in groups:
            for n,h in files.items():
                q=base/n;assert sha(q)==h,str(q);rel=str(q.relative_to(ROOT))
                assert rel not in result or result[rel]==h;result[rel]=h
    return result

verify()
spec=importlib.util.spec_from_file_location('pinned_species_mass_source',SPECIES/'source.py')
species=importlib.util.module_from_spec(spec);spec.loader.exec_module(species)
native=species.native;words=species.field_words;g=words.g
matter=g.source.matter
MASSES=tuple(map(F,read(SPECIES/'REGISTER.json')['mass_midpoints_exact']))

def parameter_macro(coexact_cap,harmonic_cap,slices,bits,h=F(1,8),mass=F(1),coupling=F(1)):
    """g=1 is literally field25; g scales currents and g² Coulomb energy.

    g=0 is the neutral/free intervention. Its field modes and lattice free
    duration remain unchanged; the physical macro time is h.
    """
    coupling=F(coupling);h=F(h);mass=F(mass)
    if coupling==1:
        yield from words.macro_word(coexact_cap,harmonic_cap,slices,bits,h,mass);return
    info=words.layout(coexact_cap,harmonic_cap)
    yield 'Z',(info['spin'][0],),native.angle(0,0,-h*mass)
    for axis in range(3):
        diagonal=list(words.spin_diagonalizer(axis,info['spin']));yield from diagonal
        if coupling:
            for cell in range(8):
                for sign in (-1,1):
                    predicate=tuple((q,(cell>>j)&1) for j,q in enumerate(info['cell']))+((info['spin'][1],int(sign<0)),)
                    for mode,(cap,wires) in enumerate(zip(info['caps'],info['modes'])):
                        gamma=coupling*g.midpoint(g.coefficient(cell,axis,sign,mode,bits))
                        yield from words.translation_word(cap,gamma,slices,bits,wires,predicate)
        target=[matter.route(cell,axis,1)[0] for cell in range(8)]
        assert target==[matter.route(cell,axis,-1)[0] for cell in range(8)]
        difference={i^j for i,j in enumerate(target)};assert len(difference)==1
        mask=difference.pop();assert mask and mask&(mask-1)==0
        yield from words.remap(native.primitive_templates()['flip'],(info['cell'][mask.bit_length()-1],))
        yield from native.negate_template(diagonal)
    yield from words.free_field_word(info,F(1),slices,bits,include_coulomb=False)
    if coupling:
        yield from words.diagonal_values_word([-coupling**2*row['coulomb'] for row in g.data()['sectors']],info['cell'])
    yield 'Z',(info['clock'],),native.angle(0,0,-h/16)

def flavored_parameter_macro(coexact_cap,harmonic_cap,slices,bits,masses=MASSES,h=F(1,8),coupling=F(1)):
    """Same task29 flavor insertion, on the complete g-parameter source."""
    info=words.layout(coexact_cap,harmonic_cap)
    flavor=(info['memory_bits'],info['memory_bits']+1)
    yield from species.flavor_mass_word(masses,h,flavor,info['spin'][0])
    remaining=iter(parameter_macro(coexact_cap,harmonic_cap,slices,bits,h,F(1),coupling))
    assert next(remaining)==('Z',(info['spin'][0],),native.angle(0,0,-F(h)))
    yield from remaining

def multi_x(values,target,scratch=()):
    if not values:return [('flip',(target,))]
    if len(values)==1:
        q,v=next(iter(values.items()))
        return ([] if v else [('flip',(q,))])+[('cx',(q,target))]+([] if v else [('flip',(q,))])
    return native.frozen.multi_control(values,target,tuple(scratch))

def signed_shift_gates(width,coordinate,rail,scratch=()):
    """Positive rail0 increments, negative rail1 decrements modulo2^width."""
    assert len(coordinate)==width
    out=[]
    for sign in (1,-1):
        for j in reversed(range(width)):
            controls={rail:int(sign<0)}
            controls.update({coordinate[k]:int(sign>0) for k in range(j)})
            out+=multi_x(controls,coordinate[j],scratch)
    return out

def gates_word(gates):
    for kind,targets in gates:yield from words.remap(native.primitive_templates()[kind],targets)

def spatial_layout(width):
    # Source graph cell order is x-major, then y, then z. Each coordinate's
    # binary bits use little endian within its block.
    coordinates=[tuple(range(2+(2-axis)*width,2+(3-axis)*width)) for axis in range(3)]
    end=2+3*width
    return {'spin':(0,1),'coordinates':coordinates,'flavor':(end,end+1),
            'scratch':tuple(range(end+2,end+2+max(0,width-2))),
            'probe':end+2+max(0,width-2),'clock':end+3+max(0,width-2),
            'bits':end+4+max(0,width-2)}

def spatial_matter_word(width,h,masses=MASSES):
    info=spatial_layout(width)
    yield from species.flavor_mass_word(masses,h,info['flavor'],info['spin'][0])
    for axis in range(3):
        diag=list(words.spin_diagonalizer(axis));yield from diag
        yield from gates_word(signed_shift_gates(width,info['coordinates'][axis],1,info['scratch']))
        yield from native.negate_template(diag)

def reduced_matter_word(k,h,masses=MASSES,flavor=(2,3)):
    """k values are rational multiples of pi; obtained from allowed modes."""
    yield from species.flavor_mass_word(masses,h,flavor,0)
    for axis,theta_pi in enumerate(k):
        label=('XX','XY','XZ')[axis]
        yield from native.pauli(label,native.angle(-F(theta_pi)))

def controlled(word,probe):
    """Control every native source primitive, retaining scalar I/J phases."""
    for kind,targets,a in word:
        if probe in targets:raise ValueError('fresh probe required')
        if kind=='H':
            yield 'Y',targets,native.angle(F(1,8))
            yield from controlled([('Z',targets,a)],probe)
            yield 'Y',targets,native.angle(F(-1,8))
        elif kind=='J':
            for label,c in (('II',F(3,4)),('ZI',F(1,4)),('IZ',F(1,4)),('ZZ',F(-1,4))):
                for ctrl,d in (('I',F(1,2)),('Z',F(-1,2))):
                    yield from words.remap(native.pauli(label+ctrl,tuple(c*d*x for x in a)),tuple(targets)+(probe,))
        else:
            assert kind in ('I','X','Y','Z')
            for ctrl,c in (('I',F(1,2)),('Z',F(-1,2))):
                yield from words.remap(native.pauli(kind+ctrl,tuple(c*x for x in a)),tuple(targets)+(probe,))

def bloch(k,h,m):
    answer=expm(-1j*float(h)*float(m)*matter.BETA)
    for axis,value in enumerate(k):answer=expm(-1j*float(value)*matter.ALPHA[axis])@answer
    return answer

@lru_cache(None)
def route(size,cell,axis,sign):
    if size==2:return matter.route(cell,axis,sign)
    graph=native.spatial.graph(size);coords={tuple(c):i for i,c in enumerate(graph.cell_coordinates)}
    c=list(graph.cell_coordinates[cell]);c[axis]=(c[axis]+sign)%size;other=coords[tuple(c)]
    if sign<0:
        destination,path,current=route(size,other,axis,1);assert destination==cell
        return other,tuple(reversed(path)),-current
    a,b=matter.spin.common.spatial.links()[axis][0]
    def internal(start,end):
        queue=deque([(start,(start,))]);seen={start}
        while queue:
            x,p=queue.popleft()
            if x==end:return p
            for y in graph.neighbours[x]:
                if y//12==start//12 and y not in seen:queue.append((y,p+(y,)));seen.add(y)
        raise AssertionError('connected cell')
    path=internal(12*cell,12*cell+a)+internal(12*other+b,12*other)
    edges={tuple(edge):j for j,edge in enumerate(graph.seams)};J=np.zeros(len(edges),np.int64)
    for x,y in zip(path,path[1:]):
        pair=tuple(sorted((x,y)));J[edges[pair]]+=1 if x==pair[0] else -1
    return other,path,J

@lru_cache(None)
def raw_generators():
    response=read(species.FIELD/'inputs/reverse-engineering-reality/code/source_selection_model/response.json')
    def q(a,b):return float(F(a))+math.sqrt(5)*float(F(b))
    out=[]
    for blocks in response['generators']:
        G=np.zeros((6,6),complex)
        for block,rows in enumerate(blocks):
            for i,row in enumerate(rows):
                for j,cell in enumerate(row):G[3*block+i,3*block+j]=q(*cell[:2])+1j*q(*cell[2:])
        out.append(G)
    return out

@lru_cache(None)
def response(controls):return expm(sum((float(F(c))*G for c,G in zip(controls,raw_generators())),np.zeros((6,6),complex)))
@lru_cache(None)
def indices(bits,targets):
    remaining=[q for q in range(bits) if q not in targets]
    return np.array([[sum(((a>>(len(targets)-1-j))&1)<<q for j,q in enumerate(targets))+
                      sum(((b>>j)&1)<<q for j,q in enumerate(remaining)) for b in range(1<<len(remaining))]
                    for a in range(1<<len(targets))])

def execute(events,bits,start):
    M=np.asarray(start,complex).reshape(1<<bits,-1).copy();columns=M.shape[1]
    initial=M.conj().T@M;failure=np.zeros((columns,columns),complex)
    bound=F()
    for event in events:
        native.verify_pulse(event);U=response(tuple(event['controls_exact']));d=1<<len(event['targets'])
        ix=indices(bits,tuple(event['targets']));old=M[ix]
        lost=np.einsum('ab,bcr->acr',U[d:,:d],old).reshape(-1,columns);failure+=lost.conj().T@lost
        M[ix]=np.einsum('ab,bcr->acr',U[:d,:d],old);bound+=F(event['full_response_operator_error_exact'])
    gram=float(np.linalg.norm(M.conj().T@M+failure-initial,2));assert gram<3e-10
    return M,failure,gram,bound

def events(word,bits=100):return list(words.pulses(word,bits))
