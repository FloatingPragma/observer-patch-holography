"""Independent six-level replay of the transferred terminal charge observable."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from collections import defaultdict
import hashlib,json,math,sys
import numpy as np
from scipy.linalg import expm

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
PINS=json.loads((HERE/'INPUTS.json').read_text())

def read(name):
    h=PINS['snapshots'][name];b=(HERE/'inputs'/h).read_bytes()
    assert hashlib.sha256(b).hexdigest()==h
    return json.loads(b)

ANNEX='openwave-tuning/next/CHAIN5/common_charge_native_review/'
BASE='openwave-tuning/next/CHAIN5/common_charge_topology/'
GEN='reverse-engineering-reality/code/source_selection_model/response.json'
source_generators=read(GEN)['generators']
GENERATORS=[]
for blocks in source_generators:
    G=np.zeros((6,6),complex)
    for b,block in enumerate(blocks):
        for i,row in enumerate(block):
            for j,entry in enumerate(row):
                a,c,d,e=map(F,entry);G[3*b+i,3*b+j]=float(a)+math.sqrt(5)*float(c)+1j*(float(d)+math.sqrt(5)*float(e))
    GENERATORS.append(G)

@lru_cache(None)
def six_response(controls):
    return expm(sum((float(F(x))*G for x,G in zip(controls,GENERATORS)),np.zeros((6,6),complex)))

@lru_cache(None)
def groups(bits,targets):
    other=[b for b in range(bits) if b not in targets]
    return np.array([[sum(((x>>(len(targets)-1-j))&1)<<b for j,b in enumerate(targets))+
                      sum(((y>>j)&1)<<b for j,b in enumerate(other)) for y in range(1<<len(other))]
                    for x in range(1<<len(targets))])

def run(tape,bits):
    M=np.eye(1<<bits,dtype=complex);failure=np.zeros_like(M);stages=[];bound=F()
    for event in tape:
        U=six_response(tuple(event['controls_exact']));d=1<<len(event['targets']);ix=groups(bits,tuple(event['targets']))
        old=M[ix].copy();lost=np.einsum('ij,jkl->ikl',U[d:,:d],old).reshape(-1,1<<bits)
        E=lost.conj().T@lost;failure+=E;stages.append(E.diagonal().real)
        M[ix]=np.einsum('ij,jkl->ikl',U[:d,:d],old);bound+=F(event['full_response_operator_error_exact'])
    gram=float(np.linalg.norm(M.conj().T@M+failure-np.eye(1<<bits),2));assert gram<2e-10
    return M,failure,stages,bound,gram

def one(M,G,bit):
    ids=[i for i in range(len(M)) if not (i>>bit)&1]
    for a in ids:
        b=a+(1<<bit);x=M[a].copy();y=M[b].copy();M[a]=G[0,0]*x+G[0,1]*y;M[b]=G[1,0]*x+G[1,1]*y

def target(row):
    r=len(row['regions']);bits=r+4;N=1<<bits;clock=bits-1;M=np.eye(N,dtype=complex)
    H=np.array([[1,1],[1,-1]],complex)/np.sqrt(2)
    if row['preparation']=='uniform_plane':
        for bit in range(3):one(M,H,bit)
    elif row['preparation']=='basis_cell2':one(M,np.array([[0,1],[1,0]]),1)
    else:
        a=math.pi/7;one(M,np.array([[math.cos(a),math.sin(a)],[-math.sin(a),math.cos(a)]]),0);one(M,H,2)
    one(M,H,clock)
    n=row['macros'];permutation=np.arange(N)^(7 if n%2 else 0);M=M[permutation]
    one(M,np.diag(np.exp([-1j*n/128,1j*n/128])),clock)
    # Independent exact parity action of H exp(i*pi*Phi) H, without Walsh gates.
    for j,region in enumerate(row['regions']):
        values=[int(0 in region)-int(12*c in region) for c in range(8)]
        assert values==row['flux_values_by_final_cell'][j]
        permutation=[i^((abs(values[i&7])%2)<<(3+j)) for i in range(N)];M=M[permutation]
    if '_Y_' in row['name']:one(M,np.diag(np.exp([1j*math.pi/4,-1j*math.pi/4])),clock)
    one(M,H,clock);return M

def pointer(tape):
    U,Fp,_,bound,gram=run(tape,3);K=U[:,[0,4]];failed=Fp[np.ix_([0,4],[0,4])];effects={}
    for law in ('C','E'):
        records=defaultdict(lambda:defaultdict(lambda:np.zeros((2,2),complex)))
        metadata={}
        for index,(u,v) in enumerate(((0,0),(0,1),(1,0),(0,2))):
            spectator=tuple([1]*(2-u-v)+[0]*(8+u+v))
            for coin in (0,1):
                delta=u-v;current=(delta-(2*coin-1)*(abs(delta)%2))//2
                label=current if law=='C' else delta;key=(coin,label)
                emitted=(u-current,v+current)+spectator
                records[key][emitted]+=K[[index,index+4]]/math.sqrt(2)
                metadata[key]={'coin':coin,'label':label,'current':current,'bit':int(2*current+2*coin-1<0)}
        effects[law]=[(metadata[k],sum((A.conj().T@A for A in blocks.values()),np.zeros((2,2),complex))) for k,blocks in sorted(records.items())]
        assert np.linalg.norm(sum(E for meta,E in effects[law])+failed-np.eye(2),2)<2e-12
    return effects,failed,bound

def tensor(es):
    out=np.ones((1,1),complex)
    for E in es:out=np.kron(E,out)
    return out

def main():
    # The independently frozen compiler checks exact rational pulse certificates.
    sys.path.insert(0,str(ROOT/'openwave-tuning/next/CHAIN4/dirac_kg_connection'))
    import connection as compiler
    compiler_pins=compiler.verify();seen=set()
    tapes=read(ANNEX+'TERMINAL_ANNEX_PROGRAMS.json');rows=read(ANNEX+'TERMINAL_ANNEX_RECORDS.json')
    failures=read(ANNEX+'TERMINAL_ANNEX_FAILURES.json');ptape=read(BASE+'NATIVE_PROGRAMS.json')['native_READY_pointer']
    for tape in [ptape,*tapes.values()]:
        for event in tape:
            key=json.dumps(event,sort_keys=True)
            if key not in seen:compiler.native.verify_pulse(event);seen.add(key)
    effects,pfail,pbound=pointer(ptape)
    maxoperator=maxgram=maxfine=maxfailure=0.;fine_count=pulses=0;worst=F();checked=[]
    for row in rows:
        name=row['name'];tape=tapes[name];bits=len(row['regions'])+4;r=bits-4
        M,failed,stages,bound,gram=run(tape,bits);ideal=target(row)
        error=float(np.linalg.norm(M-ideal,2));assert error<2e-10
        maxoperator=max(maxoperator,error);maxgram=max(maxgram,gram);pulses+=len(tape)
        for got,arch in zip(stages,failures[name]['first_export_records']):
            delta=float(np.max(abs(got-np.array(arch['basis_column_probabilities']))));assert delta<1e-22;maxfailure=max(maxfailure,delta)
        rho=sum((np.outer(M[c::8,0],M[c::8,0].conj()) for c in range(8)),np.zeros((1<<(r+1),1<<(r+1)),complex))
        for laws,archive in row['laws'].items():
            coarse=np.zeros(1<<(r+1));actual=[]
            for combination in product(*(effects[law] for law in laws)):
                records=[meta for meta,E in combination];prob=float(np.trace(tensor([E for meta,E in combination])@rho).real)
                actual.append((records,prob));coarse[sum(meta['bit']<<j for j,meta in enumerate(records))]+=prob
            assert len(actual)==len(archive['fine'])
            for (records,prob),arch in zip(actual,archive['fine']):
                assert records==arch['records'];delta=abs(prob-arch['probability']);assert delta<2e-11;maxfine=max(maxfine,delta);fine_count+=1
            fail=[]
            for j in range(r+1):
                es=[sum(E for meta,E in effects[laws[k]]) if k<j else pfail if k==j else np.eye(2) for k in range(r+1)]
                fail.append(float(np.trace(tensor(es)@rho).real))
            assert max(abs(np.array(fail)-archive['first_pointer_exports']))<2e-20
            assert max(abs(coarse-archive['coarse']))<2e-11
            assert abs(coarse.sum()+sum(fail)+failed[0,0].real-1)<2e-10
        timeout=F(6+3*r)*F(35,2**31)**2
        exact=bound+(r+1)*pbound+timeout;assert exact==F(row['complete_terminal_record_TV_bound_exact']);worst=max(worst,exact)
        checked.append({'name':name,'all_columns':1<<bits,'operator_residual':error,'Gram_residual':gram})
    result={'schema':'oph.common_charge.terminal_opposing_review.v1','experiments':len(rows),'native_pulses_executed':pulses+len(ptape),'unique_exact_pulse_certificates':len(seen),
            'fine_C_E_probabilities_recomputed':fine_count,'maximum_operator_residual':maxoperator,'maximum_Gram_residual':maxgram,'maximum_fine_probability_difference':maxfine,
            'maximum_first_export_probability_difference':maxfailure,'complete_terminal_TV_exact':str(worst),'complete_terminal_TV_display':float(worst),
            'all_reset_timeout_bounds_checked':True,'reset_fine_law_proof_review':'Separate exact 78-state reset review consumed; not a second execution of its 45385 records.',
            'original_huge_field_tape_executed':False,'terminal_observable_annex_pass':True,'cases':checked,'official_rows_changed':[]}
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    (HERE/'COMPILER_PINS.json').write_text(json.dumps(compiler_pins,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'}),flush=True)

if __name__=='__main__':main()
