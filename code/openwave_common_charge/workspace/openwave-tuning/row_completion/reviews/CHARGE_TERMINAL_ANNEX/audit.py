"""Producer-independent terminal charge/clock annex review.

No producer, first-reviewer, or shared prior-review executor is imported.
The exact rational native pulse compiler and literal reset controller are
source definitions; the M6 executor, pointer reduction and targets are fresh.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from collections import defaultdict
import hashlib, importlib.util, json, math, sys
import numpy as np
from scipy.linalg import expm

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
PRODUCER=ROOT/'openwave-tuning/next/CHAIN5/common_charge_native_review'
SHARED=ROOT/'openwave-tuning/next/CHAIN5/dirac_kg_opposing_review'
EXPECTED='7104fb54312af99ba577ed6e2a69c3bc6390954f2e4e329db3627285750a1a2b'
read=lambda p:json.loads(p.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
save=lambda n,v:(HERE/n).write_text(json.dumps(v,indent=2,default=str)+'\n')
PINS={}


def pin(p,expected=None):
    actual=sha(p)
    if expected is not None:assert actual==expected,str(p)
    PINS[str(p.relative_to(ROOT))]=actual


@lru_cache(None)
def frozen(name):
    p=PRODUCER/'inputs'/name
    if p.exists():return p
    p=SHARED/'inputs'/name
    if p.exists():return p
    # Some source definitions are nested in the consumed field provider's
    # immutable snapshot, rather than copied at the outer snapshot root.
    matches=[p for p in (SHARED/'inputs').rglob(Path(name).name) if str(p).endswith('/'+name)]
    assert matches,name
    assert len({sha(p) for p in matches})==1,('ambiguous frozen definition',name)
    return min(matches,key=lambda p:len(p.parts))


def custody():
    pin(PRODUCER/'MANIFEST.json',EXPECTED)
    m=read(PRODUCER/'MANIFEST.json')
    for name,digest in m['files'].items():pin(ROOT/name,digest)
    inputs=read(PRODUCER/'INPUTS.json')
    pin(PRODUCER/'INPUTS.json')
    pin(SHARED/'MANIFEST.json',inputs['shared_review_manifest'])
    for group in ('copied_files','shared_frozen_files'):
        for name,digest in inputs[group].items():pin(frozen(name),digest)
    path=frozen('openwave-tuning/next/CHAIN2/local12_limit_audit/refinement_adapter.py')
    spec=importlib.util.spec_from_file_location('terminal_review_source_compiler',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


class RawM6:
    def __init__(self,compiler):
        self.compiler=compiler;self.count=0;self.verified=set()
        p=frozen('reverse-engineering-reality/code/source_selection_model/response.json')
        self.G=[]
        for two_blocks in read(p)['generators']:
            matrix=np.zeros((6,6),complex)
            for b,rows in enumerate(two_blocks):
                for i,row in enumerate(rows):
                    for j,q in enumerate(row):
                        matrix[3*b+i,3*b+j]=float(F(q[0]))+math.sqrt(5)*float(F(q[1]))+1j*(float(F(q[2]))+math.sqrt(5)*float(F(q[3])))
            assert np.linalg.norm(matrix+matrix.conj().T)<1e-12
            self.G.append(matrix)

    @lru_cache(None)
    def response(self,controls):
        return expm(sum((float(F(c))*g for c,g in zip(controls,self.G)),np.zeros((6,6),complex)))

    @lru_cache(None)
    def blocks(self,bits,targets):
        remaining=tuple(b for b in range(bits) if b not in targets)
        rows=[]
        for rest in range(1<<len(remaining)):
            base=sum(((rest>>j)&1)<<b for j,b in enumerate(remaining))
            rows.append([base+sum(((word>>(len(targets)-1-j))&1)<<b for j,b in enumerate(targets)) for word in range(1<<len(targets))])
        return np.array(rows)

    def run(self,tape,bits,archived=None):
        size=1<<bits;A=np.eye(size,dtype=complex);failed=np.zeros_like(A)
        budget=F();maxstage=0.
        if archived is not None:assert len(archived['first_export_records'])==len(tape)
        for j,event in enumerate(tape):
            signature=json.dumps(event,sort_keys=True)
            if signature not in self.verified:
                assert self.compiler.verify_pulse(event)
                self.verified.add(signature)
            U=self.response(tuple(event['controls_exact']))
            groups=self.blocks(bits,tuple(event['targets']));d=groups.shape[1]
            old=A[groups].copy()
            lost=np.matmul(U[d:,:d],old).reshape(-1,size)
            stage=lost.conj().T@lost;failed+=stage
            A[groups]=np.matmul(U[:d,:d],old)
            if archived is not None:
                record=archived['first_export_records'][j]
                error=float(np.max(abs(np.diag(stage).real-record['basis_column_probabilities'])))
                maxstage=max(maxstage,error);assert error<2e-13
                assert record['pulse']==j and record['targets']==event['targets']
                assert abs(float(np.trace(stage).real)-record['trace'])<2e-12
            budget+=F(event['full_response_operator_error_exact']);self.count+=1
        gram=float(np.linalg.norm(A.conj().T@A+failed-np.eye(size),2))
        assert gram<2e-10
        if archived is not None:
            target=np.array(archived['failure_effect_real'])+1j*np.array(archived['failure_effect_imag'])
            assert np.linalg.norm(failed-target,2)<2e-12
        return A,failed,budget,gram,maxstage


def one_gate(A,G,bit):
    for low in range(A.shape[0]):
        if low&(1<<bit):continue
        indices=[low,low+(1<<bit)]
        A[indices]=G@A[indices].copy()


def ideal(row):
    r=len(row['regions']);bits=4+r;clock=bits-1
    H=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
    X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]])
    M=np.eye(1<<bits,dtype=complex)
    if row['preparation']=='uniform_plane':
        for bit in range(3):one_gate(M,H,bit)
    elif row['preparation']=='basis_cell2':one_gate(M,X,1)
    else:
        assert row['preparation']=='nonuniform_entanglement_sensitive'
        one_gate(M,expm(1j*math.pi*Y/7),0);one_gate(M,H,2)
    one_gate(M,H,clock)
    n=row['macros']
    if n%2:
        permutation=np.array([code^7 for code in range(1<<bits)])
        M=M[permutation]
    one_gate(M,np.diag([np.exp(-1j*n/128),np.exp(1j*n/128)]),clock)
    values=[]
    for j,region in enumerate(row['regions']):
        v=[int(0 in region)-int(12*c in region) for c in range(8)]
        values.append(v);one_gate(M,H,3+j)
        # Integer flux means exp(i*pi*Phi)=(-1)^Phi exactly; avoid a
        # floating phase calculation for this source's complete code.
        for code in range(1<<bits):
            if code&(1<<(3+j)):M[code]*=(-1)**v[code&7]
        one_gate(M,H,3+j)
    assert values==row['flux_values_by_final_cell']
    setting=row['name'].split('_n')[0].split('_')[-1]
    if setting=='Y':one_gate(M,np.diag([np.exp(1j*math.pi/4),np.exp(-1j*math.pi/4)]),clock)
    else:assert setting=='X'
    one_gate(M,H,clock)
    return M


def pointer_effects(engine,tape):
    A,failed,bound,gram,_=engine.run(tape,3)
    result={}
    for law in ('C','E'):
        amplitudes=defaultdict(lambda:defaultdict(lambda:np.zeros((2,2),complex)))
        metadata={}
        # Canonical pair code states retain their total-two spectator shell.
        for code,(u,v) in enumerate(((0,0),(0,1),(1,0),(0,2))):
            shell=(2-u-v)
            K=A[np.ix_([code,code+4],[0,4])]
            for coin in (0,1):
                difference=u-v
                current=(difference-(2*coin-1)*(abs(difference)%2))//2
                label=current if law=='C' else difference
                key=(coin,label);post=(u-current,v+current,shell)
                amplitudes[key][post]+=K/math.sqrt(2)
                metadata[key]={'coin':coin,'label':label,'current':current,'bit':int(2*current+2*coin-1<0)}
        rows=[]
        for key in sorted(amplitudes):
            effect=sum((K.conj().T@K for K in amplitudes[key].values()),np.zeros((2,2),complex))
            rows.append((metadata[key],effect))
        result[law]=rows
        assert np.linalg.norm(sum(e for _,e in rows)+failed[np.ix_([0,4],[0,4])]-np.eye(2),2)<1e-12
    return result,failed[np.ix_([0,4],[0,4])],bound


def tensor_little(factors):
    result=np.ones((1,1),complex)
    for e in reversed(factors):result=np.kron(result,e)
    return result


def check_laws(A,failed,row,pointers,pointer_failure):
    vector=A[:,0];r=len(row['regions']);rho=np.zeros((1<<(r+1),)*2,complex)
    for cell in range(8):
        branch=vector[cell::8];rho+=np.outer(branch,branch.conj())
    count=0;error=0.
    for lawstring,stored in row['laws'].items():
        choices=[pointers[x] for x in lawstring]
        coarse=np.zeros(1<<(r+1));fine=list(product(*choices))
        assert len(fine)==len(stored['fine'])
        for combination,old in zip(fine,stored['fine']):
            records=[metadata for metadata,_ in combination]
            prob=float(np.trace(tensor_little([E for _,E in combination])@rho).real)
            assert records==old['records']
            delta=abs(prob-old['probability']);error=max(error,delta);assert delta<2e-12
            assert prob>-2e-13
            code=sum(metadata['bit']<<j for j,(metadata,_) in enumerate(combination))
            coarse[code]+=prob;count+=1
        first=[]
        for j in range(r+1):
            factors=[]
            for k,choices_k in enumerate(choices):
                factors.append(sum(E for _,E in choices_k) if k<j else pointer_failure if k==j else np.eye(2))
            first.append(float(np.trace(tensor_little(factors)@rho).real))
        assert max(abs(coarse-stored['coarse']))<2e-12
        assert max(abs(np.array(first)-stored['first_pointer_exports']))<2e-12
        assert abs(sum(coarse)+sum(first)+failed[0,0].real-1)<2e-10
    return count,error


def reset_checks():
    from model.native_clock import configuration as c
    target=(0,0,1,1,0,0,0,0,0,0,0,0)
    probabilities=list(c.success_probabilities(full=True,mode='reset',target=target))
    archived=read(PRODUCER/'FRESH_RESET_RECORDS.json')
    assert len(probabilities)==78
    assert probabilities==list(map(F,archived['exact_success_probabilities_all78']))
    assert max(1-p for p in probabilities)<=F(35,2**31)
    count=0
    for row in archived['successful_tapes']+archived['timeout_controls']:
        state=list(row['input'])
        for seam,coin,label in row['records']:
            a,b=seam;d=state[a]-state[b];j=(d-(2*coin-1)*(abs(d)%2))//2
            assert label==(j if row['law']=='C' else d)
            state[a]-=j;state[b]+=j;count+=1
        if 'support_columns' in row:assert tuple(state)==target
        else:assert state==row['output']
    # Check every Kraus partial injection on the native domain, not just the
    # records sampled for the numerical experiments.
    domain=tuple(c.native.states());seam_checks=0
    for seam in c.quantum.seams():
        for law in ('C','E'):
            blocks=defaultdict(dict)
            for state in domain:
                for coin in (0,1):
                    a,b=seam;d=state[a]-state[b];j=(d-(2*coin-1)*(abs(d)%2))//2
                    out=list(state);out[a]-=j;out[b]+=j
                    blocks[coin,j if law=='C' else d][state]=tuple(out);seam_checks+=1
            assert all(len(v)==len(set(v.values())) for v in blocks.values())
    return {'all78_exact_success_probabilities_match':True,
            'one_block_timeout_upper_exact':str(F(35,2**31)),
            'fine_reset_records_replayed':count,'all_seam_coin_column_checks':seam_checks,
            'two_block_uniform_timeout_upper_exact':str(F(35,2**31)**2),
            'comparison_marginalizes_private_success_tape':True}


def main():
    assert not (HERE/'MANIFEST.json').exists()
    compiler=custody();engine=RawM6(compiler)
    programs=read(PRODUCER/'TERMINAL_ANNEX_PROGRAMS.json')
    failures=read(PRODUCER/'TERMINAL_ANNEX_FAILURES.json')
    rows=read(PRODUCER/'TERMINAL_ANNEX_RECORDS.json')
    original=read(frozen('openwave-tuning/next/CHAIN5/common_charge_topology/NATIVE_PROGRAMS.json'))
    effects,fail,pointer_bound=pointer_effects(engine,original['native_READY_pointer'])
    assert len(programs)==len(rows)==36
    results=[];totalfine=0;worst=F()
    for row in rows:
        name=row['name'];bits=row['all_input_columns'].bit_length()-1
        A,B,bound,gram,stage=engine.run(programs[name],bits,failures[name])
        reference=ideal(row);operator_error=float(np.linalg.norm(A-reference,2));assert operator_error<2e-10
        fine,error=check_laws(A,B,row,effects,fail);totalfine+=fine
        r=len(row['regions']);registers=6+3*r;timeout=F(registers)*F(35,2**31)**2
        assert registers==row['reset_registers'] and row['reset_blocks_per_register']==2
        assert F(row['reset_timeout_union_bound_exact'])==timeout
        pulse=bound+(r+1)*pointer_bound
        assert pulse==F(row['all_native_pulse_error_bound_exact'])
        complete=pulse+timeout
        assert complete==F(row['complete_terminal_record_TV_bound_exact'])
        assert F(row['clock_phase_exact'])==F(row['macros'],64)
        assert F(row['physical_time_exact'])==F(row['macros'],8)
        worst=max(worst,complete)
        results.append({'name':name,'pulses':len(programs[name]),'columns':1<<bits,
            'all_column_operator_error':operator_error,'Gram_residual':gram,
            'all_column_first_failure_difference':stage,'fine_probability_difference':error,
            'fine_probabilities':fine,'complete_TV_exact':str(complete)})
        print('REPLAY',name,'OK',flush=True)
    reset=reset_checks()
    report={'producer_independent':True,'producer_executor_or_review_executor_imported':False,
        'source_definition_dependency':'Frozen refinement_adapter exact pulse certificate and literal native reset controller only.',
        'experiments':results,'pointer_word_bound_exact':str(pointer_bound),
        'complete_terminal_TV_upper_exact':str(worst),'complete_terminal_TV_display':float(worst),
        'raw_M6_pulses_replayed_including_pointer':engine.count,
        'complete_fine_probabilities_recomputed':totalfine,'reset':reset,
        'floating_residuals_are_diagnostics_not_exact_error_certificates':True}
    save('NATIVE_AUDIT.json',report);save('SOURCE_PINS.json',PINS)
    print('ALL36 TERMINAL EXPERIMENTS PASS',engine.count,totalfine,float(worst),flush=True)


if __name__=='__main__':main()
