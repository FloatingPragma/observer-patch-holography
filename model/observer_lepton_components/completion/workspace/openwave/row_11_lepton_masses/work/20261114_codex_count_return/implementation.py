"""N=1 retained-count experiment; compiled native fibre, actual source words.

The tensor factorization is only a floating evaluation device. Its discarded
norms are reported; uniform full-map bounds concern untruncated source words.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import sys, json, math, time, hashlib, importlib.util, argparse
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
FAMILY=HERE.parent/'20261114_codex_family_filter'
spec=importlib.util.spec_from_file_location('sealed_family',FAMILY/'implementation.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
np=m.np;c=m.c;iv=m.iv;synth=m.synth;interpreter=m.interpreter
from model.lepton_register_round import carrier,walk as wk
from model.observer.native import native
GATE=3e-8
SVD_CUTOFF=1e-12
TRUNCATIONS=[]
ALL_WORDS=[]
BLOCKS=[]
ARRAYS={}
EFFECTS={}
NORM_CACHE={}

def sha(b):return hashlib.sha256(b).hexdigest()
def exact_l1(a):
    return sum((F(float(x.real)).__abs__()+F(float(x.imag)).__abs__() for x in np.asarray(a).ravel()),F())
def frobenius_upper(a):
    a=np.asarray(a)
    if not a.size:return F()
    peak=max(float(np.max(np.abs(a.real))),float(np.max(np.abs(a.imag))))
    return 2*F(peak)*math.isqrt(a.size+2*math.isqrt(a.size)+1)
def represented_norm_bound(a):
    """1+eta upper bound from an exact rational Gram of finite binary64 data.

    All denominators are powers of two, so one common integer scale makes
    the Gram products exact without repeated rational matrix arithmetic.
    """
    a=np.asarray(a,complex)
    if not a.shape[1]:return F()
    key=(a.shape,sha(a.tobytes()))
    if key in NORM_CACHE:return NORM_CACHE[key]
    pairs=[(float(z.real).as_integer_ratio(),float(z.imag).as_integer_ratio()) for z in a.ravel()]
    d=max(max(r[1],i[1]) for r,i in pairs)
    real=np.array([r[0]*(d//r[1]) for r,i in pairs],object).reshape(a.shape)
    imag=np.array([i[0]*(d//i[1]) for r,i in pairs],object).reshape(a.shape)
    gr=real.T@real+imag.T@imag
    gi=real.T@imag-imag.T@real
    for j in range(a.shape[1]):gr[j,j]-=d*d
    eta=F(max(sum(abs(int(gr[j,k]))+abs(int(gi[j,k])) for k in range(a.shape[1])) for j in range(a.shape[1])),d*d)
    bound=1+eta;NORM_CACHE[key]=bound
    return bound
def remember_word(name,word,units):
    interpreter.source_word(word,9)
    ALL_WORDS.append((name,word));BLOCKS.append(dict(name=name,**m.pulse_bound(len(word),units)))
    return word

def native_census():
    w=carrier.carrier_walk();cfg=wk.config(2,2);records=[];calls=[];checks={}
    exact=True;native_ok=True
    for cell in ('X','Y'):
        rows=[]
        for k,(a,b) in enumerate(w.addresses):
            image=wk.images(cfg,0,a,b)
            exact &=len(image)==1 and image[0][2]==1
            outcfg,rec,amp=image[0]
            exact &=wk.images(outcfg,rec,a,b)==((cfg,0,F(1)),)
            rows.append(rec)
            for coin in (0,1):
                moved,(_,current)=native(cfg,(a,b),coin)
                want=0 if rec==0 else 1 if rec==2 else -1
                native_ok &=moved==outcfg and current==want
                calls.append(dict(cell=cell,address=k,seam=[a,b],coin=coin,
                                  configuration=list(moved),record=rec,current=current))
        records.append(rows)
    fx=np.array([r==2 for r in records[0]],int)
    fy=np.array([r in (2,3) for r in records[1]],int)
    pairs=[dict(x=k,y=l,record_x=records[0][k],record_y=records[1][l],
                xor=int(fx[k]^fy[l]),and_flag=int(fx[k]*fy[l]))
           for k in range(60) for l in range(60)]
    checks['actual_240_pinned_calls']=bool(native_ok and len(calls)==240)
    checks['actual_120_completion_columns_and_inverse']=bool(exact)
    checks['all_3600_original_labels_retained']=len(pairs)==3600
    checks['native_count_is_2_over_9']=sum(r['xor'] for r in pairs)==800
    checks['native_AND_control_is_1_over_72']=sum(r['and_flag'] for r in pairs)==50
    # All 60 columns of each restricted one-cell projector, with actual
    # coherent source step/inverse, not merely uniform-state probabilities.
    R=np.full((60,60),F(1,30),object)-np.eye(60,dtype=int)
    for name,flags in zip(('X','Y'),(fx,fy)):
        wanted=R@np.diag(list(map(int,flags)))@R
        actual=np.empty((60,60),object);actual.fill(F())
        outside=False
        for k in range(60):
            ticked=w.step({(cfg,0,k):F(1)})
            chosen=wk.restrict(ticked,lambda key: key[1]==2 if name=='X' else key[1] in (2,3))
            back=w.inverse_step(chosen)
            for (config,rec,j),a in back.items():
                outside |= config!=cfg or rec!=0
                if config==cfg and rec==0:actual[j,k]=a
        checks['full_source_R_D_R_intertwiner_'+name]=not outside and np.array_equal(actual,wanted)
        checks['wrong_diagonal_only_projector_detected_'+name]=not np.array_equal(wanted,np.diag(flags))
    data=dict(port=2,addresses=[list(x) for x in w.addresses],records=records,
              calls=calls,pairs=pairs,native_calls=240,completion_columns=120,
              compiled_restricted_count_oracle=True,source_fibre_dimension=3600,
              fresh_full_count_repair_word_executions=0,
              compression='all count repair ticks replaced by proved R D R fibre intertwiners; 100 percent restricted-oracle compilation')
    return fx,fy,data,checks

def chart(flags):return tuple(k+64*int(flags[k]) for k in range(60))

@lru_cache(None)
def prep_word(codes,padweight):
    weights=[F()]*512
    for k in codes:weights[k]=(1-padweight)/60
    weights[511]=padweight
    return synth._state_word(weights)

def ray_word(codes,padweight):
    p=prep_word(codes,padweight);phases=[0]*512;phases[0]=iv.pi_grid()
    return c.inverse(p)+synth._phase_word(phases)+p

def reflection_word(codes):
    return ray_word(codes,F())+[c.pulse('I',math.pi,[0])]

def phase_word(and_control=False):
    beta=iv.split_angle(F(9,25),F(16,25))
    phases=[]
    for i in range(16):
        o=i&1;r=(i>>1)&1;x=(i>>2)&1;y=(i>>3)&1
        q=x*y if and_control else x^y
        phases.append(float(F(beta,iv.AQ))*r*(1-2*o)*q)
    return c.compile_diagonal(phases),beta

def compress_basis(basis,core,axis,label,other_basis):
    """Orthonormalize and report every discarded floating contraction norm."""
    u,s,vh=np.linalg.svd(basis,full_matrices=False)
    keep=max(1,int(sum(s>SVD_CUTOFF)))
    transfer=s[:,None]*vh
    moved=np.moveaxis(core,axis,0)
    discarded=np.tensordot(transfer[keep:],moved,axes=(1,0))
    loss=float(np.linalg.norm(discarded))
    represented_factor=represented_norm_bound(u[:,keep:])*represented_norm_bound(other_basis)
    bound=represented_factor*exact_l1(transfer[keep:])*frobenius_upper(moved)
    residual=float(np.linalg.norm(basis-u@transfer)*np.linalg.norm(moved)*float(represented_norm_bound(other_basis)))
    TRUNCATIONS.append(dict(label=label,before=basis.shape[1],after=keep,discarded_amplitude_norm=loss,
        discarded_amplitude_bound=str(bound),represented_basis_factor=str(represented_factor),factorization_residual=residual))
    new=np.tensordot(transfer[:keep],moved,axes=(1,0))
    return u[:,:keep],np.moveaxis(new,0,axis)

class State:
    def __init__(self,x,y,core):self.x=x;self.y=y;self.core=core
    def copy(self):return State(self.x.copy(),self.y.copy(),self.core.copy())
    def orthogonalize(self,axis,label):
        basis=self.x if axis==0 else self.y
        basis,self.core=compress_basis(basis,self.core,axis,label,self.y if axis==0 else self.x)
        if axis==0:self.x=basis
        else:self.y=basis
        return self
    def effect(self):
        t=self.core.reshape(self.core.shape[:2]+(-1,3))
        return np.einsum('ai,bj,abkc,ijkd->cd',self.x.conj().T@self.x,
                         self.y.conj().T@self.y,t.conj(),t,optimize=True)
    def store(self,label):
        ARRAYS[label+'__x']=self.x;ARRAYS[label+'__y']=self.y;ARRAYS[label+'__core']=self.core
        EFFECTS[label]=self.effect()
    def count_project(self,axis,addresses,label):
        out=self.copy();b=out.x if axis==0 else out.y
        masked=np.zeros_like(b);masked[list(addresses)]=b[list(addresses)]
        if axis==0:out.x=masked
        else:out.y=masked
        return out.orthogonalize(axis,label)
    def count_word(self,axis,word,codes,label,destructive=False):
        out=self.copy();b=out.x if axis==0 else out.y
        b=c.execute(word,9,b)
        if axis==0:out.x=b
        else:out.y=b
        good=out.count_project(axis,codes,label+'_code')
        bad=out.count_project(axis,tuple(sorted(set(range(512))-set(codes))),label+'_padding')
        if destructive:
            good.store(label+'_nonreturn');return bad
        bad.store(label+'_padding');return good
    def family_columns(self):
        rx,ry=self.core.shape[:2]
        return self.core.transpose(2,3,0,1,4).reshape(512,4,rx*ry*3)
    def set_family(self,a):
        rx,ry=self.core.shape[:2]
        self.core=a.reshape(a.shape[:2]+(rx,ry,3)).transpose(2,3,0,1,4)
        return self
    def family_word(self,word,label,destructive=False):
        out=self.copy();good,bad=m.word_export(word,self.family_columns())
        terminal=self.copy().set_family(good if destructive else bad)
        terminal.store(label+'_nonreturn' if destructive else label)
        out.set_family(bad if destructive else good)
        return out
    def inverse_tick(self,base):
        out=self.copy();value,fail=m.apply_tick(self.family_columns(),base['object_inverse'],
                    base['inverse_shift_clean'],base['inverse_shift_dirty'],True)
        rx,ry=self.core.shape[:2]
        terminal=self.copy();terminal.core=fail[0].reshape(-1,rx,ry,3).transpose(1,2,0,3)
        terminal.store('inverse_family_scratch')
        out.set_family(value);return out

def flag_transition(basis,outbit,inbit):
    result=np.zeros_like(basis)
    source=[k for k in range(512) if ((k>>6)&1)==inbit]
    target=[(k&~64)+64*outbit for k in source]
    result[target]=basis[source]
    return result

def coherent_phase(state,matrix,label,codesx,codesy):
    """Execute the full actual 16x16 source block, including offdiagonals."""
    expanded=[];transfers=[]
    for axis,basis in enumerate((state.x,state.y)):
        parts=[flag_transition(basis,o,i) for o in (0,1) for i in (0,1)]
        stack=np.concatenate(parts,axis=1)
        u,s,vh=np.linalg.svd(stack,full_matrices=False)
        # No truncation in this operator-basis expansion. Later state-specific
        # chart exports carry the separately recorded floating truncations.
        b=u
        reconstruction=float(np.linalg.norm(stack-b@(b.conj().T@stack)))
        factor=4*np.linalg.norm(state.core)*np.linalg.norm(matrix)*float(represented_norm_bound(state.y if axis==0 else state.x))
        TRUNCATIONS.append(dict(label=label+'_operator_basis_'+str(axis),before=stack.shape[1],after=b.shape[1],
            discarded_amplitude_norm=0.,discarded_amplitude_bound='0',factorization_residual=reconstruction*factor,
            raw_operator_reconstruction_residual=reconstruction))
        expanded.append(b)
        transfers.append(np.array([b.conj().T@p for p in parts]).reshape(2,2,b.shape[1],basis.shape[1]))
    bx,by=expanded;tx,ty=transfers
    old=state.core.reshape(state.core.shape[:2]+(128,4,4,3))
    new=np.zeros((bx.shape[1],by.shape[1],128,4,4,3),complex)
    for xo in (0,1):
        for xi in (0,1):
            for yo in (0,1):
                for yi in (0,1):
                    block=matrix[4*xo+8*yo:4*xo+8*yo+4,4*xi+8*yi:4*xi+8*yi+4]
                    new +=np.einsum('ak,bl,uv,klhvfc->abhufc',tx[xo,xi],ty[yo,yi],block,old,optimize=True)
    out=State(bx,by,new.reshape(bx.shape[1],by.shape[1],512,4,3))
    # Explicit full chart exports after the coupling: first X, then Y.
    bad=out.count_project(0,tuple(sorted(set(range(512))-set(codesx))),label+'_xpad');bad.store(label+'_xpad')
    out=out.count_project(0,codesx,label+'_xcode')
    bad=out.count_project(1,tuple(sorted(set(range(512))-set(codesy))),label+'_ypad');bad.store(label+'_ypad')
    out=out.count_project(1,codesy,label+'_ycode')
    # Compress the carried state, not merely the basis. The unused directions
    # from the complete operator expansion otherwise inflate reference width.
    for axis in (0,1):
        moved=np.moveaxis(out.core,axis,0)
        u,s,vh=np.linalg.svd(moved.reshape(moved.shape[0],-1),full_matrices=False)
        keep=max(1,int(sum(s>SVD_CUTOFF)))
        loss=float(np.linalg.norm(s[keep:]))
        represented_factor=represented_norm_bound(out.x)*represented_norm_bound(out.y)*represented_norm_bound(u[:,keep:])
        bound=represented_factor*exact_l1(s[keep:])*frobenius_upper(vh[keep:])
        residual=float(np.linalg.norm(moved.reshape(moved.shape[0],-1)-u@(s[:,None]*vh))*
                       float(represented_norm_bound(out.x)*represented_norm_bound(out.y)))
        TRUNCATIONS.append(dict(label=label+'_state_rank'+str(axis),before=len(s),after=keep,discarded_amplitude_norm=loss,
            discarded_amplitude_bound=str(bound),represented_basis_factor=str(represented_factor),factorization_residual=residual))
        basis=out.x if axis==0 else out.y
        basis=basis@u[:,:keep]
        newcore=(s[:keep,None]*vh[:keep]).reshape((keep,)+moved.shape[1:])
        out.core=np.moveaxis(newcore,0,axis)
        if axis==0:out.x=basis
        else:out.y=basis
    return out

def exact_effect(z):
    """Independent full432 native family contraction, with literal phase J."""
    p=m.embedding(m.rays()['prep'])[list(m.CODES)]
    phi=m.embedding(m.rays()['phi'])[list(m.CODES)]
    raw=(m.reg.full_tick()/4)@np.kron(p[:,None],np.eye(3))
    ticked=raw.reshape(144,3,3).copy()
    for k in range(144):ticked[k]*=z**m.reg.step_of(k//12)
    M=np.einsum('k,kij->ij',phi.conj(),ticked)
    return M.conj().T@M

def source_bindings():
    result={}
    for mod in list(sys.modules.values()):
        path=getattr(mod,'__file__',None)
        if path:
            p=Path(path).resolve()
            if p.suffix=='.py' and p.is_relative_to(m.SOURCE):result[str(p.relative_to(ROOT))]=sha(p.read_bytes())
    for name in ('implementation.py','SOURCE_ARRAYS.npz','RECEIPT.json','MANIFEST.sha256','PROOF.md'):
        p=FAMILY/name;result[str(p.relative_to(ROOT))]=sha(p.read_bytes())
    return dict(sorted(result.items()))

def calculate():
    fx,fy,census,checks=native_census()
    (HERE/'NATIVE_CENSUS.json').write_text(json.dumps(census,sort_keys=True,indent=2)+'\n')
    cx,cy=chart(fx),chart(fy)
    base=np.load(FAMILY/'SOURCE_ARRAYS.npz')
    oldreceipt=json.loads((FAMILY/'RECEIPT.json').read_text())
    checks['sealed_family_manifest']=sha((FAMILY/'MANIFEST.sha256').read_bytes())=='4714358385d7e0655197fe8b0653e4495de9df752ae4dc083b1d019c8c792fef'
    blank=np.zeros(512,complex);blank[0]=1
    gx=[];prep_mass=1.
    for name,codes in (('X',cx),('Y',cy)):
        word=remember_word('prepare_'+name,prep_word(codes,F()),9)
        a=c.execute(word,9,blank);target=np.zeros(512);target[list(codes)]=1/math.sqrt(60)
        checks['actual_uniform60_preparation_'+name]=np.linalg.norm(a-target)<GATE
        ARRAYS['prepared_'+name]=a
        p=np.zeros_like(a);p[list(codes)]=a[list(codes)]
        bad=a-p
        EFFECTS['prepare_'+name+'_padding']=prep_mass*float(np.vdot(bad,bad).real)*np.eye(3)
        prep_mass*=float(np.vdot(p,p).real);gx.append(p)
    # These actual family native source arrays are reused, with exact bindings.
    # They include the original preparation and full forward tick, not ideal
    # substitute matrices. The new count-dependent stages below are fresh.
    EFFECTS['family_preparation_padding']=prep_mass*float(np.vdot(base['preparation_padding'],base['preparation_padding']).real)*np.eye(3)
    EFFECTS['forward_family_scratch']=prep_mass*m.effect(base['forward_scratch_failure'])
    core=base['family_after_tick'][None,None,:,:,:].copy()
    state=State(gx[0][:,None],gx[1][:,None],core)
    state.orthogonalize(0,'initial_x').orthogonalize(1,'initial_y')
    BLOCKS.extend([dict(name='shared_family_preparation',**m.pulse_bound(oldreceipt['preparation_pulses'],9)),
                   dict(name='shared_family_forward_tick',**m.pulse_bound(2*oldreceipt['object_word_pulses']+oldreceipt['controlled_shift']['native_pulses'],2))])
    rwx=reflection_word(cx);rwy=reflection_word(cy)
    # Both sides of B† phase B are executed, even though R60 fixes the ideal
    # initial uniform vector. The source word phases are retained.
    for axis,(name,codes,word) in enumerate((('X',cx,rwx),('Y',cy,rwy))):
        state=state.count_word(axis,remember_word('forward_Rin_'+name,word,20),codes,'forward_Rin_'+name)
    pw,beta=phase_word();pi=c.inverse(pw)
    remember_word('forward_count_phase',pw,1)
    phase=c.execute(pw,4);inverse=c.execute(pi,4)
    ARRAYS['count_phase']=phase;ARRAYS['inverse_count_phase']=inverse
    target=np.diag([complex(3/5,4/5)**(((i>>1)&1)*(1-2*(i&1))*(((i>>2)&1)^((i>>3)&1))) for i in range(16)])
    checks['actual_full16_count_phase']=np.linalg.norm(phase-target)<GATE
    checks['actual_separately_reversed_phase']=np.linalg.norm(inverse@phase-np.eye(16))<GATE
    state=coherent_phase(state,phase,'forward_phase',cx,cy)
    for axis,(name,codes,word) in enumerate((('X',cx,rwx),('Y',cy,rwy))):
        state=state.count_word(axis,remember_word('forward_Rout_'+name,word,20),codes,'forward_Rout_'+name)
    weak,comp=m.filter_words('phi')
    state=state.family_word(remember_word('weak_Phi',weak,20),'first_click')
    state=state.family_word(remember_word('compensate_Phi',comp,19),'correction_padding')
    # Inverse uses the same unmeasured count archive and negative source words.
    for axis,(name,codes,word) in enumerate((('X',cx,rwx),('Y',cy,rwy))):
        state=state.count_word(axis,remember_word('inverse_Rin_'+name,c.inverse(word),20),codes,'inverse_Rin_'+name)
    remember_word('inverse_count_phase',pi,1)
    state=coherent_phase(state,inverse,'inverse_phase',cx,cy)
    for axis,(name,codes,word) in enumerate((('X',cx,rwx),('Y',cy,rwy))):
        state=state.count_word(axis,remember_word('inverse_Rout_'+name,c.inverse(word),20),codes,'inverse_Rout_'+name)
    state=state.inverse_tick(base)
    BLOCKS.append(dict(name='shared_source_inverse_family_tick',**m.pulse_bound(2*oldreceipt['object_word_pulses']+oldreceipt['inverse_shift_native_pulses'],2)))
    for axis,(name,codes) in enumerate((('X',cx),('Y',cy))):
        state=state.count_word(axis,remember_word('return_'+name,ray_word(codes,F(1,2)),20),codes,'return_'+name,True)
    state=state.family_word(remember_word('return_object',m.return_word(),20),'return_object',True)
    state.store('returned')
    # Returned padding registers remain environment coordinates. Extract only
    # the ideal selected511,511,511 branch for its Kraus-map comparison.
    piobj=m.PADDING.index(511)
    k=np.einsum('a,b,abij->ij',state.x[511],state.y[511],state.core[:,:,piobj,:3,:])
    ARRAYS['selected_return_map']=k
    A0=exact_effect(1);Ab=exact_effect(complex(3/5,4/5));A=7/9*A0+2/9*Ab
    AND=71/72*A0+1/72*Ab
    ARRAYS['A0']=A0;ARRAYS['Abeta']=Ab;ARRAYS['A1']=A;ARRAYS['AND_control']=AND
    checks['actual_N1_selected_return_map']=np.linalg.norm(k+np.eye(3)-2/5*A)<GATE
    checks['actual_N1_full_return_effect']=np.linalg.norm(EFFECTS['returned']-(np.eye(3)-2/5*A)@(np.eye(3)-2/5*A))<GATE
    checks['actual_N1_click_effect']=np.linalg.norm(EFFECTS['first_click']-16/25*A)<GATE
    nonreturn=sum(EFFECTS[n] for n in ('return_X_nonreturn','return_Y_nonreturn','return_object_nonreturn'))
    checks['actual_full_count_and_object_nonreturn']=np.linalg.norm(nonreturn-4/25*(A-A@A))<GATE
    checks['changed_count_AND_detected']=np.linalg.norm(AND-A)>1e-3
    direct=sum(Ab if r['xor'] else A0 for r in census['pairs'])/3600
    checks['all_source_labels_contract_to_same_A1']=np.linalg.norm(direct-A)<GATE
    total=sum(EFFECTS.values())
    checks['all_retained_terminal_effects_sum']=np.linalg.norm(total-np.eye(3))<GATE
    trunc_bound=sum((F(x['discarded_amplitude_bound']) for x in TRUNCATIONS),F())
    checks['floating_truncation_below_separate_gate']=trunc_bound<F(1,10**9)
    checks['factorization_reconstruction_residual_small']=sum(x['factorization_residual'] for x in TRUNCATIONS)<1e-9
    # New unique requests, including reversed and zero pulses, receive the
    # owning interval certificate. Shared fixed family blocks retain their
    # original certificate, bindings and primitive failure contracts.
    unique={(e['kind'],e['angle']) for _,w in ALL_WORDS for e in w}
    for kind,angle in sorted(unique):synth.certify_primitive(kind,angle)
    checks['all_new_unique_primitive_intervals']=True
    delta=F();gain=F(1)
    for block in BLOCKS:
        delta+=gain*F(block['complete_halfdiamond']);gain*=F(block['gain'])
    checks['complete_finite_source_bound_below_1e4']=delta<F(1,10000)
    queried=[]
    for name,word in ALL_WORDS:
        if name not in ('prepare_X','return_Y'):continue
        for pulse in (0,len(word)-1):
            d=1<<len(word[pulse]['targets'])
            a=m.primitive_failures(word,pulse,d,blank)
            queried.append(dict(word=name,pulse=pulse,complement=d,mass=float(np.linalg.norm(a)**2)))
    checks['lazy_count_processor_failure_queries']=len(queried)==4
    for name,value in EFFECTS.items():ARRAYS['effect__'+name]=value
    result=dict(task='LOCAL-M10-R11-COUNT-RETURN-20261114',candidate='lepton-encoded-count-filter',category='F',
        checks={k:bool(v) for k,v in checks.items()},all_checks_pass=bool(all(checks.values())),
        N=1,beta_witness='exp(i beta)=(3+4i)/5',beta_grid_integer=beta,angle_denominator=iv.AQ,
        count_charts={'X':list(cx),'Y':list(cy),'padding_count_each':452,'destructive_padding':511},
        actual_pinned_native_calls=240,actual_coherent_completion_columns=120,
        full_one_cell_projector_columns_checked=120,original_address_pair_records=3600,
        compiled_restricted_count_oracle=True,count_repair_compilation_share=1,
        fresh_count_repair_pulse_executions=0,shared_family_pulse_blocks=['preparation','forward_tick','inverse_tick'],
        fresh_source_words=[dict(name=n,pulses=len(w),sha256=sha(json.dumps(w,sort_keys=True).encode())) for n,w in ALL_WORDS],
        main_branch_pulses=sum(b['pulses'] for b in BLOCKS),unique_fresh_primitives=len(unique),
        complete_error_blocks=BLOCKS,complete_halfdiamond_bound=str(delta),complete_gain_bound=str(gain),
        floating_truncations=TRUNCATIONS,floating_truncation_amplitude_sum=sum(x['discarded_amplitude_norm'] for x in TRUNCATIONS),
        floating_truncation_amplitude_bound=str(trunc_bound),
        diagnostic_gate=GATE,terminal_effect_sum_error=float(np.linalg.norm(total-np.eye(3))),
        returned_map_error=float(np.linalg.norm(k+np.eye(3)-2/5*A)),
        queried_processor_first_failures=queried,source_hashes=source_bindings(),
        census_sha256=sha((HERE/'NATIVE_CENSUS.json').read_bytes()),
        output_contract='all source-code padding, first count nonreturn, object nonreturn, scratch and processor first failures retained; no success renormalization',
        supplied_services=['coherent native-fibre code/load/export','blank-double and record-zero source fibre','coherent proper-code load/gluing','fresh blank processors','analog phase control','complete trace semantics'],
        physical_row_established=False,fixed_precision_N_limit_claimed=False)
    return result

def main():
    started=time.perf_counter();result=calculate()
    (HERE/'RECEIPT.json').write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(HERE/'SOURCE_ARRAYS.npz',**ARRAYS)
    bindings=dict(archive_sha256=sha((HERE/'SOURCE_ARRAYS.npz').read_bytes()),
       receipt_sha256=sha((HERE/'RECEIPT.json').read_bytes()),implementation_sha256=sha(Path(__file__).read_bytes()),
       arrays={k:dict(shape=list(v.shape),dtype=str(v.dtype),raw_sha256=sha(v.tobytes())) for k,v in ARRAYS.items()})
    (HERE/'ARRAY_BINDINGS.json').write_text(json.dumps(bindings,sort_keys=True,indent=2)+'\n')
    summary=dict(status='pass' if result['all_checks_pass'] else 'fail',checks=len(result['checks']),
        failed=[k for k,v in result['checks'].items() if not v],seconds=time.perf_counter()-started,
        main_branch_pulses=result['main_branch_pulses'],halfdiamond=float(F(result['complete_halfdiamond_bound'])),
        returned_map_error=result['returned_map_error'],truncation_sum=result['floating_truncation_amplitude_sum'])
    (HERE/'RUN.json').write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n')
    print(json.dumps(summary,sort_keys=True,indent=2),flush=True)
    assert result['all_checks_pass'],summary['failed']

if __name__=='__main__':main()
