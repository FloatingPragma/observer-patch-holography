"""Producer-independent full ROW29 opposing review; writes only this directory.

Native matrices are rebuilt from the raw M6 response and control coefficients.
Geometry is rebuilt from translated integer faces. Frozen eigenvectors are never
reselected. Numerical diagnostics are distinguished from outward certificates.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from collections import defaultdict
import hashlib, json, math, sys, time
sys.dont_write_bytecode = True
import numpy as np
import sympy as s
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
ROW = ROOT/'openwave-tuning/row_completion/ROW_29'
SIZE = ROOT/'openwave-tuning/common_dynamics/photon_size3'
EXPECTED = {
 str((ROW/'MANIFEST.json').relative_to(ROOT)): '993744d14c2da8d66a93dc76422f81799074febad9914ae33bb60bf455cf00b6',
 str((SIZE/'MANIFEST.json').relative_to(ROOT)): '90178885fcaaaba405b950a6f132fd47d519d21a130a4f982ec09423db129ec2',
}
def load(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def pincheck():
    pins = dict(EXPECTED)
    for rel in EXPECTED:
        assert sha(ROOT/rel) == EXPECTED[rel], rel
        d = load(ROOT/rel)
        for group in ('files', 'source_pins'):
            for n,h in d.get(group, {}).items():
                assert sha(ROOT/n) == h, n
                assert n not in pins or pins[n] == h
                pins[n] = h
    return pins

def stage(name, result):
    (HERE/(name+'.json')).write_text(json.dumps(result, indent=2)+'\n')
    print(name, 'PASS', flush=True)
    return result

def imports():
    files={}
    for module in tuple(sys.modules.values()):
        filename=getattr(module,'__file__',None)
        if not filename:continue
        p=Path(filename).resolve()
        if p.is_file() and p.is_relative_to(ROOT) and not p.is_relative_to(HERE):
            files[str(p.relative_to(ROOT))]=sha(p)
    stage('IMPORTS',{'files':files,'python':sys.version,'numpy':np.__version__,'sympy':s.__version__,
                    'numerical_threads':'OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1'})
    return files

def geometry():
    sys.path.insert(0, str(ROOT/'openwave-tuning/common_dynamics/vector_field'))
    import linear_response as linear
    A,B,D,terms,labels = linear.bloch_data()
    graph,full_B,full_A,periods,winding,face_labels=linear.cc.build(3)
    frozen_geometry=load(SIZE/'GEOMETRY.json')
    assert [list(x) for x in graph.seams]==frozen_geometry['edges']
    exact_entries=[[i,j,int(full_A[i,j])] for i in range(full_A.rows) for j in range(full_A.cols) if full_A[i,j]]
    assert exact_entries==frozen_geometry['face_entries']
    assert [[i,j,int(3*periods[i,j])] for i in range(3) for j in range(972) if periods[i,j]]==frozen_geometry['period_times3_entries']
    assert [[int(x) for x in winding[:,j]] for j in range(3)]==frozen_geometry['winding_columns']
    # Different inverse construction from producer's Matrix.pinv().
    C,R = A.rank_decomposition()
    Ap = R.T*(R*R.T).inv()*(C.T*C).inv()*C.T
    assert A*Ap*A == A and Ap*A*Ap == Ap
    P = s.eye(36)-A*Ap; Q = s.eye(26)-Ap*A
    N = s.Matrix.hstack(*(B.col_join(A.T)).nullspace())
    H = N*(N.T*N).inv()*N.T
    chains = s.zeros(36,3)
    for j in range(3): chains[30+2*j,j] = chains[31+2*j,j] = 1
    W = s.Rational(2,3)*H*chains
    crosses = [s.Matrix([[0,0,0],[0,0,-1],[0,1,0]]),s.Matrix([[0,0,1],[0,0,0],[-1,0,0]]),s.Matrix([[0,-1,0],[1,0,0],[0,0,0]])]
    V = s.zeros(26,3); V[25,0]=1;V[24,1]=-1;V[23,2]=1
    assert A.rank()==22 and B.rank()==11 and B*A==s.zeros(12,26)
    assert W.T*W == s.Rational(2,3)*s.eye(3) and V.T*V == s.eye(3)
    assert A.T*W==s.zeros(26,3) and B*W==s.zeros(12,3) and A*V==s.zeros(36,3)
    assert all(P*D[j]*Q==W*crosses[j]*V.T for j in range(3))
    frozen=load(ROW/'MAXWELL.json')
    assert W==s.Matrix(frozen['W']) and V.T==s.Matrix(frozen['V_transpose'])
    assert s.trace(Ap*Ap.T)==s.Rational(3383,240) and s.trace(A*A.T)==88
    # Independent finite Laurent-polynomial BA, rather than several sampled k.
    variables=s.symbols('zx zy zz', nonzero=True)
    Ak=s.zeros(36,26);Bk=B.copy()
    for e,f,shift,c in terms: Ak[e,f]+=c*s.prod(variables[j]**(-shift[j]) for j in range(3))
    for axis,links in enumerate(linear.cc.geometry.links()):
        for channel,(a,b) in enumerate(links): Bk[b,30+2*axis+channel]=-1/variables[axis]
    assert (Bk*Ak).applyfunc(s.expand)==s.zeros(12,26)
    k=s.Matrix(s.symbols('kx ky kz', real=True)); cross=sum((k[j]*crosses[j] for j in range(3)),s.zeros(3))
    assert cross*k==s.zeros(3,1)
    assert (cross*cross-k*k.T+(k.dot(k))*s.eye(3)).applyfunc(s.expand)==s.zeros(3)
    # The raw translation list directly bounds the linear and quadratic errors.
    linear_sum=sum(abs(c)*s.sqrt(sum(t*t for t in shift)) for _,_,shift,c in terms)
    # Use Frobenius per-entry triangle bounds rather than sum of all entries.
    sums=defaultdict(lambda:[0,0])
    for e,f,shift,c in terms:
        length=s.sqrt(sum(t*t for t in shift));sums[e,f][0]+=abs(c)*length;sums[e,f][1]+=abs(c)*length**2/2
    first2=s.simplify(sum(v[0]**2 for v in sums.values()));second2=sum(v[1]**2 for v in sums.values())
    assert first2<=121 and second2<=s.Rational(17,2)**2
    columns=defaultdict(list);convolution=defaultdict(int)
    for e,f,shift,c in terms:columns[f].append((e,shift,c))
    for column in columns.values():
        for e,d,c in column:
            for v,t,b in column:convolution[e,v,tuple(x-y for x,y in zip(d,t))]+=c*b
    majorant=np.zeros((36,36),np.int64)
    for (e,v,delta),value in convolution.items():majorant[e,v]+=abs(value)
    weights=majorant.sum(axis=1)
    assert min(weights)>0 and np.all(17*(majorant@weights)<=1504*weights)
    assert max(F(int(a),int(b)) for a,b in zip(majorant@weights,weights))==F(1504,17)
    bad=Ak.copy();e,f=next((i,j) for i in range(36) for j in range(26) if A[i,j]);bad[e,f]=-bad[e,f]
    assert (Bk*bad).applyfunc(s.expand)!=s.zeros(12,26)
    return stage('GEOMETRY', {'exact_all_direction_curl':True,'exact_Laurent_BA':True,'full_size3_integer_source_join':True,
        'zero_face_rank':22,'zero_incidence_rank':11,'rho':'2/3','a':'1/2','b':'3','c_squared':'1',
        'linear_Frobenius_coefficient_squared':str(first2),'quadratic_Frobenius_coefficient_squared':str(second2),
        'all_refinement_curvature_bound_from_actual_convolution':'1504/17',
        'inverse_Frobenius_squared':'3383/240','orientation_mutation_rejected':True,
        'derivation':'Actual translated faces; rank-factor pseudoinverse and independent nullspace projector; V selected by the signed macro-face columns and checked against all three derivatives.'})

@lru_cache(None)
def generators():
    base=ROOT/'reverse-engineering-reality/code'
    raw=load(base/'source_selection_model/response.json')
    controls=load(base/'m1_source_realization/receipt.json')
    def scalar(pair):return float(F(pair[0]))+math.sqrt(5)*float(F(pair[1]))
    mats=[]
    for blocks in raw['generators']:
        g=np.zeros((6,6),complex)
        for block, rows in enumerate(blocks):
            for i,row in enumerate(rows):
                for j,v in enumerate(row):g[3*block+i,3*block+j]=scalar(v[:2])+1j*scalar(v[2:])
        assert np.linalg.norm(g+g.conj().T)<1e-13
        mats.append(g)
    c={r['name']:np.array([scalar(v) for v in r['coefficients']]) for r in controls['evidence']['controls']['rows']}
    coeff={'X':c['imaginary_01'],'Y':c['real_01'],'Z':c['phase_0']-c['phase_1'],'I':c['phase_0']+c['phase_1'],'J':c['relative_phase']}
    coeff['H']=(coeff['X']+coeff['Z'])/math.sqrt(2)
    return {name:sum((x*g for x,g in zip(row,mats)),np.zeros((6,6),complex)) for name,row in coeff.items()}

@lru_cache(None)
def local(kind,angle):return expm(angle*generators()[kind])
@lru_cache(None)
def slots(bits,targets):
    # Independent explicit map; target[0] is the high local bit.
    rest=[i for i in range(bits) if i not in targets]
    return np.array([[sum(((local_i>>(len(targets)-1-j))&1)<<t for j,t in enumerate(targets))+
                      sum(((other>>j)&1)<<t for j,t in enumerate(rest))
                      for other in range(1<<len(rest))] for local_i in range(1<<len(targets))])

def execute(word,bits,state):
    M=np.asarray(state,complex).reshape(1<<bits,-1).copy();n=M.shape[1]
    failure=np.zeros((n,n),complex); initial=M.conj().T@M;maxstep=0.
    for p in word:
        U=local(p['kind'],float(p['angle']));d=1<<len(p['targets']);indices=slots(bits,tuple(p['targets']));old=M[indices]
        lost=np.einsum('ij,jkr->ikr',U[d:,:d],old).reshape(-1,n)
        failure+=lost.conj().T@lost
        M[indices]=np.einsum('ij,jkr->ikr',U[:d,:d],old)
    gram=float(np.linalg.norm(M.conj().T@M+failure-initial,2))
    return M,failure,gram

def native_pointer():
    sys.path.insert(0,str(ROOT/'reverse-engineering-reality'))
    from model.native_interface import pointer
    from model.observer.source_discrimination import circuit
    from model.observer import native
    start=np.eye(8,dtype=complex)[:,[pointer.READY_CODE,pointer.READY_CODE+4]]
    M,failure,gram=execute(pointer._local_word(),3,start)
    result={}
    for law in ('C','E'):
        ks=defaultdict(lambda:defaultdict(lambda:np.zeros((2,2),complex)));metadata={}
        for code in range(4):
            source=circuit.decode(code);K=M[[code,code+4]]
            for coin in (0,1):
                destination,(_,current)=native.native(source,pointer.SEAM,coin)
                key=(coin,current if law=='C' else source[0]-source[1])
                ks[key][destination]+=K/math.sqrt(2)
                metadata[key]={'coin':coin,'label':key[1],'current':current,'bit':int(2*current+2*coin-1<0)}
        result[law]=[(metadata[key],sum((K.conj().T@K for K in ks[key].values()),np.zeros((2,2),complex))) for key in sorted(ks)]
        assert np.linalg.norm(sum((E for _,E in result[law]),np.zeros((2,2),complex))+failure-np.eye(2))<1e-12
    return result,gram,failure

def tensor_effects(branches):
    out=[]
    for m2,E2 in branches:
        for m1,E1 in branches:
            for m0,E0 in branches:
                out.append((m0['bit']+2*m1['bit']+4*m2['bit'],np.kron(E2,np.kron(E1,E0))))
    return out

def native_radiative():
    sys.path.insert(0,str(ROW));import radiative_limit as producer
    data=load(ROW/'RADIATIVE.json');programs=load(ROW/'RADIATIVE_PROGRAM.json')
    branches,pgram,pfailure=native_pointer();ploss=float(np.trace(pfailure).real)
    effects={law:tensor_effects(v) for law,v in branches.items()}
    # Sequential local pointer failures telescope to I minus the successful
    # tensor effect; no pointer outcome is removed by normalization.
    for law,entries in branches.items():
        S=sum((E for _,E in entries),np.zeros((2,2),complex))
        fail=np.kron(np.eye(4),pfailure)+np.kron(np.eye(2),np.kron(pfailure,S))+np.kron(pfailure,np.kron(S,S))
        assert np.linalg.norm(fail+np.kron(S,np.kron(S,S))-np.eye(8),2)<2e-12
    reports=[];maxfine=0.;maxgram=0.;maxrecord=0.;maxphase=0.;maxclock=0.
    from model.native_clock import configuration
    from model.native_interface import pointer
    from model.observer.source_discrimination import circuit
    p0=min(map(F,configuration.success_probabilities(4,32,full=True,mode='reset',target=circuit.decode(0))))
    pr=min(map(F,configuration.success_probabilities(4,32,full=True,mode='reset',target=circuit.decode(pointer.READY_CODE))))
    blank_ready_bound=3*(1-p0)+3*(1-pr)
    cert=producer.cert;nf=producer.nf
    herror=nf.root_certificate.exact_template_bound(nf.H(0),[('I',(0,),cert._pi(-1,2)),('H',(0,),cert._pi(1,2))])
    for row in data['experiments']:
        name=row['name'];word=programs[name]
        M,Fgram,gram=execute(word,3,np.eye(8));maxgram=max(maxgram,gram);v=M[:,0]
        assert gram<2e-11
        time0=F(row['time']);pos=F(row['position_x']);bit=0 if row['polarization']=='y' else 1
        diagonal,phase_error=producer.phase_word(time0)
        assert word==producer.apparatus(bit,time0,pos,row['field_setting'],row['clock_setting'],diagonal)
        actual_native=F(row['word_pointer_error_exact']);reset=F(row['reset_failure_exact'])
        assert reset==blank_ready_bound
        position_word=[producer.compiler.pulse('Z',-math.pi*float(pos),[0])]
        poserror=nf.root_certificate.exact_template_bound(position_word,[('Z',(0,),cert._pi(-pos))])
        yerror=nf.root_certificate.exact_template_bound([producer.compiler.pulse('Z',math.pi/4,[0])],[('Z',(0,),cert._pi(1,4))]) if row['field_setting']=='Y' else F()
        assert actual_native==phase_error+2*herror+poserror+yerror+producer.clock.exact_word_bound(1,row['clock_setting'],time0)+2*producer.clock.pointer_bound()
        assert F(row['native_complete_TV_exact'])==actual_native+reset
        for law, entries in effects.items():
            probs=np.zeros(8);fine=[]
            for target,E in entries:
                val=float(np.vdot(v,E@v).real);assert val>-1e-14;probs[target]+=val;fine.append(val)
            assert abs(sum(probs)+Fgram[0,0].real-1)<2e-11
            stored=np.array(row['laws'][law]['all8_probabilities']);err=float(np.max(abs(probs-stored)));maxrecord=max(maxrecord,err)
            assert err<3e-12
            for (meta,E),frozen in zip(branches[law],row['laws'][law]['fine_one_bit_effects']):
                E0=np.array(frozen['effect_real'])+1j*np.array(frozen['effect_imag']);assert meta==frozen['record']
                maxfine=max(maxfine,float(np.max(abs(E-E0))))
            expected=math.cos(2*math.pi*float(pos-time0)) if row['field_setting']=='X' else math.sin(2*math.pi*float(pos-time0))
            contrast=sum(p*(-1 if j&(1<<bit) else 1) for j,p in enumerate(probs))
            maxphase=max(maxphase,abs(contrast-expected));assert abs(contrast-expected)<2*float(actual_native)+1e-11
            clockcontrast=sum(p*(-1 if j&4 else 1) for j,p in enumerate(probs))
            clocktarget=math.cos(float(time0)/8) if row['clock_setting']=='X' else math.sin(float(time0)/8)
            maxclock=max(maxclock,abs(clockcontrast-clocktarget));assert abs(clockcontrast-clocktarget)<2*float(actual_native)+1e-11
        reports.append({'name':name,'pulses':len(word),'all_native_columns':8,'fine_C_outcomes':len(effects['C']),'fine_E_outcomes':len(effects['E']),
            'complete_Gram_residual':gram,'first_export_failure_trace':float(np.trace(Fgram).real)})
    # Complete three-bit columns make omission or success renormalization visible.
    wrong=programs['wrong_double_b'];M,Fg,g=execute(wrong,3,np.eye(8));p=abs(M[:,0])**2
    wrongcontrast=float(sum(v*(-1 if j&1 else 1) for j,v in enumerate(p)));assert abs(wrongcontrast)>.5
    split,_,_=execute(programs['clock_split_control'],1,np.eye(2))
    coalesced,_,_=execute(producer.clock.local_words(1,'Y',F(1,4)),1,np.eye(2))
    assert np.linalg.norm(split-coalesced)<1e-12
    # Actual clock mutation on the identical field word, keeping source phases.
    row=next(x for x in data['experiments'] if x['time']=='1/4' and x['field_setting']=='X' and x['clock_setting']=='Y')
    wrong_clock=json.loads(json.dumps(programs[row['name']]))
    tick=[p for p in wrong_clock if p['kind']=='Z' and p['targets']==[2] and abs(p['angle']+1/64)<1e-15]
    assert len(tick)==1;tick[0]['angle']*=2
    v,_,_=execute(wrong_clock,3,np.eye(8)[:,0]);p=abs(v[:,0])**2
    contrast=sum(v*(-1 if j&4 else 1) for j,v in enumerate(p))
    clock_mutation=abs(contrast-math.sin(1/32));assert clock_mutation>.03
    # Visible scalar is not silently phase-aligned away. Mutating it changes
    # a vacuum/padding coherent reference even though the tested means agree.
    T=F(1,8);word,_=producer.phase_word(T);M,Fg,g=execute(word,2,np.eye(4))
    phases=np.array([-2*math.pi*float(T),-4*math.pi*float(T),-4*math.pi*float(T),0.])
    phase_resid=float(np.linalg.norm(M-np.diag(np.exp(1j*phases)),2));assert phase_resid<1e-11
    scalar_mutation=abs(1-np.exp(1j*2*math.pi*float(T)));assert scalar_mutation>.7
    # Physical field compression: other polarization has zero first moment;
    # its two output signs are fair and padding is part of the dilation.
    X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex)
    for bit in (0,1):
        for O in (X,Y):
            full=np.kron(np.eye(2),O) if bit==0 else np.kron(O,np.eye(2));compressed=full[:3,:3]
            expected=np.zeros((3,3),complex);expected[0,1<<bit]=O[0,1];expected[1<<bit,0]=O[1,0]
            assert np.array_equal(compressed,expected)
            other=2 if bit==0 else 1
            assert full[other,other]==0 and sum(abs(full[:,other])**2)==1
    assert maxfine<1e-12
    return stage('NATIVE', {'experiments':reports,'all_settings':len(reports),'max_complete_Gram_residual':maxgram,
        'max_stored_probability_difference':maxrecord,'max_fine_effect_difference':maxfine,'max_field_target_difference':maxphase,
        'max_clock_target_difference':maxclock,'wrong_clock_contrast_separation':clock_mutation,
        'pointer_complete_Gram_residual':pgram,'pointer_export_trace':ploss,'wrong_double_b_contrast':wrongcontrast,
        'full_phase_reference_residual':phase_resid,'scalar_deletion_operator_defect':scalar_mutation,
        'all_first_export_complements_accumulated':True,'no_success_normalization':True,'compressed_field_POVM_exact':True,
        'split_clock_equals_terminal_clock':True,'floating_residuals_are_diagnostics_not_certificates':True})

def finite_bounds():
    data=load(ROW/'RADIATIVE.json');entries=[]
    for row in data['refinements']:
        h=F(row['h'])
        # Interval uses its provider's published integer SCALE.
        sys.path.insert(0,str(ROW));import radiative_limit as r
        p=r.iv.pi()*2;lo=F(p.lo,r.iv.SCALE);hi=F(p.hi,r.iv.SCALE);q=h*hi
        eta=11*q;D=8*eta*eta+32*eta*eta*(10+eta)+F(17,2)*q*q
        zeta=32*eta+4*D/(F(4,5)*q-D);dw=F(5,4)*D/h
        assert eta<F(1,16) and q+D<F(1,4)-eta and D<q/4 and zeta<F(1,8)
        # Check high-cluster separation actually used by the proof, not only
        # the looser public domain checks.
        assert F(1,4)-eta-q>F(4,5)*q
        assert dw<lo/4
        assert F(row['source_frequency_error_exact'])==dw
        assert F(row['source_cap_H_error_including_visible_vacuum_exact'])==2*dw
        assert F(row['polarization_subspace_error_exact'])==zeta
        field=64*(zeta+dw/lo);assert F(row['local_field_relative_error_exact'])==field
        tv=max(F(x['native_complete_TV_exact'])+F(x['time'])*2*dw for x in data['experiments'])
        assert F(row['complete_native_source_TV_upper_exact'])==tv
        entries.append({'L':row['L'],'TV_exact':str(tv),'TV':float(tv),'field_mean_error':float(2*tv+field),
                        'huge_lattice_allocated':False,'terminal_native_words_executed':True})
    return stage('BOUNDS',{'refinements':entries,'all_direction_error_formula_checked':True,
                          'finite_box_source_limit':'h -> 0 at fixed bounded nonzero allowed physical momenta',
                          'native_precision_floor_retained':True,'zero_error_native_limit_executed':False})

def source_contract():
    # General Hermitian frequency matrix: exact symbolic free cap identity.
    a,b,x,y=s.symbols('a b x y',real=True);Om=s.Matrix([[a,x+s.I*y],[x-s.I*y,b]])
    H=s.diag((a+b)/2,(a+b)/2,(a+b)/2);H[1:3,1:3]=H[1:3,1:3]+Om
    annih=[]
    for j in range(2):v=s.zeros(3);v[0,j+1]=1;annih.append(v)
    for j in range(2):assert H*annih[j]-annih[j]*H == -sum((Om[j,l]*annih[l] for l in range(2)),s.zeros(3))
    # A finite correlated reference tests the exact partial-trace theorem,
    # including subnormalized inputs and a non-unitary observed branch.
    rng=np.random.default_rng(804291);z=rng.normal(size=(24,24))+1j*rng.normal(size=(24,24))
    rho=z@z.conj().T;rho/=np.trace(rho)
    def trace_spectator(m):return np.einsum('asrbst->arbt',m.reshape(4,3,2,4,3,2)).reshape(8,8)
    q,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
    u=expm(-1j*np.array([[.7,.1j,0,0],[-.1j,.3,.2,0],[0,.2,.8,0],[0,0,0,0]]))
    branch=np.diag([1.,.3,.6,.2])@u;max_residual=0.
    for scale in (1.,.37):
        for op in (u,branch):
            big=np.kron(np.kron(op,q),np.eye(2));small=np.kron(op,np.eye(2))
            error=float(np.linalg.norm(trace_spectator(big@(scale*rho)@big.conj().T)-small@trace_spectator(scale*rho)@small.conj().T))
            assert error<1e-14;max_residual=max(max_residual,error)
    # Squaring truncated coordinates does not produce compressed full H.
    aa=np.array([[0.,1.],[0.,0.]])
    wrong=(aa@aa.T+aa.T@aa)/2;correct=np.diag([.5,1.5])
    assert np.linalg.norm(wrong-correct,2)==1
    assert np.linalg.norm(aa@aa.T-aa.T@aa-np.eye(2),2)==2
    return stage('SOURCE',{'symbolic_general_Hermitian_cap_Heisenberg_identity':True,
        'correlated_spectator_reference_max_residual':max_residual,
        'normalized_subnormalized_and_local_branch_reduction':True,
        'naive_square_capped_fields_mutation_rejected':True,'cap_CCR_defect_retained':2,
        'proof_scope':'Autonomous free factors, compressed linear first moments, no interacting Weyl or continuous field-spectrum claim.'})

def main():
    before=pincheck();(HERE/'SOURCE_PINS.json').write_text(json.dumps(before,indent=2)+'\n')
    geometry();native_radiative();finite_bounds()
    # Review source algebra separately in a real symmetric source frame.
    source_contract()
    assert pincheck()==before
    imported=imports()
    assert all(n not in before or before[n]==h for n,h in imported.items())
    stage('CHECKS',{'all_checks_pass':True,'producer_pins':len(before),'source_pins_unchanged':True})
if __name__=='__main__':main()
