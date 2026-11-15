"""Independent opposing audit. No producer verifier or full pair semigroup."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import ast
import hashlib
import importlib.util
import json
import math
import sys

import numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CHAIN = HERE.parent
BRIDGE = ROOT / 'openwave-tuning/auditor/native_field_bridge'
PINS = {}
CLAIMS = []
DETAIL = {}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def pin(p, expected=None):
    p = Path(p)
    actual = sha(p)
    if expected is not None:
        assert actual == expected, (str(p), actual, expected)
    PINS[str(p.relative_to(ROOT))] = actual
    return actual


def claim(name, condition, evidence):
    CLAIMS.append({'claim': name, 'status': 'PASS' if condition else 'FAIL', 'evidence': evidence})
    if not condition:
        print('FAIL', name, evidence, flush=True)


def load(name, p):
    spec = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def err(a, b):
    return float(np.max(abs(a-b)))


def parity_transport(mask, src, dst):
    if not (mask >> src & 1) or (mask >> dst & 1):
        return None
    phase = (-1) ** ((mask & ((1 << src)-1)).bit_count())
    cleared = mask ^ (1 << src)
    phase *= (-1) ** ((cleared & ((1 << dst)-1)).bit_count())
    return cleared | (1 << dst), phase


def independent_clock_operators(states, law):
    ix = {s: i for i, s in enumerate(states)}
    groups = {}
    for col, (mask, record) in enumerate(states):
        for seam, (a, b) in enumerate(((0,1),(1,2),(2,0))):
            d = (mask >> a & 1) - (mask >> b & 1)
            for coin in (0, 1):
                current = (d-(2*coin-1)*(abs(d) % 2)) // 2
                moved = parity_transport(mask, a, b) if current == 1 else parity_transport(mask, b, a) if current == -1 else (mask, 1)
                r = list(record)
                r[2*seam] = (r[2*seam]+current) % 6
                r[2*seam+1] = (r[2*seam+1]-current) % 6
                target = ix[(moved[0], tuple(r))]
                for chosen in (('C','E') if law == 'MIX' else (law,)):
                    key = (seam, chosen, coin, current if chosen == 'C' else d)
                    if key not in groups:
                        groups[key] = np.zeros((384,384))
                    groups[key][target,col] = moved[1]/math.sqrt(12 if law == 'MIX' else 6)
    return groups


def audit_clock():
    contract = json.loads((CHAIN/'clock/CLOCK_CONTRACT.json').read_text())
    pin(CHAIN/'clock/common_clock.py', 'eee29493c0cc8a2bfdaf2868fdce82e10d055535f8ee0e6a0b4b818ef9fb38db')
    pin(CHAIN/'clock/CLOCK_CONTRACT.json')
    for rel, expected in contract['source_sha256'].items():
        pin(ROOT/rel, expected)
    for key in ('receipt', 'state_archive'):
        pin(CHAIN/'clock'/contract[key], contract[key+'_sha256'])
    mod = load('audited_common_clock', CHAIN/'clock/common_clock.py')
    clock = mod.CommonClock()
    states = clock.model.states
    ix = {s: i for i, s in enumerate(states)}
    H = np.zeros((384,384))
    delta = np.array([1,-1]*3)
    for col, (mask, rec) in enumerate(states):
        H[col,col] = .23*sum((1 if p<3 else -1)*(mask>>p&1) for p in range(6))
        H[col,col] += .37/2*sum(((v+3)%6-3)**2 for v in rec)+2*.61
        for sign in (-1,1):
            target = (mask, tuple((np.array(rec)+sign*delta)%6))
            H[ix[target],col] -= .61
        for e, ends in enumerate(((0,1),(1,2),(2,0))):
            for side, port in enumerate(ends):
                for src, dst, increment in ((port,3+e,1),(3+e,port,-1)):
                    moved = parity_transport(mask,src,dst)
                    if moved is None:
                        continue
                    r = list(rec); r[2*e+side] = (r[2*e+side]+increment)%6
                    H[ix[(moved[0],tuple(r))],col] += (1 if side==0 else -1)*moved[1]
    claim('clock.actual_H_reconstructed', err(H,clock.h)<2e-15, {'max_error':err(H,clock.h)})
    assert hashlib.sha256(np.asarray(clock.h,dtype='<f8').tobytes()).hexdigest() == contract['source_H_array_sha256_le_float64']
    numbers = np.array([mask.bit_count() for mask,_ in states])
    spec = {q: np.linalg.eigh(H[np.ix_(numbers==q,numbers==q)]) for q in range(7)}
    ground = min(float(v[0][0]) for v in spec.values())
    reference = np.zeros((384,384)); multiplicity = 0
    for q,(ev,vec) in spec.items():
        ids=np.flatnonzero(numbers==q); selected=abs(ev-ground)<1e-10
        reference[np.ix_(ids,ids)] = vec[:,selected]@vec[:,selected].T
        multiplicity += sum(selected)
    reference /= multiplicity
    reservoir = []
    for initial in (0,1):
        for final in (0,1):
            op = np.zeros((384,384))
            for col,(mask,rec) in enumerate(states):
                full = mask | (initial << 6)
                occupation = mask & 1
                if occupation == initial:
                    if final == initial:
                        op[col,col] = (1-2*initial)/math.sqrt(2)
                else:
                    moved = parity_transport(full,0,6) if occupation else parity_transport(full,6,0)
                    if (moved[0] >> 6) == final:
                        op[ix[(moved[0]&63,rec)],col] = moved[1]/math.sqrt(2)
            reservoir.append(op)
    rho = sum(k@reference@k.T for k in reservoir).astype(complex)
    actual,branches = clock.preparation()
    claim('clock.actual_charged_preparation', err(rho,actual)<2e-14 and err(sum(k.T@k for k in reservoir),np.eye(384))<1e-14,
          {'max_error':err(rho,actual),'ground_multiplicity':int(multiplicity),'branch_probabilities':[float(np.trace(k@reference@k.T)) for k in reservoir]})
    assert hashlib.sha256(np.asarray(actual,dtype='<c16').tobytes()).hexdigest() == contract['preparation_array_sha256_le_complex128']
    all_ops = {}
    for law in ('C','E','MIX'):
        ops = independent_clock_operators(states,law)
        actual_ops={(m['seam'],m['law'],m['coin'],m['record_value']):k/math.sqrt(3)
                    for seam in clock.model.field.operators_by_law[law] for m,k in seam}
        discrepancy=max(err(k,actual_ops[key]) for key,k in ops.items())
        complete=err(sum(k.T@k for k in ops.values()),np.eye(384))
        claim('clock.actual_native_K_'+law, discrepancy<2e-15 and complete<2e-15,
              {'marks':len(ops),'operator_error':discrepancy,'completeness_error':complete,
               'basis_coin_seam_events':384*6})
        all_ops[law]=ops
    ops=all_ops['MIX']
    operator_sha=hashlib.sha256(b''.join(np.asarray(ops[k],dtype='<f8').tobytes() for k in sorted(ops))).hexdigest()
    T=.37; rate=.8; archive=np.load(CHAIN/'clock/common_clock_states.npz')
    # Ordered one-event quadrature and an independent full density ODE.
    # Neither uses the producer's vectorized counting generator.
    one=np.zeros((384,384),complex); zero=np.zeros_like(one); total_reference=np.zeros_like(one)
    for q,(ev,v) in spec.items():
        ids=np.flatnonzero(numbers==q); state=rho[np.ix_(ids,ids)]
        if np.linalg.norm(state)<1e-15:
            continue
        r=v.T@state@v; omega=ev[:,None]-ev[None,:]
        zero[np.ix_(ids,ids)] = math.exp(-rate*T)*(v@(r*np.exp(-1j*T*omega))@v.T)
        # Separate two-frequency contractions via Gauss-Legendre quadrature;
        # 18 and 24 nodes bound the reported numerical comparison empirically.
        def integrate(nodes):
            grid,weights=np.polynomial.legendre.leggauss(nodes); out=np.zeros_like(r)
            klocal=[sparse.csr_matrix(k[np.ix_(ids,ids)]) for k in ops.values()]
            for t,w in zip((grid+1)*T/2,weights*T/2):
                before=v@(r*np.exp(-1j*t*omega))@v.T
                jump=sum((k@before)@k.T for k in klocal)
                out+=w*np.exp(-1j*(T-t)*omega)*(v.T@jump@v)
            return rate*math.exp(-rate*T)*(v@out@v.T)
        q18=integrate(18); q24=integrate(24)
        assert err(q18,q24)<3e-14
        one[np.ix_(ids,ids)] = q24
        klocal=[sparse.csr_matrix(k[np.ix_(ids,ids)]) for k in ops.values()]
        h=H[np.ix_(ids,ids)]
        def actual_rhs(t,y):
            x=y.reshape(state.shape)
            return (-1j*(h@x-x@h)+rate*(sum((k@x)@k.T for k in klocal)-x)).ravel()
        solved=solve_ivp(actual_rhs,(0,T),state.ravel(),method='DOP853',rtol=3e-11,atol=3e-14)
        assert solved.success
        total_reference[np.ix_(ids,ids)]=solved.y[:,-1].reshape(state.shape)
    def archive_sum(suffix):
        return sum(archive[row['complete_final_count_state_archive_keys'][suffix]] for row in contract['preparation_labels'])
    zold,o_old,told=[archive_sum(i) for i in range(3)]
    count_errors=[err(zero,zold),err(one,o_old)]
    count_traces=[float(np.trace(v).real) for v in (zold,o_old,told)]
    probabilities=[math.exp(-rate*T),rate*T*math.exp(-rate*T),1-(1+rate*T)*math.exp(-rate*T)]
    claim('clock.complete_count_0_1_tail',max(count_errors)<3e-14 and max(abs(np.array(count_traces)-probabilities))<3e-14,
          {'independent_zero_one_errors':count_errors,'traces':count_traces,'infinite_tail_retained':True})
    psd_min=min(float(np.linalg.eigvalsh(archive[key]).min()) for row in contract['preparation_labels'] for key in row['complete_final_count_state_archive_keys'])
    claim('clock.all_reservoir_count_states',psd_min>-2e-13 and err(zold+o_old+told,archive['final_state'])<2e-13,
          {'retained_joint_states':12,'minimum_eigenvalue':psd_min})
    full_error=err(total_reference,archive['final_state'])
    tail_error=err(total_reference-zero-one,told)
    claim('clock.actual_prepared_complete_output_and_tail',max(full_error,tail_error)<2e-11,
          {'full_preparation_direct_ODE_max_error':full_error,'tail_state_max_error':tail_error})
    # Actual q=1 source sector, direct density ODE using independently rebuilt K.
    ids=np.flatnonzero(numbers==1); h=H[np.ix_(ids,ids)]
    klocal=[sparse.csr_matrix(k[np.ix_(ids,ids)]) for k in ops.values()]
    v=np.zeros(len(ids),complex);v[0]=1/math.sqrt(2);v[4]=1j/math.sqrt(2)
    r=np.outer(v,v.conj()); full=np.zeros((384,384),complex);full[np.ix_(ids,ids)]=r
    def rhs(t,y):
        x=y.reshape(r.shape); jump=sum((k@x)@k.T for k in klocal)
        return (-1j*(h@x-x@h)+rate*(jump-x)).ravel()
    independent=solve_ivp(rhs,(0,T),r.ravel(),method='DOP853',rtol=2e-11,atol=2e-13).y[:,-1].reshape(r.shape)
    got=clock.semigroup(full,T,rate)[np.ix_(ids,ids)]
    scaled=clock.semigroup(full,T,rate/2,h_scale=.5)
    half=clock.semigroup(full,T/2,rate)
    rate_only=clock.semigroup(full,T,rate/2)
    pgf=clock.semigroup(full,T,rate,attempt_mark=.43,idle_rate=.8,idle_mark=.61)
    pgf_expected=math.exp(T*(rate*(.43-1)+.8*(.61-1)))
    claim('clock.generator_and_joint_scaling',err(independent,got)<1e-10 and err(scaled,half)<1e-13,
          {'direct_density_ODE_error':err(independent,got),'joint_scaling_error':err(scaled,half),
           'rate_only_vs_half_time_frobenius':float(np.linalg.norm(rate_only-half))})
    claim('clock.independent_Poisson_attempt_idle_law',abs(np.trace(pgf)-pgf_expected)<2e-13,
          {'PGF_trace_error':float(abs(np.trace(pgf)-pgf_expected)),'expected_attempts':rate*T,'expected_services':1+1.6*T})
    _,rows=clock.native(actual,True)
    noop_mass=sum(r['probability'] for r in rows if r['applied_current']==0)
    evolved=clock.evolve(actual,T)
    _,after_rows=clock.native(evolved,True)
    readout_difference=max(abs(a['probability']-b['probability']) for a,b in zip(rows,after_rows))
    # Conditional path with one actual zero-current mark and one identity idle.
    _,early=clock.native(clock.evolve(actual,.07),True)
    mark=max((i for i,row in enumerate(early) if row['applied_current']==0),key=lambda i:early[i]['probability'])
    state,ledger=clock.trajectory(actual,T,rate,1.6,[(.07,mark),(.2,'idle')])
    branch=clock.evolve(actual,.07); k=ops[tuple(clock.ops[mark][0][x] for x in ('seam','law','coin','record_value'))]
    branch=k@branch@k.T; probability=float(np.trace(branch).real)
    conditional=clock.evolve(branch/probability,T-.07)
    density=math.exp(-1.6*T)*rate*.8*probability
    claim('clock.conditional_history_and_resources',err(state,conditional)<1e-13 and abs(density-ledger['path_density_wrt_ordered_times'])<1e-13 and ledger['native_commit_pulses']==1 and ledger['coherent_service_calls']==3,
          {'state_error':err(state,conditional),'density_error':abs(density-ledger['path_density_wrt_ordered_times']),'ledger':{k:ledger[k] for k in ('native_zero_current_attempts','uniformization_idles','native_commit_pulses','coherent_service_calls','supplied_coherent_duration')}})
    claim('clock.drop_noop_drop_tail_wrong_clock_counterexamples',noop_mass>.5 and probabilities[2]>.03 and readout_difference>.001 and np.linalg.norm(rate_only-half)>.01,
          {'deleted_noop_mass_one_native_attempt':noop_mass,'deleted_count_tail':probabilities[2],'zero_attempt_native_readout_difference':readout_difference})
    malformed={}
    def rejected(name,fn):
        try:
            value=fn();malformed[name]={'rejected':False,'nonfinite_output':bool(not np.isfinite(value).all()) if isinstance(value,np.ndarray) else False}
        except (ValueError,TypeError,IndexError,AssertionError) as e:
            malformed[name]={'rejected':True,'exception':type(e).__name__}
    cross=full.copy(); a=np.flatnonzero(numbers==0)[0];b=ids[0];cross[a,b]=cross[b,a]=.1
    rejected('cross_number_coherence',lambda:clock.semigroup(cross,0,.8))
    rejected('negative_duration',lambda:clock.semigroup(full,-1,.8))
    rejected('invalid_uniform_rate',lambda:clock.trajectory(actual,T,.8,.4,[]))
    rejected('invalid_event_mark',lambda:clock.trajectory(actual,T,.8,1.6,[(.1,-1)]))
    bad=full.copy();bad[ids[0],ids[0]]=2;bad[ids[4],ids[4]]=-1
    rejected('negative_density',lambda:clock.semigroup(bad,0,.8))
    rejected('negative_density_trajectory',lambda:clock.trajectory(bad,0,.8,1.6,[])[0])
    nan=full.copy();nan[ids[0],ids[0]]=float('nan')
    rejected('nan_density',lambda:clock.semigroup(nan,0,.8))
    claim('clock.declared_parameter_and_number_admission',all(malformed[x]['rejected'] for x in ('cross_number_coherence','negative_duration','invalid_uniform_rate','invalid_event_mark')),malformed)
    claim('clock.general_density_input_validation',all(malformed[x]['rejected'] for x in ('negative_density','negative_density_trajectory','nan_density')),
          {'classification':'API hardening gap; positivity and finiteness are mathematical input preconditions, not a false finite-member theorem','tests':{x:malformed[x] for x in ('negative_density','negative_density_trajectory','nan_density')}})
    DETAIL['clock']={'operator_array_sha256_sorted_le_float64':operator_sha,'malformed_inputs':malformed,'contract_source':contract['source_id']}
    DETAIL['clock']['malformed_fixture_basis_indices']=[int(ids[0]),int(ids[4])]


def audit_force():
    manifest_path=CHAIN/'force/MANIFEST.json'
    pin(manifest_path,'f4f2f6b490aef1a87eb0fd6b5c07c429d5f394baf4414e0f39b7dc309babba1e')
    manifest=json.loads(manifest_path.read_text())
    for rel,expected in manifest['files'].items():
        pin(CHAIN/'force'/rel,expected)
    contract=json.loads((CHAIN/'force/FORCE_CONTRACT.json').read_text())
    for directory,name,expected in ((BRIDGE,'FINAL_MANIFEST.json',contract['source_manifest_sha256']),
                                   (CHAIN/'field_limit','MANIFEST.json',contract['field_manifest_sha256'])):
        pin(directory/name,expected)
        for rel,expected in json.loads((directory/name).read_text())['files'].items():
            pin(directory/rel,expected)
    result=json.loads((CHAIN/'force/RESULT.json').read_text())
    for rel,expected in result['source_pins'].items():
        pin(ROOT/rel,expected)
    fixed=load('field_feedback',BRIDGE/'field_feedback.py')
    reciprocal=load('audited_reciprocal',BRIDGE/'reciprocal_feedback.py')
    subject=load('audited_force_measurement',CHAIN/'force/force_measurement.py')
    edges,neighbours=fixed.graph()
    # Independent Green polynomial evaluation from the 96x96 graph spectrum.
    lap=6*np.eye(96)
    for a,b in edges:lap[a,b]-=1;lap[b,a]-=1
    eigen,vec=np.linalg.eigh(lap)
    weights=np.zeros(96);weights[1:]=-np.expm1(512*np.log1p(-eigen[1:]/12))/eigen[1:]
    green=(vec*weights)@vec.T
    table=(np.diag(green)[:,None]+np.diag(green)[None,:]-2*green)/2
    frozen_table=reciprocal.pair_table(neighbours)
    claim('force.actual_96_port_field_operator',len(edges)==288 and eigen[1]>.5 and err(table,frozen_table)<3e-14,
          {'ports':96,'seams':len(edges),'spectral_gap':float(eigen[1]),'spectral_vs_stencil_table_error':err(table,frozen_table)})
    native=fixed.literal_native()
    def independent_rates(state,kappa=1,nu=1):
        rates=np.ones((2,288),float)*nu
        for side,p in enumerate(state):
            for e,(a,b) in enumerate(edges):
                if p in (a,b):
                    target=list(state);target[side]=a+b-p
                    rates[side,e]=nu*math.exp(-2*kappa*(table[tuple(target)]-table[state]))
        return rates
    def expected_law(state,kappa=1,nu=1):
        rr=independent_rates(state,kappa,nu);total=rr.sum();g=-math.expm1(-.001*total)/total
        law=[[math.exp(-.001*total),-1,-1,-1,0,0,*state,0]]
        for side,p in enumerate(state):
            for seam,(a,b) in enumerate(edges):
                for coin in (0,1):
                    target=list(state);target[side]=(a if coin else b) if p in (a,b) else p
                    j=int(p==a)-int(target[side]==a)
                    law.append([g*rr[side,seam]/2,side,seam,coin,j,j if side==0 and seam==0 else 0,*target,1])
        return np.array(law),g,total
    archive=np.load(CHAIN/'force/FIRST_EVENT_OUTCOMES.npz')
    discrepancies=[];native_events=0;cap_rows=[]
    for measurement in result['measurements']:
        laws=archive[measurement['archive_key']];kappa=measurement['kappa'];drift=0
        for x in (0,1):
            expected,g,total=expected_law((x,9),kappa)
            discrepancies.append(err(expected,laws[x]))
            for row in expected[1:]:
                _,side,seam,coin,j,_,xx,yy,_=row;side,seam,coin=map(int,(side,seam,coin))
                input_load=tuple(int(p==(x,9)[side]) for p in range(96))
                out,record=native(input_load,edges[seam],coin)
                assert sum(out)==1 and max(out)==1 and record==(coin,int(j)) and out[int((xx,yy)[side])]==1
                native_events+=1
            drift += expected[:,0]@expected[:,5]/g/2
        force=kappa*(table[0,9]-table[1,9]);z=2*force
        assert abs(drift-.5*math.sinh(z))<3e-14
        assert abs(drift-force)<=abs(z)**3*math.cosh(abs(z))/12+3e-14
        cap_rows.append({'kappa':kappa,'force':force,'drift':float(drift)})
    claim('force.full_native_cap_one_laws',max(discrepancies)<3e-14,
          {'actual_literal_native_events_checked':native_events,'outcomes_per_preparation':1153,'max_archive_error':max(discrepancies)})
    claim('force.measured_force_drift_and_Einstein_remainder',True,cap_rows)
    # Exact finite support of the first jump; no 9216-state matrix.
    def kernel(state,kappa=1):
        rr=independent_rates(state,kappa); A=defaultdict(float)
        for side,p in enumerate(state):
            for seam,(a,b) in enumerate(edges):
                for coin in (0,1):
                    target=list(state);target[side]=(a if coin else b) if p in (a,b) else p
                    A[tuple(target)]+=rr[side,seam]/2
        return A,float(rr.sum())
    cap2=[]
    for kappa in (0,1):
        for state in ((0,9),(1,9)):
            A,L=kernel(state,kappa);T=.001;g=-math.expm1(-L*T)/L
            p0=math.exp(-L*T);p1=0;p2=0
            grid,w=np.polynomial.legendre.leggauss(24);t=(grid+1)*T/2;w=w*T/2
            integral_one=0
            for dest,a in A.items():
                _,Ly=kernel(dest,kappa);d=Ly-L
                h=T*math.exp(-L*T) if abs(d)<1e-10 else math.exp(-L*T)*(-math.expm1(-d*T))/d
                p1+=a*h;p2+=a*(g-h)
                integral_one+=a*float(w@(np.exp(-L*t)*np.exp(-Ly*(T-t))))
            assert abs(p0+p1+p2-1)<3e-14 and abs(p1-integral_one)<3e-14
            if kappa==0:assert abs(p1-L*T*math.exp(-L*T))<3e-14
            cap2.append({'kappa':kappa,'initial':state,'first_step_position_support':len(A),'no_event':p0,'one_event_no_cap':p1,'cap_two':p2,'mass':p0+p1+p2})
    # Independently replay actual cap-one/two paths, retaining every hold.
    paths=[]
    for cap in (1,2):
        path=reciprocal.execute((0,9),edges,frozen_table,1.,89,horizon=.5,cap=cap)
        state=(0,9);logp=0;time=0
        for row in path['records']:
            rr=independent_rates(state);side=row['copy'];a,b=edges[row['seam']];p=state[side]
            dest=(a if row['coin'] else b) if p in (a,b) else p
            j=int(p==a)-int(dest==a);assert j==row['current']==row['observer_current']
            logp += -float(rr.sum())*row['wait']+math.log(rr[side,row['seam']]/2)
            state=list(state);state[side]=dest;state=tuple(state);time+=row['wait']
        assert path['terminal']['kind']=='cap_exhausted' and len(path['records'])==cap
        assert abs(logp-path['log_path_density'])<1e-11
        paths.append({'cap':cap,'records':cap,'holds':sum(r['current']==0 for r in path['records']),'path_density_error':abs(logp-path['log_path_density'])})
    claim('force.complete_cap_two_partition_and_executor',True,{'partitions':cap2,'actual_paths':paths})
    charge=np.eye(96)[0]-np.eye(96)[9];B=np.eye(96)[0]-np.eye(96)[1]
    bare=float(B@green@charge);selfterm=float(B@green@B)/2;force=table[0,9]-table[1,9]
    primary=archive['kappa_1'];null_mass=primary[:,0,0];holds=np.array([p[(p[:,8]==1)&(p[:,4]==0),0].sum() for p in primary])
    _,_,L=expected_law((0,9));_,_,Ltwice=expected_law((0,9),nu=2)
    wrong_clock_survival=math.exp(-Ltwice*.001)
    claim('force.self_field_required',abs(force-(bare-selfterm))<3e-14 and selfterm>.1,
          {'bare_field':bare,'self_term':selfterm,'force':float(force),'identity_error':abs(force-(bare-selfterm))})
    claim('force.drop_hold_drop_cap_wrong_clock_counterexamples',holds.min()>.4 and (1-null_mass).min()>.4 and abs(wrong_clock_survival-null_mass[0])>.1,
          {'mass_after_cap_deletion':null_mass.tolist(),'deleted_hold_mass':holds.tolist(),'wrong_clock_survival':wrong_clock_survival,'declared_survival':float(null_mass[0]),'conclusive_only_drift':0})
    # Independent integer recurrence and exact rational transfer bounds.
    transfer=contract['controlled_field_to_rate_transfer'];S=1<<40
    energy={};max_rounding_error=0
    for row in transfer['compiled_pairs']:
        x,y=row['state'];q=[int(i==x)-int(i==y) for i in range(96)];z=[0]*96
        for _ in range(512):
            numerators=[6*z[i]+sum(z[j] for j in neighbours[i])+S*q[i] for i in range(96)]
            rounded=[]
            for numerator in numerators:
                quotient,remainder=divmod(numerator,12)
                rounded.append(quotient+int(remainder>=6))
            z=rounded
        val=F(sum(a*b for a,b in zip(q,z)),2*S);assert val==F(row['exact_energy'])
        residual=[S*q[i]-6*z[i]+sum(z[j] for j in neighbours[i]) for i in range(96)]
        cert=F(2*sum(v*v for v in residual),S*S)
        assert cert==F(row['residual_field_error_squared_upper'])
        energy[(x,y)]=val
    theta=F(23,24)**512;eta=2*theta+F(408,S)*(1-theta);delta=eta+2*theta;ell=4*delta
    tv=F(96,1000)*ell/(1-ell);outward=F(transfer['complete_marked_capped_path_TV_upper'])
    assert tv<=outward and outward<F(671,10**12)
    assert F(contract['rational_force_score_transfer']['frozen_score_error_upper'])==2*outward/F(67,100000)
    # The rate envelope follows from |Delta U_exact|<=kappa/2 and frozen error.
    eupper=sum(F(1,math.factorial(k)) for k in range(9))+F(1,math.factorial(9))/(1-F(1,10))
    assert eupper<F(11,4) and F(11,4)**2/(1-16*theta)<8
    claim('force.fixedpoint40_all_compiled_pairs_and_TV_bound',len(energy)==22,
          {'pairs_checked':len(energy),'TV_upper':str(outward),'score_error_upper':contract['rational_force_score_transfer']['frozen_score_error_upper'],
           'comparison':'exact-real exponentials and ideal frozen stencil; floating library evaluation is excluded',
           'finite_cap_scope':'couple common marked rates and coins until first rate disagreement, stop at common cap'})
    malformed={}
    # Actual low-level callable accepts some out-of-domain settings.
    for name,kwargs in [('negative_horizon',{'horizon':-.001}),('negative_rate',{})]:
        rr=independent_rates((0,9));
        if name=='negative_rate':rr[0,0]=-1
        try:
            got,_,_=subject.first_event((0,9),edges,rr,native,**kwargs)
            malformed[name]={'rejected':False,'minimum_probability':float(got[:,0].min()),'sum':float(got[:,0].sum())}
        except (ValueError,AssertionError) as e:
            malformed[name]={'rejected':True,'exception':type(e).__name__}
    claim('force.general_rate_horizon_input_validation',all(x['rejected'] for x in malformed.values()),
          {'classification':'API hardening gap outside frozen beta/kappa/nu/T member','tests':malformed})
    DETAIL['force']={'malformed_inputs':malformed,'source_identity':contract['source_identity'],
                     'preparation_array_sha256_le_float64':hashlib.sha256(np.asarray([[0,9,.5],[1,9,.5]],dtype='<f8').tobytes()).hexdigest(),
                     'independent_GK_array_sha256_le_float64':hashlib.sha256(np.asarray(green,dtype='<f8').tobytes()).hexdigest()}


def main():
    audit_clock()
    print('Clock audit complete',flush=True)
    audit_force()
    claim('source_identity_boundary',True,{'quantum_C10_to_classical_force_join':False,'common_name_implies_transfer':False,
          'physical_seconds':False,'coherent_service_native_compiler':False,'dynamic_continuum_force_limit':False})
    for rel,expected in PINS.items():assert sha(ROOT/rel)==expected
    scientific_failures=[c['claim'] for c in CLAIMS if c['status']=='FAIL' and not c['claim'].endswith('input_validation')]
    result={'schema':'oph.chain.independent_opposing_clock_force_audit.v1','reviewer':'worker7 independent opposing audit',
            'status':'PASS_SCOPED_WITH_ADMISSION_GAPS' if not scientific_failures else 'FAIL',
            'scientific_blockers':scientific_failures,'claims':CLAIMS,'details':DETAIL,'input_pins':PINS,
            'reviewer_script_sha256':sha(__file__),'producer_verifiers_executed':False,'full_9216_semigroup_constructed':False,
            'official_rows_changed':[],'producer_files_modified':[]}
    (HERE/'AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'claims':len(CLAIMS),'failures':[c['claim'] for c in CLAIMS if c['status']=='FAIL'],'pins':len(PINS)},indent=2))


if __name__=='__main__':main()
