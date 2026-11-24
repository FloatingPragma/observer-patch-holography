"""Predict the collective branch intercept and gradient coefficient from two local meters."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_30_autonomous.experiment import PatchClock, SHOTS as HAND_TRIALS
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
TRIALS = 2**17
SEED = 301541
SAMPLES = 161
WORDS = np.array(list(itertools.product((0,1),repeat=3)))
I, X, Y, Z = (np.eye(2), np.array([[0,1],[1,0]]), np.array([[0,-1j],[1j,0]]),
              np.diag([1.,-1.]))
NAMES = ('exchange','drawn_link','empty_link','vacant_neighbour')


def restrict(generator, columns):
    block = columns.conj().T @ generator @ columns
    err = np.linalg.norm(generator @ columns-columns @ block)
    assert err < 1e-12
    return block, float(err)


def local_generators():
    law=candidate.native_law
    exch, contrast=law.local_term_parts()
    full=exch+candidate.KAPPA*contrast
    eye=np.eye(law.DIM)
    def blank(pair):
        return eye[:,law.SRC.index(pair)]
    def memory(write,out):
        return (law.echo() @ write).reshape(-1,law.N_SRC)[:,law.SRC.index(out)]
    drawn=memory((blank((2,0))+blank((0,2)))/np.sqrt(2),(1,1))
    empty=memory(blank((0,0)),(0,0))
    errors={}
    dimer=np.column_stack([blank((0,1)),blank((1,0))])
    hx,errors['blank_exchange']=restrict(full,dimer)
    # Independent endpoint phase service. The other physical probe port
    # holds the complementary load throughout this single-endpoint pulse.
    endpoint=np.column_stack([blank((0,0)),blank((1,0))])
    hz,errors['phase_service']=restrict(candidate.KAPPA*contrast,endpoint)
    out={'exchange':hx}
    for name,mark,neighbour in [('drawn_link',drawn,1),('empty_link',empty,1),
                                ('vacant_neighbour',drawn,0)]:
        # The complementary probe port distinguishes and conserves the two
        # equal-total-load global branches. Its addressed interaction is idle.
        cols=np.column_stack([np.kron(mark,np.eye(law.N_SRC)[:,law.SRC.index((p,neighbour))])
                              for p in (0,1)])
        np.testing.assert_allclose(cols.T@cols,I,atol=1e-12)
        out[name],errors[name]=restrict(full,cols)
    return out,hz,errors


def first_crossing(f):
    lo=0.
    for hi in np.arange(.025,16.025,.025):
        if f(lo)*f(hi)<=0:
            return brentq(f,lo,hi)
        lo=hi
    raise ArithmeticError('Unable to balance local analyzer')


def analyzers(hx,hz):
    init=np.array([1.,0.])
    tx=first_crossing(lambda t: abs((expm(-1j*hx*t)@init)[1])**2-.5)
    ux=expm(-1j*hx*tx)
    probe=ux@init
    # Cancel Ramsey contrast with a degree-phase pulse. Only equal splits
    # enter the preparation of these lab axes.
    tz=first_crossing(lambda t: abs((ux@expm(-1j*hz*t)@probe)[0])**2-.5)
    uz=expm(-1j*hz*tz)
    ops={'x':ux,'z':uz}
    analy={}
    for label,axis in [('X',X),('Y',Y),('Z',Z)]:
        target=(I+axis)/2
        for depth in range(5):
            for seq in itertools.product(ops,repeat=depth):
                u=I
                for key in seq:
                    u=ops[key] @ u
                for bit in (0,1):
                    projector=np.outer(u[bit].conj(),u[bit])
                    if np.linalg.norm(projector-target)<1e-9:
                        analy[label]=(projector, ''.join(seq),bit)
                        break
                if label in analy:
                    break
            if label in analy:
                break
        assert label in analy
    return ini_probes(probe,analy),dict(balanced_exchange_hold=tx,balanced_phase_hold=tz,
        terminal_axes={k:dict(pulses=v[1],read_port=v[2]) for k,v in analy.items()})


def ini_probes(probe,analy):
    # Under exchange from |0>, Z+iY advances. The drawn-link phase,
    # starting at +Y, advances Y−iX.
    z,y,x=[analy[n][0] for n in ('Z','Y','X')]
    return [(np.array([1.,0.]),[z,y,I-z,I-y])] + [(probe,[y,I-x,I-y,x])]*3


def addresses():
    car=federation(0)
    a,b=next((u,v) for u,v in car.seams if u//12==v//12==0 and car.degree[u]==car.degree[v])
    c=next(p for p in sorted(car.neighbours[a]) if p not in (a,b) and p>a and p//12==0)
    unused=sorted((u,v) for u,v in car.seams if u//12==v//12==1)
    clock=next((a,b) for a,b in itertools.combinations(unused,2) if set(a).isdisjoint(b))
    service=next(p for p in sorted(car.neighbours[a]) if p not in (b,c) and p>a and p//12==0)
    return dict(probe_dimer=[a,b],condition_link=[a,c],probe_phase_service=[a,service],hand_seams=clock)


def measure(counts, calibration):
    shots=counts.sum(axis=-1,keepdims=True)
    rates=(counts @ WORDS)/shots
    offset=np.mean(calibration)
    a,b=(calibration[:2]-calibration[2:])/2
    hand=rates[:,:,:,0:2].mean(axis=(0,1))
    cs=np.linalg.solve([[a,b],[b,-a]],(hand-offset).T)
    t=np.unwrap(np.arctan2(cs[1],cs[0]))
    t-=t[0]
    y=rates[:,:,:,2]
    angle=np.unwrap(np.angle((y[:,0]-y[:,2])+1j*(y[:,1]-y[:,3])),axis=1)
    tc=t-t.mean()
    slope=angle @ tc / (tc@tc)
    if not np.max(abs(np.diff(angle,axis=1)))<np.pi*.9:
        raise ArithmeticError('Unresolved phase step')
    h=slope[0]/2
    j=slope[1]/4
    return np.array([h,j,2*abs(j-h),4*j*h,*slope]),float(np.max(abs(np.diff(angle,axis=1))))


def run():
    mats,hz,errors=local_generators()
    instruments,setup=analyzers(mats['exchange'],hz)
    rng=np.random.default_rng(SEED)
    hand=PatchClock(rng)
    # Unequal predeclared sample steps. Neither commanded exposure nor the
    # eigenvalues of either probe enter the estimator.
    times=np.r_[0.,np.cumsum(np.resize([.037,.026,.031],SAMPLES-1))]
    probs=np.empty((len(NAMES),4,len(times),8))
    for index,(name,(probe,axes)) in enumerate(zip(NAMES,instruments)):
        for i,t in enumerate(times):
            psi=expm(-1j*mats[name]*t)@probe
            for k,axis in enumerate(axes):
                p=np.clip([hand.prob(t),hand.prob(t+np.pi/2),np.vdot(psi,axis@psi).real],0.,1.)
                probs[index,k,i]=np.prod(np.where(WORDS,p,1-p),axis=1)
    counts=rng.multinomial(TRIALS,probs)
    # PatchClock uses 2^22 shots per calibration point. Reuse this single
    # calibration across both meters and both controls in every resample.
    cal=np.array(hand.calibration)/HAND_TRIALS
    v,step=measure(counts,cal)
    p=counts/TRIALS
    bootstrap=np.stack([measure(rng.multinomial(TRIALS,p),rng.binomial(HAND_TRIALS,cal)/HAND_TRIALS)[0]
                        for _ in range(400)])
    lo,hi=np.quantile(bootstrap,[.025,.975],axis=0)
    names=['h','J','predicted_m','predicted_c_squared',*NAMES]
    return dict(addresses=addresses(),setup=setup,columns=errors,hand_column_residual=hand.leakage,
                exposure_control='Addressed single-seam holds; spectator couplings idle during pulses',
                inference='Open-chain Ising reduction predicts m=2|J-h| and c²=4Jh in dimer-spacing units',
                uncertainty='Percentile bootstrap of all joint counts and a shared hand calibration; simulated shot noise only',
                orientations=['Z','Y','minus Z','minus Y'],link_orientations=['Y','minus X','minus Y','X'],
                trials=TRIALS,calibration_trials=HAND_TRIALS,hand_calibration=hand.calibration,
                results={k:dict(value=v[i],bootstrap95=[lo[i],hi[i]]) for i,k in enumerate(names)},
                bootstrap_replicates=400,maximum_observed_phase_step=step,
                readout_order=['hand1','hand2','probe'],joint_words=WORDS.tolist(),
                protocol_names=NAMES,all_counts=counts.tolist())


def canonical(obj):
    if isinstance(obj,dict):
        return {k:canonical(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):
        return [canonical(v) for v in obj]
    if isinstance(obj,(float,np.floating)):
        return round(float(obj),10)
    return obj


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--output',type=Path)
    opt=cli.parse_args()
    report=canonical(run())
    if opt.output:
        with opt.output.open('x') as f:
            json.dump(report,f,separators=(',',':'),allow_nan=False)
            f.write('\n')
    elif report!=json.loads((HERE/'coupling_meter.json').read_text()):
        raise SystemExit('Local coupling-meter replay differs')
    print(report['results'])
