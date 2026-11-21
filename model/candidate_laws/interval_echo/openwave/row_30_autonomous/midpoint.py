"""A local stationary matter reference and separate two-coupling velocity calibration."""
import argparse
import gzip
import json

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

from model.candidate_laws.interval_echo import candidate
from .channel import K, FreeChannel, apply, matrices, preparation, standing
from .experiment import HERE, TIMES, SHOTS, SQ, PatchClock, canonical, same_receipt


def physical_operations():
    law = candidate.native_law
    h0,h1 = law.local_term_parts()
    h = h0 + K*h1
    eye = np.eye(law.DIM)
    def state(loads,c=None,d=0):
        offset = 0 if c is None else 1+law.LABELS.index((c,d))
        return eye[:,offset*law.N_SRC + law.SRC.index(loads)]

    def bright(loads,c):
        return (state(loads,c,1)+state(loads,c,-1))/SQ
    light = np.column_stack((bright((1,0),1),bright((0,1),0)))
    spectator = bright((0,0),0)
    swap = np.array([[0,-1],[-1,0]])
    residuals = {}
    for theta in (.43,np.pi/2,np.pi):
        U = expm(-1j*h*theta/K)*np.exp(1j*theta)
        mixer = np.eye(2)+np.expm1(1j*theta)*np.ones((2,2))/2
        residuals[f'marked_mixer_{theta:.3f}'] = float(max(np.linalg.norm(U @ light-light @ mixer),
                                                          np.linalg.norm(U @ spectator - spectator)))
    swap_q = swap.T @ np.diag([0,1]) @ swap
    residuals['gateway_swap_reference_undo'] = float(np.linalg.norm(swap_q-np.diag([1,0])))
    B = np.column_stack((state((1,0)),state((0,1))))
    bridge = expm(-1j*h*np.pi/2)
    mix = np.array([[1,1j],[1j,1]])/SQ
    residuals['blank_bridge_mixer'] = float(np.linalg.norm(bridge @ B - np.exp(-1j*(K+.5)*np.pi/2)*B @ mix))
    dark = (state((0,0),0,0)+state((0,0),1,0))/SQ
    residuals['echo_opens_dark_vacant_bridge'] = float(np.linalg.norm(law.echo() @ dark-state((0,0))))
    residuals['echo_closes_blank_vacant_bridge'] = float(np.linalg.norm(law.echo() @ state((0,0))-dark))
    # A balanced temporary load on a blank seam; after the write a local endpoint
    # occupation selects one sharp coin with equal amplitudes on its two signs.
    seed = (B[:,0]+B[:,1])/SQ
    written = law.echo() @ seed
    coin_preps = []
    for c in (0,1):
        # Each coin has one deterministic outgoing endpoint occupation in the total-one sector.
        matches = [(abs(np.vdot(bright(pair,c),written))**2,pair) for pair in ((1,0),(0,1))]
        prob,pair = max(matches)
        address = np.tile(np.array([x == pair for x in law.SRC]),1+len(law.LABELS))
        conditional = written*address/np.sqrt(prob)
        residuals[f'sharp_coin_{c}_sign_even_after_native_write'] = float(np.linalg.norm(conditional-bright(pair,c)))
        coin_preps.append(dict(coin=c,occupied_pair=list(pair),probability=float(prob)))
    assert max(residuals.values())<1e-11, residuals
    degree = np.column_stack((state((0,0)),state((1,0))))
    degree_h = degree.T @ (K*h1) @ degree
    residuals['phase_pulse_degree_restriction'] = float(np.linalg.norm((K*h1) @ degree-degree @ degree_h))
    return residuals,coin_preps,light.T @ h @ light,degree_h


def instrument(counts,rng,clock,t,prob):
    """Local clock words times all matter outcomes, retaining failures."""
    p1,p2 = clock.prob(t),clock.prob(t+np.pi/2)
    words = np.array([(1-p1)*(1-p2),(1-p1)*p2,p1*(1-p2),p1*p2])
    counts[:] = rng.multinomial(SHOTS,(words[:,None]*prob).ravel()).reshape(4,-1)


def clock_times(counts,clock):
    rates = counts.sum(axis=-1).mean(axis=0)/SHOTS
    return clock.read(np.column_stack((rates[:,2]+rates[:,3],rates[:,1]+rates[:,3])))


def calibrate(rng,clock,hopping_h,degree_h):
    counts = np.zeros((3,len(TIMES),4,2),int)
    for i,t in enumerate(TIMES):
        hop_read = abs(expm(-1j*hopping_h*t)[1,0])**2
        p = [hop_read]
        for phase in (0,np.pi/2):
            state = np.array([1,np.exp(1j*phase)])/SQ
            # Physical ordinary (load/reference=0) service on one of the two branches.
            state = expm(-1j*degree_h*t) @ state
            state = (np.eye(2)+(1j-1)*np.ones((2,2))/2) @ state
            p.append(abs(state[0])**2)
        for arm,prob in enumerate(p):
            instrument(counts[arm,i],rng,clock,t,np.array([1-prob,prob]))
    t_measured = clock_times(counts,clock)
    read = counts[:,:,:,1].sum(axis=2)/SHOTS
    def hop_error(a):
        return np.sum((read[0]-np.sin(a*t_measured)**2)**2)
    grid = np.linspace(0,16,801)
    peak = int(np.argmin([hop_error(x) for x in grid]))
    a = minimize_scalar(hop_error,bounds=(grid[peak-1],grid[peak+1]),method='bounded').x
    angle = np.unwrap(np.angle((read[1]-.5)+1j*(read[2]-.5)))
    delta = abs(np.polyfit(t_measured,angle,1)[0])
    v = delta/2
    c2 = a**3/(4*np.hypot(a,v))
    return dict(joint_counts=counts.tolist(),hopping=a,degree_contrast=delta,coefficient=np.sqrt(c2))


def direct_line(rng,clock,main,tri,mode,reference):
    phases = [1,-1,1j,-1j]
    counts = np.zeros((4,len(TIMES),4,3),int)
    U = np.array([[1,1j],[1j,1]])/SQ
    psi = main.launch(mode)
    phi = tri.launch(reference)
    # Balanced blank bridge at launch; a known degree-five reference rotation sets each arm phase.
    src = U @ np.array([0,1])
    main_spec = main.v.T @ psi
    tri_spec = tri.v.T @ phi
    for i,t in enumerate(TIMES):
        a = main.v[:2] @ (main_spec*np.exp(-1j*main.e*t))
        b = tri.v[:2] @ (tri_spec*np.exp(-1j*tri.e*t))
        for s,phase in enumerate(phases):
            # Native swaps on each arm carry their unique changed collar coin; multiplying
            # their relative gate matrices by the same scalar includes the remote spectator.
            A,B = a*src[0]*phase, b*src[1]
            failure = abs(A[1])**2+abs(B[1])**2
            # In the accepted arm the endpoints are parked at 1 while the vacant
            # bridge is un-written, then unparked. All remote configurations yield quiet.
            parked_A,parked_B = -A[0],-B[0]
            active_after_return = np.array([-parked_A,-parked_B])
            click = abs((U @ active_after_return)[1])**2
            instrument(counts[s,i],rng,clock,t,np.array([failure,1-failure-click,click]))
    t_measured = clock_times(counts,clock)
    rates = counts[:,:,:,2].sum(axis=2)/SHOTS
    trace = (rates[0]-rates[1])+1j*(rates[2]-rates[3])
    phase = np.unwrap(np.angle(trace))
    slope,intercept = np.polyfit(t_measured,phase,1)
    rms = np.sqrt(np.mean((phase-(slope*t_measured+intercept))**2))
    return dict(signed_frequency=slope,phase_residual_rms=rms,joint_counts=counts.tolist(),
                mean_parking_failure=counts[:,:,:,0].sum(axis=2).mean()/SHOTS)


def run():
    with gzip.open(HERE/'receipt.json.gz','rt') as fp:
        prior = json.load(fp)
    rng = np.random.default_rng(202611213003)
    car,cells,sites,_ = preparation()
    M = len(cells)
    d,o = matrices(car,sites)
    main = FreeChannel(d,o)
    trimer = [2,11,0]
    hands = ((6,7),(8,9))
    assert all(car.degree[t]==5 for t in trimer)
    assert trimer[0] in car.neighbours[sites[0]] and sites[0]==3
    assert not set(trimer).intersection(sites)
    assert all(b in car.neighbours[a] for a,b in zip(trimer[:-1],trimer[1:]))
    assert all(a in car.neighbours[b] for a,b in hands)
    assert not set(sum(hands,())).intersection(set(trimer)|set(sites))
    tri_d,tri_off = matrices(car,trimer)
    tr = FreeChannel(tri_d,tri_off)
    reference = np.array([1,2,1])/np.sqrt(6)
    midpoint = -K/2
    residual = np.linalg.norm(apply(tri_d,tri_off,reference)-midpoint*reference)
    clock = PatchClock(rng)
    gates,coin_preps,h_active,h_degree = physical_operations()
    local_cal = calibrate(rng,clock,h_active,h_degree)
    lines = []
    for old in prior['lines']:
        j = old['j_source']
        q = j*np.pi/M
        wave,_ = standing(M,q,old['band'],d,o)
        line = direct_line(rng,clock,main,tr,wave,reference)
        lines.append(dict(q=q,band=old['band'],target=old['target'],direct=line['signed_frequency'],
                          two_line=old['beat']-prior['midpoint_relative'],**line))
    mass = abs(lines[0]['direct'])
    held = []
    for r in lines:
        if r['band']==1 and r['target']>.30:
            w = abs(r['direct'])
            kg = np.hypot(mass,local_cal['coefficient']*r['q'])
            packet_kg = np.hypot(mass,prior['inferred_c']*r['q'])
            held.append(dict(q=r['q'],direct_rest_reference=w,local_cal_kg=kg,
                             local_cal_KG_error=abs(kg/w-1),packet_KG_error=abs(packet_kg/w-1)))
    diffs = [abs(r['direct']-r['two_line']) for r in lines]
    speed_error = abs(local_cal['coefficient']/prior['inferred_c']-1)
    checks = dict(max_all_frequency_difference=max(diffs),speed_relative_difference=speed_error,
                  max_held_local_coefficient_KG_error=max(x['local_cal_KG_error'] for x in held[:-1]),
                  max_local_column_residual=float(max(gates.values())),
                  max_source_error=float(max(main.max_source_error,tr.max_source_error)))
    return dict(trimer=trimer,physical_bridge=[int(sites[0]),trimer[0]],parking=[int(sites[1]),trimer[1]],
                hand_seams=hands,reference_residual=float(residual),clock_calibration=clock.calibration,
                shots=SHOTS,operations=gates,coin_preps=coin_preps,calibration=local_cal,rest=mass,
                held_and_challenge=held,checks=checks,lines=lines,
                passes=bool(max(diffs)<.02 and speed_error<.02 and checks['max_held_local_coefficient_KG_error']<.03))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=canonical(run())
    filename=HERE/'midpoint.json.gz'
    if args.write:
        with filename.open('xb') as raw, gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as fp:
            fp.write(json.dumps(result,separators=(',',':'),allow_nan=False).encode())
    else:
        with gzip.open(filename,'rt') as fp:
            assert same_receipt(result,json.load(fp))
    print(result['checks'],result['passes'],result['held_and_challenge'])
    print({k:v for k,v in result['calibration'].items() if k!='joint_counts'})
