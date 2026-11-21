"""Free middle-band spectroscopy and packet time of flight on the registered screen."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo.openwave.row_30_arrows.clock import instrument
from .channel import (K, FreeChannel, bloch, column_checks, compile_source, matrices,
                      packet, preparation, standing)


HERE = Path(__file__).resolve().parent
SHOTS = 4_194_304
PACKET_SHOTS = 16_777_216
TIMES = np.linspace(0,5.80,81)
WORDS = np.array(list(itertools.product(range(2),repeat=3)))
SQ = np.sqrt(2)


class PatchClock:
    def __init__(self,rng):
        self.effect,self.levels,self.leakage,self.blind = instrument()
        self.calibration = [int(rng.binomial(SHOTS,self.prob(t))) for t in np.arange(4)*np.pi/2]
        clicks = np.array(self.calibration)/SHOTS
        self.offset = clicks.mean()
        a,b = (clicks[:2]-clicks[2:])/2
        self.design = np.array([[a,b],[b,-a]])

    def prob(self,t):
        state = np.array([0,1,1])/SQ
        return float(abs(self.effect @ (state*np.exp(-1j*self.levels*t)))**2)

    def read(self,rates):
        quadratures = np.linalg.solve(self.design,(rates-self.offset).T).T
        times = np.unwrap(np.arctan2(quadratures[:,1],quadratures[:,0]))
        return times-times[0]


def spatial_index(counts,M):
    """Estimate sinusoidal cell frequency, allowing any of four port polarizations."""
    contrast = ((counts[0]-counts[1])/SHOTS).reshape(M,4)
    n = np.arange(M)
    errors = []
    for j in range(M):
        phase = (np.pi*j/M)*n
        design = np.column_stack((np.cos(phase),np.sin(phase)))
        model = design @ np.linalg.lstsq(design,contrast,rcond=None)[0]
        errors.append(np.linalg.norm(contrast-model))
    return int(np.argmin(errors))


def record_line(rng, clock, channel, M, ref, v, address, time):
    launches = [channel.launch((ref+s*v)/SQ) for s in (1,-1,1j,-1j)]
    space = [rng.multinomial(SHOTS,abs(arm)**2) for arm in launches[:2]]
    j = spatial_index(space,M)
    p = np.array([channel.clicks(arm,address,time) for arm in launches])
    counts = np.zeros((4,len(time),8),int)
    for i,t in enumerate(time):
        clock_prob = (clock.prob(t),clock.prob(t+np.pi/2))
        for arm in range(4):
            probs = np.array([*clock_prob,p[arm,i]])
            joint = np.prod(np.where(WORDS,probs,1-probs),axis=1)
            counts[arm,i] = rng.multinomial(SHOTS,joint)
    rates = (counts @ WORDS)/SHOTS
    observed_time = clock.read(rates[:,:,:2].mean(axis=0))
    y = rates[:,:,2]
    complex_signal = (y[0]-y[1]) + 1j*(y[2]-y[3])
    phase = np.unwrap(np.angle(complex_signal))
    # A regression on phase versus observed hand angle needs no generator eigenvalues.
    slope,intercept = np.polyfit(observed_time,phase,1)
    residual = phase-np.polyval([slope,intercept],observed_time)
    return dict(spatial_counts=np.array(space).tolist(), joint_counts=counts.tolist(),
                j_observed=j, q_observed=np.pi*j/M, beat=float(slope),
                beat_fit_rms=float(np.sqrt(np.mean(residual**2))),
                smallest_interference_amplitude=float(abs(complex_signal).min()),
                max_clock_time_error=float(abs(observed_time-time).max()),
                beat_se=float(np.sqrt(np.sum(residual**2)/(len(time)-2)
                                      /np.sum((observed_time-observed_time.mean())**2))))


def packet_flight(rng,clock,channel,M,q):
    # Position projectors commute. Each repetition has one address and two local fuel bits.
    coord = np.repeat(np.arange(M),4)
    times = np.array([0,5.5])
    cases = []
    for sign in (1,-1):
        initial = channel.launch(packet(M,sign*q))
        runs = []
        clock_rates = []
        for t in times:
            prof = channel.profile(initial,t)
            probs = [clock.prob(t),clock.prob(t+np.pi/2)]
            outcomes = np.array(list(itertools.product((0,1), repeat=2)))
            word_probs = np.prod(np.where(outcomes,probs,1-np.array(probs)),axis=1)
            counts = rng.multinomial(PACKET_SHOTS,(word_probs[:,None]*prof).ravel())
            counts = counts.reshape(4,-1)
            address = counts.sum(axis=0)
            mean = np.dot(address,coord)/PACKET_SHOTS
            variance = np.dot(address,(coord-mean)**2)/PACKET_SHOTS
            clock_rates.append(counts.sum(axis=1) @ outcomes/PACKET_SHOTS)
            runs.append(dict(mean=mean, mean_se=np.sqrt(variance/PACKET_SHOTS),
                             end_cells_probability=float(prof[:20].sum()+prof[-20:].sum()),
                             joint_counts=counts.tolist()))
        # Read the later hand angle in the declared positive one-turn window.
        quadrature = np.linalg.solve(clock.design,(np.array(clock_rates)-clock.offset).T)
        hand_angle = np.arctan2(quadrature[1],quadrature[0])
        dt = (hand_angle[1]-hand_angle[0])%(2*np.pi)
        speed = (runs[1]["mean"]-runs[0]["mean"])/dt
        cases.append(dict(sign=sign,clock_elapsed=dt,signed_speed=speed,
                          counting_se=np.hypot(*(a["mean_se"] for a in runs))/dt,
                          runs=runs))
    return cases


def run():
    rng = np.random.default_rng(202611210430)
    car,cells,sites,hands = preparation()
    M = len(cells)
    columns = column_checks()
    diag,off = matrices(car,sites)
    channel = FreeChannel(diag,off)
    clock = PatchClock(rng)
    ref,residual = standing(M,0.,3,diag,off)
    cases = []
    for target in (0.,.30,.45,.65,.85,1.05,2.50):
        j = round(M*target/np.pi)
        for band in ((1,) if j == 0 else (1,2)):
            v,res = standing(M,j*np.pi/M,band,diag,off)
            residual = max(residual,res)
            # Fix the receiver to cell zero for every run; choose the port for visibility.
            address = int(np.argmax(np.abs(ref[:4]*v[:4])))
            line = record_line(rng,clock,channel,M,ref,v,address,TIMES)
            cases.append(dict(target=target,j_source=j,band=band,address_index=address,
                              receiver_port=int(sites[address]), **line))
    calibration = [r for r in cases if r['target']==.30]
    midpoint = np.mean([r['beat'] for r in calibration])
    cal_w = abs(np.diff([r['beat'] for r in calibration])[0])/2
    mass = midpoint-cases[0]['beat']
    qcal = calibration[0]['q_observed']
    flights = packet_flight(rng,clock,channel,M,qcal)
    speed = (flights[1]['signed_speed']-flights[0]['signed_speed'])/2
    c2 = speed*cal_w/qcal
    summary = []
    for lower,upper in zip(cases[3::2],cases[4::2]):
        q = lower['q_observed']
        w = midpoint-lower['beat']
        kg = np.sqrt(mass*mass+c2*q*q)
        para = mass+c2*q*q/(2*mass)
        symmetric = (lower['beat']+upper['beat'])/2-midpoint
        summary.append(dict(q=q,observed=w,kg=kg,kg_error=abs(kg/w-1),
                            parabola=para,parabola_error=abs(para/w-1),
                            midpoint_shift=symmetric))
    uniform_gap = np.diff(bloch(0,degree=False)[0][1:3])[0]/2
    alt_d,alt_off = matrices(car,sites,degree=False)
    control = FreeChannel(alt_d,alt_off)
    rest,_ = standing(M,0,1,diag,off)
    prepared = channel.launch(rest)
    overlap = abs((control.v.T @ prepared))**2
    control_survival = float(abs(np.dot(overlap,np.exp(-1j*control.e*.3)))**2)
    held = summary[:-1]
    checks = dict(max_held_kg_error=max(r['kg_error'] for r in held),
                  largest_held_parabola_error=max(r['parabola_error'] for r in held),
                  max_held_midpoint_drift=max(abs(r['midpoint_shift']) for r in held),
                  all_spatial_modes_recovered=all(r['j_source']==r['j_observed'] for r in cases),
                  max_source_vector_error=channel.max_source_error, max_mode_residual=residual,
                  max_clock_time_error=max(r['max_clock_time_error'] for r in cases))
    gates = dict(held_kg=checks['max_held_kg_error']<=.03,
                 held_parabola=checks['largest_held_parabola_error']>.05,
                 held_midpoint=checks['max_held_midpoint_drift']<=.03,
                 source=channel.max_source_error**2<1e-9,
                 spatial=checks['all_spatial_modes_recovered'],control=abs(uniform_gap)<1e-11)
    return dict(cells=cells,ports=sites.tolist(),hand_seams=hands,active_seams=len(car.seams),
                full_federation_ports=car.n,closure=columns,shots=SHOTS,packet_shots=PACKET_SHOTS,
                joint_word_order=WORDS.tolist(),ensemble_times=TIMES.tolist(),
                clock_calibration=clock.calibration,clock_leakage=float(clock.leakage),
                lines=cases,flights=flights,rest_frequency=mass,midpoint_relative=midpoint,
                calibration_q=qcal,calibration_frequency=cal_w,packet_speed=speed,
                inferred_c=float(np.sqrt(c2)),held_and_challenge=summary,
                control_half_gap=uniform_gap,control_rest_survival_at_point_three=control_survival,
                checks=checks,gates={key:bool(value) for key,value in gates.items()})


def canonical(obj):
    if isinstance(obj,dict):
        return {k:canonical(v) for k,v in obj.items()}
    if isinstance(obj,(tuple,list)):
        return [canonical(v) for v in obj]
    if isinstance(obj,(float,np.floating)):
        return float(format(float(obj), ".3g")) if abs(obj)<1e-9 else round(float(obj),10)
    return obj


def same_receipt(a,b):
    """Counts and identifiers are exact. Float checks allow native BLAS rounding."""
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same_receipt(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(same_receipt(x,y) for x,y in zip(a,b))
    if isinstance(a,float):
        return bool(np.isclose(a,b,rtol=0,atol=1e-8))
    return a==b


if __name__=='__main__':
    cli = argparse.ArgumentParser()
    cli.add_argument('--write',action='store_true')
    args = cli.parse_args()
    result = canonical(run())
    import gzip
    file = HERE/'receipt.json.gz'
    if args.write:
        with file.open('xb') as raw, gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as fp:
            fp.write(json.dumps(result,separators=(',',':'),allow_nan=False).encode())
    else:
        with gzip.open(file,'rt') as fp:
            assert same_receipt(result,json.load(fp))
    for key in ('rest_frequency','packet_speed','inferred_c','checks','gates','held_and_challenge',
                'control_rest_survival_at_point_three'):
        print(key,result[key])
