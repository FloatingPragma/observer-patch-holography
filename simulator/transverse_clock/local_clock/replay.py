"""Read the frozen local winding calibration, moving aperture and unwritten-tag control."""
from pathlib import Path

import numpy as np

from tools.strict_json import load
from .. import clock
from ..receipts import HERE as PARENT, archive, receipt
from ..equipped.chart import addressed, station_times
from ..kernel import require
from .source import observe

HERE = Path(__file__).resolve().parent


def report(size,x,times,rows):
    home = (0,0,0)
    samples = []
    for sent in (0,11,39):
        at,d,request = addressed(home,[(0,x)],size,sent)
        remote = 2*(sent+d)+23
        end,back,echo = addressed(at,[(0,-x)],size,(request,remote))
        require(end==home and echo[0]==sent,'missing clock calibration exchange')
        samples.append((sent,echo[1],sent+d+back))
    cal = np.asarray(samples,dtype=np.int64)
    stamp = np.r_[2*np.asarray(times)+23].astype(np.int64)
    end,_,payload = addressed(at,[(0,-x)],size,(cal,stamp,rows.copy()))
    require(end==home,'missing local return')
    return payload


def local_times(cal,stamp):
    _, times = station_times(cal, stamp)
    return np.asarray(times, dtype=float)


def peak(rows,times,shots):
    require(rows.dtype==np.int64 and rows.shape==(len(times),9) and len(times)>=3 and
            0<shots<=np.iinfo(np.int64).max//9 and np.all((rows>=0)&(rows<=shots)) and
            np.all(rows.sum(axis=1)==shots) and np.all(np.isfinite(times)) and
            np.all(np.diff(times)>0),'invalid local aperture record')
    i = int(np.argmax(rows[:,1]))
    if i==0 or i==len(times)-1 or rows[i,1]<1000 or np.count_nonzero(rows[:,1]==rows[i,1])!=1:
        return None
    return dict(time=float(times[i]),counts=int(rows[i,1]),fraction=float(rows[i,1]/shots))


def calibrate(protocol,base,root):
    menu = {**base,'sizes':protocol['sizes'],'panels':{'reference':[16,0],**protocol['center']}}
    timing = archive(root/'timing.npz',clock.acquire(menu,'timing'))
    raw,targets,forecasts = {},{},{}
    for size in protocol['sizes']:
        chart = clock.chart_for(timing,size)
        signal = clock.read(timing[f'{size}_reference'],chart,base['shots'],base['clock_gate'])['speed']
        shots = protocol['shots_per_four_volume']*size**4
        for label,center in protocol['center'].items():
            key = f'{size}_{label}'
            xy = clock.center_for(size,center)
            schedule = [0,*[round(j*size/32) for j in range(40,97)]]
            cal,stamp,counts = report(size,0,schedule,observe(size,(0,xy[1]),0,schedule,protocol))
            raw[key+'_cal'],raw[key+'_stamp'],raw[key] = cal,stamp,counts
            rest = peak(counts,local_times(cal,stamp),shots)
            require(rest is not None,'rest winding unresolved')
            vx = clock.read(timing[key],chart,base['shots'],base['clock_gate'])['speed']
            require(abs(vx)<signal,'no subluminal scalar forecast')
            ratio = np.sqrt(1-(vx/signal)**2)
            expected = rest['time']/ratio
            position = round(vx*expected)%size
            times = [0,*[round(j*expected/32) for j in range(20,53)]]
            targets[key] = observe(size,xy,position,times,protocol,analytic=True)
            forecasts[key] = dict(rest=rest,vx=vx,signal=signal,ratio=ratio,
                                  predicted_lap=expected,station=position,ticks=times)
    archive(root/'rest.npz',raw)
    targets = archive(root/'finite_forecast.npz',targets)
    receipt(root/'forecast.json',forecasts)
    return forecasts,targets


def replay():
    protocol, base = load(HERE/'protocol.json'), load(PARENT/'protocol.json')
    root = HERE/'records'
    forecasts, targets = calibrate(protocol, base, root)
    raw,result = {},{}
    for size in protocol['sizes']:
        shots = protocol['shots_per_four_volume']*size**4
        for label,center in protocol['center'].items():
            key = f'{size}_{label}'
            f = forecasts[key]
            xy = clock.center_for(size,center)
            pos,ticks = f['station'],f['ticks']
            cal,stamp,rows = report(size,pos,ticks,observe(size,xy,pos,ticks,protocol))
            raw[key+'_cal'],raw[key+'_stamp'],raw[key] = cal,stamp,rows
            # Removing the tag changes the coherent Hilbert space. Resolve this control on its own.
            plain = observe(size,xy,pos,ticks,{**protocol,'tag_modulus':1})
            raw[key+'_unwritten'] = plain
            obs = peak(rows,local_times(cal,stamp),shots)
            decision = 'unresolved' if obs is None else (
                'within tolerance' if abs(obs['time']/f['predicted_lap']-1)<=.05 else 'outside tolerance')
            result[key] = dict(observed=obs,decision=decision,
                    maximum_count_difference=int(abs(rows-targets[key]).max()),
                    max_tag_write_effect_counts=int(abs(plain[:,0]-rows[:,:-1].sum(axis=1)).max()))
    archive(root/'opened.npz',raw)
    receipt(root/'result.json',result)
    return result
