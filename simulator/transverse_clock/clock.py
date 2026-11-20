"""Acquire independent cohorts and read the winding of a transverse pulse."""
import math

import numpy as np

from .equipped.chart import stations, deliver, reconstruct
from .kernel import packet, step, marginal, counters, require


def center_for(size,values):
    return tuple({'L/4':size//4,'L/8':size//8,'0':0}[c] if isinstance(c,str) else c for c in values)


def acquire(spec,cohort):
    require(cohort in ('timing','rest','moving'),'unregistered detector cohort')
    axis = 0 if cohort=='timing' else 1
    raw = {}
    for size in spec['sizes']:
        times = np.arange(33)*size//8
        order,cal,stamps,lap = stations(size,axis,times)
        raw[f'{size}_cal'] = cal
        raw[f'{size}_stamps'] = stamps
        raw[f'{size}_lap'] = np.asarray(lap)
        for label,center in spec['panels'].items():
            if axis==1 and label=='reference':
                continue
            if axis==0 and label in ('rest','short_rest'):
                continue
            if axis==1 and (cohort=='rest')!=(label in ('rest','short_rest')):
                continue
            modes,shift,state,_ = packet(size,center_for(size,center),spec['halfwidth'])
            previous = 0
            rows = []
            for t in times:
                for _ in range(int(t-previous)):
                    state = step(state,shift)
                previous = t
                p = marginal(state,modes,size,axis)
                if label=='erased' and axis==1:
                    # Exact uniform random translation in y of the source ensemble.
                    p = np.full(size,p.sum()/size)
                rows.append(counters(p,spec['shots']))
            raw[f'{size}_{label}'] = deliver(size,axis,order,np.asarray(rows))
    return raw


def read(rows,chart,shots,gate):
    coords,times,lap = chart
    require(rows.dtype==np.int32 and rows.shape==(len(times),len(coords)) and
            np.all(rows>=0) and np.all(rows.sum(axis=1)==shots),'invalid pulse counts')
    z = rows @ np.exp(2j*np.pi*coords/lap)/shots
    # The sum of all roots of unity on the checked complete chart vanishes exactly.
    z[np.ptp(rows,axis=1)==0] = 0
    visibility = min(abs(z))
    result = {'minimum_visibility':float(visibility)}
    if visibility<gate['min_visibility']:
        return result
    angles = np.unwrap(np.angle(z))
    if max(abs(np.diff(angles)))>=gate['max_increment']:
        return result
    duration = times[-1]-times[0]
    scale = lap/(2*np.pi*duration)
    # Cumulative midpoint counts have discrepancy <= half a count at each cut.
    # Summation by parts on the circle bounds the complex moment error by pi/shots.
    radius = 2*math.asin(min(1.,np.pi/(shots*visibility)))*scale
    result.update(speed=float((angles[-1]-angles[0])*scale),radius=float(radius),
                  turns=float((angles[-1]-angles[0])/(2*np.pi)))
    return result


def chart_for(raw,size):
    return reconstruct(*(raw[f'{size}_{field}'] for field in ('cal','stamps','lap')))


def predict(spec,x,rests):
    forecast = {}
    for size in spec['sizes']:
        cx = chart_for(x,size)
        cy = chart_for(rests,size)
        def reading(raw,chart,label):
            return read(raw[f'{size}_{label}'],chart,spec['shots'],spec['clock_gate'])
        signal = reading(x,cx,'reference')
        for label,center in spec['panels'].items():
            if label in ('reference','rest','short_rest'):
                continue
            vx = reading(x,cx,label)
            rest = reading(rests,cy,'short_rest' if label=='short' else 'rest')
            v,c = abs(vx['speed']),signal['speed']
            lo = max(0.,v-vx['radius'])/(c+signal['radius'])
            hi = (v+vx['radius'])/(c-signal['radius'])
            require(0<=lo<=hi<1 and rest['speed']>rest['radius'],'unresolved calibration')
            key = f'{size}_{label}'
            forecast[key] = dict(speed=vx['speed'],reference=c,rest=rest['speed'],
                                 rest_radius=rest['radius'],ratio=math.sqrt(1-(v/c)**2),
                                 interval=[math.sqrt(1-hi**2),math.sqrt(1-lo**2)])
    return forecast


def finite_forecasts(spec,rest):
    targets = {}
    for size in spec['sizes']:
        coords,times,_ = chart_for(rest,size)
        for label,center in spec['panels'].items():
            if label in ('reference','rest','short_rest'):
                continue
            modes,_,initial,energy = packet(size,center_for(size,center),spec['halfwidth'])
            rows = []
            for t in times:
                wave = initial*np.exp(-1j*t*energy)[...,None]
                rows.append(counters(marginal(wave,modes,size,1),spec['shots'])[coords.astype(int)])
            targets[f'{size}_{label}'] = np.asarray(rows)
    return targets


def open_pulse(spec,forecast,targets,y,rest):
    result = {}
    for size in spec['sizes']:
        chart = chart_for(y,size)
        require(all(np.array_equal(a,b) for a,b in zip(chart,chart_for(rest,size))),
                'moving detector differs from enrolled clock or chart')
        for label,center in spec['panels'].items():
            key = f'{size}_{label}'
            if key not in forecast:
                continue
            f = forecast[key]
            counts = y[key]
            observation = read(counts,chart,spec['shots'],spec['clock_gate'])
            distances = abs(counts.astype(np.int64)-targets[key]).sum(axis=1)/(2*spec['shots'])
            observation['max_coherent_forecast_tv'] = float(max(distances))
            if 'speed' in observation:
                rate = observation['speed']
                lo = (rate-observation['radius'])/(f['rest']+f['rest_radius'])
                hi = (rate+observation['radius'])/(f['rest']-f['rest_radius'])
                predlo,predhi = f['interval']
                tolerance = spec['clock_gate']['longwave_tolerance']
                observation.update(ratio=rate/f['rest'],interval=[lo,hi])
                if lo>=predhi-tolerance and hi<=predlo+tolerance:
                    status = 'within tolerance'
                elif hi<predlo-tolerance or lo>predhi+tolerance:
                    status = 'outside tolerance'
                else:
                    status = 'unresolved'
            else:
                status = 'unresolved'
            observation['decision'] = status
            result[key] = observation
    return result
