"""Independent repeat of the return experiment with an equal-window beat measurement."""
import argparse
from pathlib import Path

import numpy as np

from .shell_source import check
from .shells import protocol, read_hand, save
from .shell_spectroscopy import response, run as acquire, times_interval


def advance(z, error):
    if np.min(abs(z))<=.5 or np.max(abs(np.angle(z[1:]*z[:-1].conj())))>=np.pi-2*error:
        return None
    return np.unwrap(np.angle(z))


def over_duration(interval, duration, eta):
    check(duration>2*eta,"clock interval unresolved")
    return times_interval(interval,(1/(duration+2*eta),1/(duration-2*eta)))


def window_read(rest, line, hand, width, eta, component_error):
    """Return a bounded squared-frequency shift, or an unresolved beat."""
    check(rest.ndim==1 and line.shape==rest.shape and np.isfinite(rest).all() and np.isfinite(line).all()
          and (hand is None or (hand.shape==rest.shape and np.isfinite(hand).all()))
          and isinstance(width,(int,np.integer)) and 0<2*width<=len(rest)
          and 0<eta<np.pi/2 and np.isfinite(component_error) and component_error>0,
          "invalid beat measurement")
    lag = len(rest)-width
    if hand is None or min(hand[lag]-hand[0],hand[-1]-hand[0])<=2*eta:
        return {"status":"unresolved hand duration"},None
    r = advance(rest,eta)
    ell = advance(line,eta)
    if r is None or ell is None:
        return {"status":"unresolved single-frequency return"},None
    beat = line*rest.conj()
    means = np.array([beat[:width].mean(),beat[-width:].mean()])
    radius = np.sqrt(2)*component_error
    if min(abs(means))<=radius:
        return {"status":"unresolved visibility"},None
    phase_budget = sum(np.arcsin(radius/abs(means)))
    raw = np.angle(means[1]*means[0].conj())
    coarse = (ell[lag]-ell[0])-(r[lag]-r[0])
    beat_advance = raw+2*np.pi*np.rint((coarse-raw)/(2*np.pi))
    if abs(beat_advance-coarse)+4*eta+phase_budget>=np.pi:
        return {"status":"unresolved winding"},None
    increment = over_duration((beat_advance-phase_budget,beat_advance+phase_budget),
                                hand[lag]-hand[0],eta)
    rest_advance = r[-1]-r[0]
    rest_interval = over_duration((rest_advance-2*eta,rest_advance+2*eta),hand[-1]-hand[0],eta)
    twice_rest_plus_increment = tuple(2*r0+k for r0,k in zip(rest_interval,increment,strict=True))
    shift = times_interval(increment,twice_rest_plus_increment)
    return {"status":"resolved","kinetic_interval":list(increment),"rest_interval":list(rest_interval),
            "squared_shift_interval":list(shift),"minimum_window_visibility":float(min(abs(means)))},shift


def run(base, spec, repeat):
    spec = {**spec,"seed":repeat["seed"]}
    data,endpoint = acquire(spec,base)
    width,alpha,shots = repeat["window"],repeat["alpha"],spec["shots"]
    settings = len(spec["loads"])*len(spec["reference_terms"])*len(spec["law_members"])
    point_marginals = 8*(spec["clock_steps"]+1)*settings
    quantization = 1e-8
    eta = np.arcsin(np.sqrt(2)*(quantization+np.sqrt(2*np.log(4*point_marginals/alpha)/shots)))
    # Each component of the product reads four independent +/-1 sample means.
    # A single trial changes a window mean by at most 2/(width*shots).
    window_marginals = 4*settings
    eps = quantization+np.sqrt(8*np.log(4*window_marginals/alpha)/(width*shots))
    out = dict(window_component_error=float(eps),point_angle_error=float(eta),members={})
    for member in spec["law_members"]:
        for reference in spec["reference_terms"]:
            tag = f"{member}_{reference}"
            clocks = sum(data[f"{tag}_{q}_{mode}"].sum(axis=(0,2),dtype=np.int64)
                         for q in spec["loads"] for mode in spec["modes"])
            hand = read_hand(clocks,allow_unresolved=True)
            readings, bounds = [], []
            for q in spec["loads"]:
                rest,line = (response(data[f"{tag}_{q}_{mode}"],shots,spec["clock_steps"]+1) for mode in spec["modes"])
                reading, bound = window_read(rest,line,hand,width,eta,eps)
                readings.append(reading)
                bounds.append(bound)
            band = None if any(b is None for b in bounds) else [bounds[1][0]-bounds[0][1],bounds[1][1]-bounds[0][0]]
            out["members"][tag] = dict(loads=readings,measured_difference_interval=band,
                common_squared_shift="excluded" if band and (band[0]>0 or band[1]<0) else "unresolved",
                endpoint_interval=endpoint["members"][tag]["measured_difference_interval"])
    return data,out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    fixed = [protocol(f) for f in ("shell_protocol.json","shell_spectroscopy_protocol.json","shell_beat_protocol.json")]
    counts,summary = run(*fixed)
    save("shell_beat",counts,summary,args.output)


if __name__=="__main__":
    main()
