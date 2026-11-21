"""A moving rail analyzer for a freely opening, ordinarily captured bright pair.

A sharp head-on approach gives a marked point birth with equal rail amplitudes.
Let n count mature updates after 128 native opening services. The enrolled
station is at z=2n/7. The point birth supplies every axial momentum. Receiving
stations test the same normalized 25-address complex port on both rails and
compare the overlapping amplitudes at four phases. Each appointment uses
fresh launches. The station positions, aperture and analyzer phases are
supplied. No forces, material filters or phase plates act after formation.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from functools import lru_cache
from pathlib import Path

import numpy as np

from .. import law
from . import clock_inference as reader
from .formed_clock_readout import ordinary_birth

MEMBER = law.LAWS["torus"]
OPENING = MEMBER.maturation[2]
MUS = np.array(MEMBER.rail_contact_phases)/2
V = 2/7
WAITS = 7*np.arange(448,1121,112)
CENTRES = 2*WAITS//7
LENGTH = 16385
OFFSETS = np.arange(-12,13)
MENUS = ("low","high","bare","dephased")
OUTCOMES = ("+X","+Y","-X","-Y","marked_miss","refusal")
PHASES = np.array([1,1j,-1,-1j])
SHOTS = 2**34
GRID = 2**40
SEED = 89140724
ARCHIVE = Path(__file__).with_name("clock_ballistic_counts.npz")
FAMILY = len(MENUS)*len(WAITS)
ALPHA = 1e-6
ANGLE_GATE = .05
AMPLITUDE_GATE = .10
MIN_ADMISSION = 1e-6


def energy(k,mu):
    return np.arccos(np.cos(mu)*(3+np.cos(k))/4)


def saddles(mu):
    """The two positive roots of E'(k)=v. The reflected roots have sign −E."""
    c = np.cos(mu)
    a, b, d = 1-V**2, -6*V**2, 16*V**2/c**2-9*V**2-1
    discriminant = b*b-4*a*d
    if discriminant <= 0:
        raise ValueError("distinct nondegenerate saddles required")
    return np.arccos((-b + np.array([1,-1])*np.sqrt(discriminant))/(2*a))


def aperture(which):
    """Normalized receiver ket q; q†u(k) is its complex momentum response."""
    if which == "bare":
        return (OFFSETS == 0).astype(complex)
    index = ("low","high").index(which)
    roots = np.array([saddles(mu) for mu in MUS])
    nuisance = np.concatenate((-roots.ravel(), roots[:,1-index]))
    steering = lambda k: np.exp(1j*np.outer(OFFSETS,k))
    reject, _ = np.linalg.qr(steering(nuisance),mode="reduced")
    goal = steering([roots[:,index].mean()])[:,0]
    q = goal-reject@(reject.conj().T@goal)
    return q/np.linalg.norm(q)


def opening(k,mu):
    """Equal-slot projections of both signed branches after native birth opening.

    In coordinates u*S+w*FS, the contact coin has eigenvalues exp(-i*mu)
    on S and -exp(i*mu) on its complement. F exchanges the coordinates.
    The two signed poles follow by interpolation from two consecutive samples.
    """
    a = (3+np.cos(k))/4
    u, w = np.ones_like(k,dtype=complex), np.zeros_like(k,dtype=complex)
    d, e = np.exp(-1j*mu), -np.exp(1j*mu)
    for j in range(OPENING):
        theta = law.hop_angle(j,MEMBER.maturation)
        stay, hop = np.exp(-1j*theta)*np.array([np.cos(theta),1j*np.sin(theta)])
        x, y = e*u+(d-e)*(u+a*w), e*w
        u, w = stay*x+hop*y, stay*y+hop*x
    sample = u+a*w
    first = (e*u+(d-e)*sample)*a+e*w
    root = np.exp(-1j*energy(k,mu))
    positive = (first-root.conj()*sample)/(root-root.conj())
    return positive, sample-positive


def forecasts():
    """Order n^(-1/2) stationary amplitudes for the two filtered menus."""
    coefficients = np.empty((2,2),complex)
    frequencies = np.empty((2,2))
    opening_weights = np.empty((2,2,2))
    for r,mu in enumerate(MUS):
        roots = saddles(mu)
        opening_weights[r] = np.abs(opening(roots,mu))
        for which,k in enumerate(roots):
            E = energy(k,mu)
            curvature = (np.cos(mu)*np.cos(k)/4 - V**2*np.cos(E))/np.sin(E)
            p = np.vdot(aperture(MENUS[which]),np.exp(1j*OFFSETS*k))
            coefficients[which,r] = (p*opening(np.array([k]),mu)[0][0]
                                      *np.exp(-.25j*np.pi*np.sign(curvature))
                                      /np.sqrt(2*np.pi*abs(curvature)))
            frequencies[which,r] = E-V*k
    leading = (coefficients[:,None,:]*np.exp(-1j*WAITS[None,:,None]*frequencies[:,None,:])
               /np.sqrt(WAITS)[None,:,None])
    return frequencies, leading, opening_weights


def amplitudes(length, waits, centres):
    """Exact finite-strip contraction for low, high and bare common receivers."""
    if length < len(OFFSETS) or len(waits) != len(centres):
        raise ValueError("distinct aperture addresses and paired appointments required")
    k = 2*np.pi*np.fft.fftfreq(length)
    taps = np.array([aperture(name).conj() for name in MENUS[:3]])
    result = np.empty((3,len(waits),2),complex)
    for r,mu in enumerate(MUS):
        positive,negative = opening(k,mu)
        E = energy(k,mu)
        for i,(n,z) in enumerate(zip(waits,centres,strict=True)):
            factor = np.exp(-1j*n*E)
            wave = np.fft.ifft(positive*factor+negative*factor.conj())
            result[:,i,r] = taps@wave[(z+OFFSETS)%length]
    return result


def birth():
    """Local one-service chart, lifted into the long strip by finite propagation."""
    refused, selected = ordinary_birth(5,0.)
    natal = selected[0,0,0,0]
    b = float(np.vdot(natal,natal).real)
    if (np.max(abs(refused-[0,0,.5])) > 1e-13 or abs(b-.5)>1e-13 or
            np.linalg.norm(natal-law.EQUAL*np.vdot(law.EQUAL,natal))>1e-13):
        raise ArithmeticError("head-on birth, natal slots or wire census changed")
    return b


@lru_cache(maxsize=1)
def _source():
    if OPENING+WAITS[-1] >= 8192 or 2*(OPENING+WAITS[-1]+max(OFFSETS)) >= LENGTH:
        raise ValueError("enrolled age or strip wraps")
    h = amplitudes(LENGTH,WAITS,CENTRES)
    b = birth()
    ports = b/8*abs(h[:,:,0,None]+PHASES*h[:,:,1,None])**2
    washed = np.repeat(ports[:1].sum(axis=-1,keepdims=True)/4,4,axis=-1)
    all_ports = np.concatenate((ports,washed),axis=0)
    words = np.concatenate((all_ports,(b-all_ports.sum(-1))[...,None],
                            np.full((*all_ports.shape[:2],1),1-b)),axis=-1)
    if np.any(words < 0) or np.max(abs(words.sum(-1)-1)) > 1e-13:
        raise ArithmeticError("incomplete original-launch outcomes")
    return h,words


def source():
    return _source()[1].copy()


def weights(words):
    result = np.zeros(words.shape,dtype=np.int64)
    result[...,:4] = np.rint(GRID*words[...,:4]).astype(np.int64)
    result[...,5] = GRID//2
    result[...,4] = GRID//2-result[...,:4].sum(-1)
    return result


def draw(weight):
    result = np.empty_like(weight)
    for key in np.ndindex(weight.shape[:-1]):
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED,*key])))
        result[key] = rng.multinomial(SHOTS,weight[key]/GRID)
    return result


def _increment(first,last):
    """Unwrap the endpoint arc near its own measured center; retain its circle."""
    if not first or not last:
        return None
    lo,hi = reader._difference(first,last)
    shift = 2*np.pi*np.floor(((lo+hi)/2+np.pi)/(2*np.pi))
    return [lo-shift,hi-shift]


def _distance_bound(arc, forecast):
    if arc is None:
        return np.pi
    half = (arc[1]-arc[0])/2
    center_error = abs(np.angle(np.exp(1j*((arc[0]+arc[1])/2-forecast))))
    return min(np.pi,half+center_error)


def score(counts):
    h,words = _source()
    if (counts.shape != words.shape or counts.dtype != np.int64 or np.any(counts<0)
            or np.any(counts>SHOTS) or np.any(counts.sum(-1)!=SHOTS)):
        raise ValueError("complete original-start count words required")
    tail = reader._settings(FAMILY,FAMILY,ALPHA)
    # Per-port difference from nearest dyadic probability is at most half a unit.
    arcs = [[reader._arc(row,tail,1/GRID) for row in rows] for rows in counts]
    frequency,leading,opening_weights = forecasts()
    rates = frequency[:,1]-frequency[:,0]
    delta = WAITS-WAITS[0]
    inc = words[...,:4].sum(-1)
    checks = []
    records = {}
    for m,name in enumerate(MENUS):
        advances = [_increment(arcs[m][0],a) for a in arcs[m][1:]]
        records[name] = {"admission":inc[m].tolist(),"advance_intervals_radians":advances,
                         "admitted_counts":counts[m,:,:4].sum(-1).tolist()}
        if m < 2:
            error = abs(h[m]/leading[m]-1)
            bounds = [_distance_bound(a,dt*rates[m]) for a,dt in zip(advances,delta[1:],strict=True)]
            checks.extend([(name,"conditional_amplitude",np.max(error)<AMPLITUDE_GATE),
                           (name,"admission",np.min(inc[m])>MIN_ADMISSION),
                           (name,"every_signed_advance",max(bounds)<ANGLE_GATE)])
            # Competing zero, reversed sign and other positive saddle-line forecast.
            wrong = np.array([0.,-rates[m],rates[1-m]])*delta[-1]
            endpoint = advances[-1]
            excludes = (endpoint is not None and all(
                abs(np.angle(np.exp(1j*((endpoint[0]+endpoint[1])/2-x))))
                >(endpoint[1]-endpoint[0])/2+ANGLE_GATE for x in wrong))
            checks.append((name,"reject_competing_endpoint",excludes))
            records[name].update(leading_relative_errors=error.tolist(),
                                 mature_waits_outside_amplitude_gate=WAITS[
                                     np.any(error>=AMPLITUDE_GATE,axis=1)].tolist(),
                                 count_forecast_errors_upper=bounds)
    checks.append(("dephased","null_all_words",all(a is None for a in arcs[3])))
    return {"velocity":V,"mature_updates":WAITS.tolist(),"receiver_centres":CENTRES.tolist(),
            "per_service_rail_rates":rates.tolist(),"conditional_opening_magnitudes":opening_weights.tolist(),
            "counted_launches":int(counts.sum()),"family_failure_probability":ALPHA,
            "records":records,"failed_gates":[[n,g] for n,g,passed in checks if not passed]}


def experiment(path=ARCHIVE, *, acquire=False, seed_replay=False):
    path = Path(path)
    expected = weights(source())
    if acquire:
        with path.open("xb") as handle:
            np.savez_compressed(handle,weights=expected,counts=draw(expected))
    blob = path.read_bytes()
    digest = hashlib.sha256(blob).hexdigest()
    if not acquire and digest != ARCHIVE.with_suffix(".sha256").read_text(encoding="ascii").split()[0]:
        raise ValueError("ballistic acquisition checksum mismatch")
    with np.load(io.BytesIO(blob),allow_pickle=False) as archive:
        if sorted(archive.files)!=["counts","weights"]:
            raise ValueError("unexpected ballistic archive keys")
        archived = archive["weights"]
        if (archived.shape!=expected.shape or archived.dtype!=np.int64 or
                np.any(archived < 0) or np.any(archived > GRID) or np.any(archived.sum(-1)!=GRID)
                or np.any(archived[...,5]!=GRID//2) or np.max(abs(archived-expected))>1
                or np.max(abs(archived[...,:4]/GRID-source()[...,:4]))>1/GRID):
            raise ValueError("different six-outcome source")
        counts = archive["counts"]
    if np.any((archived==0) & (counts!=0)):
        raise ValueError("nonzero count on an impossible outcome")
    if seed_replay and not np.array_equal(counts,draw(archived)):
        raise ValueError("different draws on this multinomial backend")
    return {"acquisition_sha256":digest,**score(counts)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--acquire",action="store_true")
    parser.add_argument("--records",type=Path,default=ARCHIVE)
    parser.add_argument("--seed-replay",action="store_true")
    args = parser.parse_args()
    report = experiment(args.records,acquire=args.acquire,seed_replay=args.seed_replay)
    print(json.dumps(report,indent=2))
