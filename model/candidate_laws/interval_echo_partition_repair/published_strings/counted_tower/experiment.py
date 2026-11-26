"""Finite endpoint counts and simultaneous ten-frame small-ball bounds."""
import argparse
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
import json
import math
from pathlib import Path

import numpy as np
import mpmath as mp

from ...proper_balls.measure import FRAMES, PI
from .. import experiment as survey, finite_rods
from ..archive import source_densities
from ..memories import THETA
from ..native_pair import NativePair, endpoint
from matched_einstein_limit.frame import decode

HERE = Path(__file__).resolve().parent
STOCKS = ("blank", "ordinary", "held", "half_defect", "blind_pulse")
SHAPES = tuple(dict.fromkeys(FRAMES.values()))
BANKS_PER_LEVEL = len(STOCKS)*len(SHAPES)
ALPHA = 0.01
SEED = 26112719
GAIN = 4096
COEFFICIENT = 4*math.pi/15


def settings(level, banks_per_level=BANKS_PER_LEVEL):
    if not isinstance(level, int) or level < 0:
        raise ValueError("expected a nonnegative enrollment level")
    d = 32*4**level
    shots = 2**48*(level+1)**2
    bits = 40 + level.bit_length()
    beta = ALPHA/((level+1)*(level+2))
    # A dyadic emitter contributes at most one probability-grid unit.
    with localcontext() as ctx:
        ctx.prec=65
        correction=D(2*banks_per_level*100*(level+1)*(level+2)).ln()/(2*shots)
        error=math.nextafter(float(correction.sqrt()+D(2)**(-bits)),math.inf)
    rounding=2.**(-45)/(level+1)
    return dict(inverse_radius=d, coherent_repeats=d**4//GAIN, shots_per_bank=shots,
                probability_bits=bits, level_failure=beta, probability_radius=error,
                outward_entropy_margin=rounding,
                entropy_difference_radius=2*(GAIN*math.sqrt(2)*math.log(3)*error+rounding))


def _decimal(x):
    return D(x.numerator)/D(x.denominator)


def _encoded_probability(rho, collar, pulse_rho, shape, design):
    d = design["inverse_radius"]
    u = tuple(F(x).limit_denominator(3) for x in survey._frames()[shape])
    r = F(1,d)
    rod = finite_rods.read_ball(finite_rods.angle_gain(r)*(rho+collar), r, u)
    stress = pulse_rho*(1+u[1]**2+u[2]**2)
    with localcontext() as ctx:
        ctx.prec = 65
        work = 4*PI*_decimal(rod.weighted_over_4pi*stress)
    return phase_numerator(work, design)


def phase_numerator(work, design):
    """Read the dyadic terminal rate on the native work-driven seam for either material."""
    m, bits = (design[k] for k in ("coherent_repeats", "probability_bits"))
    with localcontext() as ctx:
        ctx.prec = 65
        t = work*(2*PI)/(D(5)/3).ln()/(D(3).sqrt()/4)
        x = 2*m*t
        if abs(x) > D(1)/4:
            raise ArithmeticError("stock exceeds the registered phase interval")
        sine = sum((-1)**k*x**(2*k+1)/D(math.factorial(2*k+1)) for k in range(17))
        cosine = 1+sum((-1)**k*x**(2*k)/D(math.factorial(2*k)) for k in range(1,17))
        p = D(3)/8 + D(3).sqrt()*sine/8 + (1-cosine)/8
        numerator = int((p*2**bits).to_integral_value())
    # The frozen dyadic table is independently checked on the actual native seam.
    native = NativePair(0, "right").state(THETA+m*float(t))
    measured = float(endpoint(native, "left")[1,1].real)
    if abs(measured-numerator/2**bits) > 2.**(-bits)+8e-15:
        raise ArithmeticError("native terminal read differs from its finite emitter")
    return numerator


def emit(level, density_records):
    cfg = settings(level)
    ordinary, collar = (F(x).limit_denominator(33) for x in density_records["held"])
    sources = dict(blank=(F(0),F(0),F(0)), ordinary=(ordinary,F(0),ordinary),
                   held=(ordinary,collar,ordinary+collar),
                   half_defect=(ordinary,collar/2,ordinary+collar),
                   blind_pulse=(ordinary,collar,ordinary))
    out = {}
    for j,shape in enumerate(SHAPES):
        out[shape] = {}
        for i,stock in enumerate(STOCKS):
            numerator = _encoded_probability(*sources[stock], shape, cfg)
            rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED, level, j, i])))
            # Aggregate independent binary banks in blocks with exactly representable counts.
            batch = 2**42
            chunks = cfg["shots_per_bank"]//batch
            hits = int(sum(map(int, rng.binomial(batch, numerator/2**cfg["probability_bits"], size=chunks))))
            out[shape][stock] = dict(ones=hits, zeros=cfg["shots_per_bank"]-hits, dyadic_one_numerator=numerator)
    return out


def _probability(shift):
    return 1/4 + math.sin(THETA+shift)**2/2


def read_count(ones, zeros, cfg):
    """Forecast the unamplified regional entropy, using only terminal counts and fixed settings."""
    if min(ones,zeros) < 0 or ones+zeros != cfg["shots_per_bank"]:
        raise ValueError("terminal record is incomplete")
    # Direct high-precision binary entropy independently implements the inverse.
    # Reserve enough digits before subtracting the two close regional entropies.
    with mp.workdps(4*len(str(cfg["inverse_radius"]))+55):
        frequency = mp.mpf(ones)/(ones+zeros)
        width = mp.mpf(cfg["probability_radius"])
        theta=mp.pi/6
        def population(x):
            return mp.mpf(1)/4+mp.sin(theta+x)**2/2
        def entropy(p):
            p=max(population(-0.125),min(population(0.125),p))
            phase=mp.asin(mp.sqrt(2*p-mp.mpf(1)/2))-theta
            q=population(phase/cfg["coherent_repeats"])
            def h(v):
                return -v*mp.log(v)-(1-v)*mp.log1p(-v)
            return (h(q)-h(mp.mpf(3)/8))*cfg["inverse_radius"]**4
        margin=cfg["outward_entropy_margin"]
        return [float(entropy(frequency)),
                math.nextafter(float(entropy(frequency-width)-margin),-math.inf),
                math.nextafter(float(entropy(frequency+width)+margin),math.inf)]


def infer(counts, cfg, rod_rows):
    readings = {shape: {stock: read_count(rec["ones"], rec["zeros"],cfg)
                         for stock,rec in banks.items()} for shape,banks in counts.items()}
    result = {}
    for label, upper, lower in [(stock,stock,"blank") for stock in STOCKS[1:]] + [("collar","held","ordinary")]:
        s,a,t,g,intervals = [],[],[],[],[]
        for shape in FRAMES.values():
            hi,lo = readings[shape][upper],readings[shape][lower]
            shat = hi[0]-lo[0]
            interval = (hi[1]-lo[2],hi[2]-lo[1])
            row = rod_rows[shape]
            def difference(key):
                return row[upper][key]-row[lower][key]
            s.append(shat)
            intervals.append(interval)
            a.append(difference("area_over_ell4"))
            t.append(difference("action_stress"))
            g.append(difference("metric_charge"))
        bulk = max(max(abs(v-2*math.pi*COEFFICIENT*t[i]) for v in ends) for i,ends in enumerate(intervals))
        geometry = max(abs(ai+COEFFICIENT*gi) for ai,gi in zip(a,g,strict=True))
        stationarity = max(max(abs(v+a[i]) for v in ends) for i,ends in enumerate(intervals))
        guaranteed_discrepancy = max(max(0.,ends[0]+a[i],-a[i]-ends[1]) for i,ends in enumerate(intervals))
        tensor = list(decode([x/(2*math.pi*COEFFICIENT) for x in s]))
        result[label] = dict(entropy_tensor_estimate=tensor, independent_action_tensor=list(decode(t)),
                             area_tensor_estimate=list(decode([-x/COEFFICIENT for x in a])),
                             entropy_tensor_sampling_radius=(17/8)*cfg["entropy_difference_radius"]/(2*math.pi*COEFFICIENT),
                             bulk_upper=bulk, geometry_upper=geometry, stationarity_upper=stationarity,
                             stationarity_lower=guaranteed_discrepancy,
                             theorem_component_bound=(17/8)*(bulk+geometry+stationarity)/COEFFICIENT,
                             ten_entropy_intervals_over_ell4=intervals)
    return result


def measure():
    densities = source_densities()
    board = survey._board(densities, finite=True)
    levels = []
    for j,d in enumerate(survey.SCALES):
        cfg = settings(j)
        raw = emit(j,densities)
        rods = {r["shape"]:r["trials"] for r in board if r["inverse_ell"] == d}
        levels.append(dict(design=cfg, banks=raw, comparisons=infer(raw,cfg,rods)))
    return dict(schema="oph.counted_string_tower.v1", simultaneous_coverage=1-ALPHA,
                seed=SEED, signed_frame_order=list(FRAMES.keys()), levels=levels)


def verify():
    data = measure()
    old = json.loads((HERE/"results.json").read_text())
    if any(data[k] != old[k] for k in ("schema","simultaneous_coverage","seed","signed_frame_order")):
        raise ValueError("counted source or frame convention changed")
    for a,b in zip(data["levels"], old["levels"], strict=True):
        if a["banks"] != b["banks"] or a["design"].keys() != b["design"].keys():
            raise ValueError("raw positive or zero counts or frozen sampling plan changed")
        for key,val in a["design"].items():
            if isinstance(val,int):
                if type(b["design"][key]) is not int or val != b["design"][key]:
                    raise ValueError("integer sampling plan changed")
            elif not math.isclose(val,b["design"][key],rel_tol=5e-15,abs_tol=0):
                raise ValueError("sampling error allocation changed")
        if a["comparisons"].keys() != b["comparisons"].keys():
            raise ValueError("counted stock comparisons changed")
        for stock,record in a["comparisons"].items():
            previous = b["comparisons"][stock]
            if record.keys() != previous.keys():
                raise ValueError("incomplete counted tensor comparison")
            for key, value in record.items():
                np.testing.assert_allclose(value,previous[key],rtol=3e-10,atol=1e-11)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    args = parser.parse_args()
    if args.write:
        (HERE/"results.json").write_text(json.dumps(measure(),indent=2,allow_nan=False)+"\n")
    else:
        verify()
        print("Replayed 105 terminal banks on the common seven-cut, ten-frame material survey.")


if __name__ == "__main__":
    main()
