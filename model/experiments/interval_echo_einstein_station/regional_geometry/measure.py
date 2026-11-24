"""Test an operational fixed-volume tangent with independent return and free-pulse rulers."""
import argparse
import hashlib
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.experiments.contact import radar
from model.experiments.interval_echo_einstein_station import refinement
from model.experiments.interval_echo_null_balance.regional import verify
from model.experiments.interval_echo_metric_response.measure import hand_columns
from model.experiments.material_inertia.dark_collar import label

HERE = Path(__file__).resolve().parent
PAIRS = list(combinations(range(6), 2))
FORWARD = np.roll(np.eye(3), 1, axis=1)
DIRECTIONS = np.stack([(FORWARD.T-FORWARD)/36, (2*np.eye(3)-FORWARD-FORWARD.T)/36])
NAMES = ("first_return", "free_pulse")
WORD_SHOTS = 2**24
EPS = .25


def hand_times():
    """Decode native endpoint probabilities to a signed, unwrapped material hand."""
    law = candidate.native_law
    exchange, agreement = law.local_term_parts()
    h = exchange + candidate.KAPPA*agreement
    pair = np.eye(law.N_SRC)[:, law.SRC.index((1, 1))]
    arms = np.column_stack([np.kron(label(2, s), pair) for s in (1, -1)])
    sharp = law.SRC.index((2, 0))
    prepared = arms.T @ law.echo()[:, sharp]
    contact = expm(-1.5j*np.pi*h)
    analyses = np.stack([(law.echo() @ arms)[sharp], (contact @ law.echo() @ arms)[sharp]])
    phases = np.exp(-1j*np.outer([32., 33.], np.arange(2*radar.TICKS+1)*radar.DT))
    yes = abs(analyses @ (prepared[:,None]*phases))**2
    xy = 2*yes-1
    angles = np.unwrap(np.arctan2(xy[1],xy[0]))
    return angles-angles[0]


def read_word_laws(level, max_error_check, *, record_independent=None):
    config = json.loads((HERE.parent/"refinement.json").read_text())
    source = dict(source_ports=config["source_ports"], hand_seam=config["hand_seam"],
                  source_word_energies=config["energies"])
    base, fixed, _ = refinement.station(level, source, tuple(config["protected_tag_seam"]))
    sites = list(fixed) + [4, fixed[2]+1]
    hand = hand_times()
    half_return = hand[np.add.outer(np.arange(1,radar.TICKS+1), np.arange(1,radar.TICKS+1))]/2
    pulse_time = hand[1:radar.TICKS+1]
    laws, alternate_laws = [], []
    discrepancy = 0.

    def moments(legs, intensity):
        joint = radar.round_trip_law(legs)
        records = []
        for i,j in PAIRS:
            np.testing.assert_allclose(joint[i,j].sum(),1,rtol=0,atol=2e-13)
            returned = joint[i,j,:-1,:-1]
            pulses = (intensity[i,j]+intensity[j,i])/(2*radar.TICKS)
            records.append([[returned.sum(), np.sum(half_return*returned)],
                            [pulses.sum(), pulse_time @ pulses]])
        return records

    for first, second in product(config["energies"], repeat=2):
        h = base.tolil()
        a, _, c = source["source_ports"]
        h[a,a] -= first
        h[c,c] -= second
        h = h.tocsr()
        legs, intensity = refinement.sparse_hit_laws(h, sites, PAIRS, return_intensity=True)
        if max_error_check:
            independent, _ = refinement.projector_check(h, sites, PAIRS)
            discrepancy = max(discrepancy,float(np.max(abs(legs-independent))))
            # With the repeated detector turned off, compute the second ruler
            # directly by powers of an independently obtained dense unitary.
            unitary = expm(-1j*radar.DT*h.toarray())
            releases = np.eye(h.shape[0])[:,sites]
            free = []
            for _ in range(radar.TICKS):
                releases = unitary @ releases
                free.append(abs(releases[sites])**2)
            free = np.stack(free, axis=-1)
            alternate_laws.append(moments(independent, free))
            if record_independent is not None:
                record_independent(first, second, independent, free, hand)
        laws.append(moments(legs,intensity))
    layout = (3,3,len(PAIRS),2,2)
    alternate = np.array(alternate_laws).reshape(layout) if alternate_laws else None
    return sites, np.array(laws).reshape(layout), discrepancy, alternate


def simplex(ranges, dr, sites):
    indices = np.ix_(sites, sites)
    sq = ranges[indices]**2
    dsq = 2*ranges[indices]*dr[indices]
    gram = (sq[0,1:,None]+sq[None,0,1:]-sq[1:,1:])/2
    dgram = (dsq[0,1:,None]+dsq[None,0,1:]-dsq[1:,1:])/2
    if np.min(np.linalg.eigvalsh(gram)) <= 0:
        raise ArithmeticError("these distances do not span a positive simplex")
    measure = np.sqrt(np.linalg.det(gram))/(6 if len(sites)==4 else 2)
    dlog = np.trace(np.linalg.solve(gram, dgram))/2
    return measure, dlog, gram


def geometry(word_law, mixture, tangent):
    mixed = np.einsum("ij,ijpk->pk", mixture, word_law)
    change = np.einsum("ij,ijpk->pk", tangent, word_law)
    values = mixed[:,1]/mixed[:,0]
    slopes = (change[:,1]-values*change[:,0])/mixed[:,0]
    d, dd = np.zeros((6,6)), np.zeros((6,6))
    for (i,j), x, dx in zip(PAIRS, values, slopes, strict=True):
        d[i,j] = d[j,i] = x
        dd[i,j] = dd[j,i] = dx
    volume, dv, gram = simplex(d, dd, [0,1,2,3])
    face = [simplex(d, dd, part)[:2] for part in combinations(range(4),3)]
    area = sum(a for a,_ in face)
    da = sum(a*diff for a,diff in face)/area
    # Locate the two held stations from their four anchor ranges. Their mutual
    # distance is withheld and the reference records determine the axes.
    x = np.linalg.cholesky(gram)
    b = (d[0,1:4,None]**2 + d[None,0,4:6]**2 - d[1:4,4:6]**2)/2
    held = np.linalg.solve(x, b)
    predicted = np.linalg.norm(held[:,0]-held[:,1])
    return dict(volume=float(volume), surface_area=float(area), dlog_volume=float(dv),
                dlog_area=float(da), held_observed=float(d[4,5]), held_predicted=float(predicted),
                held_relative_error=float(abs(predicted/d[4,5]-1)),
                anchor_radius_relative_error=float(np.max(abs(np.linalg.norm(held,axis=0)/d[0,4:6]-1))))


def evaluate(laws, ref, coefficients):
    tangent = np.einsum("i,ijk->jk", coefficients, DIRECTIONS)
    result = {}
    for axis,name in enumerate(NAMES):
        data = laws[:,:,:,axis,:]
        d = geometry(data, ref, tangent)
        step = .02
        lo = geometry(data, ref-step*tangent, tangent)
        hi = geometry(data, ref+step*tangent, tangent)
        d["central_dlog_volume"] = float(np.log(hi["volume"]/lo["volume"])/(2*step))
        d["central_dlog_area"] = float(np.log(hi["surface_area"]/lo["surface_area"])/(2*step))
        result[name] = d
    return result


def run():
    # This also checks the competing restored-energy contact on all 99-state load columns.
    hand_columns()
    turns = json.loads((HERE.parent/"refinement.json").read_text())["source_units"]
    ref = np.array(turns[1])/36
    rng = np.random.default_rng(26112342)
    screens = []
    transferred = None
    for level in (0,1,2):
        sites, laws, error, alternate = read_word_laws(level, level==0)
        slopes = np.array([geometry(laws[:,:,:,0,:],ref,d)["dlog_volume"] for d in DIRECTIONS])
        if np.linalg.norm(slopes) < 1e-11:
            raise ArithmeticError("calibration has no informative first-return volume response")
        coeff = np.array([slopes[1], -slopes[0]])/np.linalg.norm(slopes)
        if coeff[0] < 0:
            coeff = -coeff
        if level==0:
            transferred = coeff
        local = evaluate(laws, ref, coeff)
        common = evaluate(laws, ref, transferred)
        tangent = np.einsum("i,ijk->jk", transferred, DIRECTIONS)
        counts = {}
        read_entropy = []
        for name,sign in (("minus",-1),("reference",0),("plus",1)):
            prior = ref + sign*EPS*tangent
            np.testing.assert_allclose(prior.sum(axis=0),1/3,atol=2e-16)
            np.testing.assert_allclose(prior.sum(axis=1),1/3,atol=2e-16)
            counts[name] = rng.multinomial(WORD_SHOTS,prior.ravel()).reshape(3,3)
        score = -np.log(counts["reference"]/WORD_SHOTS)
        measured = (counts["plus"]-counts["minus"])/(2*EPS*WORD_SHOTS)
        for name in ("minus","plus"):
            p = counts[name]/WORD_SHOTS
            read_entropy.append(np.sum(p*(score-np.sum(score*p))**2)/(4*EPS**2*WORD_SHOTS))
        pref = counts["reference"]/WORD_SHOTS
        read_entropy.append((np.sum(measured**2/pref)-measured.sum()**2)/WORD_SHOTS)
        # Without the resident, or with both occupied energies restored, every
        # word has the same flight law (first=second=0 up to scalar phases).
        blind = evaluate(np.broadcast_to(laws[0,0],laws.shape), ref, transferred)
        slope_disagreement = None
        if alternate is not None:
            alt = evaluate(alternate, ref, transferred)
            slope_disagreement = max(abs(common[n][k]-alt[n][k])
                                     for n in NAMES for k in ("dlog_volume","dlog_area"))
            if slope_disagreement >= 1e-12:
                raise ArithmeticError("independent area or volume slopes disagree")
        screens.append(dict(cells=20*4**level, stations=sites, word_moments=laws.tolist(),
                            moments_layout="word x word x unordered station pair x ruler x [success, hand-weighted success]",
                            projector_renewal_max_error=error if level==0 else None,
                            independent_area_volume_slope_disagreement=slope_disagreement,
                            calibration_gradient=slopes.tolist(), local_coefficients=coeff.tolist(),
                            local=local, transferred=common,
                            direct_entropy_slope=float(-np.sum(tangent*np.log(ref))),
                            measured_source_charge_slope=float(np.sum(score*measured)),
                            source_slope_standard_error=float(np.sqrt(sum(read_entropy))),
                            word_counts={n:v.tolist() for n,v in counts.items()},
                            blind_control_maximum=max(abs(d[k]) for d in blind.values()
                                                      for k in ("dlog_volume","dlog_area"))))
    return dict(design_sha256=hashlib.sha256((HERE/"DESIGN.md").read_bytes()).hexdigest(),
                shots_per_source=WORD_SHOTS, seed=26112342, pairs=[list(p) for p in PAIRS],
                reference=ref.tolist(), directions=DIRECTIONS.tolist(),
                transferred_coefficients=transferred.tolist(), screens=screens)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    else:
        verify(json.loads((HERE/"results.json").read_text()),result)
    for s in result["screens"]:
        print(s["cells"], s["local_coefficients"], s["direct_entropy_slope"],
              {n:(v["dlog_volume"], v["dlog_area"], v["held_relative_error"])
               for n,v in s["transferred"].items()})


if __name__=="__main__":
    main()
