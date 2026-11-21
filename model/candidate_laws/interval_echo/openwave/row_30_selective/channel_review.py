"""Independent analysis of the autonomous-channel raw receipt and competing dispersions."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from .source import HERE

WAVES = ("kg_velocity", "kg_frequency", "parabola_velocity", "parabola_frequency",
         "cosine_velocity", "lattice_kg_velocity", "native_band")


def predict(mass, freq, speed, anchor, q):
    return np.stack(np.broadcast_arrays(
        np.sqrt(mass**2 + speed*freq/anchor*q**2),
        np.sqrt(mass**2 + (freq**2-mass**2)*(q/anchor)**2),
        mass + speed*q**2/(2*anchor),
        mass + (freq-mass)*(q/anchor)**2,
        mass + speed*(1-np.cos(q))/np.sin(anchor),
        np.sqrt(mass**2 + 4*speed*freq/np.sin(anchor)*np.sin(q/2)**2),
        4*np.sqrt(3-np.sqrt(6+2*np.cos(q))),
    ), axis=-1)


def spatial(r, cells):
    counts = np.asarray(r["spatial_counts"], dtype=float)
    counts /= counts.sum(axis=1)[:, None]
    signal = np.diff(counts, axis=0).reshape(cells, 4)
    n = np.arange(cells)
    def fit(q):
        X = np.column_stack([np.cos(q*n), np.sin(q*n)])
        return float(np.linalg.norm(signal - X @ np.linalg.lstsq(X, signal, rcond=None)[0]))
    step = np.pi/cells
    errors = np.array([fit(j*step) for j in range(cells)])
    order = np.argsort(errors)
    j, rival = order[:2]
    lo, hi = max(0, (j-.51)*step), min(np.pi, (j+.51)*step)
    optimum = minimize_scalar(fit, bounds=(lo,hi), method="bounded", options=dict(xatol=1e-12))
    measured = min([lo, optimum.x, hi, j*step], key=fit)
    return dict(target=r["target"], band=r["band"], best=int(j), runner_up=int(rival),
                best_residual=errors[j], rival_residual=errors[rival],
                q_grid=j*step, q_continuous=measured)


def draw_line(r, shots, word, offset, design, rng):
    counts = np.asarray(r["joint_counts"])
    rates = counts @ word/shots
    # Sample the marginal count rates jointly using the observed multinomial covariance.
    second = np.einsum("ati,ij,ik->atjk", counts/shots, word, word)
    cov = (second - rates[:,:,:,None]*rates[:,:,None,:])/shots
    eigenvalues, vectors = np.linalg.eigh(cov)
    if eigenvalues.min() < -1e-16:
        raise ArithmeticError("Measured joint-count covariances are inconsistent")
    chol = vectors * np.sqrt(np.maximum(eigenvalues, 0))[:, :, None, :]
    rates = rates[None] + np.einsum("atij,datj->dati", chol,
                                      rng.normal(size=(len(offset),*rates.shape)))
    clock_rates = rates[:,:,:,:2].mean(axis=1)
    time_quad = np.einsum("dij,dtj->dti", design, clock_rates-offset[:,None,None])
    times = np.unwrap(np.arctan2(time_quad[:,:,1],time_quad[:,:,0]))
    y = rates[:,:,:,2]
    phase = np.unwrap(np.angle((y[:,0]-y[:,1])+1j*(y[:,2]-y[:,3])))
    t = times-times.mean(axis=1,keepdims=True)
    return (t*phase).sum(axis=1)/(t*t).sum(axis=1)


def draw_flight(row, shots, offset, inverse, rng):
    words = np.array([[0,0],[0,1],[1,0],[1,1]])
    positions = []
    angles = []
    for run in row["runs"]:
        counts = np.asarray(run["joint_counts"])
        coords = np.tile(np.repeat(np.arange(counts.shape[1]//4),4),4)
        features = np.column_stack((coords, np.repeat(words,counts.shape[1],axis=0)))
        mean = counts.ravel() @ features / shots
        centered = features-mean
        covariance = (centered.T @ ((counts.ravel()/shots)[:,None]*centered))/shots
        sims = rng.multivariate_normal(mean,covariance,size=len(offset))
        positions.append(sims[:,0])
        phase = np.einsum("dij,dj->di", inverse, sims[:,1:]-offset[:,None])
        angles.append(np.arctan2(phase[:,1], phase[:,0]))
    return (positions[1]-positions[0])/((angles[1]-angles[0])%(2*np.pi))


def continuous_q_error(line, cells, shots):
    """Delta-method frequency error with the two quadratures of each port as nuisances."""
    p = np.asarray(line["spatial_counts"]) / shots
    signal = (p[0] - p[1]).reshape(cells, 4)
    n = np.arange(cells)
    phase = n * line["q_observed"]
    design = np.column_stack((np.cos(phase), np.sin(phase)))
    derivative = np.column_stack((-n*np.sin(phase), n*np.cos(phase)))
    score = derivative @ np.linalg.lstsq(design, signal, rcond=None)[0]
    score -= design @ np.linalg.lstsq(design, score, rcond=None)[0]
    score = score.ravel()
    variance = (p.sum(axis=0) @ score**2 - np.sum((p @ score)**2)) / shots
    return np.sqrt(variance) / (score @ score)


def run():
    receipt = HERE.parent/"row_30_autonomous"/"receipt.json.gz"
    with gzip.open(receipt) as f:
        r = json.load(f)
    spatial_modes = [spatial(line,len(r["cells"])) for line in r["lines"]]
    q_errors = [continuous_q_error(line, len(r["cells"]), r["shots"]) for line in r["lines"]
                if line["q_observed"] > 0]
    print("Nonzero continuous-q local standard error range", min(q_errors), max(q_errors))
    # Parametric Gaussian approximation to the multinomial, with one shared clock calibration.
    B = 1500
    rng = np.random.default_rng(212122)
    shots = r["shots"]
    clock = rng.binomial(shots, np.asarray(r["clock_calibration"])/shots, size=(B,4))/shots
    offset = clock.mean(axis=1)
    a, b = ((clock[:,:2]-clock[:,2:])/2).T
    inverse = np.linalg.inv(np.stack((a,b,b,-a),axis=1).reshape(-1,2,2))
    draws = np.array([draw_line(line,shots,np.array(r["joint_word_order"]),offset,inverse,rng)
                      for line in r["lines"]])
    midpoint = (draws[1]+draws[2])/2
    mass = midpoint-draws[0]
    freq = (draws[2]-draws[1])/2
    right,left = [draw_flight(arm,r["packet_shots"],offset,inverse,rng) for arm in r["flights"]]
    speed = (left-right)/2
    anchor = r["calibration_q"]
    c_frequency = np.sqrt(freq**2-mass**2)/anchor
    c_speed = np.sqrt(speed*freq/anchor)
    def ci(values):
        return np.percentile(values,[2.5,97.5]).tolist()
    rows = []
    for index, row in enumerate(r["held_and_challenge"]):
        mean = predict(r["rest_frequency"],r["calibration_frequency"],r["packet_speed"],anchor,row["q"])
        model = predict(mass,freq,speed,anchor,row["q"])
        ratios = 100*(model/(midpoint-draws[3+2*index])[:,None]-1)
        point_error = 100*(mean/row["observed"]-1)
        rows.append(dict(q=row["q"], frequency=row["observed"],
                         models={name:dict(percent_error=point_error[i],
                                           count_only_95_percent_error=ci(ratios[:,i]))
                                 for i,name in enumerate(WAVES)}))
    summary = dict(spatial_fits=spatial_modes, c_from_calibration_line=np.sqrt(
        r["calibration_frequency"]**2-r["rest_frequency"]**2)/anchor,
        c_from_calibration_line_95=ci(c_frequency), c_from_packet_95=ci(c_speed),
        count_only_mass_95=ci(mass), packet_speed_95=ci(speed), rows=rows)
    def native(v):
        if isinstance(v,dict):
            return {k:native(x) for k,x in v.items()}
        if isinstance(v,list):
            return [native(x) for x in v]
        if isinstance(v,(float,np.floating)):
            return round(float(v),9)
        return v
    return native(summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    report = run()
    if args.output:
        with args.output.open("x") as f:
            json.dump(report,f,indent=2,allow_nan=False)
            f.write("\n")
    elif report != json.loads((HERE/"channel_review.json").read_text()):
        raise SystemExit("Independent channel audit diverged from archive")
    print("c from rest and one line", report["c_from_calibration_line"],report["c_from_calibration_line_95"])
    print("c from packet", report["c_from_packet_95"])
    print("last held velocity-matched parabola", report["rows"][-2]["models"]["parabola_velocity"])
