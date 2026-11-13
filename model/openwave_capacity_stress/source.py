"""Exact finite work records and a correlated, fixed-anchor expansion test."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent


def decimal(x, places=12, upper=False):
    """Directed decimal enclosure of an exact rational."""
    x = F(x)
    scale = 10 ** places
    n = -((-x.numerator * scale) // x.denominator) if upper else x.numerator * scale // x.denominator
    sign = '-' if n < 0 else ''
    q, r = divmod(abs(n), scale)
    return f'{sign}{q}.{r:0{places}d}'


def sqrt_interval(x, places=16):
    x = F(x)
    assert x >= 0
    scale = 10 ** places
    q = isqrt(x.numerator * scale * scale // x.denominator)
    lo = F(q, scale)
    return lo, lo if lo * lo == x else F(q + 1, scale)


def enclosure(pair):
    return [decimal(pair[0]), decimal(pair[1], upper=True)]


def energy(volume, k, epsilon=F(1), volume_coupled=True):
    assert k > 0 and volume > 0 and epsilon > 0
    return epsilon * (volume if volume_coupled else 1) / k


def cycle(k, v0=F(1), v1=F(3, 2)):
    """All controller work, including changing the active capacity."""
    states = [(v0, k), (v1, k), (v1, k+1), (v0, k+1), (v0, k)]
    es = [energy(v, n) for v, n in states]
    ws = [b-a for a, b in zip(es, es[1:])]
    assert sum(ws) == 0
    return {'k': k, 'energies': list(map(str, es)),
            'work_strokes': list(map(str, ws)), 'total_work': str(sum(ws)),
            'spurious_work_if_capacity_strokes_omitted': str(ws[0]+ws[2]),
            'complete_appended_reference_records': [
                {'word': word, 'reference_bit': bit, 'probability': str(F(1,2**(k+1)))}
                for word in range(2**k) for bit in (0,1)],
            'capacity_reference': 'a degenerate reference bit in the maximally mixed state; both bit outcomes retained'}


def phase_records(k, phase_pi, cosine):
    """Hadamard preparation, native free phase, complete X and word records.

    A phase is specified as an exact multiple of pi. Time is phase/deltaE;
    the experiment measures the declared Hamiltonian gap, not an inferred mass.
    """
    delta = energy(F(3, 2), k) - energy(F(1), k)
    records = [{'word': word, 'X': sign, 'probability': str(F(1+sign*cosine, 2 * 2**k))}
               for word in range(2**k) for sign in (-1, 1)]
    assert sum(F(r['probability']) for r in records) == 1
    return {'k': k, 'active_dimension': 2**k, 'capacity_bits': k,
            'capacity_nats_exact': f'{k}*log(2)', 'phase_pi': str(phase_pi),
            'time_over_pi': str(F(phase_pi)/delta), 'energy_gap': str(delta),
            'preparation_work': str(delta/2), 'mean_final_energy': str(energy(F(1),k)+delta/2),
            'mean_X_readout_work': '0', 'records': records}


def interpolate(z, rows, column):
    if z <= rows[0][0]:
        return rows[0][column]
    if z >= rows[-1][0]:
        return rows[-1][column]
    for a, b in zip(rows, rows[1:]):
        if a[0] <= z <= b[0]:
            return a[column] + (b[column]-a[column]) * (z-a[0]) / (b[0]-a[0])
    raise AssertionError('interpolation domain')


def load_data(name):
    return [tuple(F(x) for x in line.split(',')[:3])
            for line in (HERE/'inputs'/name).read_text().splitlines()
            if line.strip() and not line.startswith('#')]


def load_systematics():
    return [[F(x) for x in line.split()] for line in
            (HERE/'inputs/data_MM20.dat').read_text().splitlines()
            if line.strip() and not line.startswith('#')]


def covariance_components(data, sps_column=4):
    rows = load_systematics()
    diag = [sig**2 for _, _, sig in data]
    u = [[h*interpolate(z, rows, c)/100 for z,h,_ in data] for c in (1,sps_column)]
    return diag, u


def inverse_apply(diag, u, residual):
    """Exact rank-two Woodbury formula, retaining every covariance entry."""
    assert len(u) == 2 and all(d > 0 for d in diag)
    b = [[F(int(i==j)) + sum(u[i][k]*u[j][k]/diag[k] for k in range(len(diag)))
          for j in range(2)] for i in range(2)]
    det = b[0][0]*b[1][1]-b[0][1]*b[1][0]
    assert det > 0
    v = [sum(ui*r/d for ui,r,d in zip(col,residual,diag)) for col in u]
    w = [(b[1][1]*v[0]-b[0][1]*v[1])/det,
         (b[0][0]*v[1]-b[1][0]*v[0])/det]
    return [(r-u[0][i]*w[0]-u[1][i]*w[1])/diag[i] for i,r in enumerate(residual)]


def chi_square_interval(diag, u, residual_intervals):
    mid = [(a+b)/2 for a,b in residual_intervals]
    rad = [(b-a)/2 for a,b in residual_intervals]
    y = inverse_apply(diag,u,mid)
    centre = sum(a*b for a,b in zip(mid,y))
    # C >= D > 0 implies C^-1 <= D^-1; no independence assumption.
    error = 2*sum(abs(a)*b for a,b in zip(y,rad)) + sum(e*e/d for e,d in zip(rad,diag))
    return max(F(0),centre-error), centre+error


def hubble_squared(z, h0=F('67.4'), om=F('.315'), radiation=F(0)):
    assert z >= 0 and h0 > 0 and 0 < om < 1 and 0 <= radiation < 1-om
    return h0*h0*(om*(1+z)**3 + radiation*(1+z)**4 + 1-om-radiation)


def comparison(name, sps_column=4):
    data = load_data(name)
    diag,u = covariance_components(data,sps_column)
    predictions = [sqrt_interval(hubble_squared(z)) for z,_,_ in data]
    residual = [(a-h,b-h) for (a,b),(_,h,_) in zip(predictions,data)]
    chi = chi_square_interval(diag,u,residual)
    # This entire box is enclosed; the endpoints are not a parameter fit.
    box = [sqrt_interval(hubble_squared(z,F('66.9'),F('.308')))[0:1] +
           sqrt_interval(hubble_squared(z,F('67.9'),F('.322'),F('.0001')))[1:2]
           for z,_,_ in data]
    radpred = [sqrt_interval(hubble_squared(z,radiation=F('.0001'))) for z,_,_ in data]
    maxrad = max(b-a for (_,b),(a,_) in zip(radpred,predictions))
    cov = [[str((diag[i] if i==j else 0)+sum(col[i]*col[j] for col in u))
            for j in range(len(data))] for i in range(len(data))]
    from scipy.special import gammaincc
    return {'table': name, 'rows': len(data), 'fitted_chronometer_parameters': 0,
            'anchors': {'H0_km_s_Mpc': '67.4', 'Omega_m': '0.315'},
            'sps_column': sps_column, 'chi_square': enclosure(chi),
            'conditional_chi_square_tail_diagnostic': float(gammaincc(len(data)/2,float(sum(chi)/4))),
            'chi_square_anchor_and_radiation_box': enclosure(chi_square_interval(
                diag,u,[(a-h,b-h) for (a,b),(_,h,_) in zip(box,data)])),
            'max_radiation_shift_km_s_Mpc_upper': decimal(maxrad,upper=True),
            'systematics_clamped_at_upper_endpoint_rows': [i for i,(z,_,_) in enumerate(data) if z>F('1.475')],
            'full_covariance_rational': cov,
            'records': [{'z':str(z),'observed_H':str(h),'quoted_sigma':str(s),
                         'predicted_H_interval':enclosure(pred),
                         'residual_H_interval':enclosure(res)}
                        for (z,h,s),pred,res in zip(data,predictions,residual)]}


def run():
    provenance = json.loads((HERE/'inputs/provenance.json').read_text())
    for name in provenance['retained_files']:
        assert hashlib.sha256((HERE/'inputs'/name).read_bytes()).hexdigest() == provenance['files'][name]['sha256']
    records = [phase_records(k,p,c) for k in range(1,7)
               for p,c in ((F(0),1),(F(1,2),0),(F(1),-1),(F(2),1))]
    return {'verdict':'PASS: finite pressure/work records and conditional correlated expansion comparison',
            'source': {'Hamiltonian':'epsilon*V/k times identity on the active k-bit register',
                       'fixed_k_pressure':'-epsilon/k', 'fixed_k_density':'epsilon/k',
                       'fixed_k_w':'-1', 'capacity_only_control_pressure':'0',
                       'capacity_only_control_density':'epsilon/(k*V)',
                       'varying_capacity_exchange':'Q=-rho*(dk/dt)/k',
                       'mechanical_pressure_during_capacity_change':'-rho',
                       'separately_conserved_effective_w':'-1+(1/3)*dlog(k)/dlog(a)',
                       'conservation_boundary':'effective w absorbs capacity work; it is not automatically mechanical pressure'},
            'complete_phase_records':records, 'complete_cycles':[cycle(k) for k in range(1,7)],
            'comparisons':[comparison('HzTable_MM_BC03.dat'),comparison('HzTable_MM_M11.dat'),
                           comparison('HzTable_MM_BC03.dat',3)],
            'calibration_boundary':'H0 and Omega_m are supplied Planck2018 rounded inputs; Lambda amplitude, horizon, and metric-source join are not source-derived',
            'openwave_row28_promotion':False,
            'theorem_reuse':['EinsteinBranch/FixedCapacityWLaw.lean','CapacityNonidentifiability.lean'],
            'numerical_scope':'rational covariance and outward sqrt/chi-square enclosures; tail probability is a floating diagnostic'}
