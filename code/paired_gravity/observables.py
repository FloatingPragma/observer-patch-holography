"""Candidate exact response algebra and interval/measurement calculations."""
from fractions import Fraction as F
import mpmath as mp
from .format import load


def exp_bounds(x, terms=20):
    if x < 0:
        a, b = exp_bounds(-x, terms)
        return 1/b, 1/a
    if not 0 <= x <= F(1, 100):
        raise ValueError('small finite exponent required')
    term = total = F(1)
    for k in range(1, terms+1):
        term *= x/k
        total += term
    tail = term*x/(terms+1)/(1-x/(terms+2))
    return total, total+tail


def interval_strings(bounds):
    scale = 10**60
    lo, hi = bounds
    return [str(F((lo*scale).__floor__(), scale)), str(F((hi*scale).__ceil__(), scale))]


def build(data_path, sources):
    data = load(data_path)
    output = []
    for gamma in (-1, 0, 1, 2):
        strength = F(1, 10000)
        delta = -strength/2
        lo, hi = exp_bounds(delta)
        rows = []
        for row in data['deflection']['rows']:
            measured = 1+F(row['gamma_minus_one'])/2
            sd = F(row['sd_gamma'])/2
            predicted = F(1+gamma, 2)
            rows.append(dict(name=row['name'], measured=str(measured), sd=str(sd),
                             prediction=str(predicted), residual=str(predicted-measured)))
        output.append(dict(gamma=gamma, stiffness=[[str(1+gamma**2), str(-gamma)], [str(-gamma), '1']],
                           schur='1', marginal_covariance_coefficient='1',
                           clock_exponent=str(delta), clock_ratio=interval_strings((lo, hi)),
                           optical_exponent=str((1+gamma)*strength),
                           asymptotic_linear_deflection=str(2*(1+gamma)*strength),
                           exterior_tensor_coefficient=str(1-gamma), comparisons=rows))
    red = data['redshift']
    residual = 1-F(red['ratio'])
    systematic = F(red['systematic_bound'])
    profile = max(F(0), abs(residual)-systematic)
    finite = []
    for item in sources:
        u, w = list(map(F, item['u'])), list(map(F, item['w']))
        eps_u, eps_w = [F(x[1]) for x in item['error_intervals']]
        delta = u[42]-u[21]
        optical = []
        for g in (-1, 0, 1, 2):
            exponent = (1+g)*u[21]+w[21]
            optical.append(dict(gamma=g, exponent=str(exponent), index=interval_strings(exp_bounds(exponent)),
                                log_stationary_error=str(abs(1+g)*eps_u+eps_w)))
        finite.append(dict(projection=item['projection'], emit=[21, 33], receive=[42, 33],
                           clock_exponent=str(delta), clock_ratio=interval_strings(exp_bounds(delta)),
                           log_stationary_clock_error=str(2*eps_u), optical=optical))
    rays = []
    for raw in ('0', '1/1000000', '1/10000', '1/1000'):
        p = F(raw)
        with mp.workdps(70):
            x = mp.mpf(p.numerator)/p.denominator
            def integrand(t):
                denominator = 2-t*t+(mp.expm1(-2*x*t*t)/(t*t) if t else -2*x)
                return 2/mp.sqrt(denominator)
            angle = 2*mp.quad(integrand, [0, 1])-mp.pi
            numerical = F(int(mp.nint(angle*10**45)), 10**45)
        rays.append(dict(p=raw, numerical_angle=str(numerical), linear_angle=str(2*p),
                         rigorous_remainder=str(11*p*p)))
    return dict(branches=output, finite_reads=finite, nonlinear_rays=rays, redshift=dict(prediction='1', residual=str(residual),
                profiled_statistical_residual=str(profile/F(red['statistical_sd'])),
                reference_fraction=red['reference_fraction']),
                outcome='SOURCE_CLASS_DOES_NOT_SELECT_PAIRED_GRAVITY',
                physical_establishment=False, raw_measurements_reanalysed=False)
