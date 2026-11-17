"""Root audit of the frozen full-source radical bounds and literal trial.

No producer code imported. Consumes its CAR compression, symmetry-average
identity and norm bounds; independently rebuilds trial moments and checks
all claimed radical inequalities by rational determinant signs.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sympy as s

HERE = Path(__file__).resolve().parent
SOURCE = HERE/'hadrons/shared_auxiliary'


def lower_bound(a, b, c2, bound):
    # The comparison 2x2 matrix minus bound is positive semidefinite.
    return a >= bound and b >= bound and (a-bound)*(b-bound) >= c2


def upper_bound(a, b, c2, bound):
    # Bound lies between the two roots, hence above the lower eigenvalue.
    return (a-bound)*(b-bound) <= c2


def frac(x):
    return F(int(x.p), int(x.q))


def main():
    paths = [SOURCE/'MANIFEST.json', SOURCE/'FESHBACH_BOUNDS.json', SOURCE/'LOCAL_ALGEBRA.json',
             HERE/'hadrons/RAW_COEFFICIENT_INPUTS.json', HERE/'hadrons/SPIN_FLAVOUR_SPECTRUM.json',
             HERE/'hadrons/spin_flavour_cache/block_3_3_3_3.json', HERE/'SPIN_FLAVOUR_REVIEW_MANIFEST.json']
    pins = {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result, local, inputs, old, block = [json.loads(p.read_text()) for p in paths[1:6]]
    norms = {r['N']: r for r in local['norms']}
    mass = F(inputs['432']['declared_light_midpoint'])
    unit = F(inputs['electron_anchor_GeV'])/mass
    G = s.Matrix(block['G'])
    M, K, E = [s.Matrix(block['components'][name]) for name in ('M', 'K', 'E')]
    Gi = G.inv(method='DM')
    current_form = s.Rational(32, 5)*G-s.Rational(2, 45)*M*Gi*M+E/180
    checks = []
    for member, ancestor in zip(result['results'], old['results']):
        kap = F(member['kappa'])
        assert kap == F(ancestor['kappa'])
        trial = member['vacuum_trial']
        v = s.Matrix(trial['integer_source_coordinates'])
        norm = frac((v.T*G*v)[0])
        H = s.Rational(mass.numerator, mass.denominator)*M+K/8+s.Rational(kap.numerator, kap.denominator)*E/6
        p = frac((v.T*H*v)[0])/norm
        d = frac((v.T*current_form*v)[0])/norm
        kinetic = frac((v.T*K*v)[0])/(8*norm)
        cross = frac((v.T*current_form*Gi*H*v)[0])/norm
        assert norm == F(trial['exact_source_norm'])
        assert p == F(trial['exact_bare_Rayleigh'])
        assert d == F(trial['exact_mixed_current_expectation'])
        assert kinetic == F(trial['exact_hopping_expectation'])
        assert cross == F(trial['exact_G2_H_expectation'])
        qtrial = 3*kap+(cross-kinetic)/d
        ctrial = d/64
        upper = F(trial['exact_full_source_variational_upper'])
        assert qtrial == F(trial['exact_auxiliary_trial_diagonal'])
        assert ctrial == F(trial['exact_offdiagonal_squared'])
        assert upper_bound(p, qtrial, ctrial, upper)
        sectors = {}
        for r in member['all_N_sector_bounds']:
            N = r['N']
            prior = next(z for z in ancestor['complete_baryon_number_ground_records'] if z['N'] == min(N, 24-N))
            assert r['P0_ground_interval_exact'] == prior['exact_interval']
            plo, phi = map(F, r['P0_ground_interval_exact'])
            d1, d2, q = map(F, (r['Q0_one_aux_lower'], r['Q0_two_aux_lower'], r['Q0_lower_exact']))
            a2, b2 = map(F, (r['P0_Q0_coupling_norm_squared_bound'], r['Q0_one_two_coupling_norm_squared_bound']))
            assert a2 == F(norms[N]['P_zero_to_one_aux_norm_squared'])
            assert b2 == F(norms[N]['Q_one_to_two_aux_norm_squared_bound'])
            count = min(N, 24-N)
            if count:
                omega_upper = (6*kap-d2)/count
                assert omega_upper > 0 and omega_upper**2 >= F(1, 64)+mass**2
                assert d1 == 3*kap-(count-2 if N == 12 else count)*omega_upper
            else:
                assert d1 == 3*kap and d2 == 6*kap
            assert lower_bound(d1, d2, b2, q)
            lo, hi = map(F, r['full_source_ground_interval_exact'])
            assert lower_bound(plo, q, a2, lo)
            assert hi == (upper if N == 12 else phi)
            sectors[N] = (lo, hi)
        assert set(sectors) == set(range(0, 25, 3))
        vlo, vhi = sectors[12]
        assert all(z[0] > vhi for N, z in sectors.items() if N != 12)
        comp = min(F(r['lowest_energy_exact_interval'][0]) for r in ancestor['all_competing_irrep_records']
                   if r['N'] == 12 and r['U4_highest_weight'] != [3, 3, 3, 3])
        vacuum_record = next(r for r in member['all_N_sector_bounds'] if r['N'] == 12)
        spin_lower = F(member['other_spin_or_isospin_lower_exact'])
        assert lower_bound(comp, F(vacuum_record['Q0_lower_exact']),
                           F(vacuum_record['P0_Q0_coupling_norm_squared_bound']), spin_lower)
        assert (spin_lower > vhi) == member['spin0_isospin0_vacuum_certified']
        for gap in member['baryon_gap_bounds']:
            blo, bhi = sectors[gap['N']]
            pair = blo-vhi, bhi-vlo
            assert pair == tuple(map(F, gap['exact_gap_interval']))
            assert all(abs(float(x*unit)-y) < 1e-12 for x, y in zip(pair, gap['conditional_GeV_interval']))
        remainder = member['H_eff2_uniform_operator_remainder']
        gamma, perturbation, energy, margin = map(F, [remainder[k] for k in
            ('H0_gap', 'perturbation_norm_upper', 'energy_absolute_upper', 'resolvent_margin')])
        assert gamma == 4*kap/3 and energy == max(abs(vlo), abs(vhi))
        assert margin == gamma-perturbation-energy
        assert remainder['bound_exists'] == (margin > 0)
        if margin > 0:
            assert F(remainder['bound_exact']) == F(47, 96)*(perturbation+energy)/(gamma*margin)
        checks.append({'kappa': str(kap), 'full_number_sector_bounds_checked': len(sectors),
                       'literal_integer_trial_moments_rebuilt': True,
                       'N12_vacuum_certified': True, 'singlet_vacuum_certified': spin_lower > vhi,
                       'singlet_vacuum_margin_source_units': float(spin_lower-vhi),
                       'lowest_baryon_gap_GeV': member['baryon_gap_bounds'][0]['conditional_GeV_interval']})
        print('verified full source bounds and trial', kap, flush=True)
    assert pins == {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    out = {'schema': 'oph.shared_auxiliary.root_bound_review.v1', 'all_checks_pass': True,
           'checks': checks, 'pins': pins,
           'consumed_not_rederived': ['CAR compression and colour/spin character',
               'SU4 averaged mixed-current identity', 'P/Q coupling norms', 'double-commutator trial formula',
               'upper-number particle-hole map'],
           'independent': ['trial rational moments', 'all 27 radical block bounds by determinant signs',
               'vacuum competition and gap subtractions', 'spin comparison and remainder arithmetic'],
           'full_source_baryon_spin_proved': False, 'physical_mass_solution': False}
    (HERE/'SHARED_AUXILIARY_BOUND_REVIEW.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__':
    main()
