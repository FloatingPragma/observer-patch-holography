"""Root checks of actual current covariances and six full-source certificates.

Consumes frozen full-CAR component matrices and Q/character bounds, but
imports no producer code. Reassembles actual current-image Gram matrices,
literal trials and rational Schur positivity with root's own inertia.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy import sparse
import sympy as s
from nuclear_threshold import inertia

HERE = Path(__file__).resolve().parent
SOURCE = HERE/'hadrons/auxiliary_spin_refinement'


def rational(q):
    return s.Rational(q.numerator, q.denominator)


def fraction(x):
    return F(int(x.p), int(x.q))


def exact_dot(A, B):
    bound = A.shape[1]*int(np.max(np.abs(A), initial=0))*int(np.max(np.abs(B), initial=0))
    if bound < 2**63:
        return A@B
    return A.astype(object)@B.astype(object)


def main():
    paths = [SOURCE/'MANIFEST.json', SOURCE/'BARYON_STATES.json', SOURCE/'MATRIX_FESHBACH.json',
             SOURCE/'Q_COLOUR_BOUNDS.json', SOURCE/'BRANCHING.json', SOURCE/'SUMMARY.json',
             HERE/'hadrons/RAW_COEFFICIENT_INPUTS.json', HERE/'nuclear_threshold.py']
    states = json.loads(paths[1].read_text())['records']
    bounds = json.loads(paths[2].read_text())['results']
    summaries = json.loads(paths[5].read_text())['results']
    inputs = json.loads(paths[6].read_text())
    mass = F(inputs['432']['declared_light_midpoint'])
    unit = F(inputs['electron_anchor_GeV'])/mass
    cached, covariance_checks = {}, []
    for state in states:
        name = state['name']
        frame_path, op_path = SOURCE/('CAR_STATE_'+name+'.npz'), SOURCE/('OPERATORS_'+name+'.npz')
        paths.extend([frame_path, op_path])
        data, ops = np.load(frame_path), np.load(op_path)
        R = data['R']
        assert R.shape[1] == 20 and all(int(w).bit_count() == 9 for w in data['words'])
        Graw = exact_dot(R.T, R)
        G = s.Matrix(Graw.tolist())
        assert G == s.Matrix(state['G'])
        assert inertia([[-F(int(v)) for v in row] for row in G.tolist()]) == 20
        Gi = G.inv(method='DM')
        components = {}
        for tag in ('M', 'K', 'E', 'G24'):
            op = sparse.csr_matrix((ops[tag+'_data'], ops[tag+'_indices'], ops[tag+'_indptr']), shape=(len(R), len(R)))
            # Protect exact integer arithmetic before multiplying the raw image.
            row_abs = int(np.max(np.asarray(abs(op).sum(axis=1)), initial=0))
            assert row_abs*int(np.max(np.abs(R), initial=0)) < 2**63
            image = op@R
            form = exact_dot(R.T, image)
            if tag == 'G24':
                assert np.array_equal(image, data['G24R'])
                assert s.Matrix(form.tolist()) == s.Matrix(state['mixed_current24'])
                D = s.Matrix(form.tolist())/24
                full_square = exact_dot(image.T, image)
                assert s.Matrix(full_square.tolist()) == s.Matrix(state['mixed_current24_square'])
                D2 = s.Matrix(full_square.tolist())/576
            else:
                assert s.Matrix(form.tolist()) == s.Matrix(state['components'][tag])
                components[tag] = s.Matrix(form.tolist())
        covariance = D2-D*Gi*D
        assert covariance == covariance.T
        # Positive semidefinite by the exact Gram projection construction above.
        rank = covariance.rank()
        covariance_checks.append({'name': name, 'full_CAR_carrier_dimension': len(R),
            'actual_current_image_reassembled': True, 'outside_low_block_current_covariance_rank': rank,
            'compressed_square_is_not_full_square': rank > 0})
        cached[name] = (G, Gi, components, D, D2, covariance)
    checks = []
    for member, summary in zip(bounds, summaries):
        kap = F(member['kappa'])
        for rec in member['records']:
            G, Gi, comp, D, D2, covariance = cached[rec['name']]
            H = rational(mass)*comp['M']+comp['K']/8+rational(kap)*comp['E']/6
            E, q, dmax, high = map(F, [rec[key] for key in
                ('full_source_lower_exact', 'Q_lower_exact', 'D_operator_norm_upper_exact', 'other_P0_ground_lower_exact')])
            assert q <= F(rec['Q_lower_exact_before_rounding']) and q > E
            assert high == min(F(r['ground_lower_exact']) for r in rec['all_other_SU4_competitors'])
            c = F(1, 64)/(q-E)
            margin = high-E-c*dmax
            assert margin > 0 and margin == F(rec['second_Feshbach_margin_exact'])
            assert c == F(rec['coefficient_c_exact'])
            B = H-rational(E)*G-rational(c)*D-rational(c*c/margin)*covariance
            count = inertia([[-fraction(x) for x in row] for row in B.tolist()])
            assert count == 20
            tr = rec['trial']
            v = s.Matrix(tr['integer_source_coordinates'])
            norm = fraction((v.T*G*v)[0])
            p = fraction((v.T*H*v)[0])/norm
            d = fraction((v.T*D*v)[0])/norm
            kin = fraction((v.T*comp['K']*v)[0])/(8*norm)
            cross = fraction((v.T*D*Gi*H*v)[0])/norm
            assert [norm, p, d, kin, cross] == [F(tr[key]) for key in
                ('source_norm_exact', 'bare_Rayleigh_exact', 'Gsum_expectation_exact',
                 'hopping_expectation_exact', 'Gsum_H_expectation_exact')]
            auxiliary = 3*kap+(cross-kin)/d
            coupling = d/64
            upper = F(tr['full_source_upper_exact'])
            assert auxiliary == F(tr['one_auxiliary_diagonal_exact'])
            assert coupling == F(tr['coupling_squared_exact'])
            assert (p-upper)*(auxiliary-upper) <= coupling
            checks.append({'kappa': str(kap), 'name': rec['name'], 'exact_negative_inertia_of_minus_B': count,
                           'literal_full_source_trial_verified': True})
        nucleon_upper = F(next(r for r in member['records'] if r['name'] == 'nucleon')['trial']['full_source_upper_exact'])
        actual = all(F(r['lower_exact']) > nucleon_upper for r in member['all_competing_P0_present_baryon_spin_sectors'])
        assert actual == member['full_source_baryon_spin_half_isospin_half_certified']
        assert actual == summary['lowest_Bplusminus1_spin_half_isospin_half_certified']
        blo, bhi = map(F, summary['baryon_energy_interval_exact'])
        vlo, vhi = map(F, summary['vacuum_energy_interval_exact'])
        interval = blo-vhi, bhi-vlo
        assert interval == tuple(map(F, summary['same_source_baryon_gap_exact_interval']))
        assert all(abs(float(x*unit)-y) < 1e-12 for x, y in zip(interval, summary['conditional_baryon_gap_GeV_interval']))
        print('verified actual-source spin certificates and gap', kap, actual, flush=True)
    pins = {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    out = {'schema': 'oph.baryon_spin.root_review.v1', 'all_checks_pass': True,
           'current_covariance_checks': covariance_checks, 'matrix_checks': checks,
           'exact_matrix_certificate_count': len(checks), 'pins': pins,
           'consumed_not_independently_rederived': ['full-CAR component operators',
                'Q-sector colour/spin character and lower bounds', 'old-source other-irrep enclosures',
                'particle-hole auxiliary intertwiner', 'double-commutator trial identity'],
           'independently_reassembled': ['actual raw-carrier frame/component forms',
                'full current square and its outside-block covariance',
                'six rational Schur-positivity certificates', 'six literal trial moments and upper bounds',
                'spin comparisons and same-vacuum mass subtraction'],
           'full_source_component_result': 'spin/isospin-half lowest baryons at kappa1 and2, conditional on consumed Q/character certificates',
           'physical_rows_solved': False}
    (HERE/'BARYON_SPIN_REVIEW.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__':
    main()
