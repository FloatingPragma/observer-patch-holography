"""Isospin operators and explicitly calibrated archived hadron fit outputs."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
UPSTREAM = ROOT / 'model/quenched_hadron_binding'


def flavour_checks():
    """Matrices label q_a anti-q_b ket amplitudes, with unit Hilbert-Schmidt norm."""
    plus = np.array([[0., 1.], [0., 0.]])
    minus = plus.T
    zero = np.diag([1., -1.]) / np.sqrt(2)
    charge = np.diag([2/3, -1/3])
    triplet = [plus, zero, minus]
    comm = lambda a, b: a @ b - b @ a
    residuals = [np.linalg.norm(comm(charge, a) - q*a)
                 for a, q in zip(triplet, [1, 0, -1])]
    residuals.extend([np.linalg.norm(comm(plus, zero) + np.sqrt(2)*plus),
                      np.linalg.norm(comm(plus, minus) - np.sqrt(2)*zero)])
    gram = [[float(np.trace(a.T @ b)) for b in triplet] for a in triplet]
    return dict(valence_order=['u anti-d', '(u anti-u - d anti-d)/sqrt(2)', 'd anti-u'],
                electric_charges=[1, 0, -1], gram=gram,
                all_traces_zero=all(np.trace(a) == 0 for a in triplet),
                max_charge_and_isospin_residual=float(max(residuals)),
                charged_conjugates_equal=bool(np.array_equal(plus.T, minus)),
                connected_correlator='delta_ab C(t), for identical u,d propagators',
                disconnected_neutral='zero: Tr_flavour(M_zero)=0 with identical u,d propagators',
                exact_domain='degenerate u,d Wilson operator, no electromagnetic field')


def multiply(value, stat, sensitivity, scale):
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('energy anchor must be finite and positive')
    return dict(value_GeV=float(value*scale), statistical_GeV=float(stat*scale),
                fit_sensitivity_GeV=float(sensitivity*scale))


def calibrate(chiral, ratio, rho_GeV, string_GeV):
    if ratio <= 0 or not np.isfinite(ratio):
        raise ValueError('positive finite light-quark tuning ratio required')
    p, syst = chiral['primary'], chiral['systematic']
    rho = dict(pion=multiply(ratio, 0, 0, rho_GeV),
               rho=multiply(1, 0, 0, rho_GeV),
               nucleon=multiply(*p['mN_over_mrho_phys'], syst['mN_over_mrho_phys'], rho_GeV))
    string = dict(rho=multiply(*p['mrho_over_sqrtsigma_phys'],
                              syst['mrho_over_sqrtsigma_phys'], string_GeV),
                  nucleon=multiply(*p['mN_over_sqrtsigma_phys'],
                                  syst['mN_over_sqrtsigma_phys'], string_GeV))
    string['pion'] = {key: ratio*value for key, value in string['rho'].items()}
    return dict(rho_anchor=rho, string_anchor=string,
                conditional_inverse_spacing_GeV=dict(rho=rho_GeV/p['mrho_phys'][0],
                                                     string=string_GeV/p['sqrt_sigma'][0]),
                pion_rho_anchor_zero_error_reason='ratio and rho mass are fixed calibration inputs; no total physical uncertainty assigned',
                statistical_domain='paired delete-one-gauge-stream jackknife of the retained fit chain',
                fit_sensitivity_domain='finite fit-choice envelope, not a confidence interval or a discretization bound',
                anchor_and_model_uncertainty_GeV=None)


def certificate(replay=False):
    upstream = json.loads((UPSTREAM/'results.json').read_text())
    tables = json.loads((UPSTREAM/'recorded_tables.json').read_text())
    if not upstream['all_checks_pass']:
        raise ValueError('upstream hadron receipt failed')
    production = upstream['production']['B']
    if replay:
        from model.quenched_hadron_binding.analysis import analyse
        actual = analyse(tables['lattices']['B'])
        for name in ['primary', 'systematic']:
            for key, value in actual['chiral'][name].items():
                if not np.allclose(value, production['chiral'][name][key], rtol=1e-10, atol=1e-12):
                    raise ValueError(f'archived fit replay differs: {name}/{key}')
    measured = upstream['measured']
    convention = upstream['convention']
    ratio = tables['lattices']['B']['ratio_phys_input']
    rho_GeV = measured['m_rho_MeV']/1000
    string_GeV = convention['sqrt_sigma_MeV']/1000
    values = calibrate(production['chiral'], ratio, rho_GeV, string_GeV)
    operators = flavour_checks()
    references = {
        'pi+': dict(value_GeV=.13957039, uncertainty_GeV=.00000018,
                    source='https://pdg.lbl.gov/2025/listings/rpp2025-list-pi-plus-minus.pdf'),
        'pi0': dict(value_GeV=.1349768, uncertainty_GeV=.0000005,
                    source='https://pdgweb.lbl.gov/2025/listings/rpp2025-list-pi-zero.pdf'),
        'pi-': dict(value_GeV=.13957039, uncertainty_GeV=.00000018,
                    source='https://pdg.lbl.gov/2025/listings/rpp2025-list-pi-plus-minus.pdf')}
    particles = []
    for label, content, charge in zip(references, operators['valence_order'], [1, 0, -1]):
        particles.append(dict(particle=label, valence_operator=content, electric_charge=charge,
                              spin_parity='0-', source='pseudoscalar Wilson quark-antiquark correlator',
                              rho_anchor=values['rho_anchor']['pion'],
                              string_anchor=values['string_anchor']['pion'],
                              experimental_reference=references[label],
                              physical_source_derived_mass_GeV=None))
    files = [HERE/'readout.py', HERE/'PROOF.md', HERE/'test_readout.py',
             UPSTREAM/'analysis.py', UPSTREAM/'literature.py', UPSTREAM/'lattice.py',
             UPSTREAM/'recorded_tables.json', UPSTREAM/'results.json']
    hashes = {str(f.relative_to(ROOT)): hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
    checks = dict(upstream_pass=upstream['all_checks_pass'],
                  isospin_and_charge=operators['max_charge_and_isospin_residual'] < 1e-14,
                  orthonormal_triplet=bool(np.allclose(operators['gram'], np.eye(3))),
                  neutral_disconnected_cancellation=operators['all_traces_zero'],
                  ratio_input_unchanged=ratio == measured['mpi_over_mrho'],
                  retained_all_four_kappas=len(production['hadrons']) == 4,
                  rho_anchor_pion_is_calibration=abs(values['rho_anchor']['pion']['value_GeV']-.13799628)<1e-14,
                  both_scale_choices_exposed=abs(values['rho_anchor']['pion']['value_GeV']-
                                                values['string_anchor']['pion']['value_GeV'])>.03)
    return dict(schema='oph.pion_mass_readout.v1', all_checks_pass=all(checks.values()), checks=checks,
                source_sha256=hashes, scope='conditional GeV conversion of archived quenched hadron fits',
                particles=particles, flavour_operators=operators, calibrations=values,
                empirical_inputs=dict(light_quark_tuning_mpi_over_mrho=ratio,
                                      rho_anchor_GeV=rho_GeV, alternative_sqrt_sigma_GeV=string_GeV,
                                      physical_point_kappa=production['chiral']['primary']['kappa_phys']),
                source_lattice={k: production[k] for k in ['lattice','beta','configurations','streams']},
                lattice_mass_rows=[{k:v for k,v in row.items() if k=='kappa' or k in
                                    ['mass_pion','mass_rho','mass_nucleon','mpi_over_mrho','mN_over_mrho']}
                                   for row in production['hadrons']],
                unresolved_nucleon_kappas=production['nucleon_not_resolved_at'],
                pion_charged_neutral_splitting=dict(source_GeV=0,
                    experimental_difference_GeV=references['pi+']['value_GeV']-references['pi0']['value_GeV'],
                    missing_terms=['electromagnetism','nondegenerate up/down masses']),
                physical_attachment_obligations=['quark loops','finite-volume control','continuum extrapolation',
                    'electromagnetic and up/down splitting','native observer/source/clock transfer'],
                archive_scope='reduced correlator and Wilson-loop jackknife tables; raw gauge configurations absent',
                detector_scope='archived correlator fits; no additional complete native observer instrument is claimed')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--replay-archive', action='store_true')
    args=parser.parse_args()
    result=certificate(args.replay_archive)
    if not result['all_checks_pass']:
        raise SystemExit('pion checks failed')
    data=json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n'
    if args.write:
        (HERE/'results.json').write_text(data)
    if args.check and data != (HERE/'results.json').read_text():
        raise SystemExit('pion receipt is stale')
    print(json.dumps(dict(all_checks_pass=True, checks=len(result['checks']),
                         rho_pion_GeV=result['calibrations']['rho_anchor']['pion']['value_GeV'],
                         string_pion_GeV=result['calibrations']['string_anchor']['pion']['value_GeV'])))


if __name__=='__main__':
    main()
