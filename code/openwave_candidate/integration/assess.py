#!/usr/bin/env python3
"""Dispatch every frozen criterion; local partial evidence never changes official status."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import re
import platform
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'candidate'))
from oph_openwave.construction import Construction
from oph_openwave.quantum_records import run as joint_probe
from oph_openwave.geometry import carrier
from oph_openwave.observer_protocol import run as observer_probe
from oph_openwave import record_preparation, record_energy, interacting_records, flux_controls, regional_charge, general_cut, modular_native, sparse_native
from oph_openwave.runner import validate_result


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def probes():
    c = carrier(mass=.173)
    rng = np.random.default_rng(9317)
    loads, records = tuple(map(int, rng.integers(-19, 20, c.ports))), (0,)*len(c.seams)
    initial = c.gauss_register(loads, records)
    for _ in range(2000):
        loads, records, _ = c.committed_event(loads, records, int(rng.integers(len(c.seams))), int(rng.integers(2)))
        require(c.gauss_register(loads, records) == initial, 'gate failed: c.gauss_register(loads, records) == initial')
        require(all(isinstance(x, int) for x in loads + records), 'gate failed: all(isinstance(x, int) for x in loads + records)')
    broken = (records[0]+1, *records[1:])
    require(c.gauss_register(loads, broken) != initial, 'gate failed: c.gauss_register(loads, broken) != initial')
    b = c.incidence
    eigenvalues, modes = np.linalg.eigh(b@b.T)
    lam, mode = eigenvalues[1], modes[:, 1]
    state = np.r_[mode, np.zeros(len(c.seams))].astype(complex)
    energy = np.sqrt(lam+c.mass**2)
    times = np.linspace(0, 7, 43)
    measured = np.array([np.sum(np.abs(c.evolve(state,t)[:c.ports])**2) for t in times])
    predicted = 1-lam/energy**2*np.sin(energy*times)**2
    clock_error = float(np.max(np.abs(measured-predicted)))
    require(clock_error < 1e-12, 'gate failed: clock_error < 1e-12')
    # Broken dispersion and broken clock controls must be rejected by the same observables.
    wrong_clock = 1-lam/energy**2*np.sin(1.2*energy*times)**2
    require(np.max(np.abs(measured-wrong_clock)) > .05, 'gate failed: np.max(np.abs(measured-wrong_clock)) > .05')
    h = c.hamiltonian
    target = np.zeros_like(h)
    target[:c.ports,:c.ports] = b@b.T+c.mass**2*np.eye(c.ports)
    target[c.ports:,c.ports:] = b.T@b+c.mass**2*np.eye(len(c.seams))
    kg_error = float(np.max(np.abs(h@h-target)))
    require(kg_error < 1e-12, 'gate failed: kg_error < 1e-12')
    wrong = h.copy(); wrong[c.ports:,c.ports:] *= -1
    require(np.max(np.abs(wrong@wrong-target)) > .1, 'gate failed: np.max(np.abs(wrong@wrong-target)) > .1')
    return {
      1: {'observable':'integer native charge and n+div(E) conservation', 'events':2000, 'gauss_error':0, 'broken_current_rejected':True,
          'remaining':'Topological identification with physical electric charge and no fractional charged excitations across the full state space.'},
      3: {'observable':'bare port/transit control on a nonzero Laplacian mode; no automatic spectral transfer to integer links', 'max_probability_error':clock_error, 'frequency':float(2*energy), 'wrong_frequency_rejected':True,
          'remaining':'Observed frequency is 2E, not the criterion rest clock; localized physical electron and native calibration remain missing.'},
      30: {'observable':'finite graph massive wave relation H²=m²+diag(BBᵀ,BᵀB)', 'max_error':kg_error, 'wrong_rest_sign_rejected':True,
           'remaining':'The bare carrier Laplacian is a component control: the integer-record lift has a different spectrum. A physical row requires its actual joint dispersion, continuum limit, common spacetime dimension and matter-state identification.'}}


def validate_evidence_links(receipt):
    """Every compact row reference resolves to retained scientific output."""
    def visit(value):
        if isinstance(value, dict):
            if 'evidence_key' in value:
                key = value['evidence_key']
                require(isinstance(key, str) and bool(key), 'invalid evidence reference')
                target = receipt
                for part in key.split('.'):
                    require(isinstance(target, dict) and part in target, 'missing evidence: '+key)
                    target = target[part]
                require(isinstance(target, dict) and bool(target), 'empty scientific evidence: '+key)
            for child in value.values(): visit(child)
        elif isinstance(value, list):
            for child in value: visit(child)
    visit(receipt)


def assess():
    path = ROOT/'platform_packet/packet.json'
    criteria = [dict(r, simplest_test=r['literal_gate'], oph_official_status='untested') for r in json.loads(path.read_text())['rows']]
    observations = probes()
    joint = joint_probe()
    require(all(joint['checks'].values()) and all(joint['mutations'].values()), 'joint probe failed')
    observations[1]['joint_record_instrument'] = joint
    observer = observer_probe()
    require(all(v is True for v in observer['checks'].values()) and all(v is True for v in observer['mutations'].values()), 'observer probe failed')
    observations[3]['joint_monitoring'] = {'evidence_key':'observer_protocol.monitoring_experiment'}
    observations[4] = {'observable': 'source support leakage under three monitored joint protocols',
        'evidence_key': 'observer_protocol.monitoring_experiment',
        'remaining': 'The prepared localized source leaks in these finite controls. No stable physical particle preparation or stability theorem is supplied.'}
    extensions = {name: module.run() for name, module in (
        ('record_preparation', record_preparation), ('record_energy', record_energy),
        ('interacting_records', interacting_records), ('flux_controls', flux_controls),
        ('regional_charge', regional_charge), ('general_cut', general_cut), ('modular_native', modular_native),
        ('sparse_native', sparse_native))}
    for name, result in extensions.items():
        validate_result(name, result)
    observations[3]['record_preparation_control'] = {'evidence_key':'record_extensions.record_preparation'}
    observations[1]['regional_apparatus'] = {'evidence_key':'record_extensions.regional_charge'}
    observations[1]['general_cut_identity'] = {'evidence_key':'record_extensions.general_cut'}
    observations[1]['local_coarse_modular_control'] = {'evidence_key':'record_extensions.modular_native'}
    observations[1]['canonical_native_observer'] = {'evidence_key':'record_extensions.sparse_native'}
    observations[1]['remaining'] = 'Physical electric-charge/apparatus dictionary and topological nonfragmentation across the admitted signed-sector/intervention domain. Modular Gauss, integer Gauss and one-load sharp outcomes retain distinct domains; this row need not await mass or spin results.'
    for ident, observable, missing in (
        (6, 'supplied static-flux rotation and antipodal-doublet controls', 'No spatial angular momentum representation for identified joint particles.'),
        (7, 'rotation rephasing and Z6 statistics controls', 'A rephasable rotation cube and special doublet do not identify physical spin-1/2 or exchange statistics.'),
        (20, 'static integer record energy and selected interacting instrument', 'No calibrated two-charge force or large-distance Coulomb force has been derived on this joint instrument.'),
        (22, 'integer string versus real-flux energy diagnostic', 'No physical confinement scale, isolated colour-charge obstruction or string-breaking observable is identified.')):
        observations[ident] = {'observable': observable, 'remaining': missing,
            'evidence_key': 'record_extensions.flux_controls' if ident in (6, 7) else 'record_extensions.record_energy'}
    observations[20]['selected_interacting_instrument'] = {'evidence_key':'record_extensions.interacting_records'}
    diagnostics = {4, 6, 7, 20, 22}
    require(sorted(r['id'] for r in criteria) == list(range(1,32)), "gate failed: sorted(r['id'] for r in criteria) == list(range(1,32))")
    rows=[]
    for criterion in criteria:
        ident=criterion['id']
        row={'id':ident, 'criterion':criterion['criterion'], 'required_test':criterion['simplest_test'],
             'official_status':criterion['oph_official_status'], 'full_criterion_met':False,
             'local_status':'diagnostic_only' if ident in diagnostics else ('partial_support' if ident in observations else 'blocked'),
             'missing_prerequisites':observations[ident]['remaining'] if ident in observations else criterion['required_before_run']}
        if ident in observations: row['measurement']=observations[ident]
        rows.append(row)
    sources=[Path(__file__).resolve(), path, *sorted((ROOT/'candidate/oph_openwave').glob('*.py'))]
    result = {'schema':'oph.integration.assessment.v1','scope':'local supporting probes; no platform submissions or official promotions',
            'seed':9317, 'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'supplied':{'mass':.173,'graph':'canonical twelve-port carrier; triangle joint-register control declared separately','time_unit':'inverse hop amplitude'},
            'sources':{str(p.relative_to(ROOT)):digest(p) for p in sources}, 'rows':rows,
            'observer_protocol':observer, 'record_extensions':extensions,
            'summary':{'dispatched':31,'partial_support':len(observations)-len(diagnostics),'diagnostic_only':len(diagnostics),'blocked':31-len(observations),'full_criteria_met':0}}
    validate_evidence_links(result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--verify',type=Path,help='recompute all probes and reject altered or stale receipts in the same runtime')
    args=parser.parse_args()
    if args.verify:
        receipt=json.loads(args.verify.read_text())
        expected = assess()
        if json.dumps(receipt, sort_keys=True, allow_nan=False) != json.dumps(expected, sort_keys=True, allow_nan=False): raise SystemExit('stale or tampered receipt: differs from source-bound recomputation')
        print('Receipt source hashes verified'); return
    result=assess()
    encoded=json.dumps(result,indent=2)+'\n'
    if args.output: args.output.write_text(encoded)
    print(json.dumps(result['summary']))

if __name__=='__main__': main()
