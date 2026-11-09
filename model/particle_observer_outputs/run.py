"""Replay complete finite observer outputs and authenticate their source inputs."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'
BASE=ROOT/'model/particle_calibration_selection/results.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    from . import higgs_records, charged_higgs_records, readout_contract, field_route, thermal_clock
    from . import holonomy_transport, consensus_records
    from . import gauge_input, bare_effective
    from . import controller_readout
    inherited=json.loads(BASE.read_text())['input_sha256']
    for path,digest in inherited.items():
        if sha(ROOT/path)!=digest:
            raise ValueError('stale inherited source: '+path)
    paths=set(inherited)|{str(BASE.relative_to(ROOT))}
    paths.update(str(p.relative_to(ROOT)) for p in HERE.iterdir()
                 if p.is_file() and p.suffix in ('.py','.md'))
    for module in (higgs_records,charged_higgs_records,readout_contract,field_route,thermal_clock,
                   holonomy_transport,consensus_records,gauge_input,bare_effective,controller_readout):
        paths.update(getattr(module,'SOURCE_INPUTS',()))
    pins={p:sha(ROOT/p) for p in sorted(paths)}
    report={'schema':'oph.complete_finite_observer_outputs.v1',
            'input_sha256':pins,
            'neutral_scalar':higgs_records.certificate(),
            'charged_scalar':charged_higgs_records.certificate(),
            'complete_apparatus':readout_contract.certificate(),
            'thermal_clock':thermal_clock.certificate(),
            'gauge_field':field_route.certificate(),
            'holonomy_transport':holonomy_transport.certificate(),
            'consensus_records':consensus_records.certificate(),
            'gauge_input':gauge_input.certificate(),
            'bare_effective':bare_effective.certificate(),
            'controller_readout':controller_readout.certificate(),
            'empirical_mass_inputs':[],
            'physical_particle_assignments':[],
            'source_predicted_GeV_masses':[],
            'empirical_calibrations_eliminated':[],
            'finite_output_connections':['same-source neutral and charged spectral response',
                'complete scalar and energy instruments with coherent-control compile bounds',
                'faithful Gibbs modular clock on the same interacting Hamiltonian',
                'preparation-weight cancellation in a complete gap-readout protocol',
                'gauge-coordinate cutoff and exact numerical complete-record error bounds',
                'finite link-holonomy birth, work exchange and regional protection',
                'complete record copying, meaning comparison and finite event-function gluing',
                'reuse of operational native/scalar couplings and the conditional gauge-ray receipt',
                'charged spectral coefficient inversion cancelling the charge-sector offset',
                'controller attenuation and phase-class reconstruction with a retained source-normalization premise']}
    if any(sha(ROOT/p)!=h for p,h in pins.items()):
        raise ValueError('input changed during complete-output replay')
    return json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',action='store_true')
    args=parser.parse_args(); result=build()
    if args.write:OUTPUT.write_text(result)
    elif not OUTPUT.is_file() or OUTPUT.read_text()!=result:
        raise ValueError('stale complete observer-output report')
    print('PASS: finite observer outputs, complete records and source bindings')


if __name__=='__main__':main()
