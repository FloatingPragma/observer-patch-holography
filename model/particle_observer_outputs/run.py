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
    from model.particle_defect_attachment import attachment
    from model.particle_field_causality import bounds as causality
    from model.particle_gravity_response import response as energy
    from model.particle_scalar_ground import ground
    from model.particle_regional_geometry import regional
    from model.particle_field_preparation import preparation
    from model.particle_ground_charged import spectrum
    from model.particle_prepared_field_readout import experiment
    from model.particle_local_gauge_state import state as local_gauge
    from model.particle_local_wilson_density import density as wilson_density
    from model.particle_local_history_family import history as local_history
    inherited=json.loads(BASE.read_text())['input_sha256']
    for path,digest in inherited.items():
        if sha(ROOT/path)!=digest:
            raise ValueError('stale inherited source: '+path)
    paths=set(inherited)|{str(BASE.relative_to(ROOT))}
    paths.update(str(p.relative_to(ROOT)) for p in HERE.iterdir()
                 if p.is_file() and p.suffix in ('.py','.md'))
    for module in (higgs_records,charged_higgs_records,readout_contract,field_route,thermal_clock,
                   holonomy_transport,consensus_records,gauge_input,bare_effective,controller_readout,
                   attachment,causality,energy,ground,regional,preparation,spectrum,experiment,
                   local_gauge,wilson_density,local_history):
        paths.update(getattr(module,'SOURCE_INPUTS',()))
        owner=Path(module.__file__).parent
        paths.update(str(p.relative_to(ROOT)) for p in owner.iterdir()
                     if p.is_file() and p.suffix in ('.py','.md'))
    paths.add('model/matched_einstein_limit/frame.py')
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
            'defect_field_attachment':attachment.certificate(),
            'field_causality':causality.certificate(),
            'energy_and_gravity_response':energy.certificate(),
            'full_scalar_ground':ground.certificate(),
            'regional_modular_response':regional.certificate(),
            'gauge_field_preparation':preparation.certificate(),
            'full_ground_charged_records':spectrum.certificate(),
            'prepared_spatial_field_records':experiment.certificate(),
            'local_gauge_density_family':local_gauge.certificate(),
            'local_wilson_spectral_density':wilson_density.certificate(),
            'fixed_time_local_polynomial_records':local_history.certificate(),
            'empirical_mass_inputs':[],
            'physical_particle_assignments':[],
            'source_predicted_GeV_masses':[],
            'empirical_calibrations_eliminated':[],
            'physical_requirement_status':[
                {'id':1,'name':'preparation and complete instruments','status':'finite scalar and composed spatial preparations and complete readouts established; physical spatial-state attachment unproved'},
                {'id':2,'name':'common operational position, momentum, time and energy','status':'same-H finite energies and state-specific Gibbs modular clock established; common physical spatial dictionary unproved'},
                {'id':3,'name':'interacting relativistic field behavior','status':'fixed-time local polynomial record limit and local Wilson spectral control established; interacting exact-sign comparison and chiral relativistic refinement unproved'},
                {'id':4,'name':'physical masses, mixing, scattering and decays','status':'finite spectral responses established; physical mass and scattering matching unproved'},
                {'id':5,'name':'defect-to-species attachment','status':'declared protected-label interaction established; physical species identification unproved'},
                {'id':6,'name':'controlled errors and same-source gravity','status':'finite cutoff, volume and energy identities established; uniform refinement and geometric Einstein premises unproved'}],
            'finite_output_connections':['same-source neutral and charged spectral response',
                'complete scalar and energy instruments with coherent-control compile bounds',
                'faithful Gibbs modular clock on the same interacting Hamiltonian',
                'preparation-weight cancellation in a complete gap-readout protocol',
                'gauge-coordinate cutoff and exact numerical complete-record error bounds',
                'finite link-holonomy birth, work exchange and regional protection',
                'complete record copying, meaning comparison and finite event-function gluing',
                'reuse of operational native/scalar couplings and the conditional gauge-ray receipt',
                'charged spectral coefficient inversion cancelling the charge-sector offset',
                'controller attenuation and phase-class reconstruction with a retained source-normalization premise',
                'protected memory and regional holonomy controlling complete scalar records, with preparation work explicit',
                'fixed-spacing graded locality and finite-volume complete-history bounds on the same gauge-coordinate field',
                'same-action energy continuity, complete current records and fixed-state action response']}
    report['finite_output_connections'] += [
        'full-carrier scalar singlet ground, certified preparation and complete same-clock occupation histories',
        'actual internal regional modular pair interaction and source-admitted entropy response',
        'finite gauge-word preparation, complete seed records and same-clock Wilson-gap accuracy bounds',
        'complete compensated charged spectral records about the full scalar ground and offset-invariant shape tests',
        'one prepared finite field experiment composing parent, volume, cutoff and complete local detector errors',
        'referenced mixed gauge-state family with volume-independent energy densities and all-time local roughness at fixed parameters',
        'all-time local Wilson low-mode weight and one-body sign-polynomial spectral error at fixed parameters',
        'exact local marginal consistency, fixed-time complete polynomial record limits and finite word-code error moduli']
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
