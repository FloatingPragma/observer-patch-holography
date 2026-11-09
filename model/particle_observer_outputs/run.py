"""Replay complete finite observer outputs and authenticate their source inputs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'
BASE=ROOT/'model/particle_calibration_selection/results.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def capture_certificate(module):
    """Optional replay timing; the certificate itself is unchanged."""
    progress = os.environ.get('OPH_REPLAY_PROGRESS') == '1'
    if progress:
        started = time.monotonic()
        print('START '+module.__name__, file=sys.stderr, flush=True)
    result = module.certificate()
    if progress:
        print('DONE '+module.__name__+' {:.3f}s'.format(time.monotonic()-started),
              file=sys.stderr, flush=True)
    return result


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
    from model.particle_sign_response_boundary import boundary as sign_boundary
    from model.particle_smooth_sign_field import smooth
    from model.particle_joint_translation_readout import readout as joint_translation
    from model.particle_physical_causality import cone as physical_cone
    from model.particle_field_state_response import state_response as field_thermal
    from model.particle_geometric_response import geometry as spatial_geometry
    from model.particle_spatial_ground_state import ground as spatial_ground
    from model.particle_ground_spectral_records import records as ground_records
    from model.particle_gauge_characteristics import characteristics as gauge_characteristics
    from model.particle_common_contract import contract as common_contract
    from model.particle_gauss_ground_state import gauss
    from model.particle_charged_residue import residue
    from model.particle_spatial_modular_response import regional as spatial_modular
    from model.particle_link_cutoff_limit import cutoff
    from model.particle_physical_boundary_state import boundary as physical_boundary
    from model.particle_primary_locality import locality as primary_locality
    from model.particle_yang_mills_gap_contract import transfer as gap_transfer
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
                   local_gauge,wilson_density,local_history,sign_boundary,smooth,joint_translation,
                   physical_cone,field_thermal,spatial_geometry,
                   spatial_ground,ground_records,gauge_characteristics,
                   common_contract,gauss,residue,spatial_modular,cutoff,physical_boundary,primary_locality,gap_transfer):
        paths.update(getattr(module,'SOURCE_INPUTS',()))
        owner=Path(module.__file__).parent
        paths.update(str(p.relative_to(ROOT)) for p in owner.iterdir()
                     if p.is_file() and p.suffix in ('.py','.md'))
    paths.add('model/matched_einstein_limit/frame.py')
    paths.add('model/particle_common_contract/contract.json')
    paths.add('model/particle_joint_reconstruction/dependencies.json')
    pins={p:sha(ROOT/p) for p in sorted(paths)}
    report={'schema':'oph.complete_finite_observer_outputs.v1',
            'input_sha256':pins,
            'neutral_scalar':capture_certificate(higgs_records),
            'charged_scalar':capture_certificate(charged_higgs_records),
            'complete_apparatus':capture_certificate(readout_contract),
            'thermal_clock':capture_certificate(thermal_clock),
            'gauge_field':capture_certificate(field_route),
            'holonomy_transport':capture_certificate(holonomy_transport),
            'consensus_records':capture_certificate(consensus_records),
            'gauge_input':capture_certificate(gauge_input),
            'bare_effective':capture_certificate(bare_effective),
            'controller_readout':capture_certificate(controller_readout),
            'defect_field_attachment':capture_certificate(attachment),
            'field_causality':capture_certificate(causality),
            'energy_and_gravity_response':capture_certificate(energy),
            'full_scalar_ground':capture_certificate(ground),
            'regional_modular_response':capture_certificate(regional),
            'gauge_field_preparation':capture_certificate(preparation),
            'full_ground_charged_records':capture_certificate(spectrum),
            'prepared_spatial_field_records':capture_certificate(experiment),
            'local_gauge_density_family':capture_certificate(local_gauge),
            'local_wilson_spectral_density':capture_certificate(wilson_density),
            'fixed_time_local_polynomial_records':capture_certificate(local_history),
            'exact_sign_response_boundary':capture_certificate(sign_boundary),
            'smooth_interacting_local_records':capture_certificate(smooth),
            'joint_position_translation_records':capture_certificate(joint_translation),
            'refinement_uniform_causality':capture_certificate(physical_cone),
            'interacting_field_thermal_records':capture_certificate(field_thermal),
            'declared_spatial_geometry_response':capture_certificate(spatial_geometry),
            'spatial_low_energy_preparation':capture_certificate(spatial_ground),
            'ground_energy_spectral_records':capture_certificate(ground_records),
            'gauge_matter_characteristics':capture_certificate(gauge_characteristics),
            'common_source_contract':capture_certificate(common_contract),
            'Gauss_constrained_preparation':gauss.certificate(include_records=True),
            'accessible_charged_threshold':capture_certificate(residue),
            'spatial_regional_modular_response':capture_certificate(spatial_modular),
            'fixed_volume_link_cutoff_limit':capture_certificate(cutoff),
            'physical_boundary_state':capture_certificate(physical_boundary),
            'primary_locality':capture_certificate(primary_locality),
            'continuum_gap_transfer':capture_certificate(gap_transfer),
            'common_construction_specification_resolved':True,
            'finite_physical_preparation_including_boundary_charges_resolved':True,
            'finite_complete_boundary_records_and_locality_resolved':True,
            'empirical_mass_inputs':[],
            'physical_particle_assignments':[],
            'source_predicted_GeV_masses':[],
            'empirical_calibrations_eliminated':[],
            'canonical_joint_dependency_map':'model/particle_joint_reconstruction/dependencies.json',
            'comparison_requirement_domains':[
                {'id':1,'name':'preparation and complete instruments','status':'selected finite Gauss-constrained polynomial preparation, physical boundary sectors and complete central-charge readout constructed; joint refinement compatibility and a physical particle detector dictionary are unproved'},
                {'id':2,'name':'common operational position, momentum, time and energy','status':'selected Hamiltonian source clock and supported polynomial modular generator fixed; separate Gibbs, position and crystal-momentum controls retain their own states; physical particle observables and proper-time dictionary unproved'},
                {'id':3,'name':'interacting relativistic field behavior','status':'fixed-positive-smoothing local records, refinement-uniform causal bounds and same-clock quadratic gauge-matter cone matching established on declared families; convergent physical states, interacting Lorentz transformations and chiral matter are unproved'},
                {'id':4,'name':'physical masses, mixing, scattering and decays','status':'finite spectral responses and a conditional same-clock local-correlation gap criterion established; interacting four-dimensional Yang-Mills decay premises, physical mass and scattering matching unproved'},
                {'id':5,'name':'defect-to-species attachment','status':'declared protected-label interaction established; physical species identification unproved'},
                {'id':6,'name':'controlled errors and same-source gravity','status':'fixed-time regulator/volume/cutoff bounds, refinement-uniform influence tails, same-action energy and declared spatial strain responses established; observer-derived proper geometry, physical state refinement and Einstein regional premises unproved'}],
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
        'exact local marginal consistency, fixed-time complete polynomial record limits and finite word-code error moduli',
        'native inverse-gap electric-response counterfamily distinguishing low-mode weight from dynamical control',
        'declared smooth gauge-coordinate action with exponentially summable interactions and complete fixed-time local quantum records',
        'exact joint source/reference translations, conserved crystal momentum and complete position-then-momentum records on all occupations',
        'physical-distance influence bound uniform under declared spacing refinement, with complete-history and implementation errors',
        'unchanged208-state quantum-flux Gibbs target, complete reference-preserving polynomial preparation and169 exact rational flux histories',
        'same-H variational witness excluding the spatial degree-zero tracial-Fock preparation as a finite ground',
        'declared spatial coframe coupling with six source strain responses, frame rotation covariance and fixed-regulator dilation identity',
        'compact constant-coframe causal bounds on the same parent Wilson action, with outgoing-link electric cross terms onsite',
        'shared208-state thermal and pressure control with complete stress instruments, retaining target Gibbs versus prepared polynomial state distinctions',
        'gap-free low-energy preparation plans for full spatial CAR and product coordinate cutoffs, with reference-preserving polynomial error',
        'conjugate-charge spectral symmetry and complete phase records on the same cooled208-state control, with a finite Lean downward-record bound',
        'native trace-dual electric metric and Wilson plaquette Hessian give a common twelve-direction quadratic gauge cone and an explicit matter-speed matching relation']
    report['finite_output_connections'] += [
        'one authenticated open Wilson source, Gauss polynomial state, retained idle reference and source clock, with typed actual-producer joins',
        'source-only physical Gauss preparation and nonzero plaquette dynamics, with a separate complete78-state directional control',
        'nonzero original charged-threshold residues and complete40-outcome collision records on the native208-state ground trial',
        'physical regional support floors on the selected polynomial state, with separately calculated warm internal modular controls',
        'fixed-volume electric cutoff removal with same-clock complete bounded-history subsequences, distinct from joint physical refinement']
    report['finite_output_connections'] += [
        'actual physical Wilson-loop boundary sectors and complete central-charge records on the selected Gauss polynomial state',
        'same-state complete boundary and reference records under a remote plaquette intervention, with a rational refinement-uniform influence bound',
        'dense local Euclidean-decay criterion for a same-clock continuum gap, with noncollapse and mixed-state boundaries explicit']
    common_contract.validate_inventory(report)
    selected=report['common_source_contract']['fixed_primary_member']
    common_contract.verify_result(selected, report['spatial_regional_modular_response']['selected_common_state_consumer'])
    physical_boundary.verify_bound(selected,report['physical_boundary_state']['selected_state_boundary'])
    primary_locality.verify_bound(selected,report['primary_locality']['selected_member'])
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
