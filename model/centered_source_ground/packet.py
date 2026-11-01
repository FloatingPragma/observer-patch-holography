"""Cold source-bound replay for the changed clock and its actual ground."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import bounds,cube,dimer,faithfulness,gap,particles

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).resolve().parent/'results.json'
PRIORS=('interacting_vacuum','thermodynamic_ground','continuum_field_net',
        'continuum_ground_constraints','spectral_readout','continuum_cpt_constraint')


def source_bindings():
    paths=integration_source_paths(ROOT,{
        'domain_and_proof':'model/centered_source_ground/README.md',
        'source_receipts':[f'model/{name}/results.json' for name in PRIORS],
        'independent_checks':[],
    })
    paths.add('paper/tex_fragments/CENTERED_SOURCE_GROUND.tex')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def validate(value):
    if value.get('schema')!='oph.centered_source_ground.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid centered-source receipt')
    if value.get('source_bindings') != source_bindings():
        raise ValueError('centered-source receipt has incomplete or stale custody')
    row=value['cube']['bounds']
    if row['centered_global_ground_number']!=8 or row['old_N4_ground_remains_centered_ground'] is not False:
        raise ValueError('centered ground sector changed')
    if row['unique_N8_ground_proved_by_imported_Lieb'] is not True:
        raise ValueError('imported finite uniqueness result lost')
    if row['neutral_gap_from_charge_bounds'] is not False or row['thermodynamic_ground_uniqueness_proved'] is not False:
        raise ValueError('finite uniqueness promoted beyond its domain')
    finite = value['finite_gap']
    if (finite['full_Fock_dimension'] != 65536 or finite['balanced_N8_dimension'] != 4900
            or finite['full_centered_gap_strict_lower'] != '1/128'
            or finite['Gram_residual_strict_upper'] != '1/512'
            or finite['integer_proposals_recomputed_and_verified'] is not True):
        raise ValueError('complete finite centered gap certificate changed')
    for key in ('numerical_eigenvalue_or_residual_used_for_admission',
                'full_source_couplings_changed', 'thermodynamic_gap_proved',
                'relativistic_mass_identified', 'localized_topological_protection_proved',
                'native_ground_preparation_executed'):
        if finite[key] is not False:
            raise ValueError('finite gap promoted beyond its proved domain')
    if (finite['normalized_trial_ground_infidelity_strict_upper'] != '1/1000000'
            or finite['ground_energy_enclosure_strict_lower'] != '-206203/50000'
            or finite['ground_energy_enclosure_strict_upper'] != '-82481/20000'):
        raise ValueError('certified finite reference approximation changed')
    charged = value['charged_particles']
    selection = charged['selection']
    repair = charged['complete_repair']
    if (selection['ground_twice_spin_twice_eta'] != [0, 0]
            or selection['charged_twice_spin_twice_eta'] != [1, 1]
            or selection['joint_singlet_isotypic_dimension'] != 588
            or selection['charged_joint_irrep_multiplicity'] != 1344
            or selection['charged_joint_isotypic_dimension'] != 5376
            or selection['full_K_cyclic_response_inclusion'] is not True
            or selection['every_isotypic_copy_reachable_claimed'] is not False
            or selection['internal_spin_is_physical_Lorentz_spin'] is not False
            or repair['full_number_sectors'] != [7, 9]
            or repair['joint_input_dimension'] != 22880
            or repair['native_raw_slots_retained'] != 120
            or repair['outside_slots_retained'] != 10
            or repair['repair_support_number'] != 6
            or repair['nontrivial_native_defect_record_obtained'] is not False
            or repair['new_native_preparation_or_detector_executed'] is not False
            or charged['new_gap_or_particle_mass_claimed'] is not False
            or charged['physical_species_identified'] is not False):
        raise ValueError('complete charged-sector scope changed')
    for key in ('all_native_maps_zero_on_entire_domain',
                'every_outside_map_identity_on_entire_domain',
                'arbitrary_reference_cross_terms_and_finite_histories_preserved'):
        if repair[key] is not True:
            raise ValueError('complete charged-sector repair identity changed')
    response = charged['response_and_current']
    if (response['oriented_bond_currents_checked'] != 12
            or response['local_charge_and_continuity_equations_checked'] != 8
            or response['all_local_charge_measures_equal'] is not True
            or response['physical_electromagnetic_normalization_or_CPT_claimed'] is not False):
        raise ValueError('charged response/current domain changed')
    if charged['onsite_mass_order']['actual_unique_ground_order_expectations'] != ['0']*4:
        raise ValueError('a finite singlet was assigned a classical mass order')
    for key in ('finite_many_body_gap_is_evidence_for_onsite_mass_order',
                'mass_order_fluctuations_claimed_zero',
                'interacting_gap_or_thermodynamic_symmetry_breaking_excluded'):
        if charged['onsite_mass_order'][key] is not False:
            raise ValueError('zero finite order expectation promoted beyond its domain')
    lieb=value['cube']['imported_Lieb_theorem']
    if lieb['proof_imported_not_reproved'] is not True:
        raise ValueError('external theorem attribution lost')
    if (lieb['real_symmetric_hopping_checked'] is not True
            or lieb['connected_vertices']!=list(range(8))
            or lieb['bipartition_sizes']!=[4,4] or lieb['strict_positive_U']!='1'
            or lieb['particle_number']!=8 or lieb['even_half_filling'] is not True
            or lieb['ground_hidden_spin']!=0 or lieb['ground_multiplicity']!=1):
        raise ValueError('imported spin theorem hypotheses or conclusion changed')
    local=value['cube_local_state']
    if (local['full_four_state_local_CAR_faithful'] is not True
            or local['double_occupation_strict_lower']!='1/15000'
            or local['double_occupation_strict_upper']!='1/8'
            or local['infinite_ground_local_faithfulness_inferred'] is not False
            or local['geometric_modular_charge_identified'] is not False):
        raise ValueError('finite centered local-state theorem changed or promoted')
    for key in ('old_ground_reused_without_check','interacting_Lorentz_continuum_proved','physical_mass_or_species_identified'):
        if value['bounds'][key] is not False:raise ValueError('unproved centered physical conclusion')


def canonical(value):return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def _produce():
    pins=source_bindings()
    for name in PRIORS:
        prior=strict_json.loads((ROOT/f'model/{name}/results.json').read_text())
        if prior.get('all_pass') is not True:raise ValueError('prior proof did not pass')
        for path,digest in prior['source_bindings'].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('prior proof source changed: '+path)
    value={'schema':'oph.centered_source_ground.v1','all_pass':True,
           'cube':cube.certificate(),'cube_local_state':faithfulness.certificate(),
           'finite_gap':gap.certificate(),'charged_particles':particles.certificate(),
           'dimer':dimer.certificate(),'bounds':bounds.certificate(),
           'source_bindings':pins,
           'proof_scope':'New charge-clock and new ground. Exact full-source symmetry, imported all-N lower bound and Lieb theorem with checked finite hypotheses, rational rank-one Gram proof of a full finite gap, certified trial-state error, complete charged-sector and current identities, full-CAR one-site faithfulness, dimer readout and centered thermodynamic construction. Floating eigensolvers only propose integer witnesses; they do not determine admission. No physical particle, topological protection or thermodynamic mass is inferred.'}
    validate(value)
    if source_bindings()!=pins:raise ValueError('source changed during replay')
    return value


def report():
    pins=source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-centered-ground-') as temp:
        path=Path(temp)/'worker.json';env=dict(os.environ)
        env['OPENBLAS_NUM_THREADS']=env['OMP_NUM_THREADS']='1'
        subprocess.run([sys.executable,'-W','error','-X','pycache_prefix='+str(Path(temp)/'pycache'),
                        '-m','model.centered_source_ground.packet','--worker',str(path)],
                       cwd=ROOT,env=env,check=True,capture_output=True,text=True)
        value=strict_json.loads(path.read_text())
    if source_bindings()!=pins or value['source_bindings']!=pins:
        raise ValueError('source changed across cold replay')
    validate(value);return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',nargs='?',const=OUTPUT,type=Path)
    group.add_argument('--check',nargs='?',const=OUTPUT,type=Path)
    group.add_argument('--worker',type=Path)
    group.add_argument('--diagnostic',action='store_true')
    args=parser.parse_args()
    if args.diagnostic:print(canonical(cube.numerical_diagnostic()));return
    if args.worker:args.worker.write_text(canonical(_produce()),encoding='utf-8',newline='\n');return
    value=report()
    if args.write:args.write.write_text(canonical(value),encoding='utf-8',newline='\n')
    else:
        retained=strict_json.loads(args.check.read_text());validate(retained)
        if canonical(value)!=canonical(retained):raise ValueError('centered-source replay changed')
    print(json.dumps({'all_pass':True,'source_bindings':len(value['source_bindings']),
                      'physical_identification':False}))


if __name__=='__main__':main()
