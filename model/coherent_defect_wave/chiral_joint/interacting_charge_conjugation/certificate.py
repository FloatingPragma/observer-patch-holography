"""Fresh exact symmetry replay with validated frozen source prerequisites."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
from . import symmetry
from model.particle_common_contract.provenance import import_closure
PARENTS=('atomic_spin_geometry','interacting_ground_preparation','interacting_charged_response')

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def raw_certificate():
    pins={};parents={}
    for name in PARENTS:
        path=ROOT/'model/coherent_defect_wave/chiral_joint'/name/'results.json'
        value=json.loads(path.read_text())
        for relative,expected in value['input_sha256'].items():
            if digest(ROOT/relative)!=expected:raise ValueError('changed inherited input: '+relative)
            if relative in pins and pins[relative]!=expected:raise ValueError('inconsistent ancestry')
            pins[relative]=expected
        relative=path.relative_to(ROOT).as_posix();pins[relative]=digest(path);parents[name]=value
    exact={'links':symmetry.link_checks(),'local':symmetry.local_checks(),
           'paths':[symmetry.path_checks(n) for n in (2,4)],
           'occupied_controls':symmetry.fixture_checks(),'ground_character':symmetry.ground_character()}
    seeds={'model/coherent_defect_wave/chiral_joint/neutral_local_spectrum/symmetries.py',
           'model/coherent_defect_wave/chiral_joint/oriented_links/rotations.py',
           'model/coherent_defect_wave/chiral_joint/peter_weyl_links/clebsch.py',
           'model/coherent_defect_wave/chiral_joint/peter_weyl_links/model.py',
           'model/coherent_defect_wave/chiral_joint/atomic_spin_geometry/momentum.py',
           'model/particle_common_contract/provenance.py'}
    for relative in import_closure(seeds):pins[relative]=digest(ROOT/relative)
    channel=parents['interacting_charged_response']['same_H_interacting_channels']['ec']['full_interacting_energy_readout']
    return {'schema':'oph.interacting_charge_conjugation.v1','exact_source_controls':exact,
      'same_source_spectral_consequences':{
        'source_member':'N=2; lambda=2^-34; full PW119; actual unique physical ground',
        'all_native_origin_addition_removal_spectral_measures_equal':True,
        'every_origin_mode_mean_occupation':'1/2',
        'all_five_ground_family_pump_failure_probabilities':'1/2',
        'ec_removal_bin_probability_inherited':channel['conditional_probabilities']['line'],
        'paired_addition_bin_has_identical_source_measure':True,
        'raw_L_transformation':'L -> 16-L','centered_raw_L':'L-8',
        'source_opposite_charge_energies_equal_without_calibration':True},
      'reference_boundary':{'old_reference_H':'4 N_ref','particle_hole_transformed_reference_H':'120-4 N_ref',
        'old_reference_clock_charge_conjugation_invariant':False,
        'same_reference_total_energy_equality_claimed':False,
        'conjugate_boundary_charges':['+1','-1'],
        'equal_total_sector_energies_require_explicit_equal_reference_cost':True,
        'new_reference_extension_implemented':False},
      'scope':{'finite_group_odd_channel':'G1','continuous_physical_spin_proved':False,
        'full_native_CPT_or_Lorentz_claimed':False,'empirical_particle_identity':False,
        'changed_H_or_transport_added':False},
      'input_sha256':dict(sorted(pins.items())),
      'owned_source_sha256':{p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}}


REFERENCE_SHA='0bbfd57371b0b00d74caacc6321d54a8ac073fa1296633541b64fa46e9225606'

def validate_pins(pins):
    for name,expected in pins.items():
        path=Path(name)
        if path.is_absolute() or '..' in path.parts:raise ValueError('repository-relative source pin required')
        if digest(ROOT/path)!=expected:raise ValueError('changed source input: '+name)

def compare_research(value):
    path=HERE/'research_reference.json'
    if digest(path)!=REFERENCE_SHA:raise ValueError('immutable research receipt changed')
    reference=json.loads(path.read_text());validate_pins(reference['input_sha256'])
    ignored={'schema','input_sha256','owned_source_sha256'}
    fresh={k:v for k,v in value.items() if k not in ignored}
    # The public spelling makes the antiunitary interpretation explicit.
    fresh=json.loads(json.dumps(fresh))
    links=fresh['exact_source_controls']['links']
    links['link_CP_permutation']=links.pop('link_charge_reflection_permutation')
    old={k:v for k,v in reference.items() if k not in ignored}
    if fresh!=old:raise ValueError('fresh exact source calculation differs from research result')
    return {'reference_sha256':REFERENCE_SHA,'all_inherited_scientific_fields_equal':True,
            'only_field_rename':'link_CP_permutation -> link_charge_reflection_permutation',
            'fresh_exact_source_identities_replayed':True}

def certificate():
    value=raw_certificate();comparison=compare_research(value)
    value.pop('owned_source_sha256')
    seeds={p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.add((HERE/'research_reference.json').relative_to(ROOT).as_posix())
    for relative in import_closure(seeds):value['input_sha256'][relative]=digest(ROOT/relative)
    value['input_sha256']=dict(sorted(value['input_sha256'].items()))
    value['independent_research_comparison']=comparison
    value['antiunitary_time_and_interpretation']={
        'source_H_preserved':True,'source_time_evolution':'Theta U(t) Theta^-1=U(-t)',
        'unitary_physical_CP_symmetry_claimed':False,'physical_CPT_or_spin_statistics_claimed':False,
        'absence_of_measured_weak_CP_violation_inferred':False,
        'complex_family_CKM_or_weak_CP_observables_derived':False,
        'Yukawa_domain':'the declared real native coefficients'}
    return value
