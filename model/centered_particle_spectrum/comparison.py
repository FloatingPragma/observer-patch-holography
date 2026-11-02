"""Retrospective fixed-fixture test of the explicit finite-gap mass dictionary."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json

ROOT=Path(__file__).resolve().parents[2]
TARGET='model/particle_species_detector/physical_targets.json'


def minimax_floor(source_ratio,target_ratio):
    if type(source_ratio) not in (int,F) or type(target_ratio) not in (int,F) or source_ratio<=0 or target_ratio<=0:
        raise ValueError('positive exact ratios required')
    r,q=F(source_ratio),F(target_ratio)
    return max(F(0),(q-r)/(q+r))


def certificate(spectrum):
    rows=spectrum['rows']
    if len(rows)!=228 or not all(r['strictly_positive_weight_proved'] for r in rows):
        raise ValueError('the complete certified centered response is required')
    lo=F(rows[0]['excitation_energy_interval'][0]);hi=F(rows[-1]['excitation_energy_interval'][1])
    if not 0<lo<hi or not hi/lo<11:raise ArithmeticError('centered visible energy envelope changed')
    targets=json.loads((ROOT/TARGET).read_text())
    if targets['edition']!='PDG 2025 fixed comparison baseline':raise ValueError('fixed comparison edition changed')
    leptons=targets['charged_leptons']
    if [r['name'] for r in leptons]!=['electron','muon','tau']:raise ValueError('charged-lepton fixture changed')
    electron=F(leptons[0]['mass_MeV_c2']);electron_error=F(leptons[0]['reported_uncertainty_MeV_c2'])
    if not 0<=electron_error<electron:raise ValueError('invalid electron uncertainty box')
    comparisons=[]
    for row in leptons[1:]:
        mass=F(row['mass_MeV_c2']);error=F(row['reported_uncertainty_MeV_c2'])
        if not 0<=error<mass:raise ValueError('invalid heavy-lepton uncertainty box')
        ratio=(mass-error)/(electron+electron_error);bound=minimax_floor(hi/lo,ratio)
        if bound<=0:raise ArithmeticError('direct finite-gap dictionary is not excluded')
        comparisons.append({'light':'electron','heavy':row['name'],'target_ratio_lower':str(ratio),
                            'source_ratio_upper':str(hi/lo),'relative_error_lower':str(bound),
                            'relative_error_diagnostic':float(bound),'exact_assignment_possible':False,
                            'sources':[leptons[0]['source'],row['source']]})
    return {'dictionary':'Assign any distinct visible centered-cube positive CAR excitation gaps to charged-lepton rest energies using one positive multiplicative scale, without an additive offset.',
            'status':'retrospective fixed-fixture conditional dictionary exclusion',
            'source_energy_envelope':[str(lo),str(hi)],'source_ratio_strict_upper':'11',
            'target_path':TARGET,'target_sha256':hashlib.sha256((ROOT/TARGET).read_bytes()).hexdigest(),
            'comparisons':comparisons,'all_three_distinct_frequency_assignments_excluded':True,
            'finite_gap_is_relativistic_mass_assumed_only_for_this_test':True,
            'physical_rest_masses_derived':False,'new_eligible_empirical_comparisons':0,
            'unobserved_states_composites_continuum_or_OPH_excluded':False}
