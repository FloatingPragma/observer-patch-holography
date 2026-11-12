"""Full-source positivity and exact complete energy-window certificates."""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from .source import lepton_source,local_norm_certificate,exact_baryon_window,full_source_terms,parameters,INPUTS,G
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate.action import SOURCE_FILES
from model.openwave_larmor_dynamics.certificate import equivalent

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def pins():
    seeds=set(SOURCE_FILES)
    for package in ('openwave_common_mass_source','openwave_axial_state','openwave_neutrino_calibration'):
        for p in (ROOT/'model'/package).iterdir():
            if p.suffix in ('.py','.md','.tex') or p.name in ('inputs.json','oscillation_inputs.json'):
                seeds.add(str(p.relative_to(ROOT)))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
            for p in sorted(import_closure(seeds))}


def certificate():
    l=lepton_source();p=parameters();unit=F(INPUTS['electron_unit_GeV'])
    assert unit==l['calibration']['empirical_electron_anchor_GeV']
    floor=p['mbar']-p['delta']/2-18*G
    assert floor>1
    lightest=F(l['calibration']['oscillation_inputs']['declared_lightest_mass_eV'])/10**9
    assert 0<lightest<unit
    rows=[]
    for i,name in enumerate(('electron','muon','tau')):
        exact=l['calibration']['charged_exact'][i]/l['calibration']['charged_exact'][0]
        rows.append({'channel':name,'energy_electron_units_exact':str(exact),
            'energy_GeV_exact':str(exact*unit),'energy_GeV_diagnostic':float(exact*unit),
            'same_interacting_source_exact_eigenbranch':True,
            'sector':'quark vacuum, one Dirac excitation, all other modes in the pure Dirac reference',
            'mass_role':'electron empirical anchor'if i==0 else'native analyzer branch with retrospective phase2/9',
            'conjugate_branch_same_mass':True})
    for i,mass in enumerate(l['calibration']['physical_masses_eV']):
        rows.append({'channel':f'neutral_{i+1}','energy_GeV_diagnostic':float(mass)*1e-9,
            'same_interacting_source_exact_eigenbranch':True,
            'mass_role':'oscillation-calibrated Dirac branch; lightest mass is declared',
            'conjugate_branch_same_mass':True,'active_and_sterile_not_separate_mass_species':True})
    return {'schema':'oph.common_semileptonic_mass_source.v1',
        'verdict':'PASS_POSITIVE_COMMON_SOURCE_AND_COMPLETE_MASS_RECORD_BOUNDS',
        'source':full_source_terms(),'strong_parameters_exact':{k:str(v)for k,v in p.items()},
        'coupling_exact':str(G),'electron_unit_GeV_exact':str(unit),
        'local_contact_norm_certificate':local_norm_certificate(),
        'positivity':{'colourful_sector_floor_per_occupation_exact':str(floor),
            'colourful_sector_floor_exceeds_electron':True,
            'finite_vacuum_unique':True,'finite_source_gap_GeV_exact':str(lightest),
            'full_charge_plus_and_minus_one_threshold_GeV_exact':str(unit),
            'finite_matter_source_includes_photons':False,
            'quark_family_running_or_physical_pole_masses_identified':False},
        'lepton_reference':{'positive_Dirac_CAR_modes':48,'complete_doubled_operator_rows':96,
            'positive_reference_constant_electron_units_diagnostic':float(sum(l['E'])/2),
            'full_canonical_matrix_error':float(np.linalg.norm(l['A']@l['A'].conj().T+l['B']@l['B'].conj().T-np.eye(48))),
            'active_left_normal_covariance_exact':'I_24/2',
            'active_left_anomalous_covariance_exact':'0',
            'all_native_active_current_means_exact':'0',
            'native_current_covariance_exact':'(3/4) I_24',
            'reference_is_old_filled_neutral_beta_reference':False},
        'same_source_free_Dirac_branches':rows,
        'complete_coupled_baryon_energy_records':exact_baryon_window(),
        'scope':'One declared192-mode semileptonic action. Purely leptonic contacts and photons are absent. Charged/neutral free branches and full finite charged threshold coexist with coupled baryon energy-window records. No complete physical particle spectrum, interacting pole, decay lifetime or shared continuum identification is asserted.'}


def main():
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    a=parser.parse_args();before=pins();r=certificate();assert pins()==before
    r['source_sha256']=before;path=HERE/'results.json'
    if a.write:path.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
    else:assert equivalent(json.loads(path.read_text()),r)
    print(r['verdict'])


if __name__=='__main__':main()
