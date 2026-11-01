"""Exact rejection of a stated finite-gap charged-lepton dictionary.

The benchmark is retrospective. It does not equate finite gaps with rest
masses, reject a continuum sector, or fit any source coefficient.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from model.retarded_matter_response import spectrum


ROOT=Path(__file__).resolve().parents[2]
TARGET_PATH='model/particle_species_detector/physical_targets.json'


def minimax_relative_error_lower(energy_lower,energy_upper,mass_ratio_lower):
    if (any(type(x) not in (int,F) for x in
            (energy_lower,energy_upper,mass_ratio_lower)) or
            not 0<energy_lower<=energy_upper or mass_ratio_lower<=0):
        raise ValueError('positive ordered exact energy and mass-ratio bounds required')
    spread=F(energy_upper)/F(energy_lower)
    ratio=F(mass_ratio_lower)
    return max(F(0),(ratio-spread)/(ratio+spread))


def certificate():
    targets=json.loads((ROOT/TARGET_PATH).read_text())
    if targets['edition']!='PDG 2025 fixed comparison baseline':
        raise ValueError('fixed retrospective comparison edition changed')
    leptons=targets['charged_leptons']
    if [x['name'] for x in leptons]!=['electron','muon','tau']:
        raise ValueError('the retained charged-lepton dictionary changed')
    local=spectrum.certificate()
    first,last=local['rows'][0],local['rows'][-1]
    lower=F(first['gap_above_ground_interval'][0])
    upper=F(last['gap_above_ground_interval'][1])
    # The other nonzero single-charge local ground response is lambda=sqrt3/16.
    # 1<sqrt3<2 puts it between1/16 and1/8, strictly inside this envelope.
    if not (F(first['gap_above_ground_interval'][1])<F(1,16) and
            F(last['gap_above_ground_interval'][0])>F(1,8)):
        raise ArithmeticError('up-channel charged line left the checked energy envelope')
    spread=upper/lower
    if not spread<62:raise ArithmeticError('the complete finite-gap spread changed')
    e=leptons[0]
    electron=F(e['mass_MeV_c2']);electron_error=F(e['reported_uncertainty_MeV_c2'])
    if not 0<=electron_error<electron:raise ValueError('invalid electron comparison interval')
    rows=[]
    for heavy in leptons[1:]:
        mass=F(heavy['mass_MeV_c2']);error=F(heavy['reported_uncertainty_MeV_c2'])
        if not 0<=error<mass:raise ValueError('invalid heavy comparison interval')
        ratio=(mass-error)/(electron+electron_error)
        bound=minimax_relative_error_lower(lower,upper,ratio)
        if not bound>0:raise ArithmeticError('finite-gap dictionary is no longer excluded')
        rows.append({'light':'electron','heavy':heavy['name'],
                     'mass_ratio_lower_on_reported_uncertainty_box':str(ratio),
                     'energy_ratio_upper':str(spread),
                     'common_positive_scale_minimax_relative_error_lower':str(bound),
                     'diagnostic_error_lower':float(bound),
                     'exact_assignment_possible':False,
                     'sources':[e['source'],heavy['source']]})
    return {'dictionary_tested':'Assign visible positive single-charge ground-probe gaps on the fixed cube to charged-lepton rest energies using one positive multiplicative unit factor and no additive shift.',
            'comparison_status':'retrospective fixed PDG2025 benchmark; central value plus/minus one reported uncertainty is a declared comparison box, not an exact physical confidence guarantee',
            'target_path':TARGET_PATH,
            'target_sha256':hashlib.sha256((ROOT/TARGET_PATH).read_bytes()).hexdigest(),
            'source_energy_envelope':[str(lower),str(upper)],
            'source_spread_upper':str(spread),
            'all_78_down_lines_and_positive_up_line_included':True,
            'pair_tests':rows,
            'inequality_proof':'If both relative errors are<=e<1, their normalized predictions obey (1-e)/(1+e)<=R/Q; hence e>=(Q-R)/(Q+R) whenQ>R.',
            'finite_gap_as_rest_energy_dictionary_rejected':True,
            'unit_factor_or_source_parameters_fitted':False,
            'finite_gaps_proved_to_be_relativistic_rest_masses':False,
            'unobserved_sectors_composites_or_continuum_excluded':False,
            'whole_OPH_model_falsified':False,
            'eligible_physical_comparisons_added':0}
