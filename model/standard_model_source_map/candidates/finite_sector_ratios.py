"""Compare actual finite sector thresholds without labelling them particle masses."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

from model.particle_common_contract import contract
from tools import strict_json

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
SOURCE_INPUTS = tuple(str((HERE / name).relative_to(ROOT)) for name in
    ('finite_sector_ratios.py', 'test_finite_sector_ratios.py', 'FINITE_SECTOR_RATIOS.md')) + (
    'model/particle_primary_sector_threshold/results.json',
    'model/particle_primary_branch_certificate/results.json',
    'model/particle_primary_charge_conjugation/results.json',
    'model/particle_primary_sector_threshold/PROOF.md',
    'model/particle_primary_branch_certificate/PROOF.md',
    'model/particle_primary_charge_conjugation/PROOF.md',
)


def ratio_enclosure(numerator, denominator, ground, same_sector=False):
    """Use one shared energy zero; endpoint extrema are exact rational bounds."""
    li, ui = map(F, numerator)
    lj, uj = map(F, denominator)
    l0, u0 = map(F, ground)
    if not l0 <= u0 < min(li, lj) or li > ui or lj > uj:
        raise ValueError('positive sector gaps on a common ground interval required')
    if same_sector:
        if (li, ui) != (lj, uj):
            raise ValueError('same sector needs the same energy enclosure')
        return F(1), F(1)
    return (min((li-e)/(uj-e) for e in (l0, u0)),
            max((ui-e)/(lj-e) for e in (l0, u0)))


def _receipt(path):
    value = strict_json.load(ROOT / path)
    for name, expected in value['source_sha256'].items():
        if sha256((ROOT / name).read_bytes()).hexdigest() != expected:
            raise ValueError('source receipt input changed: '+name)
    bounds = ([value['selected_member_bound']] if 'selected_member_bound' in value
              else value.get('members', []))
    for bound in bounds:
        if contract.digest(bound['result']) != bound['result_sha256']:
            raise ValueError('source result digest changed')
    return value


def calculate(sectors, branch):
    keys = ('member_id', 'source_id', 'state_id', 'reference_id', 'clock_id')
    if any(sectors[k] != branch[k] for k in keys):
        raise ValueError('sector comparison changed source, state, reference or clock')
    l0 = max(F(sectors['global_ground_lower']), F(branch['full_source_ground_interval'][0]))
    u0 = min(F(sectors['global_ground_upper']), F(branch['full_source_ground_interval'][1]))
    rows = []
    for group in sectors['sector_results']:
        signed = group['signed_sectors']
        first = signed[0]
        lower = F(first['physical_sector_energy_lower'])
        upper = F(first['sector_bottom_excitation_interval'][1]) + F(sectors['global_ground_lower'])
        if first['delta'] == branch['physical_native_increment']:
            lower = max(lower, F(branch['full_odd_sector_bottom_interval'][0]))
            upper = min(upper, F(branch['full_odd_sector_bottom_interval'][1]))
        if not l0 <= u0 < lower <= upper:
            raise ValueError('no positive finite threshold enclosure')
        if len(signed) != 2 or signed[1]['delta'] != [-x for x in first['delta']]:
            raise ValueError('conjugate occupation pair changed')
        rows.append({'source_id': group['label'], 'native_increment': first['delta'],
                     'conjugate_increment': signed[1]['delta'],
                     'CAR_parity': sum(abs(x) for x in first['delta']) % 2,
                     'total_native_Q6': sum(q*x for q,x in zip(sectors['native_charge_order'],first['delta'],strict=True)),
                     'sector_bottom_energy_interval': [lower, upper],
                     'sector_bottom_excitation_interval': [lower-u0, upper-l0]})
    ratios = []
    for i in rows:
        for j in rows:
            interval = ratio_enclosure(i['sector_bottom_energy_interval'],
                j['sector_bottom_energy_interval'], (l0,u0), i['source_id']==j['source_id'])
            ratios.append({'numerator_id':i['source_id'], 'denominator_id':j['source_id'],
                           'finite_excitation_ratio_interval':list(interval)})
    upper = max(x['finite_excitation_ratio_interval'][1] for x in ratios)
    odd_ids = {x['source_id'] for x in rows if x['CAR_parity']==1}
    odd_upper = max(x['finite_excitation_ratio_interval'][1] for x in ratios
                   if x['numerator_id'] in odd_ids and x['denominator_id'] in odd_ids)
    if not upper < F(529,25) or not odd_upper < F(1531,100):
        raise ArithmeticError('advertised outward ratio bounds changed')
    return contract.serial({
        'schema':'oph.finite_sector_ratio_control.v1', **{k:sectors[k] for k in keys},
        'ground_energy_interval':[l0,u0], 'source_sector_rows':rows, 'ordered_ratio_rows':ratios,
        'all_three_family_ratio_upper':upper, 'all_three_family_simple_ratio_upper':F(529,25),
        'odd_family_ratio_upper':odd_upper, 'odd_family_simple_ratio_upper':F(1531,100),
        'same_sector_and_charge_conjugate_threshold_ratio':1,
        'one_positive_common_unit_conversion_cancels':True,
        'observational_comparison':{
            'source':'https://physics.nist.gov/cuu/Constants/Table/allascii.txt',
            'quantity':'2022 CODATA muon-electron mass ratio',
            'central_value':'206.7682827', 'standard_uncertainty':'0.0000046',
            'used_to_construct_source_or_bounds':False,
            'comparison_is_to_reported_value_with_uncertainty_not_exact_truth':True},
        'direct_common_scale_identification_with_charged_lepton_hierarchy_excluded_at_member_zero':True,
        'exclusion_applies_to_other_capacities_or_other_source_excitations':False,
        'exclusion_is_a_no_go_for_OPH':False,
        'ratios_are_point_eigenvalue_calculations':False,
        'physical_dispersion_or_particle_species_identified':False,
        'finite_threshold_intervals_are_physical_rest_masses':False,
        'physical_mass_table_GeV':[]})


def certificate():
    sector = _receipt(SOURCE_INPUTS[3])['selected_member_bound']['result']
    branch = _receipt(SOURCE_INPUTS[4])['selected_member_bound']['result']
    conjugation = _receipt(SOURCE_INPUTS[5])['members'][0]['result']
    member = contract.member(0)
    if any(sector[k] != member['ids'][k] for k in ('source_id','state_id','reference_id','clock_id')):
        raise ValueError('selected source identity changed')
    if any(sector[k] != conjugation[k] for k in
           ('member_id','source_id','state_id','reference_id','clock_id')):
        raise ValueError('charge conjugation uses a different source domain')
    if not (conjugation['opposite_occupation_sector_spectra_equal']
            and conjugation['original_H_and_energy_zero_preserved']):
        raise ValueError('conjugate threshold equality is not certified')
    result = calculate(sector,branch)
    return {'schema':'oph.finite_sector_ratio_certificate.v1', 'result':result,
            'result_sha256':contract.digest(result),
            'source_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_INPUTS}}
