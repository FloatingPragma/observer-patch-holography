"""Exact unit transport without conflating different particle constructions."""
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[2]
PARENTS = ('openwave_single_valley_leptons', 'openwave_neutrino_calibration',
           'particle_joint_mass_identity', 'openwave_vectorlike_flavour',
           'openwave_common_mass_source')

@dataclass(frozen=True)
class UnitMap:
    source: str
    reference: str
    clock: str
    GeV_per_source: F
    def __post_init__(self):
        if not self.source or not self.reference or not self.clock or self.GeV_per_source <= 0:
            raise ValueError('a named source/reference/clock and positive energy unit are required')
    def interval(self, endpoints):
        lo, hi = map(F, endpoints)
        if lo > hi: raise ValueError('reversed energy interval')
        return [str(lo*self.GeV_per_source), str(hi*self.GeV_per_source)]
    def inverse_time(self, source_time):
        """Return GeV^-1 time; seconds additionally require hbar in GeV s."""
        return F(source_time)/self.GeV_per_source
    def assert_same_source(self, other):
        if (self.source,self.reference,self.clock)!=(other.source,other.reference,other.clock):
            raise ValueError('unit conversion does not transfer a source, reference or clock')


def anchored_interval(energy, anchor_energy, anchor_GeV):
    """Positive exact interval arithmetic; the caller states correlations."""
    e0,e1=map(F,energy);a0,a1=map(F,anchor_energy);g0,g1=map(F,anchor_GeV)
    if not (0<=e0<=e1 and 0<a0<=a1 and 0<g0<=g1):
        raise ValueError('nonnegative energy and positive ordered anchor intervals required')
    return [str(g0*e0/a1),str(g1*e1/a0)]


def transport_complete_measure(unit, labels, energies, weights, failure, source_time):
    """Push forward every spectral outcome, retaining original-input weights."""
    if len(labels)!=len(energies) or len(labels)!=len(weights) or len(set(labels))!=len(labels):
        raise ValueError('complete distinct labels, energies and weights must have equal lengths')
    w=list(map(F,weights));f=F(failure)
    if min(w+[f])<0 or sum(w,F(0))+f!=1:
        raise ValueError('the original input, including failure, must have total probability one')
    e=list(map(F,energies));t=F(source_time);tg=unit.inverse_time(t)
    eg=[x*unit.GeV_per_source for x in e]
    assert all(x*t==y*tg for x,y in zip(e,eg))
    return {'labels':list(labels),'energy_GeV_exact':list(map(str,eg)),
            'original_input_weights_exact':list(map(str,w)), 'failure_probability_exact':str(f),
            'time_GeV_inverse_exact':str(tg),'total_probability_exact':'1',
            'all_phase_products_preserved_exactly':True}


def read_parent(name, root=ROOT):
    p=root/'model'/name/'results.json';d=json.loads(p.read_text())
    pins=d.get('source_sha256',d.get('source_bindings'))
    if not isinstance(pins,dict) or not pins: raise ValueError(f'{name}: no source binding')
    for rel,sha in pins.items():
        target=(root/rel).resolve()
        if not target.is_relative_to(root.resolve()) or not target.is_file():
            raise ValueError(f'{name}: invalid source path {rel}')
        if hashlib.sha256(target.read_bytes()).hexdigest()!=sha:
            raise ValueError(f'{name}: changed source {rel}')
    if 'all_checks_pass' in d and d['all_checks_pass'] is not True:
        raise ValueError(f'{name}: failed producer')
    if 'verdict' in d and not str(d['verdict']).startswith('PASS'):
        raise ValueError(f'{name}: failed producer verdict')
    return d,pins


def exact_unit_execution():
    # A noncommuting four-outcome instrument tests phase and probability
    # preservation, rather than comparing a diagonal matrix to itself.
    H=np.array([[0,1,0,0],[1,2,1j,0],[0,-1j,3,2],[0,0,2,4]],complex)
    psi=np.array([1,1j,2,-1],complex)/np.sqrt(7)
    unit=UnitMap('literal-four-level-control','declared coherent state','source t',F('0.00051099895069'))
    t=F(7,5);reference=F(13,7)
    source=expm(-1j*H*float(t))@psi
    physical=expm(-1j*(float(unit.GeV_per_source)*(H-float(reference)*np.eye(4)))*float(unit.inverse_time(t)))@psi
    probabilities=np.abs(source)**2;other=np.abs(physical)**2
    error=float(np.max(abs(probabilities-other)))
    assert error<2e-14 and abs(sum(probabilities)-1)<2e-14
    complete=transport_complete_measure(unit,['lower','middle','upper'],['-2','1/3','7'],['1/8','1/4','1/8'],'1/2',t)
    return {'four_level_probability_error':error,'complete_probabilities':probabilities.tolist(),
            'reference_shift_is_only_a_global_phase':True,'exact_spectral_measure':complete,
            'theorem':'H_GeV=u(H-E_ref I), t_GeV^-1=t/u; all complete POVM probabilities are invariant',
            'physical_clock_identification_from_unit_conversion':False}


def base_catalogue(parents):
    single=parents['openwave_single_valley_leptons'];neutral=parents['openwave_neutrino_calibration']
    joint=parents['particle_joint_mass_identity'];vector=parents['openwave_vectorlike_flavour'];common=parents['openwave_common_mass_source']
    rows=[]
    for r in single['charged_lepton_GeV']:
        rows.append({'particle':r['species'],'source_member':'native Higgs-Wilson charged family',
            'source_receipt':'model/openwave_single_valley_leptons/results.json',
            'observable':'unique free colourless charged band minimum',
            'mass_GeV':r['GeV'],'mass_GeV_exact':r['GeV_exact'],
            'numerical_enclosure_GeV':r['analyzer_enclosure_GeV'],
            'enclosure_role':'analyzer numerical enclosure; no experimental or model uncertainty included',
            'input_ancestry':['electron mass anchor','retrospectively selected analyzer phase 2/9','declared Higgs-Wilson coefficient r=1'],
            'physical_interacting_pole_identified':False,'antiparticle_same_mass':True,
            'openwave_rows':[3,11]})
    cat=neutral['calibrated_neutral_catalogue']
    for r in cat['rows']:
        rows.append({'particle':f"nu_{r['mass_eigenstate']}", 'source_member':'calibrated native Dirac neutral family',
            'source_receipt':'model/openwave_neutrino_calibration/results.json',
            'observable':'Dirac mass eigenvalue; flavour states are superpositions',
            'mass_GeV':r['conditional_mass_GeV'],'numerical_enclosure_GeV':r['mass_GeV_interval'],
            'enclosure_role':'numerical square-root enclosure at fixed supplied oscillation inputs',
            'input_ancestry':['NuFIT 6.0 IC24/SK six oscillation inputs','normal ordering','declared m1=0.001 eV'],
            'physical_interacting_pole_identified':False,'antiparticle_same_mass':True,
            'sterile_right_chirality_is_same_mass_species':True,'openwave_rows':[9,10]})
    reference=[]
    for r in joint['conditional_GeV_targets']:
        box=list(map(F,r['conditional_rest_coefficient_GeV_interval']))
        assert box[0]>=0 and box[0]**2<=F(r['mass_squared_GeV2_exact'])<=box[1]**2
        reference.append({'particle':r['species'],'mass_GeV':float(sum(box)/2),
            'numerical_enclosure_GeV':list(map(str,box)), 'source_member':'calibrated quadratic comparison',
            'source_receipt':'model/particle_joint_mass_identity/results.json',
            'observable':r['input_convention'],'input_ancestry':r['empirical_inputs'],
            'physical_matching_required':r['physical_matching'],'native_prediction':False,
            'is_measured_pole_mass':False})
    # The optional Weinberg comparison has a different m1 from the Dirac
    # source. It is retained, not silently used to replace the Dirac masses.
    assert len(reference)==17
    extra=[]
    for r in vector['quark_branches']:
        extra.append({'channel':r['channel'],'family':r['family'],
            'light_GeV':float(r['conditional_light_GeV']),'heavy_GeV':float(r['conditional_heavy_GeV']),
            'source_member':'declared vectorlike flavour control G=M=1',
            'observed_quark_mass_identification':False,'new_particle_prediction':False,
            'retained_original_occupation_outcomes':len(r['complete_channel_records']['full_CAR_addition_records']['complete_unconditional_native_occupation_probabilities']),
            'creation_failure_probability':r['complete_channel_records']['full_CAR_addition_records']['creation_failure_probability']})
    common_rows=common['same_source_free_Dirac_branches']
    for r in common_rows[:3]:
        item=next(x for x in rows if x['particle']==r['channel'])
        assert F(item['mass_GeV_exact'])==F(r['energy_GeV_exact'])
    baryon=[]
    for r in common['complete_coupled_baryon_energy_records']['complete_energy_window_instrument_bounds']:
        if r['half_width_electron_units']!='1/2':continue
        mean=F(r['mean_energy_GeV']);width=F(r['half_width_GeV'])
        baryon.append({'particle':r['assigned_channel'],'source_member':'four-input axial matter plus full weak contact',
            'source_receipt':'model/openwave_common_mass_source/results.json',
            'observable':'prepared-state spectral mean and complete energy-window bound',
            'mean_energy_GeV':float(mean),'energy_window_GeV':[float(mean-width),float(mean+width)],
            'conditional_window_probability_lower':float(F(r['conditional_inside_probability_lower'])),
            'original_window_probability_lower':float(F(r['original_inside_probability_lower'])),
            'preparation_failure_probability':float(F(r['preparation_failure_probability'])),
            'physical_nucleon_pole_identified':False,'input_ancestry':['electron anchor','proton mass','neutron mass','charged pion mass','nucleon axial coupling'],
            'openwave_rows':[13,14,15,25]})
    return rows,reference,extra,baryon


def certificate():
    parents={name:read_parent(name)[0]for name in PARENTS}
    rows,reference,extra,baryon=base_catalogue(parents)
    return {'schema':'oph.openwave_mass_catalogue.v1','verdict':'PASS: exact GeV transport with source-qualified particle accounting',
        'unit_transport':exact_unit_execution(),'native_and_calibrated_Dirac_masses':rows,
        'interacting_baryon_energy_windows':baryon,'complete_calibrated_comparison_catalogue':reference,
        'vectorlike_branches_retained':extra,
        'regulator_modes':parents['openwave_single_valley_leptons']['former_corner_modes_retained'],
        'neutral_Higgs_Wilson_zero_modes_retained':48,
        'source_joins':{'charged_rest_numbers_equal_in_90_and_192_mode_members':True,
            'all_sources_identified_as_one_interacting_theory':False,
            'unit_conversion_preserves_complete_instruments':True,
            'unit_conversion_proves_physical_particle_identity':False},
        'spectrum_accounting':{'Standard_Model_comparison_species':17,'charged_native_families':3,
            'calibrated_Dirac_neutral_masses':3,'vectorlike_light_and_heavy_branches':12,
            'physical_quark_poles_supplied':False,'transverse_photons_on_common_192_source_supplied':False,
            'weak_vector_and_Higgs_interacting_poles_on_common_192_source_supplied':False,
            'all_19_OpenWave_partials_validated':False},
        'scope':'Exact positive unit transport and numerical source accounting; finite spectra, calibrated coefficients and physical poles retain their different meanings.'}
