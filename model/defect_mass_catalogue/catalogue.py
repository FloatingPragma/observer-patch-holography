"""Every Standard Model target has an explicit source attachment status."""
from pathlib import Path
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
import sympy as s
from .calibration import Anchor, common_scale, convert_energy, ratio, strings

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def read(name):
    return json.loads((ROOT/name).read_text())


def build():
    lepton=read('model/lepton_cycle/physical_comparison.json')
    matter=read('model/common_source/calibrated_matter_results.json')
    quark=read('model/native_quark_response/results.json')
    scalar=read('model/defect_scalar_modes/results.json')
    joined=read('model/defect_mass_join/results.json')
    family=read('model/defect_family_vertices/results.json')
    neutral_records=family['neutral_channels']['actual_gauge_singlet_neutral_composite']['complete_records']
    neutral_lines=[row for setting in neutral_records['recorded_settings'] for row in setting['complete_target_energy_records']]
    neutral_min=min(neutral_lines,key=lambda row:F(row['excitation_energy_interval'][0]))['excitation_energy_interval']
    with localcontext() as context:
        context.prec=50
        e=Decimal(lepton['experimental_inputs']['electron_MeV']['value'])/1000
        de=Decimal(lepton['experimental_inputs']['electron_MeV']['standard_uncertainty'])/1000
        ratios=[Decimal(1),Decimal(lepton['architecture']['muon_electron_ratio']),Decimal(lepton['architecture']['tau_electron_ratio'])]
        lepton_rows=[{'species':name,'value_GeV':str(e*r),
                     'electron_anchor_standard_uncertainty_only_GeV':str(de*r),
                     'role':'empirical_anchor' if i==0 else 'conditional_fixed_effect_ratio',
                     'source_to_relativistic_rest_map_proved':False}
                    for i,(name,r) in enumerate(zip(['electron','muon','tau'],ratios))]
    native={}
    for family in ['uc','dc','Q']:
        records=quark['full_interacting_prepared_records'][family]['records']
        native[family]=[{'index':row['index'],
            'source_energy_interval':row['source_excitation_relative_actual_ground'],
            'original_attempt_probability_lower':row['original_attempt_probability_lower']}
            for row in records]
    # Test a specified dictionary: record zero in uc and dc -> the u and d
    # running masses at the same frozen SMDR scale. No record is omitted
    # from the stored native catalogue; this dictionary is a control only.
    clock='native-eight-cell:H0+2^-34T;hbar=1'
    target={name:F(str(matter['running_mass_parameters_numeric_GeV'][group][0]))
            for name,group in [('u','up'),('d','down')]}
    anchors=[Anchor(clock,tuple(native[f][0]['source_energy_interval']),
                    (target[name],target[name]),'hypothetical record-zero to MSbar mass parameter',
                    'pure MSbar Standard Model; mu=91.1876 GeV')
             for name,f in [('u','uc'),('d','dc')]]
    scale=common_scale(anchors)
    if scale is not None:
        raise ArithmeticError('the declared u/d dictionary unexpectedly admits one scale')
    u_scale=anchors[0].scale_interval()
    calibrated_native={f:[dict(row,energy_GeV_interval=strings(convert_energy(clock,clock,
         row['source_energy_interval'],u_scale))) for row in rows] for f,rows in native.items()}
    target_rows=[]
    charged=[('up','up',0),('charm','up',1),('top','up',2),
             ('down','down',0),('strange','down',1),('bottom','down',2),
             ('electron','charged_lepton',0),('muon','charged_lepton',1),('tau','charged_lepton',2)]
    for name,group,index in charged:
        target_rows.append({'species':name,'native_channel':group,
            'physical_source_rest_mass_GeV':None,
            'empirically_calibrated_action_parameter_GeV':matter['running_mass_parameters_numeric_GeV'][group][index],
            'parameter_convention':'pure MSbar Standard Model at mu=91.1876 GeV; y*v/sqrt(2)',
            'attachment':'native contact tensor exists; three source-selected families and physical pole map unproved',
            'antiparticle_mass_equality':'conditional conjugate calibrated action; not a second fit',
            'native_source_prediction':False})
    neutral=matter['optional_phenomenological_neutrino_extension']
    for i,value in enumerate(neutral['neutrino_mass_parameters_eV_numeric']):
        target_rows.append({'species':f'neutrino_mass_eigenstate_{i+1}',
            'physical_source_rest_mass_GeV':None,
            'empirically_calibrated_action_parameter_GeV':value/1e9,
            'parameter_convention':'optional phenomenological Weinberg operator; normal ordering and lightest mass zero supplied',
            'attachment':'minimal R15 has no right-handed neutrino; complete local neutral-composite energy records exist, physical neutrino dispersion and masses unproved',
            'native_source_prediction':False})
    bosonic=matter['bosonic_tree_control']
    for name,key in [('W','charged_W_mass_squared'),('Z','Z_mass_squared'),('Higgs','Higgs_radial_mass_squared')]:
        target_rows.append({'species':name,'physical_source_rest_mass_GeV':None,
            'empirically_calibrated_action_parameter_GeV':float(s.sqrt(s.sympify(bosonic[key]))),
            'parameter_convention':'tree comparison coefficient at matched v; not a measured pole or the loop action Hessian',
            'attachment':'dynamical gauge and broken-vacuum response not supplied by the common defect source',
            'native_source_prediction':False})
    for name in ['photon','gluon']:
        target_rows.append({'species':name,'physical_source_rest_mass_GeV':None,
            'empirically_calibrated_action_parameter_GeV':0,
            'parameter_convention':'unbroken calibrated gauge-action quadratic direction',
            'attachment':'no photon or gluon propagation operator in the common defect source',
            'native_source_prediction':False})
    if len(target_rows)!=17 or len({r['species'] for r in target_rows})!=17:
        raise ArithmeticError('complete target catalogue failed')
    return {'schema':'oph.defect_mass_catalogue.v1',
        'conditional_electron_anchored_leptons':lepton_rows,
        'lepton_comparison':{'muon_residual_ppm':lepton['muon']['relative_residual_ppm'],
            'tau_residual_ppm':lepton['tau']['relative_residual_ppm'],
            'uncertainty':'anchor-only errors exclude model/attachment error; no confidence interval for the theory'},
        'single_clock_quark_test':{'source_clock':clock,'specified_dictionary':'record zero uc -> u, record zero dc -> d',
            'source_ratio_interval':strings(ratio(anchors[0].source_energy,anchors[1].source_energy)),
            'target_fixed_point_ratio':str(target['u']/target['d']),
            'scale_intervals_GeV_per_source_energy':[strings(a.scale_interval()) for a in anchors],
            'common_scale_interval':strings(scale),
            'one_scale_matches_both':False,
            'interpretation':'conditional dictionary excluded even before the absent renormalization/physical rest map',
            'target_uncertainty_propagated':False,
            'u_anchor_calibrated_all_bins':calibrated_native,
            'failure_reference_energy':0,
            'failure_total_source_energy_GeV':None,
            'original_probability_accounting':{f:{k:quark['full_interacting_prepared_records'][f][k]
                for k in ['complete_law','preparation_success_lower','prepared_pump_failure_interval',
                          'prepared_outside_all_bins','success_reference_energy','failure_reference_energy']}
                for f in ['uc','dc','Q']},
            'energy_subtraction':'source_excitation_relative_actual_ground is consumed unchanged; reference cost is not subtracted again',
            'complement_outcome_retained':'source receipt contains its probability bound; its energy is not a single listed bin',
            'probabilities_changed_by_unit_conversion':False},
        'target_catalogue':target_rows,
        'native_neutral_composite':{'complete_original_record_count':neutral_records['complete_joint_record_count'],
            'complete_success_energy_records':len(neutral_lines),'recorded_settings':len(neutral_records['recorded_settings']),
            'minimum_source_excitation_interval':neutral_min,
            'all_success_excitation_energies_positive':neutral_records['all_success_excitation_energies_strictly_positive'],
            'physical_neutrino_assignment':False,'GeV_calibration':None},
        'additional_physical_targets':[{'sector':'gravitational radiation and possible spin-two quanta',
            'physical_source_rest_mass_GeV':None,'source_operator_in_common_defect_model':None,
            'attachment':'gravity has a separate source and continuum premise contract; no graviton spectrum is derived here'}],
        'extra_source_modes':[
            {'sector':'compact spatial fermions','count':scalar['counting']['complete_spatial_band_rank'],'meaning':'complete isolated one-body modes, not four assigned particle species','source':'model/defect_mode_overlap/results.json'},
            {'sector':'scalar spatial profiles','count':scalar['complete_scalar_modes_count'],'meaning':'each has two doublet and two conjugate components under global occupation cap one','source':'model/defect_scalar_modes/results.json'},
            {'sector':'common CAR carrier','count':joined['joined_source']['CAR_modes'],'meaning':'30 modes per spatial site; complete Fock spectrum is not executed','source':'model/defect_mass_join/results.json'},
            {'sector':'sterile spectator','count':2,'meaning':'optional two-spin R16 spectator modes, not part of minimal R15','source':'model/defect_family_vertices/results.json'},
            {'sector':'calibrated sine-regulator valleys','count':8,'meaning':'all valleys retained in their separate quadratic source; not identified with antiparticles','source':'model/charged_lepton_response/relativistic_mass.py'}],
        'scope':{'all_target_species_accounted_for':True,'all_physical_masses_derived':False,
            'all_native_manybody_eigenvalues_computed':False,'per_species_GeV_fit':False,
            'single_source_SM_spectrum_established':False},
        'source_bindings':bindings()}


def bindings():
    paths=['model/defect_mass_catalogue/calibration.py','model/defect_mass_catalogue/catalogue.py',
           'model/defect_mass_catalogue/PROOF.md','model/lepton_cycle/physical_comparison.json',
           'model/common_source/calibrated_matter_results.json','model/native_quark_response/results.json',
           'model/defect_mass_join/results.json','model/defect_scalar_modes/results.json','model/defect_family_vertices/results.json']
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def main():
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');a=p.parse_args()
    result=build();dest=HERE/'results.json'
    if a.check:
        if result!=json.loads(dest.read_text()):raise SystemExit('mass catalogue differs')
    else:dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('17 physical targets accounted for; conditional GeV values; no complete physical mass derivation')

if __name__=='__main__':main()
