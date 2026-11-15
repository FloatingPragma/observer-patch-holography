"""Exact coherent-source interfaces and their remaining physical joins."""
from copy import deepcopy
NAMES=('coherent_native_architecture','coherent_native_propagation','active_native_charge','native_charge_readout','native_coherent_clock','native_calibrated_energy','native_euclidean_transfer')
TITLES=('Complete coherent event architecture','Coherent native propagation','Active protected memory charge','Complete charge readout and compensated recovery','Repeated coherent source clock','Same-source interaction and calibrated energies','Finite Euclidean transfer and gap bounds')

def project(data):
 selected={n:{k:deepcopy(v)for k,v in data[n].items()if k!='source_sha256'}for n in NAMES}
 for n in NAMES:
  if not selected[n]['verdict'].startswith('PASS'):raise ValueError(('missing source acceptance',n))
 arch=selected['coherent_native_architecture'];walk=selected['coherent_native_propagation'];clock=selected['native_coherent_clock'];energy=selected['native_calibrated_energy'];charge=selected['active_native_charge'];readout=selected['native_charge_readout'];euclidean=selected['native_euclidean_transfer']
 if walk['construction']!=arch['construction'] or clock['shared_charge_clock_program']['common_hardware']!=arch['construction'] or readout['shared_construction']!=arch['construction']:raise ValueError('different hardware descriptors')
 if arch['official_OpenWave_row_award'] or walk['physical_OpenWave_row_award']:raise ValueError('finite architecture promoted to official physics')
 if clock['scope']['literal_OpenWave_row3_complete'] or clock['scope']['three_spectator_values_are_three_lepton_species']:raise ValueError('prepared clock promoted to physical oscillator or leptons')
 if energy['energy_attachment']['physical_rest_pole_or_particle_identification'] or energy['energy_attachment']['new_independent_GeV_prediction']:raise ValueError('calibrated transitions promoted to physical masses')
 if charge['scope']['physical_electric_charge_identified'] or charge['scope']['physical_quark_identified']:raise ValueError('memory charge promoted to physical charge')
 if readout['matching_clock_recovery_join']['uncompensated_clock_theorem_transferred']:raise ValueError('charged clock counterexample erased')
 for k in ('four_dimensional_field_limit','uniform_physical_continuum_gap','pure_yang_mills_identification','OpenWave_criterion_awarded'):
  if euclidean['scope'][k]:raise ValueError(('finite transfer promoted',k))
 return {'source_packages':selected,'shared_declared_hardware':deepcopy(arch['construction']),
  'named_programs':{'committed_repair':'Native original-record commitment, with exact dissipative marginal.',
   'spatial_walk':'Coherent record recoding and uniform-event reflection on a one-load federation.',
   'active_charge':'Two-load native transport plus declared controlled memory feedback and phase-sensitive readout.',
   'compensated_clock':'Inverse active feedback before original weighted recovery; two-hop coherent phase and complete quadrature readout.'},
  'exact_interfaces':{'record_commit_recovers_native_instrument':True,'walk_compression_uses_native_population_transfer':True,'compensated_charge_clock_branch_identity':readout['matching_clock_recovery_join']['compensated_identity'],'clock_and_entangling_interaction_share_compressed_generator':True,'finite_Euclidean_transfer_from_actual_population':True},
  'uncompensated_clock_counterexample':deepcopy(readout['matching_clock_recovery_join']['nonzero_charge_witness']),
  'joint_spatial_charge_clock_interacting_field_implemented':False,'all_physical_observables_reproduced':False,
  'physical_premises_closed':[],'new_independent_GeV_rest_masses':[],
  'official_OpenWave_validated_rows':0,'official_OpenWave_untested_rows':31,'independent_Nature_matches_added':0}

def render(value):
 s=value['source_packages'];energy=s['native_calibrated_energy']['energy_attachment'];rows=energy['rows']
 lines=['### Coherent construction and exact interfaces','',
 'One declared native event architecture supports named committed, coherent-walk, active-charge and compensated-clock programs. Their exact interfaces are stated below; the programs are not treated as one interchangeable evolution.','']
 summaries=[
 'The complete local coherent event reproduces the original C/E instrument when its original records are committed. Reversible record recoding changes the uncommitted observation policy and retains the inverse readout.',
 'An exact cyclic reduction propagates complete quantum port probabilities on three source sizes. The finite dispersion follows cos(omega)=1-lambda/(2M), with a zero constant mode and declared reflection. Neither a massive physical dispersion nor a laboratory speed follows.',
 'A six-state active memory extends the source algebra and supplies protected character projectors. Its transport holonomy is nontrivial. Physical electric normalization, spatial quark winding and independent protected daughter charges are not established.',
 'The actual compiled charge feedback and phase readout retain all source/memory/padding outcomes. Uncorrected charge character one damps a matching coherence by 85/96. In the exact declared service contract, inverse feedback before source recovery restores the complete original clock channel branchwise, with explicit added gates and retained failure complements. Stored primitive errors and physical implementation precision are additional requirements.',
 'All three signed gaps of the corrected 512-address source have consecutive complete quadrature records and finite-step probability bounds. These are one binary mode with three spectator values. Coherent preparation, phase and readout are declared; no self-starting physical oscillator is inferred.',
 'The same compressed generator has a certified binary–spectator entangling interaction. Equating its smallest absolute transition gap to the electron anchor gives ratios 1:3:12, which fail the charged-lepton hierarchy. The calibrated transitions have no stable physical rest-pole identification.',
 'The actual native population transfer gives a finite reflection-positive stationary path measure, a logarithmic transfer Hamiltonian and positive finite gap bounds. A nonzero occupation fourth cumulant does not supply a four-dimensional interacting field or the Clay construction.']
 for name,title,body in zip(NAMES,TITLES,summaries,strict=True):lines += [f'- [{title}](../model/openwave_{name}/PROOF.md). {body}']
 lines += ['', 'These seven component contracts have the stated separate schedules. The joint charged-field and tracked-clock implementation is specified in the [joint observer register](../claims/joint_observer.json); physical particle identification, photon/metric observations and a complete interacting continuum require further proofs. These results add zero independent Nature validations and zero official OpenWave awards. Calibrated source values retain their empirical input ancestry.','']
 return lines
