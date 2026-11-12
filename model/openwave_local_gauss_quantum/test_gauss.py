from .gauss import exact_charge_tree,finite_CAR_gauss_control,charged_block

def test_native_reference_and_integral_Gauss_units():
 d=exact_charge_tree();assert sorted(d['native_pair_charge_units_one_third'])==[-3,-1,-1,-1,2,2,2]
 assert d['native_reference_charge_per_orbital']=='0'

def test_complete_CAR_gauge_dressing_rejects_bare_hop():
 d=finite_CAR_gauss_control();assert d['complete_Fock_words']==16 and d['undressed_hopping_Gauss_defect']>0
 for r in d['every_complete_charge_flux_reference_record']:assert sum(r['two_integer_link_fluxes'])+r['reference_integer_charge']==0

def test_complete_charged_block_keeps_field_cost_and_reference():
 d=charged_block();assert d['star_exact_electric_energy_increment']=='1/16'
 assert d['source_plus_compensator_ground_energy_interval'][1]<16.5
 assert d['conditional_GeV_using_unchanged_electron_unit']>0.00051099895069
 assert abs(sum(x['probability']for x in d['complete_native_orbital_charge_link_flux_reference_records'])-1)<1e-12
