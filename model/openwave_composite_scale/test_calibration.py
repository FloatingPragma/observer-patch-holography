from decimal import Decimal
from .calibration import certificate

def test_single_scale_and_rejected_physical_assignments():
 r=certificate(); assert r['physical_particle_rest_mass_claim'] is False
 scale=Decimal(r['calibration']['energy_unit_GeV'])
 for x in r['energy_rows']:
  assert abs(Decimal(x['conditional_GeV'])-scale*Decimal(x['source_energy']))<Decimal('1e-25')
 assert all(Decimal(x['ratio_to_comparison'])<Decimal('.5') for x in r['physical_assignment_tests'])
 assert r['calibration']['parameters_fitted_to_hadrons'] is False

def test_family_binding_cancellation():
 r=certificate();rows={x['channel']:Decimal(x['conditional_GeV']) for x in r['energy_rows']}
 delta=rows['light/second-family spin-zero meson']-rows['charged spin-zero meson']
 lepton=rows['second charged lepton']-rows['light charged lepton']
 assert abs(delta-lepton)<Decimal('1e-25')
 assert rows['charge-zero spin-half baryon']==rows['charge-one spin-half baryon']
