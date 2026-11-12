from fractions import Fraction as F
from .levels import local_census,certificate

def test_exact_native_local_character_and_axial_multiplicity():
 assert local_census(3)['exact_colour_singlet_dimension']==20
 assert local_census(6)['exact_colour_singlet_dimension']==50
 assert local_census(3)['complete_spin_isospin_multiplets'][0]['multiplicity']==1

def test_complete_low_bands_and_retained_unwanted_state():
 d=certificate();assert len(d['same_source_graphs'])==5
 for row in d['same_source_graphs']:
  L=row['vertices'];assert row['nuclear_total_low_dimension_including_magnetic_components']==50*L+200*L*(L-1)
  assert all(x['strictly_bound']for x in row['binding_channels']if x['isospin']==0)
  for sector in row['six_quark_complete_low_bands']:
   law=sector['complete_original_preparation_and_spectral_records']
   assert F(law['preparation_success'])+F(law['preparation_failure'])==1
   assert len(sector['complete_low_energy_intervals_relative_to_6mbar_minus_delta_Iz'])==sector['per_magnetic_component_low_dimension']
