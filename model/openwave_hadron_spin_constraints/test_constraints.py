from .constraints import extension,selected_energies,extended_matrix
from .hyperfine import native_hyperfine_check
from model.openwave_baryon_weak_formfactors.native import formfactors

def test_state_preserving_positive_controls():
 for k in(0.,191.,622.):
  p=extension(k,0.);e=selected_energies(p);f=formfactors(extended_matrix('baryon',p));assert abs(f['gA']-1.2753)<2e-8
  assert abs(e['proton']-.93827208943)<2e-10 and abs(e['charged_pion']-.13957039)<2e-10

def test_native_hyperfine_not_composite_projector():
 for kind in('baryon','meson'):assert native_hyperfine_check(kind)['full_native_intertwiner_error']<2e-11
