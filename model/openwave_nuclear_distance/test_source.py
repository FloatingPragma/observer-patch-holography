from .nuclear import carrier,channel,fitted_eta
from .spatial import carrier as spatial_carrier,graph
from model.openwave_baryon_weak_formfactors.native import formfactors
from model.openwave_axial_state.source import matrix

def test_exact_native_parity_channels():
 d=carrier();assert d['colour_weight_dimension']==240
 assert sum(channel(1,0,p)['rank']for p in(-1,1))==27
 assert channel(0,0,-1)['rank']>0

def test_state_preserving_eta_shortcut():
 p=fitted_eta(4.);assert abs(formfactors(matrix('baryon',p))['gA']-1.2753)<2e-8

def test_actual_graph_native_colour_isometry():
 d=spatial_carrier(3);assert len(d['labels'])==364 and d['colour_error']<2e-14
