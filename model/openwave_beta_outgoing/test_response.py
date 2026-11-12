from fractions import Fraction as F
from .response import source,pair_response,wick_column_norm,parity_bound
from model.openwave_dressed_baryon_lines.source import construction,parity_basis
import numpy as np

def test_full_pair_column_matches_nonlocal_Gaussian():
 r=pair_response(8,F(1,10000));assert len(r['all_first_pair_channels'])==24
 assert abs(r['full_native_first_column_squared_norm_floating']-wick_column_norm(8))<1e-9

def test_same_native_current_and_strict_parity():
 d=source();assert 1.2752<d['ga']<1.2754
 for L in(8,16,32,64):assert F(parity_bound(L,F(1,10000))['strict_negative_exact_parity_difference_upper'])<0
 r=pair_response(8,F(1,10000));assert r['selected_actual_proton_light_electron_all_neutrals']['odd_difference']<0

def test_scalar_control_keeps_exact_reflection():
 d=construction('axial',0);B=parity_basis(d);assert np.max(abs((d['H']@B-B@(B.T@d['H']@B)).data),initial=0)<3e-10
