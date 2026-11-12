from fractions import Fraction as F
import numpy as np
from .window import compare_native,certified_universal,universal,taylor_first,first_light
from .certificate import finite_control,family_bound

def test_complete_native_factorization():
 for L in(2,4):assert compare_native(L)['all_nine_blocks_retained']

def test_all_channel_feedback_bound():
 for L in(2,4):assert finite_control(L,F(1,10))['all_channels_in_comparison']

def test_exact_surrogate_and_family():
 u=certified_universal(4);v,c=taylor_first(u,F(1,10));w=first_light(universal(4),F(1,10))
 assert np.linalg.norm(v-w)<float(F(c['total_numeric_vector_error_upper_exact']))
 assert family_bound()['absolute_C01_lower']=='19/100'
