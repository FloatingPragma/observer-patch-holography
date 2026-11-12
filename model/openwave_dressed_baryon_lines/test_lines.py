import numpy as np
from fractions import Fraction as F
from .source import construction,parity_basis,exact_parity_matrix

def test_full_old_source_transfer_and_proton_complement():
 from model.openwave_baryon_beta.decay import construction as old
 d=construction('mass_only',0);a=old();assert np.array_equal(d['keep'],a['keep']);assert np.max(abs((d['H']-a['H']).toarray()))<3e-11
 assert(d['charge_dimension'],d['spin_dimension'])==(1400,398)
 p=construction('mass_only',1);assert(p['charge_dimension'],p['spin_dimension'])==(360,102)
 assert len(set(d['ne'][d['keep']%70]))==3

def test_exact_reflection_and_current_shifts():
 for charge in(0,1):
  d=construction('axial',charge);B=parity_basis(d);h=B.T@d['H']@B
  assert np.max(abs((d['H']@B-B@h).data),initial=0)<3e-10
  assert abs(np.vdot(d['initial'],d['V']@d['initial']).real-(2*charge-1))<1e-11
  assert B.shape[1]==(199 if charge==0 else 51)

def test_exact_source_is_the_full_contact():
 d=construction('axial',1);H,S=exact_parity_matrix(d);B=parity_basis(d)
 assert H==H.T
 assert np.max(abs(np.asarray(H.evalf(17),float)-(B.T@d['H']@B).toarray().real))<3e-10
