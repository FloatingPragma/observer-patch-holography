import sympy as s
from .charge import native_data,four_mode_controls,neutral_floor,block_attachment

def test_actual_native_reference_and_complete_charge_pairs():
 d=native_data();assert len(d['pairs'])==7 and all(d['reference_checks'])
 assert sum(s.Rational(x['charge'])for x in d['pairs'])==0
 assert len(d['S'])==7 and d['Q'][d['neutral'][0],d['neutral'][0]]==0

def test_full_CAR_chirality_spin_and_native_charge_identities():
 rows=four_mode_controls();assert len(rows)==6
 assert max(r['charge_identity_norm']for r in rows)<1e-14
 assert max(max(r['three_spin_identity_norms'])for r in rows)==0

def test_neutral_complement_cannot_lower_charged_band():
 d=neutral_floor();assert all(s.Rational(x)>0 for x in d['exact_LDL_diagonal'])
 r=block_attachment();assert r['dressed_native_energy']<16.5
 assert r['isolated_charge_plus_or_minus_one_band_dimension']==4

def test_omitting_quark_singlets_loses_neutral_reference():
 d=native_data();electron=next(x for x in d['pairs']if x['channel']=='charged_lepton')
 assert 2*d['Q'][electron['right'],electron['right']]==2
 assert 2*sum(d['Q'][j,j]for j in d['S'])==0
