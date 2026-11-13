from fractions import Fraction as F
import numpy as np
from . import source as s
from .certificate import certificate

def test_positive_native_gradient_has_the_declared_nontrivial_action():
    d=s.components('baryon');G=d['hparts']['G4']
    assert d['positive_gradient_identity_exact']
    assert np.linalg.norm(G-np.trace(G)/len(G)*np.eye(len(G)))>1
    p=s.parameters();assert p['spin_gradient']==F(5,2)and p['kappa_colour']==4
    assert p['mbar']>p['delta']/2>0

def test_complete_spectra_exclude_both_odd_channels_and_keep_the_extra_sector():
    r=certificate();n=r['complete_nuclear_spectrum'];assert n['complete_singlet_dimension']==2520
    assert sum(x['fixed_Ms0_Mi0_dimension']for x in n['channels'])==240
    rows={(x['spin'],x['isospin'],x['cell_parity']):x for x in n['channels']}
    assert rows[1,0,1]['strictly_bound_in_this_native_channel']
    assert rows[0,0,-1]['above_this_channel_breakup_threshold']
    assert rows[1,1,-1]['above_this_channel_breakup_threshold']
    assert any(x['spin']==3 and x['isospin']==0 for x in r['unassigned_bound_native_channels'])
    assert all(F(x['margin_above_trial'])>0 for x in r['all_fragment_comparisons']['global_colour_controls'][-1]['all_non_3plus3_partitions'])

def test_calibrated_current_and_complete_failed_architecture_controls():
    r=certificate();a,b=map(F,r['native_axial_ground']['gA_interval'])
    assert F('1.2752')<a<b<F('1.2754')
    controls={x['gradient']:x for x in r['retrospective_gradient_controls']}
    assert controls['0']['odd_binding_to_pn_GeV']>0
    assert controls['4']['odd_S1I1_binding_to_pn_GeV']>0
    assert controls['5/2']['odd_binding_to_pn_GeV']<0 and controls['5/2']['odd_S1I1_binding_to_pn_GeV']<0
    rec=r['selected_nuclear_record']['complete_native_instrument']
    assert F(rec['rank_one_filter_success'])+F(rec['filter_failure'])==1
