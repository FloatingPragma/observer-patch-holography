from fractions import Fraction as F
from .source import native_current,packet_inner,translate,cosine_mean,KAPPA

def test_every_translation_return_path_enters_packet_moments():
    state=((0,0,0,0,0),)
    vector={state:F(1),translate(state,1):F(-1)}
    interval,words=packet_inner(vector,vector)
    c=cosine_mean(KAPPA)
    assert words=={0:F(2),-1:F(-1),1:F(-1)}
    assert interval==(2-2*c[1],2-2*c[0])
    assert interval[0]>0

def test_full_native_current_energy_and_complete_outputs():
    d=native_current()
    assert len(d['u'])==522 and len(d['h'])==19140
    assert F('2.418')<d['mean'][0]*d['unit']<F('2.419')
    assert F('0.085')<d['variance'][0]*d['unit']**2<F('0.087')
    assert 0<d['success']<1

def test_resolvent_to_window_inequality_and_density_countercontrol():
    # Exact two-point spectral measure: positive finite-eta response at E=0
    # coexists with zero weight in every window of half-width below one.
    energies=[F(-2),F(2)];weights=[F(1,2)]*2;eta=F(1,8);width=F(1)
    effect=sum(p*eta**2/(e**2+eta**2)for e,p in zip(energies,weights))
    outside=eta**2/(width**2+eta**2)
    assert effect>0 and max(F(0),(effect-outside)/(1-outside))==0
    assert all(abs(e)>width for e in energies)
