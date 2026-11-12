"""Independent native Gauss columns, changed source action and Green controls."""
from collections import defaultdict
import random
import sympy as s
import numpy as np
import pytest
from . import source,result
from model.coherent_defect_wave.chiral_joint.peter_weyl_links import model as pw
from model.coherent_defect_wave.chiral_joint.periodic_family.model import Configuration,kinetic_terms


def test_actual_native_full_factor_current_and_Gauss_projection():
    d,g=source.cochains(2);rng=random.Random(327)
    for _ in range(16):
        key=(rng.getrandbits(240),tuple(rng.randrange(119) for _ in range(24)),
             tuple(rng.randrange(5) for _ in range(8)))
        reference=rng.getrandbits(30)
        rho=source.native_density(key,reference)
        bare=source.native_density(key,0)
        b=source.native_flux(key[1])
        for x in range(8):
            native=pw.gauss_column(key,x,11)
            assert set(native)<={key}
            ref=(rho[0]-bare[0]) if x==0 else 0
            assert native.get(key,0)+ref==(rho-d.T*b)[x]
    independent=np.linalg.pinv(np.asarray(d.T*d,dtype=float))
    assert np.max(abs(independent-np.asarray(g,dtype=float)))<1e-15
    with pytest.raises(ValueError,match='zero total'):source.gauss_minimum([1]*8)


def physical_basis_input():
    path=source.physical_path(2,[(0,0,0),(1,0,0)])
    # Both weak components occupied at the arrival cell form the actual
    # colour/weak singlet of Y=-1; reference ec has Y=+1.
    y=4;mask=(1<<(30*y+10))|(1<<(30*y+11))
    return (mask,tuple(path['full_native_link_labels']),(0,)*8),1<<9


def joint_gauss(state,cell,current):
    out=defaultdict(lambda:s.Integer(0))
    for key,amplitude in state.items():
        for target,coefficient in pw.gauss_column(key,cell,current).items():out[target]+=amplitude*coefficient
        if cell==0 and current==11:out[key]+=amplitude
    return {k:s.simplify(v) for k,v in out.items() if s.simplify(v)!=0}


def test_actual_changed_full_source_column_preserves_every_Gauss_current():
    key,reference=physical_basis_input()
    assert source.native_density(key,reference)==s.Matrix([1,0,0,0,-1,0,0,0])
    assert source.coulomb_diagonal(key,reference)==s.Rational(7,48)
    value=source.changed_column(key,reference)
    assert value
    for cell in range(8):
        for current in range(12):
            assert not joint_gauss({key:1},cell,current)
            assert not joint_gauss(value,cell,current)
    two=source.changed_column(key,reference,kappa=s.Integer(2))
    difference={k:s.expand(two.get(k,0)-value.get(k,0)) for k in set(two)|set(value)}
    assert {k:v for k,v in difference.items() if v}=={key:s.Rational(7,48)}
    # The exact changed-current commutator equals the original one because
    # every source difference is diagonal in total local Y, including scalars.
    old=pw.column(key,Configuration())
    hopping=pw.hopping_column(key,kinetic_terms(Configuration()))
    lam=s.Rational(1,2**34)
    for k,v in hopping.items():old[k]=old.get(k,0)+(lam-1)*v
    old[key]=old.get(key,0)+4
    diff={k:s.expand(value.get(k,0)-old.get(k,0)) for k in set(value)|set(old)}
    assert {k:v for k,v in diff.items() if v}=={key:s.Rational(7,48)-1}


def test_path_dependent_native_flux_and_path_independent_Green_energy():
    projection=result.current_projection()
    short,long=projection['same_endpoint_different_actual_path_control']
    assert short['length']==1 and long['length']==3
    assert short['unit_kappa_Coulomb_energy']==long['unit_kappa_Coulomb_energy']
    assert short['native_flux_squared_norm']!=long['native_flux_squared_norm']
    assert short['native_oriented_flux']!=long['native_oriented_flux']
    assert short['minimum_Gauss_field']==long['minimum_Gauss_field']
    assert s.Rational(short['transverse_native_flux_squared_norm'])>0


def test_positive_additive_term_retains_all_unscreened_binding_lower_bounds():
    data=result.old_source_obstruction()
    for row in data['all_N2_constrained_rows']:
        old=list(map(s.Rational,row['actual_unchanged_source_binding_interval']))
        changed=list(map(s.Rational,row['additive_extension_binding_interval_at_kappa_one']))
        assert changed[0]>=old[0] and changed[1]>=old[1]
        if row['family']=='L':assert old==[s.Integer(row['distance'])]*2
    assert not data['unrestricted_screened_cloud_confinement_or_all_source_Coulomb_exclusion']


def test_complete_existing_readout_and_independent_coupling():
    data=result.complete_existing_readout()
    rows=data['records'];assert len(rows)==10
    arrivals=[r for r in rows if 'unit_kappa_quadratic_Gauss_flux_observable' in r]
    assert len(arrivals)==6 and all(r['unit_kappa_quadratic_Gauss_flux_observable']=='7/48' for r in arrivals)
    assert all(s.Rational(r['probability_interval'][0])>0 for r in arrivals)
    assert all('quadratic_Gauss_flux_observable' in r for r in rows if r not in arrivals)
    key,ref=physical_basis_input()
    assert source.coulomb_diagonal(key,ref,kappa=2)==2*source.coulomb_diagonal(key,ref,kappa=1)
    with pytest.raises(ValueError):source.coulomb_diagonal(key,ref,kappa=-1)


def test_missing_source_and_wrong_pin_rejected(monkeypatch):
    from . import certificate as custody
    monkeypatch.setattr(custody,'DEPENDENCIES',(*custody.DEPENDENCIES,'missing-native-Coulomb.json'))
    with pytest.raises(FileNotFoundError):custody.input_hashes()
    with pytest.raises(ValueError):custody.validate_pins({'input_sha256':{'model/openwave_current_coulomb/source.py':'0'*64}})
