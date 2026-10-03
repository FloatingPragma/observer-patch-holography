"""Independent whole-state certification and false-promotion mutations."""
from __future__ import annotations
import builtins,copy,importlib.util,sys
from fractions import Fraction
from pathlib import Path
import pytest
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import verify_native_record_renewal as verifier

@pytest.fixture(scope='module')
def receipt():return verifier.strict_load(verifier.OUTPUT)

def test_full_native_replay_without_producer_imports(receipt,monkeypatch):
 original=builtins.__import__
 forbidden={'native_record_renewal','native_repair_flavor_constraints','record_counting_mechanism_certificate',
            'source_repair_generator_certificate','conditional_quark_mass_replay','paper_math'}
 def guarded(name,*args,**kwargs):
  if name.split('.')[-1] in forbidden:raise AssertionError('forbidden producer or target import: '+name)
  return original(name,*args,**kwargs)
 monkeypatch.setattr(builtins,'__import__',guarded)
 verifier.independent_data.cache_clear();verifier.independently_verified_payload.cache_clear()
 result=verifier.verify_payload(receipt)
 assert result['independent_states']==6077 and result['independent_A5_orbits']==136
 assert result['period']==2 and result['detailed_balance'] is False
 assert result['mean_accepted_repairs']=='11/24'
 assert result['source_selected_physical_vacuum'] is False
 assert len(receipt['stationary_Q_spectral_weights'])==31
 assert Fraction(receipt['fourth_moment']['isotropic_Gaussian_radial_defect'])!=0
 assert Fraction(receipt['stationary_current_witness']['net_current'])!=0

def test_exact_producer_rebuild(receipt):
 spec=importlib.util.spec_from_file_location('renewal_producer_test',HERE/'native_record_renewal.py')
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 assert verifier.serialized(module.build_payload())==verifier.serialized(receipt)

@pytest.mark.parametrize('path,value',[
 (('scope','source_selected_vacuum'),True),
 (('scope','conditional_driven_relative_state_selected'),False),
 (('scope','stationary_total_distribution'),True),
 (('scope','aperiodic'),True),
 (('scope','physical_clock_derived'),True),
 (('scope','quark_Higgs_response_derived'),True),
 (('scope','quark_masses_emitted'),True),
 (('scope','quark_targets_loaded'),True),
 (('scope','quark_fit_performed'),True),
 (('scope','physical_mass_scheme_attached'),True),
 (('scope','thermal_equilibrium_derived'),True),
 (('scope','Gaussian_stationary_law'),True),
 (('scope','stationary_law_is_empirical_prior'),True),
 (('specification','state','P_used'),True),
 (('specification','state','projection_boundary'),'Different protected totals are the same physical vacuum'),
 (('specification','state','mass_normalization_boundary'),'The quotient selects absolute total-dependent Higgs normalization'),
 (('specification','preparation','temporal'),'Arbitrarily correlated drives give the same law'),
 (('specification','preparation','probability'),'1/12'),
 (('specification','repair','scheduler'),'Uniform choice includes inadmissible signed updates'),
 (('specification','repair','physical_clock'),'One step is one physical Planck time'),
 (('specification','stationarity','convergence'),'Every initial law converges stepwise'),
 (('specification','readout','averaging_boundary'),'Averaged and sampled spectra are identical'),
 (('specification','readout','non_Gaussian_test'),'This record law certifies all Gaussian quark response histories'),
 (('states',),6076),
 (('A5_orbits',),135),
 (('period',),1),
 (('irreducible',),1),
 (('max_accepted_moves_per_preparation',),2),
 (('complete_height_sha256',),'0'*64),
 (('full_kernel_nonzero_entries',),1),
 (('quotient_nonzero_entries',),999),
 (('detailed_balance',),True),
 (('stationary_current_witness','net_current'),'0'),
 (('stationary_current_witness','reverse_probability'),'0'),
 (('stationary_mean_accepted_moves_per_preparation',),'1/2'),
 (('stationary_mean_Q',),'nonzero family selector'),
 (('mean_trace_Q3',),'1/1000000000000000000000000000000000000000000'),
 (('stationary_residue_weights','0'),'1/6'),
 (('stationary_port_covariance','P5'),'1/5'),
 (('stationary_quadrupole_covariance_scalar',),'1/5'),
 (('fourth_moment','isotropic_Gaussian_radial_defect'),'0'),
 (('fourth_moment','mean_squared_trace_Q2'),'0'),
 (('orbits',0,'per_state_weight'),'1'),
 (('orbits',0,'orbit_size'),60),
 (('orbits',0,'mean_accepted_moves'),'3'),
 (('numeric','mean_trace_Q2'),1.5),
])
def test_scientific_mutations_rejected(receipt,path,value):
 mutant=copy.deepcopy(receipt);cursor=mutant
 for key in path[:-1]:cursor=cursor[key]
 cursor[path[-1]]=value
 with pytest.raises(ValueError):verifier.verify_payload(mutant)

def test_normalized_stationary_weights_cannot_be_reselected(receipt):
 mutant=copy.deepcopy(receipt)
 # Preserve positivity and total normalization while disturbing the balance.
 a,b=mutant['orbits'][:2];pa=Fraction(a['stationary_weight']);pb=Fraction(b['stationary_weight']);delta=min(pa,pb)/1000
 a['stationary_weight']=str(pa+delta);b['stationary_weight']=str(pb-delta)
 with pytest.raises(ValueError,match='stationarity'):verifier.verify_payload(mutant)

def test_transition_corruption_rejected(receipt):
 mutant=copy.deepcopy(receipt);key=next(iter(mutant['orbit_kernel'][0]));mutant['orbit_kernel'][0][key]='0'
 with pytest.raises(ValueError):verifier.verify_payload(mutant)

def test_extra_mass_claim_rejected(receipt):
 mutant=copy.deepcopy(receipt);mutant['physical_quark_masses']=[1,2,3,4,5,6]
 with pytest.raises(ValueError):verifier.verify_payload(mutant)

def test_cached_verification_checks_live_pins(receipt,monkeypatch):
 original=Path.read_bytes;target=verifier.ROOT/verifier.PINS[0]
 def changed(path):return original(path)+(b' ' if path==target else b'')
 monkeypatch.setattr(Path,'read_bytes',changed)
 with pytest.raises(ValueError,match='source pin'):verifier.verify_payload(receipt)

@pytest.mark.parametrize('text',['{"x":0,"x":1}','{"x":NaN}','{"x":Infinity}'])
def test_invalid_json_rejected(tmp_path,text):
 p=tmp_path/'receipt.json';p.write_text(text)
 with pytest.raises(ValueError):verifier.strict_load(p)
