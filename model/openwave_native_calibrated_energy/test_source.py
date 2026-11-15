import numpy as np
from fractions import Fraction as F
from .source import partial_transpose_binary

def test_full_syndrome_partial_transpose_detects_entanglement():
 # Binary/spectator Bell pair in the complete36-coordinate chart.
 v=np.zeros(36);v[0]=v[18]=1/np.sqrt(2);rho=np.outer(v,v)
 assert abs(sum(max(0.,-x)for x in np.linalg.eigvalsh(partial_transpose_binary(rho)))-.5)<1e-14
 assert np.array_equal(partial_transpose_binary(partial_transpose_binary(rho)),rho)

def test_calibrated_ratios_and_positive_full_source_witness():
 import json
 from pathlib import Path
 p=Path(__file__).with_name('results.json')
 if not p.exists():return
 d=json.loads(p.read_text());assert d['energy_attachment']['absolute_gap_ratios_exact']==['1','3','12']
 assert all(F(x['complete_binary_vs_spectator_syndrome_negativity_lower_exact'])>0 for x in d['interaction']['finite_entanglement_controls'])
 assert not d['energy_attachment']['literal_three_lepton_ratio_match']
 assert not d['energy_attachment']['physical_rest_pole_or_particle_identification']
