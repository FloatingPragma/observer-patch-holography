from fractions import Fraction as F
import numpy as np
from .source import parameters,sector,INPUTS
from .calculation import certificate

def test_positive_parameters_and_exact_same_source_mass_gap():
 p=parameters();assert p['mbar']>p['delta']/2>0 and p['A']>0
 assert F(INPUTS['electron_unit_GeV'])*(p['delta']-1)==F('0.00078233355931')

def test_full_colour_complements_and_native_intertwiner():
 for kind,full,small in (('baryon',2024,120),('meson',576,64)):
  d=sector(kind);assert d['H'].shape==(full,full)and d['h'].shape==(small,small)
  assert np.linalg.norm(d['H']@d['B']-d['B']@d['h'])<2e-8

def test_explicit_fits_withheld_failures_and_complete_records():
 d=certificate();assert len(d['fits'])==3 and all(abs(x['residual_GeV'])<3e-11 for x in d['fits'])
 assert all(not x['physical_assignment_supported']for x in d['held_out_comparisons'])
 for row in d['complete_particle_records']:
  assert abs(sum(row['complete_success_native_occupation_probabilities'])+sum(row['complete_failure_native_occupation_probabilities'])-1)<2e-12
