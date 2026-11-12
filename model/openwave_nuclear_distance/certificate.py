from functools import lru_cache
from fractions import Fraction as F
import json
from pathlib import Path
from .nuclear import nuclear_certificate
from .spatial import spatial_certificate

@lru_cache(None)
def certificate():
 n=nuclear_certificate();s=spatial_certificate();root=Path(__file__).resolve().parents[2];a=json.loads((root/'model/openwave_axial_state/results.json').read_text());b=[]
 for q in(0,1):b.append(next(row for row in a['complete_colour_spectra']['baryon']if row['charge']==q)['certificate']['all_eigenvalue_intervals'][0])
 threshold=[sum(F(row[i])for row in b)for i in(0,1)];cert=[]
 for name,r in n['exact_ground_channel_certificates'].items():
  e=list(map(F,r['all_eigenvalue_intervals'][0]));bound=[threshold[0]-e[1],threshold[1]-e[0]];assert bound[1]<0
  cert.append({'channel':name,'binding_interval_rational':list(map(str,bound)),'strictly_above_two_baryon_copy_threshold':True})
 return {'schema':'oph.native_nuclear_parity_and_distance.v1','verdict':'PASS_COMPLETE_NUCLEAR_PARITY_CENSUS_AND_NATIVE_SPATIAL_CONTROL','nuclear':n,'certified_four_input_unbound_channels':cert,'spatial':s,
 'no_silent_source_transfer':'The eta1 four-input source, eta-fitted fifth-input nuclear control, and fixed-parameter graph controls are separate declared constructions. Their numerical successes do not combine into one physical spectrum.',
 'scope':'Exact finite source/census/positivity identities and two certified nuclear channel ground bounds; complete numerical nuclear and baryon graph spectra and original records. No physical fm calibration, continuum orbital identification or graph-uniform nuclear binding theorem.'}
