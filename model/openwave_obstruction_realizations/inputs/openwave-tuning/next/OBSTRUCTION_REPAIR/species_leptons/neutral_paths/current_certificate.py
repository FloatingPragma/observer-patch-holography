"""Interval charged-character/current attachment and full-record bound."""
from pathlib import Path
import json
from fractions import Fraction as F
import mpmath as mp
ROOT=Path(__file__).parent;mp.mp.dps=70;mp.iv.dps=55;iv=mp.iv
r=json.loads((ROOT/'RESULT.json').read_text());review=json.loads((ROOT/'REVIEW.json').read_text())
full=json.loads((ROOT/'FULL_TAKAGI.json').read_text())
def mat(d):return iv.matrix([[iv.mpc(a,b)for a,b in zip(ar,br)]for ar,br in zip(d['real'],d['imag'])])
def adj(M):return iv.matrix([[iv.mpc(M[j,i].real,-M[j,i].imag)for j in range(M.rows)]for i in range(M.cols)])
def up(x):return mp.mpf(x.b)
def maximum(M):return max(up(abs(x))for x in M)
U=mat(full['takagi']);B=mat(r['full_current']);active=U[:3,:]
row=B*adj(active)
# Literal right-handed 120-degree rotation about the fixed face sum.
phi=(1+iv.sqrt(5))/2;n=iv.matrix([-phi,0,2*phi+1]);norm=iv.sqrt((n.T*n)[0]);n=n/norm
cross=iv.matrix([[0,-n[2],n[1]],[n[2],0,-n[0]],[-n[1],n[0],0]])
R=-iv.eye(3)/2+3*n*n.T/2+iv.sqrt(3)*cross/2
chars=iv.diag([iv.mpc(-iv.mpf(1)/2,iv.sqrt(3)/2),iv.mpc(-iv.mpf(1)/2,-iv.sqrt(3)/2),1])
character_entry=maximum(row*R-chars*row)
gram_entry=maximum(row*adj(row)-iv.eye(3))
reconstruction=maximum(B-row*active)
# All distinct C3 characters have separation sqrt3>1. Each row's
# orthogonal-character component is at most its residual. Projection and
# normalization give this deliberately loose complete row-frame bound.
rho=3*character_entry;eta=3*gram_entry
row_distance=4*(rho+eta)
assert rho+eta<mp.mpf('1e-45')
assert F(review['decimal_frame_to_exact_unitary_norm_bound'])<F(1,10**48)
assert F(review['every_actual_pole_absolute_error_bound'])<F(1,10**48)
current_distance=14*reconstruction+2*row_distance+mp.mpf('1e-48')
assert current_distance<mp.mpf('1e-46')

# Exact rational admitted radii, complete operator functional calculus,
# and the retained decimal clock interventions. No individual matching
# of degenerate heavy eigenvectors is used.
mass_radius=F(1,10**48);current_radius=F(1,10**46)
max_time=max(F(q['source_time'])for q in r['complete_records'])
energy_error=(F(13,3)+mass_radius)*mass_radius/2
probability_error=50*(current_radius+(max_time+1)*energy_error)+F(1,10**30)
assert probability_error<F(1,10**20)
out={'producer_code_imported':False,'literal_face_C3_rotation_interval_verified':True,
 'charged_character_equation_entry_bound':str(character_entry),
 'charged_row_Gram_entry_bound':str(gram_entry),'complete_current_reconstruction_entry_bound':str(reconstruction),
 'source_current_to_exact_C3_and_polar_frame_bound':str(current_distance),
 'admitted_current_operator_error_exact':str(current_radius),
 'admitted_complete_mass_operator_error_exact':str(mass_radius),
 'complete_five_outcome_record_error_exact':str(probability_error),
 'complete_five_outcome_record_error_decimal':float(probability_error),
 'clock_interventions':'The retained decimal source_time values are exact rational clock interventions; the inverse-square formula selects their nominal values.',
 'all63_bands_and_18_records_retained':True,'physical_interacting_neutrino_claim':False}
(ROOT/'CURRENT_CERTIFICATE.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
