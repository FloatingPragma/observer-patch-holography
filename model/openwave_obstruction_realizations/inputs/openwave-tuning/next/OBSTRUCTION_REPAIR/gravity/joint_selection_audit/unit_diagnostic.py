"""Conditional common-clock gravity scale comparison; no source transfer."""
from pathlib import Path
import json,hashlib,math

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
clock_path=ROOT/'openwave-tuning/next/OBSTRUCTION_REPAIR/species_leptons/COMMON_CLOCK.json'
gravity_path=ROOT/'reverse-engineering-reality/model/openwave_gravity_unit_constraint/results.json'
curve_path=HERE.parent/'curved_ball/PLANAR_NORMALIZATION.json'
quartic_path=HERE.parent/'CANDIDATE_CHECKS.json'
clock=json.loads(clock_path.read_text());references=json.loads(gravity_path.read_text())['external_reference_inputs']
U=clock['register_432_source_unit_GeV']
hc=float(references['hbar_c_GeV_fm']);G=float(references['Newton_G_natural_GeV_minus2'])
eta1=json.loads(curve_path.read_text())['entropy_per_coordinate_face']
eta4=json.loads(quartic_path.read_text())['area_quadrature_refinement']['entropy_per_face']
c4=16*eta4
rows=[]
for R,eta in [(1,eta1),(4,eta4)]:
    predicted=1/(4*eta*U*U)
    rows.append(dict(regulator_in_coordinate_cells=R,eta_per_coordinate_area=eta,
                     conditional_coordinate_length_fm=hc/U,conditional_proper_regulator_length_fm=R*hc/U,
                     conditional_G_GeV_minus2=predicted,ratio_to_repository_G_reference=predicted/G))
out=dict(scope='Counterfactual unit attachment only; no clock, Hamiltonian, or native-source transfer is established',
         common_source_transfer_proved=False,register_energy_unit_GeV=U,
         references=references,conditional_members=rows,
         measured_area_constant_c_at_R4=c4,
         hierarchy_required_if_c4_used=dict(
             Rphys_over_common_coordinate_length=math.sqrt(4*c4*G*U*U),
             inverse_Rphys_GeV=1/math.sqrt(4*c4*G),
             meaning='Required by comparison, NOT a derived source value and NOT selected for the candidate'),
         length_only_repair_allowed=False,
         reason='With c_source=1 and common hbar, energy U fixes length hbar*c/U. Changing length alone changes the limiting speed.',
         inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [clock_path,gravity_path,curve_path,quartic_path]})
(HERE/'UNIT_DIAGNOSTIC.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
