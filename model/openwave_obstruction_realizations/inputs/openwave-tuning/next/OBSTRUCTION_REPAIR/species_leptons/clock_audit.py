"""Exact audit of the two lepton analyzer normalizations and shared clock."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import mpmath as mp

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
shared=ROOT/'openwave-tuning/next/CHAIN3/species_reuse_connection'
oldpath=ROOT/'reverse-engineering-reality/model/openwave_single_valley_leptons/results.json'
inputs=[shared/'REGISTER.json',shared/'RESULT.json',shared/'PROOF.md',oldpath]
pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
r=json.loads((shared/'REGISTER.json').read_text());c=json.loads((shared/'RESULT.json').read_text());old=json.loads(oldpath.read_text())
mass=list(map(F,r['mass_midpoints_exact']));electron=mass[1]
anchor=F('0.00051099895069');unit=anchor/electron
assert unit==F(c['source_energy_unit_GeV_exact'])
assert r['selected_rounds']==2**20 and min(mass)>0
mp.mp.dps=60
theta=mp.mpf(2)/9
raw=[(1+mp.sqrt(2)*mp.cos(theta+2*mp.pi*k/3))**2 for k in range(3)]
register_limit=[x/27 for x in raw];old_limit=[x*5/108 for x in raw]
assert max(abs(a/b-mp.mpf(4)/5) for a,b in zip(register_limit,old_limit))<mp.mpf('1e-55')
old_e=F(sum(map(F,old['charged_lepton_GeV'][0]['source_effect_interval']))/2)
proton=F('0.93827208943')
spin_floor=F(3,16)*unit
assert spin_floor>proton
out=dict(schema='oph.lepton_common_clock_audit.v1',input_sha256=pins,
    register_432_midpoints_exact=list(map(str,mass)),register_432_midpoints=list(map(float,mass)),
    register_432_rounds=r['selected_rounds'],register_432_source_unit_GeV_exact=str(unit),
    register_432_source_unit_GeV=float(unit),original_3600_electron_coefficient_exact=str(old_e),
    original_3600_source_unit_GeV=float(anchor/old_e),
    exact_asymptotic_register_to_old_effect_factor='4/5',
    register_432_limit_coefficients=[mp.nstr(x,40) for x in register_limit],
    mass_clock_phase_ratio_exact=list(map(str,[16*m for m in mass])),
    phase_convention=dict(mass_interference='2 h m',clock_interference='h/8',
        ratio='16m',interpretation='Compiled mass-slot coefficient, not an interacting rest-energy pole or physical elapsed-clock derivation.'),
    raw_colour_spin_floor_GeV_exact=str(spin_floor),raw_colour_spin_floor_GeV=float(spin_floor),
    proton_comparison_source_units_exact=str(proton/unit),
    raw_colour_typed_scope='Conditional unit comparison only: same raw H clock, J=1/4, nonnegative number/kinetic/colour terms, three quarks. Sum local S_x² >=3/4, so energy>=3/16 even if all bare quark masses vanish. No electromagnetic/colour native transfer is claimed.',
    light_electron_coefficient_alone_repairs_proton=False)
(HERE/'COMMON_CLOCK.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ('register_432_midpoints','register_432_source_unit_GeV','original_3600_source_unit_GeV','raw_colour_spin_floor_GeV','light_electron_coefficient_alone_repairs_proton')},indent=2))
