"""Exact common-unit compatibility, with the source stiffness fixed."""
from fractions import Fraction as F
from pathlib import Path
import json
from model.openwave_capacity_stress.source import sqrt_interval, enclosure
import model.openwave_transverse_gravity.units as parent
ROOT=Path(parent.__file__).resolve().parents[2]

def pi_interval(terms=60):
    """Machin identity with alternating rational arctangent remainders."""
    def atan(q):
        total=sum((-1)**n*F(1,(2*n+1)*q**(2*n+1)) for n in range(terms))
        nxt=(-1)**terms*F(1,(2*terms+1)*q**(2*terms+1))
        return min(total,total+nxt),max(total,total+nxt)
    a,b=atan(5),atan(239)
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]

def solve(g, u, a, hc, target):
    g,u,a,hc,target=map(F,(g,u,a,hc,target))
    if min(g,u,a,hc,target)<=0:
        raise ValueError('positive source and unit data required')
    speed=u*a/hc
    strength=g*a/(u*hc)
    invariant=g/u**2
    assert strength/speed==invariant
    light_length=hc/u
    gravity_length=target*u*hc/g
    needed_g=target*u*u
    return dict(speed=str(speed), strength=str(strength), invariant=str(invariant),
                light_matched_length_fm=str(light_length),
                gravity_matched_length_fm=str(gravity_length),
                gravity_after_light_matching=str(invariant),
                speed_after_gravity_matching=str(target/invariant),
                required_dimensionless_g=str(needed_g),
                required_stiffness_ratio=str(g/needed_g),
                electron_anchor_preserved=True,
                length_alone_can_match_both=(g==needed_g))

def run():
    receipt=json.loads((ROOT/'model/openwave_transverse_gravity/results.json').read_text())
    # Find the actual parent's unit comparison without maintaining a duplicate table.
    def find(x):
        if isinstance(x,dict):
            if 'dimensionless_Newton_response' in x:return x
            for v in x.values():
                r=find(v)
                if r is not None:return r
        if isinstance(x,list):
            for v in x:
                r=find(v)
                if r is not None:return r
    p=find(receipt); assert p is not None
    inputs=p['external_reference_inputs']
    kappa=F(parent.profile()['kappa_exact'])
    plo,phi=pi_interval()
    glo,ghi=1/(4*phi*kappa),1/(4*plo*kappa)
    g=(glo+ghi)/2
    u=F(p['energy_unit_GeV_exact']); a=F(p['length_unit_fm_exact'])
    hc=F(inputs['hbar_c_GeV_fm']); target=F(inputs['Newton_G_natural_GeV_minus2'])
    out=solve(g,u,a,hc,target)
    scans=[{'length_multiplier':str(s),**solve(g,u,a*s,hc,target)}
           for s in (F(1,1000),F(1),F(1000))]
    energy_for_both=(sqrt_interval(glo/target)[0],sqrt_interval(ghi/target)[1])
    assert glo/u**2>target
    return {'schema':'oph.common_unit_compatibility.v1',
            'verdict':'PASS_EXACT_UNIT_ELIMINATION_FIXED_SOURCE_INCOMPATIBLE',
            'source_parent':'model/openwave_transverse_gravity/results.json',
            'formula':'G/c_ratio = g/u^2; simultaneous agreement requires g=G_observed*u^2',
            'source_stiffness_exact':str(kappa),
            'dimensionless_response_interval_exact':[str(glo),str(ghi)],
            'pi_interval_exact':[str(plo),str(phi)],
            'coefficient_scope':'parent exact dyadic stiffness and Machin rational pi enclosure; midpoint diagnostic plus enclosing intervals',
            'light_matched_G_interval':enclosure((glo/u**2,ghi/u**2)),
            'required_stiffness_ratio_interval_exact':[str(glo/(target*u*u)),str(ghi/(target*u*u))],
            'fixed_source':out,'length_controls':scans,
            'energy_unit_required_if_electron_anchor_abandoned_GeV':enclosure(energy_for_both),
            'required_energy_over_electron_unit':enclosure(tuple(x/u for x in energy_for_both)),
            'diagnostic_float':{k:float(F(v)) for k,v in out.items() if isinstance(v,str)},
            'source_stiffness_calibrated_to_G':False,
            'scope':{'physical_G_and_c_match':False,'common_strong_geometry_transfer':False,
                     'hardware_uniqueness_required':False},
            'external_reference_inputs':inputs}
