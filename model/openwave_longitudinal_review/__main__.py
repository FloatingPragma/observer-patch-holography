"""Independently reconstruct Gauss labels, rational moments and source floors."""
import argparse,hashlib,json
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from model.openwave_common_electromagnetism import source as owner
from model.openwave_common_mass_source.source import exact_baryon_window,G,INPUTS
from model.openwave_axial_state.source import sector,parameters,weights

ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent

def pins():
    from model.openwave_common_electromagnetism.__main__ import pins as parent_pins
    target=ROOT/'model/openwave_common_electromagnetism/results.json';saved=json.loads(target.read_text());out=parent_pins()
    assert saved['source_sha256']==out
    for f in(target,HERE/'__init__.py',HERE/'__main__.py'):
        out[str(f.relative_to(ROOT))]=hashlib.sha256(f.read_bytes()).hexdigest()
    return dict(sorted(out.items()))

def calculate():
    target=ROOT/'model/openwave_common_electromagnetism/results.json';saved=json.loads(target.read_text());bindings=pins()
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==v for p,v in bindings.items())
    cases=0;negative=0
    for V in(2,3,5):
        for numerators in product(range(-2,3),repeat=V):
            Q=list(map(lambda n:F(n,3),numerators));R=-sum(numerators);E=[V*n+R for n in numerators]
            assert sum(E)==0 and all(E[x]-3*V*Q[x]-R==0 for x in range(V))
            assert owner.gauss_label(Q)==(tuple(E),R)
            if any(E[x]-3*Q[x]-R for x in range(V)):negative+=1
            cases+=1
    assert negative>0
    assert owner.exact_ward()==saved['native_Ward_certificate']
    prior=exact_baryon_window();d=sector('baryon');p=parameters();unit=F(INPUTS['electron_unit_GeV']);delta=F(1,64);comparisons=[]
    for name in('neutron','proton'):
        v=list(map(int,prior[name+'_integer_trial']));N=sum(a*a for a in v);q=[];dq=[]
        for mask in d['basis']:
            charge=[F(0),F(0)];selfenergy=F(0)
            for j,(kind,x,spin,weak,colour)in enumerate(d['modes']):
                if mask>>j&1:
                    z=F(2 if weak==0 else -1,3)*(1 if kind==0 else -1);charge[x]+=z;selfenergy+=z*z
            q.append(charge);dq.append(delta*((charge[0]-charge[1])**2-selfenergy))
        h=[F(0)for _ in v]
        for part,w in weights(p).items():
            A=d['parts'][part].tocsr()
            for i in range(len(v)):
                h[i]+=w*sum(int(A.data[k])*v[int(A.indices[k])]for k in range(A.indptr[i],A.indptr[i+1]))
        old=sum(x*y for x,y in zip(v,h))/N
        h=[a+z*x for a,z,x in zip(h,dq,v)];mean=sum(x*y for x,y in zip(v,h))/N
        variance=sum((a-mean*x)**2 for a,x in zip(h,v))/N+G*G*F(prior['contact_second_moment_exact'])
        seed=prior['native_occupation_seeds'][name]['index'];work=F(2*v[seed],N)*(v[seed]*mean-h[seed])
        row=next(x for x in saved['changed_baryon_energy_records']['rows']if x['assigned_channel']==name)
        assert mean*unit==F(row['mean_energy_GeV_exact']) and (mean-old)*unit==F(row['electric_mean_shift_GeV_exact'])
        assert variance==F(row['full_source_variance_electron_units_squared_exact'])
        assert work==F(row['preparation_mean_work_electron_units_exact'])
        assert F(row['conditional_inside_probability_lower_exact'])==1-4*variance
        assert F(row['original_inside_probability_lower_exact'])+F(row['original_outside_probability_upper_exact'])+F(row['preparation_failure_probability_exact'])==1
        comparisons.append({'channel':name,'exact_mean_variance_work_and_complete_bounds_equal':True})
    floor=p['mbar']-p['delta']/2-18*G-4*delta/9;assert floor>1 and 2*(1-delta)>1
    assert F(saved['positivity']['colourful_floor_per_particle_exact'])==floor
    assert F(saved['positivity']['exact_charge_plus_and_minus_one_threshold_GeV'])==unit
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==v for p,v in bindings.items())
    return {'verdict':'PASS_INDEPENDENT_LONGITUDINAL_SOURCE_REVIEW','source_sha256':bindings,
      'independent_rational_Gauss_records':cases,'incorrect_flux_factor_negative_controls':negative,
      'full_192_mode_native_Ward_replay':True,'complete_native_baryon_moment_reconstructions':comparisons,
      'normal_ordered_full_source_threshold_bound_verified':True,
      'review_scope':'Exact native charge and Gauss/background algebra, rational full-source trial moments, complete bounds and positivity. The source is the four-input axial/common-Dirac member; no transfer to the five-input spin-gradient source or transverse field is inferred.'}

def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();r=calculate();out=HERE/'results.json'
    if a.write:out.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    else:assert json.loads(out.read_text())==r
    print(r['verdict'])
if __name__=='__main__':main()
