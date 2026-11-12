"""Exact local current algebra and a finite-hopping axial bound."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json
import numpy as np
from model.openwave_calibrated_hadrons.source import sector,parameters,angular4,integral
from .native import currents,formfactors
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate.action import SOURCE_FILES
from model.openwave_larmor_dynamics.certificate import equivalent

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def pins():
    seeds=set(SOURCE_FILES)
    for name in ('openwave_baryon_weak_formfactors','openwave_calibrated_hadrons'):
        seeds.update(str(p.relative_to(ROOT))for p in (ROOT/'model'/name).iterdir()if p.suffix in ('.py','.md','.tex')or p.name=='inputs.json')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))}

def exact_local():
    c=currents();d=c['native'];labels=d['old']['labels']
    seed=next(i for i,row in enumerate(labels)if all(x==0 for x,s,w in row)and sum(s==0 for x,s,w in row)==2 and sum(w==0 for x,s,w in row)==1)
    v=d['Binteger'][:,seed].toarray().ravel().astype(np.int64)
    for spin in (True,False):v=15*v-angular4(d['modes'],d['basis'],(0,1),spin)@v
    v=v+integral(c['full_parity'])@v
    assert np.any(v)
    jp=integral(sum(c['full'][x][0][0]for x in range(2)))
    jz=integral(sum(c['full'][x][3][0]for x in range(2)))
    p=jp@v;norm=int(v@v);assert int(p@p)==norm
    gv=F(int(p@(jp@v)),norm);ga=F(int(p@(jz@v)),norm)
    leak=3*(jz@v)-5*p;leak2=F(int(leak@leak),9*norm)
    assert gv==1 and ga==F(5,3)and leak2==F(8,9)
    return {'gV':str(gv),'gA':str(ga),'axial_leakage_norm_squared':str(leak2),
            'integer_neutron_vector':[str(int(x))for x in v],'integer_vector_norm_squared':str(norm)}

def calculate():
    p=parameters();f=formfactors(sector('baryon')['h']);local=exact_local()
    gamma=F(8,3)*p['U']-F(4,9)*p['A']-F(3,4)*(p['Js']+p['Jw'])-3*p['t']
    assert gamma>0
    sine=3*p['t']/gamma;error=12*sine;interval=[F(5,3)-error,F(5,3)+error]
    comparison=json.loads((HERE/'inputs.json').read_text());observed=F(comparison['axial_ratio_magnitude'])
    assert interval[0]>observed+10*F(comparison['quoted_uncertainty'])
    return {'schema':'native-baryon-weak-formfactors-v1','source_sha256':pins(),
        'source':'openwave_calibrated_hadrons; same complete two-cell native colour action',
        'exact_local_current_algebra':local,'finite_hopping_certificate':{
            'complement_gap_lower_rational':str(gamma),'state_sine_error_upper_rational':str(sine),
            'axial_ratio_interval_rational':[str(x)for x in interval],'axial_ratio_interval':[float(np.nextafter(float(interval[0]),-np.inf)),float(np.nextafter(float(interval[1]),np.inf))]},
        'floating_evaluation':{k:f[k]for k in ('gV','gA','current_leakage_operator_norm')},
        'comparison':comparison,'observed_axial_ratio_reproduced':False,
        'native_vector_charge_normalization_proved':True,
        'scope':'Exact native internal-spin current form factors at zero graph momentum. No physical weak lifetime, renormalized QCD axial current or relativistic spacetime dictionary.'}

def main():
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group();g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=parser.parse_args()
    before=pins();out=calculate();assert before==pins();path=HERE/'results.json'
    if args.check:assert equivalent(json.loads(path.read_text()),out)
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k in ('floating_evaluation','finite_hopping_certificate','observed_axial_ratio_reproduced')},indent=2))

if __name__=='__main__':main()
