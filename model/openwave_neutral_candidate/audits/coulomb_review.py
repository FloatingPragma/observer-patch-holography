"""Independent native Gauss, projected flux and replacement-energy controls."""
from collections import defaultdict
from itertools import product
from pathlib import Path
import hashlib,json,random
import numpy as np
import sympy as s

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=ROOT/'model/openwave_current_coulomb'


def incidence(n):
    points=list(product(range(n),repeat=3));index={p:i for i,p in enumerate(points)}
    d=s.zeros(3*n**3,n**3)
    for p,i in index.items():
        for j in range(3):
            q=list(p);q[j]=(q[j]+1)%n
            d[3*i+j,i]=-1;d[3*i+j,index[tuple(q)]]=1
    return d


def cochain_review():
    from model.openwave_current_coulomb import source
    errors=[]
    for n in (2,4):
        d=incidence(n);provided,g=source.cochains(n);v=n**3;p0=s.ones(v)/v
        assert d==provided and d.T*d*g==s.eye(v)-p0 and g*s.ones(v,1)==s.zeros(v,1)
        if n==2:assert g==(d.T*d+p0).inv()-p0
        rho=s.zeros(v,1);rho[0]=1;rho[n*n]=-1
        field=d*g*rho
        assert d.T*field==rho
        energy=(field.T*field)[0]/2
        assert energy==(rho.T*g*rho)[0]/2==s.Rational(v-1,6*v)
        # A curl around one actual four-edge plaquette is divergence-free.
        loop=s.zeros(3*v,1);loop[0]=1;loop[3*n*n+1]=1;loop[3*n]= -1;loop[1]=-1
        assert d.T*loop==s.zeros(v,1)
        assert (field.T*loop)[0]==0
        assert ((field+loop).T*(field+loop))[0]==(field.T*field)[0]+(loop.T*loop)[0]
        errors.append({'side':n,'nearest_unit_pair_energy':str(energy),'projected_dimension':v-1})
    small=source.physical_path(4,[(0,0,0),(1,0,0)])
    detour=source.physical_path(4,[(0,0,0),(0,1,0),(1,1,0),(1,0,0)])
    d=incidence(4);delta=s.Matrix(detour['native_oriented_flux'])-s.Matrix(small['native_oriented_flux'])
    assert d.T*delta==s.zeros(64,1) and delta!=s.zeros(192,1)
    assert small['unit_kappa_Coulomb_energy']==detour['unit_kappa_Coulomb_energy']
    assert detour['length']-small['length']==2
    # Mean-zero pseudoinversion cannot manufacture neutrality.
    d=incidence(2);g=source.cochains(2)[1];charged=s.ones(8,1)
    assert d.T*d*g*charged==s.zeros(8,1) and charged!=s.zeros(8,1)
    return {'exact_cochain_controls':errors,'nonzero_transverse_path_difference_retained':True,
        'unit_link_cost_difference':2,'new_quadratic_cost_difference':'0',
        'nonneutral_full_carrier_is_not_a_physical_Gauss_solution':True}


def native_operator_review():
    from model.openwave_current_coulomb import source
    from model.coherent_defect_wave.chiral_joint.peter_weyl_links import model as pw
    from model.coherent_defect_wave.chiral_joint.periodic_family.model import Configuration,kinetic_terms
    d=incidence(2);rng=random.Random(60173);checks=0
    for _ in range(12):
        key=(rng.getrandbits(240),tuple(rng.randrange(119) for _ in range(24)),tuple(rng.randrange(5) for _ in range(8)))
        ref=rng.getrandbits(30);rho=source.native_density(key,ref);bare=source.native_density(key)
        b=source.native_flux(key[1]);expected=rho-d.T*b
        for cell in range(8):
            gauss=pw.gauss_column(key,cell,11)
            reference=(rho-bare)[cell]
            assert gauss.get(key,0)+reference==expected[cell];checks+=1
    # A different axis and spin from the producer's physical column test.
    path=source.physical_path(2,[(0,0,0),(0,0,1)])
    key=((1<<(30+25))|(1<<(30+26)),tuple(path['full_native_link_labels']),(0,)*8)
    reference=1<<24;lam=s.Rational(1,2**34)
    old=defaultdict(lambda:s.Integer(0),pw.column(key,Configuration()))
    for out,value in pw.hopping_column(key,kinetic_terms(Configuration())).items():old[out]+=(lam-1)*value
    old[key]+=4
    rho=source.native_density(key,reference);g=source.cochains(2)[1]
    quadratic=(rho.T*g*rho)[0]/2
    assert quadratic==s.Rational(7,48)
    manual=dict(old);manual[key]+=quadratic-sum(i!=0 for i in key[1])
    manual={k:s.expand(v) for k,v in manual.items() if s.expand(v)!=0}
    actual=source.changed_column(key,reference)
    assert manual==actual
    for cell in range(8):
        for current in range(12):
            for state in ({key:s.Integer(1)},actual):
                output=defaultdict(lambda:s.Integer(0))
                for inp,amplitude in state.items():
                    for out,value in pw.gauss_column(inp,cell,current).items():output[out]+=amplitude*value
                    if cell==0 and current==11:output[inp]+=amplitude
                assert all(s.expand(value)==0 for value in output.values())
    for out in set(old)|set(actual):
        change=source.native_density(out,reference)-rho
        assert (actual.get(out,0)-old.get(out,0))*change==s.zeros(8,1)
    doubled=source.changed_column(key,reference,kappa=2)
    for out in set(actual)|set(doubled):
        assert s.expand(doubled.get(out,0)-actual.get(out,0))==(quadratic if out==key else 0)
    return {'independent_full_current_Gauss_columns':checks,
        'different_spin_axis_physical_column_constraints':96,
        'changed_full_H_reconstructed_from_original_column':True,
        'local_density_commutator_unchanged':True,'independent_kappa_scaling_checked':True}


def green_profile_review():
    from model.openwave_field_observables.lattice import green
    n=32;angles=2*np.pi*np.arange(n)/n
    lam=2*(3-np.cos(angles)[:,None,None]-np.cos(angles)[None,:,None]-np.cos(angles)[None,None,:])
    inverse=np.zeros_like(lam);np.divide(1.,lam,out=inverse,where=lam>0)
    def value(r):return float(np.sum(np.cos(r*angles)[:,None,None]*inverse)/n**3)
    direct=[value(r) for r in (0,3,4)]
    fft=green(n);error=max(abs(x-fft[r,0,0]) for x,r in zip(direct,(0,3,4)))
    assert error<2e-15
    pair=direct[0]-direct[1];force=direct[1]-direct[2]
    assert pair>0 and force>0
    ratio=force/((1/3-1/4)/(4*np.pi))
    assert .9<ratio<1.2
    return {'independent_real_Fourier_grid':n,'largest_Fourier_vs_fft_error':float(error),
        'unit_kappa_pair_energy_at_three':pair,'static_step_force_three_to_four':force,
        'ratio_to_continuum_step_force':ratio,'alpha_or_common_physical_length_inferred':False}


def report():
    paths=[TARGET/name for name in ('source.py','result.py','PROOF.md','results.json')]+[Path(__file__)]
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    receipt=json.loads((TARGET/'results.json').read_text())
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in receipt['input_sha256'].items())
    result={'schema':'oph.independent_native_current_coulomb_review.v1','input_sha256':pins,
        'cochains':cochain_review(),'native_operator':native_operator_review(),'Green_profile':green_profile_review(),
        'scope':'The native-current Gauss quadratic is a valid original-source observable. The original unscreened domain has the retained linear link-cost obstruction. Removing that penalty defines a distinct source with Coulomb static trial energies, without transferring the original interacting ground or transport probabilities; kappa and physical unit conversion are not derived.',
        'verdict':'PASS_INDEPENDENT_NATIVE_CURRENT_COULOMB_REVIEW'}
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in pins.items())
    return result


if __name__=='__main__':
    result=report();(HERE/'coulomb_review.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(result['verdict'])
