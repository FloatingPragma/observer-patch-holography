#!/usr/bin/env python3
"""Exact optional Weinberg interaction on the same conditional source module.

One existing lepton doublet and the same declared Higgs suffice algebraically.
The dimension-five coefficient, scale and validity domain remain new inputs.
"""
from __future__ import annotations
from evidence_contract import public_receipt
import hashlib,json,sys
from pathlib import Path
sys.dont_write_bytecode=True
import sympy as s
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PRIOR=HERE
sys.path.insert(0,str(PRIOR))
import representation_identification as rep


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def mat(p):return s.SparseMatrix(*p['shape'],{(i,j):s.sympify(v) for i,j,v in p['entries']})
def zero(a):return all(s.simplify(v)==0 for v in a)


def main():
    source=PRIOR/'results/representation_identification_operators.json'
    raw=json.loads(source.read_text());_,dg=rep.fock_data();ix=raw['even_indices']
    one=[mat(x) for x in raw['spin_lifted_one_body_generators']]
    currents=[dg(g).extract(ix,ix) for g in one];higgs=[g[3:,3:] for g in one]
    q=mat(raw['even_identification_operators']['q']);t3=mat(raw['even_identification_operators']['T3'])
    indices=[i for i in range(15) if q[i,i]==-3];assert len(indices)==2
    equations=s.Matrix.vstack(*[s.kronecker_product(g.extract(indices,indices),s.eye(2))+
        s.kronecker_product(s.eye(2),h) for g,h in zip(currents,higgs)])
    kernel=equations.nullspace();assert len(kernel)==1
    coupling=s.zeros(15,2)
    for a,i in enumerate(indices):
        for b in range(2):coupling[i,b]=s.conjugate(kernel[0][2*a+b])
    scale=s.sqrt((coupling[:,1].H*coupling[:,1])[0]);coupling=s.simplify(coupling/scale)
    checks=0
    x,y=s.symbols('H0 H1');h=s.Matrix([x,y]);ell=coupling*h;mass=ell*ell.T
    assert mass.T==mass
    for g,hg in zip(currents,higgs):
        derivative=coupling*hg*h
        assert zero(g.T*ell+derivative)
        assert zero(g.T*mass+mass*g+derivative*ell.T+ell*derivative.T)
        checks+=1
    neutral=mass.subs({x:0,y:1})
    assert neutral.rank()==1
    charge=t3+q/6;assert zero(charge.T*neutral+neutral*charge)
    for i in range(3):
        for j in range(i,3):
            k=s.zeros(3);k[i,j]=1;k[j,i]=1
            complete=s.kronecker_product(neutral,k)
            assert complete.T==complete
    eps=s.Matrix([[0,1],[-1,0]])
    physical=s.kronecker_product(eps,neutral)
    assert physical.T==-physical
    artifact=HERE/'results/neutrino_action_operators.json'
    artifact.write_text(json.dumps({'LH_invariant':rep.sparse(coupling),
        'neutral_mass_kernel':rep.sparse(neutral)},indent=2)+'\n', encoding='utf-8', newline='\n')
    paths=[source,PRIOR/'representation_identification.py',
        ROOT/'Lean/Screen/WeylYukawaConventions.lean']
    result={'status':'PASS_OPTIONAL_SAME_HIGGS_NEUTRINO_OPERATOR_WITH_EXPLICIT_NEW_COEFFICIENTS',
        'producer_sha256':sha(__file__),'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
        'artifacts':{artifact.name:sha(artifact)},'empirical_inputs':[],
        'same_M15_current_count_checked':checks,'same_Higgs_Y':'1/2',
        'invariant':'(L H)(L H), with the derived weak contraction and two-component Lorentz epsilon',
        'invariant_LH_tensor_dimension':1,'neutral_single_family_mass_rank':1,
        'three_family_symmetric_complex_coefficients':6,
        'mass_formula':'M_nu=kappa*v^2/Lambda in the declared neutral-Higgs-norm convention',
        'dimension_five_operator_is_in_minimal_renormalizable_SM':False,
        'scope':'An optional effective interaction on the same internal module and declared Higgs; neither kappa, Lambda, a neutrino mechanism nor the selected physical vacuum is derived. This is not a massive-neutrino prediction.',
        'source_selected_neutrino_masses':False}
    (HERE/'results/neutrino_action_extension.json').write_text(json.dumps(public_receipt(result),indent=2)+'\n', encoding='utf-8', newline='\n')
    print(result['status'])


if __name__=='__main__':main()
