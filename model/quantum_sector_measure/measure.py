"""Exact complete native Cayley/position channel and electric measures."""
from __future__ import annotations
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sympy as s
from model.native_quantum_gauss import algebra as native
from model.common_source import native_matter_attachment as matter
from model.particle_common_contract.provenance import import_closure

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
DURATION=s.Rational(1,8)


def require(condition,message):
    if not condition:raise ValueError(message)


def clean(a):return a.applyfunc(s.simplify)
def zero(a):return all(s.simplify(x)==0 for x in a)
def record(a):return [[str(x)for x in row]for row in a.tolist()]


def pulse(d,duration=DURATION):
    require(type(d)is int and d>0,'positive integer multiplet dimension required')
    duration=native.exact_real(duration,nonnegative=True)
    h=s.Matrix([[0,s.sqrt(d)],[s.sqrt(d),1]])
    u=clean((s.eye(2)+s.I*duration*h/2).inv()*(s.eye(2)-s.I*duration*h/2))
    require(zero(u.H*u-s.eye(2)),'Cayley unitarity failed')
    return h,u


@lru_cache(maxsize=1)
def source():
    blocks=native.native_blocks();parts=tuple(indices for _,indices in blocks)
    representation=matter.representation();selected=list(representation['selected'])
    weak=representation['operators']['weak_C2'].extract(selected,selected)
    result=[]
    for label,indices in blocks:
        d=len(indices);a,b=native.reference_states(15,indices)
        require(native.hamiltonian_action(a,15,blocks=parts)==b,'ambient A action failed')
        require(native.hamiltonian_action(b,15,blocks=parts)==native.sum_states(native.scale(a,d),b),
                'ambient B action failed')
        require(sum(abs(v)**2 for v in a.values())==d and sum(abs(v)**2 for v in b.values())==d*d,
                'source normalization failed')
        for g in native.native_currents():
            for vertex in (0,1):
                require(not native.gauss_action(g,a,vertex)and not native.gauss_action(g,b,vertex),
                        'ambient Gauss constraint failed')
        c=weak[indices[0],indices[0]]
        require(c in (0,s.Rational(3,4))and all(weak[i,i]==c for i in indices),'actual weak sector failed')
        h,u=pulse(d);p=s.simplify(abs(u[1,0])**2)
        result.append(dict(label=label,d=d,weak_C2=c,h=h,u=u,p=p))
    return result


def channel(rho):
    rows=source();u=s.diag(*(row['u']for row in rows))
    p0=s.diag(*([1,0]*5));p1=s.eye(10)-p0
    require(rho.shape==(10,10),'the admitted carrier has dimension ten')
    return clean(sum((p*u*rho*u.H*p for p in (p0,p1)),s.zeros(10)))


def adjoint(a):
    rows=source();u=s.diag(*(row['u']for row in rows))
    p0=s.diag(*([1,0]*5));p1=s.eye(10)-p0
    return clean(sum((u.H*p*a*p*u for p in (p0,p1)),s.zeros(10)))


def electric_projector():
    return s.diag(*[x for row in source()for x in (0,int(row['weak_C2']!=0))])


def closure_certificate():
    rows=source();pj=electric_projector();effect=adjoint(pj)
    r=next(i for i,row in enumerate(rows)if row['d']==2);u=rows[r]['u']
    z=s.simplify((u.H*s.diag(0,1)*u)[0,1])
    plus=s.zeros(10,1);minus=s.zeros(10,1)
    plus[2*r]=minus[2*r]=1/s.sqrt(2)
    plus[2*r+1]=1/s.sqrt(2);minus[2*r+1]=-1/s.sqrt(2)
    before=[s.simplify((v.H*pj*v)[0])for v in (plus,minus)]
    after=[s.simplify((v.H*effect*v)[0])for v in (plus,minus)]
    require(before==[s.Rational(1,2)]*2 and after[0]!=after[1],'coherent nonclosure witness failed')
    # Even diagonal states in the same electric j=0 fibre have different exits.
    singlet=next(i for i,row in enumerate(rows)if row['weak_C2']==0)
    exits=[effect[2*r,2*r],effect[2*singlet,2*singlet]]
    require(exits[0]>0 and exits[1]==0,'within-sector diagonal witness failed')
    return dict(electric_projector_half=record(pj),adjoint_half=record(effect),
        scalar_sector_adjoint_closure=False,coherent_witness=dict(block_label=rows[r]['label'],
        states='(a_r+b_r)/sqrt(2), (a_r-b_r)/sqrt(2)',initial_half_probabilities=list(map(str,before)),
        next_half_probabilities=list(map(str,after)),adjoint_offdiagonal=str(z)),
        within_j_zero_diagonal_exit_probabilities=list(map(str,exits)),
        restricted_closed_family='rho=sum_r (x_r |a_r><a_r|+y_r |b_r><b_r|), x_r,y_r>=0; sum(x_r+y_r)=1',
        family_kernels=[dict(label=row['label'],dimension=row['d'],p=str(row['p']),
          kernel=record(s.Matrix([[1-row['p'],row['p']],[row['p'],1-row['p']]])))for row in rows])


def stationary_certificate():
    rows=source();pairs=[]
    # Every output is endpoint-diagonal. Its (r,s) fixed equation has two unknowns.
    for i,left in enumerate(rows):
        for j,right in enumerate(rows):
            u,v=left['u'],right['u']
            m=s.Matrix([[u[0,0]*s.conjugate(v[0,0]),u[0,1]*s.conjugate(v[0,1])],
                        [u[1,0]*s.conjugate(v[1,0]),u[1,1]*s.conjugate(v[1,1])]])
            m=clean(m-s.eye(2));same=left['d']==right['d']
            require(m.rank()==(1 if same else 2),'full fixed-operator pair rank failed')
            if same:require(zero(m*s.ones(2,1)),'equal-copy stationary operator failed')
            pairs.append(dict(r=i,s=j,rank=m.rank(),fixed_equal_endpoints=same))
    require(sum(2-x['rank']for x in pairs)==7,'full fixed-operator dimension failed')
    weak_weight=sum(row['weak_C2']!=0 for row in rows)
    tracial=s.eye(10)/10
    require(channel(tracial)==tracial,'full tracial state not fixed')
    centre=s.eye(10)-2*electric_projector()
    balance=s.trace(tracial*centre)
    require(balance==s.Rational(3,5)and weak_weight==2,'centre trace failed')
    return dict(pair_fixed_equations=pairs,complex_fixed_operator_dimension=7,
        density_family_real_affine_dimension=6,
        full_stationary_family='direct sum over d=1,2,3,6 of M_d tensor I_endpoint/2; M_1,M_2,M_6 are nonnegative scalars, M_3 is a positive 2x2 matrix on the two distinct d=3 copies; sum_d Tr M_d=1',
        dephased_stationary_family='x_r=y_r=w_r/2; w_r>=0, sum_r w_r=1',
        family_weights_selected_by_channel=False,
        endpoint_convergence_factors={row['label']:str(1-2*row['p'])for row in rows},
        centre_balance='1-w_(d=2)-w_(d=6)',centre_balance_range=['0','1'],
        tracial_ten_dimensional_preparation=dict(weights=['1/5']*5,half_sector='1/5',centre_balance=str(balance)),
        native_default_d6_preparation=dict(initial='a_(d=6)',limiting_half_sector='1/2',limiting_centre_balance='0'),
        imported_chain_output=dict(contract='v=E_cell*b*<A>, A=(-1)^(2j)',
            dimensionless_v_over_Ecell_b_from_tracial_preparation=str(balance),
            dimensionless_v_over_Ecell_b_from_native_d6_preparation='0',
            physical_electroweak_scale_selected=False,physical_gauge_coupling_selected=False))


def histories(rounds=3):
    require(type(rounds)is int and 1<=rounds<=8,'one to eight complete rounds required')
    rows=[]
    for row in source():
        branches={'':s.Matrix([1,0])}
        for _ in range(rounds):
            branches={word+str(z):clean(p*row['u']*v)for word,v in branches.items()
                      for z,p in enumerate((s.diag(1,0),s.diag(0,1)))}
        probs={word:str(s.simplify((v.H*v)[0]))for word,v in branches.items()}
        require(sum(s.Rational(p)for p in probs.values())==1,'complete history probability lost')
        # Literal parent replay uses the unnormalized Gram basis and R/L labels.
        parent=native.complete_history(row['d'],duration=DURATION,rounds=rounds)
        replay={x['word'].translate(str.maketrans('RL','01')):x['probability']for x in parent['outcomes']}
        require(replay==probs,'parent Gram-basis complete histories differ')
        rows.append(dict(label=row['label'],complete_records=probs,parent_complete_history=parent))
    return rows


def bindings():
    names={p.relative_to(ROOT).as_posix()for p in HERE.iterdir()if p.suffix in ('.py','.md')}
    names.update(('model/native_quantum_gauss/PROOF.md','model/seam_class_balance/PROOF.md',
                  'model/repair_sector_bridge/PROOF.md','model/source_strong_coupling/PROOF.md'))
    # The imported source resolves its geometric carrier through this retained tree.
    retained=ROOT/'evidence/observer_dynamics_20260925/oph-physics-sim'
    for package in ('oph_exact','oph_fpe'):
        names.update(p.relative_to(ROOT).as_posix()for p in (retained/package).rglob('*.py'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(names))}


def certificate():
    before=bindings();rows=source()
    result=dict(schema='oph.quantum_sector_measure.v1',source_bindings=before,
        source=dict(parent='native_quantum_gauss.algebra',link_dimension=60,
          ambient_one_matter_reference_dimension=27000,admitted_invariant_Gauss_subspace_dimension=10,
          declared_Cayley_duration=str(DURATION),J='1',epsilon='1',
          complete_one_round_records=['endpoint A','endpoint B'],
          external_reference_preserved=True,finite_carrier_exact=True,
          infinite_SU2_rotor_approximation_claimed=False,cutoff_tail_bound_required=False,
          blocks=[dict(label=r['label'],dimension=r['d'],weak_C2=str(r['weak_C2']),
                       H=record(r['h']),U=record(r['u']))for r in rows]),
        closure=closure_certificate(),stationarity=stationary_certificate(),complete_histories=histories(),
        physical_promotions=[])
    require(bindings()==before,'source changed during computation')
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--check',action='store_true')
    args=p.parse_args();result=certificate();data=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if args.write:(HERE/'results.json').write_text(data)
    if args.check:require((HERE/'results.json').read_text()==data,'sector receipt differs')
    print(json.dumps(dict(sector_closure=False,stationary_operator_dimension=7,
                         tracial_centre_balance='3/5',bindings=len(result['source_bindings']))))


if __name__=='__main__':main()
