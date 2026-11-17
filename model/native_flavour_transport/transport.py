"""Source-derived finite-word CP with unchanged native singular values."""
from __future__ import annotations
import argparse
from functools import lru_cache
import hashlib
from itertools import product
import json
from pathlib import Path
import numpy as np
import sympy as s
from model.common_source import native_repair_mass as repair
from model.common_source import native_family_record_memory as memory
from model.common_source import calibrated_matter as matter
from model.common_source.chiral_candidate import action
from model.particle_common_contract.provenance import import_closure
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
ALPHABET=(2,3,10)
PRIMARY_WORD=(2,3,10)
PHASE=-s.I


def require(c,m):
    if not c:raise ValueError(m)


def clean(a):return a.applyfunc(s.expand)
def zero(a):return all(s.expand(v)==0 for v in a)
def record(a):return [[str(s.expand(v))for v in row]for row in a.tolist()]
def numeric(a):return np.array(a.evalf(25),complex)


@lru_cache(maxsize=1)
def source():
    data=repair.architecture();projectors=memory._projectors(False)
    y=data['Y'];h=clean(y*y);x=s.Symbol('x')
    poly=s.Poly(h.charpoly(x).as_expr(),x,extension=s.sqrt(5));disc=s.expand(s.discriminant(poly.as_expr(),x))
    require(disc.is_positive is True,'native squared spectrum is not simple')
    require(y.det()!=0 and s.trace(y)==1,'native positive coefficient changed')
    return dict(y=y,h=h,projectors=projectors,discriminant=disc,
                characteristic_polynomial=poly.as_expr())


def word_unitary(word,phase=PHASE):
    require(all(type(k)is int and k in ALPHABET for k in word),'word outside predeclared alphabet')
    require(s.simplify(phase*s.conjugate(phase))==1,'unit phase required')
    u=s.eye(3)
    for k in word:u=clean((s.eye(3)+(phase-1)*source()['projectors'][k])*u)
    require(zero(u.H*u-s.eye(3)),'native geometric pulse word lost unitarity')
    return u


def invariant(word,phase=PHASE):
    v=word_unitary(word,phase);h=source()['h'];hd=clean(v*h*v.H)
    comm=clean(h*hd-hd*h);cube=s.expand(s.trace(clean(comm*comm*comm)))
    j=s.radsimp(s.simplify(cube/(6*s.I*source()['discriminant'])))
    require(s.simplify(s.im(j))==0,'CP invariant is not real')
    return v,cube,s.expand(j)


@lru_cache(maxsize=1)
def controls():
    rows=[]
    for length in range(4):
        for word in product(ALPHABET,repeat=length):
            _,cube,j=invariant(word)
            rows.append(dict(word=list(word),CP_trace_cube_exact=str(cube),J_exact=str(j),
                             J_numeric=float(s.N(j,18))))
    require(len(rows)==40,'predeclared short-word set incomplete')
    return rows


def full_source_records():
    data=source();v,cube,j=invariant(PRIMARY_WORD);y=data['y']
    require(j!=0,'predeclared primary word has no CP contrast')
    eig,f=np.linalg.eigh(numeric(y))
    require(np.all(eig>0)and np.min(np.diff(eig))>1e-8,'native positive simple spectrum lost')
    for col in range(3):
        pivot=np.argmax(abs(f[:,col]));f[:,col]*=np.conjugate(f[pivot,col])/abs(f[pivot,col])
    mix=f.conj().T@numeric(v)@f
    jnum=float(np.imag(mix[0,0]*mix[1,1]*np.conj(mix[0,1]*mix[1,0])))
    require(abs(jnum-float(s.N(j,18)))<2e-13,'invariant/mass-frame CP sign or value differs')
    # Same native full45 invariant tensors, with unchanged positive Y in every channel.
    ys={'up':y,'down':clean(v*y),'charged_lepton':y}
    k=s.zeros(45);labels=action.ordering()['all_left_components']
    frames=[f.copy()for _ in labels]
    ward_checks=0
    for name,yy in ys.items():
        cs=matter._oriented_components(name)
        for primitive,g in zip(action._representation()[1],action.generators(families=False)):
            hh=primitive[3:,3:]
            if action._representation()[4][name]['conjugate_Higgs']:hh=-hh.T
            for component in range(2):
                require(zero(g.T*cs[component]+cs[component]*g+sum((hh[b,component]*cs[b]for b in range(2)),s.zeros(15))),
                        'actual two-component Higgs Ward identity failed')
                ward_checks+=1
        for internal in [matter._oriented_components(name)[1]]:
            directed=s.kronecker_product(internal,yy.conjugate());k+=directed+directed.T
        for (left,right),sign in matter._oriented_components(name)[1].todok().items():
            frames[right]=float(sign)*f
            if name=='down':frames[left]=numeric(v)@f
    frame=np.zeros((45,45),complex)
    for i,b in enumerate(frames):frame[3*i:3*i+3,3*i:3*i+3]=b
    charge=action.electric_charge()
    require(zero(charge.T*k+k*charge),'full native electric-neutral pairing failed')
    delta=s.kronecker_product(s.Matrix([[0,1],[-1,0]]),k.conjugate())
    require(zero(delta+delta.T),'full90 CAR pairing is not antisymmetric')
    diagonal=frame.T@numeric(k)@frame
    errors=[]
    for name in ys:
        for(left,right),sign in matter._oriented_components(name)[1].todok().items():
            errors.append(float(np.max(abs(diagonal[3*left:3*left+3,3*right:3*right+3]-np.diag(eig)))))
    current=frame.conj().T@numeric(matter.weak_raising_current())@frame
    quarkpairs=[(i,k)for(i,k),z in matter.weak_raising_current(families=False).todok().items()
                if labels[i]['internal_type']=='Q']
    require(len(quarkpairs)==3,'three colour weak-current pairs required')
    errors.extend(float(np.max(abs(current[3*i:3*i+3,3*k:3*k+3]-mix)))for i,k in quarkpairs)
    require(max(errors)<2e-12,'actual full45 mass/current frame differs')
    rows=[];effect=np.zeros((2,2),complex)
    for r in range(3):
        for t in range(3):
            a=mix[r,0]*mix[t,1];b=mix[r,1]*mix[t,0]
            for label,sign in [('plus',1),('minus',-1)]:
                kraus=np.array([[a,-sign*1j*b]])/np.sqrt(2)
                effect+=kraus.conj().T@kraus
                orientation=0 if r==t else (1 if (r,t)in((0,1),(1,2),(2,0))else -1)
                rows.append(dict(pair=[r,t],Y=label,probability=float(abs((a-sign*1j*b)/2)**2),
                    exact_probability_formula=f'(B_{r}0*B_{t}1+B_{r}1*B_{t}0)/4-({sign*orientation})*J/2'))
    target={r['Y']:r['probability']for r in rows if r['pair']==[0,1]}
    require(np.max(abs(effect-np.eye(2)))<2e-12 and abs(sum(r['probability']for r in rows)-1)<2e-12,
            'complete18 current records differ')
    require(abs(target['minus']-target['plus']-jnum)<2e-13,'original CP contrast differs')
    conjugate=invariant(PRIMARY_WORD,s.I)[2]
    require(s.simplify(conjugate+j)==0,'complex-conjugate word control failed')
    return dict(primary_word=list(PRIMARY_WORD),pulse_phase='-i',transport=record(v),
        unchanged_native_mass_coefficient=record(y),native_squared_characteristic_polynomial=str(data['characteristic_polynomial']),
        native_squared_discriminant_exact=str(data['discriminant']),
        CP_trace_cube_exact=str(cube),J_exact=str(j),J_numeric=jnum,
        conjugate_word_J_exact=str(conjugate),singular_values_numeric=list(map(float,eig)),
        full45_Higgs_kernel=record(k),
        exact_projector_probability_definition='B_ri=Tr(P_r(Y) V P_i(Y) V-dagger), P_i=product_(j!=i)(Y-lambda_j I)/(lambda_i-lambda_j), ordered simple roots lambda_i of the source Y polynomial',full45_mass_current_max_numeric_residual=max(errors),
        all_three_colour_current_blocks_checked=True,both_Higgs_component_internal_Ward_checks=ward_checks,
        full45_electric_neutrality_and_full90_pairing_antisymmetry_exact=True,
        mass_frame_mixing_real=mix.real.tolist(),mass_frame_mixing_imag=mix.imag.tolist(),
        complete18_original_records=rows,complete_effect_identity_analytic=True,
        exact_record_definition='V_mass=F^T V F with ordered exact orthonormal eigenframe F of positive Y; K_rs,+/-=(V_r0 V_s1, -/+ i V_r1 V_s0)/sqrt(2)',
        signed_contrast_exact=str(j),numeric_records_are_floating_point_diagnostics=True,
        physical_family_assignment_or_CKM_fit=False)



def minimal_pulse_certificate():
    a,b,cost,sint=s.symbols('a b cost sint',real=True)
    shift=cost-s.I*sint-1
    quartet=(1+a*shift)*(1+b*shift)*s.conjugate(shift)**2*a*b
    proposed=-2*a*b*(1-a-b)*sint*(1-cost)
    error=s.Poly(s.expand(s.im(quartet)-proposed),sint)
    remainder=error.rem(s.Poly(sint**2+cost**2-1,sint)).as_expr()
    require(s.expand(remainder)==0,'generic rank-one pulse CP formula failed')
    y=source()['y'];p=source()['projectors'][2];x=s.Symbol('x')
    gram=s.Matrix(3,3,lambda i,j:s.expand(s.trace(p*clean(y**(i+j)))))
    disc=s.expand(s.discriminant(y.charpoly(x).as_expr(),x))
    product_weight=s.radsimp(s.simplify(gram.det()/disc))
    require(product_weight.is_positive is True,'actual port2 lacks three nonzero mass-frame weights')
    native_j=invariant((2,))[2]
    require(s.simplify(native_j+2*product_weight)==0,'actual port2 rank-one shortcut differs')
    return dict(generic_formula='J=-2*a*b*c*sin(phi)*(1-cos(phi)); a=p0^2,b=p1^2,c=p2^2; a+b+c=1',
        generic_unit_circle_polynomial_remainder='0',source_port=2,pulse_phase='-i',
        actual_Krylov_Gram=record(gram),native_Y_discriminant_exact=str(disc),
        product_of_three_mass_frame_weights_exact=str(product_weight),
        source_identity='a*b*c=det[Tr(P2 Y^(i+j))]_(i,j=0..2)/disc(det(xI-Y))',
        all_three_weights_nonzero_exact=True,single_pulse_J_exact=str(native_j),
        single_pulse_J_numeric=float(s.N(native_j,18)),
        zero_length_word_J='0',minimum_nonzero_pulse_count_in_predeclared_alphabet=1,
        physical_CKM_value_or_unique_architecture_selected=False)


def bindings():
    names={p.relative_to(ROOT).as_posix()for p in HERE.iterdir()if p.suffix in('.py','.md')}
    names.update(repair.source_bindings())
    names.update(('model/native_flavor_orbit/PROOF.md','model/source_flavour_identifiability/PROOF.md'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(names))}


def certificate():
    pins=bindings();out=dict(schema='oph.native_flavour_transport.v1',source_bindings=pins,all_checks_pass=True,
        source_contract=dict(family_factor=3,full_native_internal_rows=45,
            empirical_mass_or_mixing_inputs=[],architecture_changed='down coefficient Yd=V_word Y, Yu=Ye=Y',
            source_forces_chosen_word=False,original_record_memory_replaced=False,
            mean_coefficient_service_inherited=True,finite_correlated_memory_identified_with_static_Y=False),
        predeclared_scan=dict(alphabet=list(ALPHABET),pulse_phase='-i',all_lengths=[0,1,2,3],
            primary_word=list(PRIMARY_WORD),selection='listed geometric ports in ascending order; primary fixed before CP evaluation'),
        short_word_controls=controls(),primary=full_source_records(),minimal_rank_one_pulse=minimal_pulse_certificate(),physical_promotions=[])
    require(bindings()==pins,'source changed during transport calculation');return out


def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--check',action='store_true');a=p.parse_args()
    r=certificate();text=json.dumps(r,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if a.write:(HERE/'results.json').write_text(text)
    if a.check:require((HERE/'results.json').read_text()==text,'transport receipt differs')
    print(json.dumps(dict(J=r['primary']['J_numeric'],words=len(r['short_word_controls']))))


if __name__=='__main__':main()
