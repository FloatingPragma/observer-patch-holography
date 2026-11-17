"""Exact transfer obstructions and full-band Dirac-current identity; source-bound rotation control."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement
import hashlib
import importlib.util
import json
import numpy as np
import sympy as s

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4];RER=ROOT/'reverse-engineering-reality'
paths=[HERE/'PLAN.md',ROOT/'openwave-tuning/next/CHAIN7/color_weak_native_join/PROOF.md',
 RER/'model/openwave_common_mass_source/source.py',RER/'model/openwave_boosted_records/source.py',
 ROOT/'openwave-tuning/common_dynamics/autonomous_dirac/PROOF.md',
 HERE.parent/'amended_current.py',HERE.parent/'AMENDED_CURRENT.json',
 HERE.parent.parent/'hadrons/audit.py',HERE.parent.parent/'hadrons/inputs/0_source.py']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();pins={str(p.relative_to(ROOT)):sha(p)for p in paths}

I=s.eye(2);sig=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
b=s.symbols('b0:4',real=True);be=b[0]*I+sum((b[j+1]*sig[j]for j in range(3)),s.zeros(2))
sol=s.solve([z for a in sig for z in be*a+a*be],b,dict=True)
assert sol==[{q:0 for q in b}]
# Single unit-charge Dirac carrier to up/down spin carrier: no nonzero
# charge-preserving injection exists because charge spectra are disjoint.
Qq=s.diag(s.Rational(2,3),s.Rational(2,3),-s.Rational(1,3),-s.Rational(1,3))
assert (s.eye(4)-Qq).det()!=0

px,py,pz,m,E=s.symbols('px py pz m E',real=True);mom=[px,py,pz]
P=sum((mom[j]*sig[j]for j in range(3)),s.zeros(2));q=E+m
relation=E**2-m*m-px*px-py*py-pz*pz
def reduce(expr):
    n,d=s.fraction(s.cancel(expr));return s.rem(s.Poly(s.expand(n),E),s.Poly(relation,E)).as_expr()
norm=((q*q)*I+P*P)/(2*E*q)
assert all(reduce(z)==0 for z in norm-I)
projected=[]
for j,a in enumerate(sig):
    actual=(q*q*a+P*a*P)/(2*E*q)
    claimed=m*a/E+mom[j]*P/(E*q)
    assert all(reduce(z)==0 for z in actual-claimed)
    projected.append([[str(s.factor(z))for z in row]for row in claimed.tolist()])
# Positive/negative spinor frames are [qI;P] and [-P;qI], each /sqrt(2Eq).
up=(q*I).col_join(P);um=(-P).col_join(q*I)
assert (up.H*um).applyfunc(s.simplify)==s.zeros(2)
assert all(reduce(z)==0 for z in (up*up.H+um*um.H)/(2*E*q)-s.eye(4))
gamma5=s.zeros(4);gamma5[:2,2:]=I;gamma5[2:,:2]=I
assert all(reduce(z)==0 for z in up.H*gamma5*up/(2*E*q)-P/E)
for j,a in enumerate(sig):
    alpha=s.zeros(4);alpha[:2,2:]=a;alpha[2:,:2]=a
    assert all(reduce(z)==0 for z in up.H*alpha*up/(2*E*q)-mom[j]*I/E)
isotropic=m/E+(px*px+py*py+pz*pz)/(3*E*q)
assert reduce(isotropic-s.Rational(1,3)-2*m/(3*E))==0

spec=importlib.util.spec_from_file_location('axial_bridge_existing_current',paths[5])
cur=importlib.util.module_from_spec(spec);spec.loader.exec_module(cur)
src=cur.src;had=cur.had
labels=list(combinations_with_replacement(range(8),3));counts=[tuple(lab.count(j)for j in range(8))for lab in labels]
index={c:i for i,c in enumerate(counts)};assert len(counts)==120
vecs=[]
for weakup in [1,2]:
    d=src.data(2,F(1),weakup);vals=had.extension_values(d,'current_square')
    H=src.dense(d)[:d['dimension'],:d['dimension']]+np.diag([float(v)for v in vals[:d['dimension']]])
    ev,U=np.linalg.eigh(H);v=np.zeros(120,complex)
    for i,c in enumerate(d['counts']):v[index[c]]=U[i,0]
    # Match the existing certificate's deterministic real eigensolver convention.
    vecs.append(v)
n,pr=vecs
J=[]
for mu in [0,3]:
    M=np.zeros((120,120),complex)
    for col,c in enumerate(counts):
        for cell in range(2):
            for spin in range(2):
                dest,z=src.act(c,4*cell+2*spin,4*cell+2*spin+1)
                if dest is not None:M[index[dest],col]+=src.value(z)*(1 if mu==0 or spin==0 else -1)
    J.append(M)
base=[np.vdot(pr,M@n)for M in J]
if base[0].real<0:pr=-pr;base=[-x for x in base]
prior=json.loads(paths[6].read_text())
assert abs(base[0].real-prior['vector_float'])<1e-10
assert abs(base[1].real-prior['axial_float'])<1e-10
# Local proper spin rotation -i sigma_x at cell1, lifted to all3 bosons in
# the exact colour-epsilon singlet representation. No occupation discarded.
Rot=np.zeros((120,120),complex)
for col,c in enumerate(counts):
    target=c[:4]+c[6:8]+c[4:6]
    Rot[index[target],col]=(-1j)**sum(c[4:])
assert np.max(abs(Rot.conj().T@Rot-np.eye(120)))==0
charge=np.array([sum(F(2,3)*c[4*x+2*sp]-F(1,3)*c[4*x+2*sp+1]for x in range(2)for sp in range(2))for c in counts],float)
assert np.max(abs(np.diag(charge)@Rot-Rot@np.diag(charge)))==0
nn,pp=Rot@n,Rot@pr
transported=[np.vdot(pp,(Rot@M@Rot.conj().T)@nn)for M in J]
untransported=[np.vdot(pp,M@nn)for M in J]
assert max(abs(x-y)for x,y in zip(base,transported))<1e-13

out=dict(schema='oph.axial_dirac_transfer_audit.v1',input_sha256=pins,
 local_massive_3plus1_Dirac_spin2_extension_exists=False,
 local_Clifford_scope='fixed cell and fixed weak/colour charge: two physical spin coordinates; all3 spatial Pauli alphas held fixed',
 beta_anticommutant_solution={str(k):str(v)for k,v in sol[0].items()},
 unit_charge_Dirac_to_fractional_quark_intertwiner_exists=False,
 charge_scope='single unit-charge one-particle carrier preserving the same charge observable; does not exclude an amended fractional-charge field source',
 exact_positive_spinor_isometry=True,exact_both_band_completeness=True,
 exact_positive_axial_projection=projected,
 exact_positive_vector_projection=['I','px I/E','py I/E','pz I/E'],
 exact_positive_axial_time_projection='(p dot sigma)/E',
 isotropic_factor='1/3+2m/(3E), conditional on an isotropic separable orbital average',
 physical_baryon_axial_factor_computed=False,
 actual_amended_source_control=dict(dimension=120,
   original_vector=float(base[0].real),original_axial=float(base[1].real),
   consistently_transported_vector=float(transported[0].real),consistently_transported_axial=float(transported[1].real),
   state_rotated_current_not_transported_vector=float(untransported[0].real),
   state_rotated_current_not_transported_axial=float(untransported[1].real),
   exact_charge_commutator_zero=True,current_transport_error=float(max(abs(x-y)for x,y in zip(base,transported)))),
 existing_bridges={
  'CHAIN7':'spin2 x weak2 quarks; zero free lepton generator; explicitly not identified with193-mode Dirac source',
  'common192_CAR':'144 colourful modes use number/scalar hopping; Bogoliubov Dirac squares are installed only on48 lepton modes',
  'boosted_records':'positive lepton kinetic extension on one-lepton/quark-vacuum reducing sector',
  'autonomous_Dirac':'four spin components, unit electric charge, actual field/clock; no colour/weak quark attachment'},
 missing_transfer=['extra Dirac/chiral quark coordinates or a proved staggered reconstruction retaining local charge',
  'same-source relativistic quark Hamiltonian and interacting reference/ground',
  'fractional electric and colour Gauss operator map',
  'charged vector/axial current intertwining and complement leakage',
  'shared kinetic/rest/field clock and units'],
 source_changed=False,empirical_axial_parameter_inserted=False,physical_axial_ratio_solved=False)
assert pins=={str(q.relative_to(ROOT)):sha(q)for q in paths}
(HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out['actual_amended_source_control'],indent=2))
