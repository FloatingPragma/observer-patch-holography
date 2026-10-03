#!/usr/bin/env python3
"""Independent 2I/McKay verifier using generator closure, not the producer.

This file intentionally has no import from mckay_golden_field_certificate and
contains no OPH dependency.  It generates 2I from two exact quaternionic
matrices, independently reconstructs classes/characters/fusion, verifies the
producer receipt, and only then compares each derived graph with affine E8.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
from collections import deque

ROOT=Path(__file__).resolve().parents[2]
RAW_DEFAULT=ROOT/"code/a5_closure/receipts/arithmon_mckay_golden_field_certificate.json"
OUT_DEFAULT=ROOT/"code/a5_closure/receipts/arithmon_mckay_golden_field_reference.json"

# This verifier uses raw rational pairs rather than the producer's Q5 class.
R=tuple[F,F]
ZERO=(F(0),F(0)); ONE=(F(1),F(0)); HALF=(F(1,2),F(0))
SQRT5=(F(0),F(1)); PHI=(F(1,2),F(1,2)); PHI_INV=(F(-1,2),F(1,2))

def radd(x,y):return (x[0]+y[0],x[1]+y[1])
def rneg(x):return (-x[0],-x[1])
def rsub(x,y):return radd(x,rneg(y))
def rmul(x,y):return (x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def rdiv(x,y):
 d=y[0]*y[0]-5*y[1]*y[1]
 if not d:raise ZeroDivisionError
 return ((x[0]*y[0]-5*x[1]*y[1])/d,(x[1]*y[0]-x[0]*y[1])/d)
def rgalois(x):return (x[0],-x[1])
def rtext(x):
 a,b=x
 def f(q):return str(q.numerator) if q.denominator==1 else f"{q.numerator}/{q.denominator}"
 if not b:return f(a)
 if not a:return f(b)+"*sqrt(5)"
 return f(a)+("+" if b>0 else "-")+f(abs(b))+"*sqrt(5)"
def quadratic_minpoly_text(value):
 a,b=value
 if b==0:raise ValueError("rational trace has no quadratic field generator")
 linear=-2*a; constant=a*a-5*b*b
 if radd(radd(rmul(value,value),rmul((linear,F(0)),value)),(constant,F(0)))!=ZERO:
  raise ValueError("computed quadratic does not annihilate exact trace")
 def term(c,var):
  if not c:return ""
  body=(var if abs(c)==1 else f"{abs(c)}*{var}") if var else str(abs(c))
  return ("+" if c>0 else "-")+body
 return "X^2"+term(linear,"X")+term(constant,"")
def rint(x,gate):
 if x[1] or x[0].denominator!=1:raise ValueError(f"{gate}: noninteger exact value {rtext(x)}")
 return x[0].numerator

Q=tuple[R,R,R,R]
ID=(ONE,ZERO,ZERO,ZERO); MINUS=(rneg(ONE),ZERO,ZERO,ZERO)
def qmul(x:Q,y:Q)->Q:
 a,b,c,d=x;e,f,g,h=y
 return (rsub(rsub(rsub(rmul(a,e),rmul(b,f)),rmul(c,g)),rmul(d,h)),
         radd(rsub(radd(rmul(a,f),rmul(b,e)),rmul(d,g)),rmul(c,h)),
         radd(rsub(radd(rmul(a,g),rmul(c,e)),rmul(b,h)),rmul(d,f)),
         radd(rsub(radd(rmul(a,h),rmul(b,g)),rmul(c,f)),rmul(d,e)))
def qbar(x):return (x[0],rneg(x[1]),rneg(x[2]),rneg(x[3]))
def qnorm(x):
 out=ZERO
 for t in x:out=radd(out,rmul(t,t))
 return out
def qkey(x):return tuple((z[0].numerator,z[0].denominator,z[1].numerator,z[1].denominator) for z in x)
def qsort(x):return sorted(x,key=qkey)

def generated_group():
 # Two exact standard H4 quaternion roots generate the carrier.  This is a
 # generator-closure algorithm, independent of direct coordinate enumeration.
 g=(ZERO,HALF,rmul(PHI,HALF),rmul(PHI_INV,HALF))
 h=(HALF,HALF,HALF,HALF)
 steps=(g,qbar(g),h,qbar(h))
 seen={ID}; todo=deque([ID])
 while todo:
  x=todo.popleft()
  for s in steps:
   y=qmul(x,s)
   if y not in seen:seen.add(y);todo.append(y)
 return seen,(g,h)

def group_center(group):return {x for x in group if all(qmul(x,y)==qmul(y,x) for y in group)}

def conjugacy(group):
 rem=set(group);out=[]
 while rem:
  x=min(rem,key=qkey)
  cl={qmul(qmul(g,x),qbar(g)) for g in group}
  if not cl<=rem:raise ValueError("independent conjugacy classes overlap")
  out.append(cl);rem-=cl
 if set().union(*out)!=group:raise ValueError("independent classes fail to partition group")
 return out

def trace_doublet(q):return rmul((F(2),F(0)),q[0])
def sigchar(chi):return {g:rgalois(v) for g,v in chi.items()}
def cin(a,b,G):
 total=ZERO
 for g in G:total=radd(total,rmul(a[g],b[g]))
 return rdiv(total,(F(len(G)),F(0)))
def ceq(a,b,G):return all(a[g]==b[g] for g in G)
def cscale_sub(a,n,b,G):return {g:rsub(a[g],rmul((F(n),F(0)),b[g])) for g in G}

def reconstruct_irreps(G,chi):
 powers=[{g:ONE for g in G},dict(chi)]
 for _ in range(len(G)-1):
  powers.append({g:rsub(rmul(chi[g],powers[-1][g]),powers[-2][g]) for g in G})
 known=[];labels=[]
 for n,cand in enumerate(powers):
  for lab,c in ((f"Sym^{n}(V)",cand),(f"sigma(Sym^{n}(V))",sigchar(cand))):
   residual=dict(c)
   for k in known:
    m=rint(cin(residual,k,G),"projection")
    if m<0:raise ValueError("negative constituent multiplicity")
    if m:residual=cscale_sub(residual,m,k,G)
   if all(v==ZERO for v in residual.values()):continue
   if rint(cin(residual,residual,G),"residual norm")!=1:continue
   if any(rint(cin(residual,k,G),"cross norm")!=0 for k in known):raise ValueError("irreducible candidates not orthogonal")
   known.append(residual);labels.append(lab+" residual")
   ds=[rint(x[ID],"degree") for x in known]
   if sum(d*d for d in ds)==len(G):
    for i,a in enumerate(known):
     for j,b in enumerate(known):
      if rint(cin(a,b,G),"orthogonality")!=(1 if i==j else 0):raise ValueError("complete system not orthonormal")
    return known,labels,ds
 raise ValueError("independent symmetric-power reconstruction incomplete")

def fusion(G,irreps,faithful):
 ds=[rint(x[ID],"dimension") for x in irreps];M=[]
 for i,a in enumerate(irreps):
  tens={g:rmul(faithful[g],a[g]) for g in G}
  row=[rint(cin(tens,b,G),"multiplicity") for b in irreps]
  if any(m<0 for m in row):raise ValueError("negative fusion multiplicity")
  if sum(x*d for x,d in zip(row,ds))!=2*ds[i]:raise ValueError("fusion dimension conservation failed")
  M.append(row)
 if any(M[i][j]!=M[j][i] for i in range(len(M)) for j in range(len(M))):raise ValueError("fusion matrix nonsymmetric")
 return M,ds

def graph_info(M,ds):
 n=len(M)
 if any(M[i][i] for i in range(n)):raise ValueError("McKay loop")
 if any(M[i][j] not in (0,1) for i in range(n) for j in range(n)):raise ValueError("non-simply-laced McKay edge")
 edges=[[i,j] for i in range(n) for j in range(i+1,n) if M[i][j]]
 seen={0};q=[0]
 while q:
  u=q.pop()
  for v in range(n):
   if M[u][v] and v not in seen:seen.add(v);q.append(v)
 eig=all(sum(M[i][j]*ds[j] for j in range(n))==2*ds[i] for i in range(n))
 return {"vertices":n,"edges":edges,"edge_count":len(edges),"connected":len(seen)==n,
         "tree":len(seen)==n and len(edges)==n-1,"simply_laced":True,
         "dimension_vector":ds,"dimension_vector_2_eigenvector":eig}

# Independent target graph: affine E8 has a trivalent node with arm lengths 1,2,5.
AFFINE_E8_EDGES=((0,1),(0,2),(2,3),(0,4),(4,5),(5,6),(6,7),(7,8))

def isomorphisms(edges_a,edges_b,n):
 A=[set() for _ in range(n)];B=[set() for _ in range(n)]
 for u,v in edges_a:A[u].add(v);A[v].add(u)
 for u,v in edges_b:B[u].add(v);B[v].add(u)
 order=sorted(range(n),key=lambda v:(-len(A[v]),v))
 candidates={v:[w for w in range(n) if len(A[v])==len(B[w])] for v in range(n)}
 result=[];mapping={};used=set()
 def rec(k):
  if k==n:result.append(tuple(mapping[i] for i in range(n)));return
  v=order[k]
  for w in candidates[v]:
   if w in used:continue
   ok=True
   for u,z in mapping.items():
    if ((u in A[v]) != (z in B[w])):ok=False;break
   if not ok:continue
   mapping[v]=w;used.add(w);rec(k+1);used.remove(w);del mapping[v]
 rec(0)
 return result

def coordinate_set_digest(group):
 rows=sorted([qkey(x) for x in group])
 return hashlib.sha256(json.dumps(rows,separators=(",",":")).encode()).hexdigest()

def verify(raw):
 G,gens=generated_group()
 if len(G)!=120:raise ValueError(f"generator closure has {len(G)} elements, expected the standard 2I order 120")
 if group_center(G)!={ID,MINUS}:raise ValueError("generator-closure center is not +/-1")
 if any(qnorm(g)!=ONE for g in G):raise ValueError("generator closure escaped unit quaternions")
 if any(qmul(a,b) not in G for a in G for b in G):raise ValueError("generated carrier not closed")
 if coordinate_set_digest(G)!=raw["source_group"]["canonical_element_set_sha256"]:
  raise ValueError("generator-closure carrier differs from direct coordinate producer set")
 chi={g:trace_doublet(g) for g in G}
 if rint(cin(chi,chi,G),"spin norm")!=1:raise ValueError("faithful doublet character not irreducible")
 if chi[ID]!=(F(2),F(0)):raise ValueError("spin dimension mismatch")
 if chi[MINUS]!=(F(-2),F(0)):raise ValueError("central -1 action mismatch")
 classes=conjugacy(G)
 sizes=sorted(len(cl) for cl in classes)
 irreps,labels,dims=reconstruct_irreps(G,chi)
 M,fdims=fusion(G,irreps,chi); gi=graph_info(M,fdims)
 chis=sigchar(chi)
 if ceq(chi,chis,G):raise ValueError("Galois doublet is not distinct")
 if rint(cin(chis,chis,G),"Galois spin norm")!=1:raise ValueError("Galois doublet is reducible")
 Ms,ds=fusion(G,irreps,chis);gs=graph_info(Ms,ds)
 # Exact golden field generation from a trace witness with nonzero sqrt(5) coefficient.
 witnesses=[(g,v) for g,v in chi.items() if v[1]!=0]
 if not witnesses:raise ValueError("trace values do not generate Q(sqrt(5))")
 g,w=witnesses[0]
 generated_sqrt5=rdiv(rsub(w,(w[0],F(0))), (w[1],F(0)))
 if generated_sqrt5!=SQRT5:raise ValueError("trace witness failed to recover sqrt(5)")
 phi_minus=rsub((F(1),F(0)),PHI)
 if rgalois(PHI)!=phi_minus:raise ValueError("sigma(phi) != 1-phi")
 iso=isomorphisms(gi["edges"],AFFINE_E8_EDGES,9)
 isos=isomorphisms(gs["edges"],AFFINE_E8_EDGES,9)
 if raw["golden_character_field"]["witness_minimal_polynomial"] != quadratic_minpoly_text(w):
  raise ValueError("producer receipt has an unverified trace minimal polynomial")
 raw_checks={
  "group_order":raw["source_group"]["order"]==len(G),
  "group_center":raw["source_group"]["center_size"]==len(group_center(G)),
  "element_set_digest":True,
  "class_count":raw["conjugacy_classes"]["count"]==len(classes),
  "class_sizes":raw["conjugacy_classes"]["sizes"]==sizes,
  "irrep_dimensions":raw["irreducibles"]["dimensions"]==dims,
  "fusion_matrix":raw["fusion"]["matrix"]==M,
  "galois_fusion_matrix":raw["galois_conjugate_fusion"]["matrix"]==Ms,
  "trace_field_containment":raw["golden_character_field"]["all_traces_in_Qsqrt5"],
  "trace_field_generation":raw["golden_character_field"]["nonrational_trace_witness"]["trace"]==rtext(w),
 }
 if not all(raw_checks.values()):raise ValueError("producer receipt disagrees with independent reconstruction: "+str([k for k,v in raw_checks.items() if not v]))
 if not iso or not isos:raise ValueError("derived McKay graph failed explicit affine-E8 graph isomorphism")
 same_graph_type=True
 verdict=("MCKAY_DERIVES_GOLDEN_CHARACTER_FIELD","GALOIS_EMBEDDING_REMAINS_UNSELECTED")
 return {
  "schema":"arithmon.mckay_golden_field.verified.v1",
  "independent_verifier":{"construction":"generator closure from two exact H4 quaternions",
      "imports_producer":False,"group_order":len(G),"group_center_size":len(group_center(G)),
      "group_element_set_sha256":coordinate_set_digest(G),"conjugacy_class_sizes":sizes,
      "irreducible_dimensions":dims,"sum_squared_dimensions":sum(d*d for d in dims),
      "raw_receipt_checks":raw_checks},
  "fusion":{"matrix":M,**gi},
  "affine_e8":{"target":"trivalent tree with arm lengths 1,2,5",
      "explicit_isomorphism":list(iso[0]),"isomorphism_count":len(iso),"isomorphic":True},
  "golden_character_field":{"containment":"all traces represented in Q(sqrt(5))",
      "witness_quaternion":[rtext(x) for x in g],"witness_trace":rtext(w),
      "recover_sqrt5_from_witness":rtext(generated_sqrt5),
      "field_exactly":"Q(sqrt(5))","reason":"trace witness has nonzero sqrt(5) coefficient and Q(sqrt(5)) is quadratic"},
  "galois":{"automorphism":"sqrt(5)->-sqrt(5); i fixed","phi":rtext(PHI),"psi":rtext(phi_minus),
      "sigma_phi_equals_psi":True,"doublet_distinct":True,"doublet_faithful":True,
      "doublet_irreducible":True,"fusion_recomputed":True,"fusion_matrix":Ms,
      "same_labeled_graph":M==Ms,"same_affine_e8_graph_type":True,
      "conjugate_graph":{"edges":gs["edges"],"explicit_isomorphism":list(isos[0]),
                         "isomorphism_count":len(isos),"dimension_vector":ds}},
  "primary_verdict":list(verdict),
  "trust_boundary":{"oph_used_as_derivation_input":False,"koide_used":False,
      "experimental_data_used":False,"physical_identification":False,
      "sl2f5_explicit_isomorphism":False}}

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--input",type=Path,default=RAW_DEFAULT);p.add_argument("--output",type=Path,default=OUT_DEFAULT);p.add_argument("--no-write",action="store_true");a=p.parse_args()
 raw=json.loads(a.input.read_text());result=verify(raw);text=json.dumps(result,indent=2,sort_keys=True)+"\n"
 if not a.no_write:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
 print(text,end="")
if __name__=="__main__":main()
