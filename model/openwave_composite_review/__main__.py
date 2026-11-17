import argparse,hashlib,json
from pathlib import Path
from itertools import combinations,product
import numpy as np
import scipy.sparse as sp
from model.openwave_native_composites.spectrum import native,one_body,selected,MASS,T,U,J
from model.openwave_native_composites.__main__ import pins as composite_pins
from model.openwave_composite_scale.__main__ import pins as scale_pins
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent

def fragment(number,conjugate=False):
 """Independently assemble the entire conserved native fragment sector."""
 q,a,_,_,_=native();g=a if conjugate else q
 modes=list(product(range(2),range(2),range(2),range(3)));index={v:i for i,v in enumerate(modes)}
 basis=[sum(1<<i for i in word)for word in combinations(range(24),number)]
 kinetic=np.zeros((24,24),complex)
 for i,(x,s,f,c)in enumerate(modes):kinetic[i,i]=T;kinetic[i,index[(1-x,s,f,c)]]=-T
 H=MASS*number*sp.eye(len(basis))+one_body(kinetic,basis)
 pauli=[np.array([[0,1],[1,0]],complex)/2,np.array([[0,-1j],[1j,0]])/2,np.diag([.5,-.5])]
 for x in range(2):
  for axis in range(8):
   current=np.zeros((24,24),complex)
   for s,f,c,d in product(range(2),range(2),range(3),range(3)):current[index[(x,s,f,c)],index[(x,s,f,d)]]=g[f][axis][c,d]
   many=one_body(current,basis);H+=U*(many@many)
  for sigma in pauli:
   current=np.zeros((24,24),complex)
   for s,z,f,c in product(range(2),range(2),range(2),range(3)):current[index[(x,s,f,c)],index[(x,z,f,c)]]=sigma[s,z]
   many=one_body(current,basis);H+=J*(many@many)
 energies=np.linalg.eigvalsh(H.toarray())
 return {'native_occupation':number,'conjugate':conjugate,'complete_dimension':len(basis),'ground_energy':float(energies[0]),'all_eigenvalues':list(map(float,energies))}

def certificate():
 q=fragment(1);a=fragment(1,True);pair=fragment(2)
 one_bound=MASS+4*U/3+3*J/4;pair_bound=2*MASS+4*U/3
 assert abs(q['ground_energy']-one_bound)<1e-12 and abs(a['ground_energy']-one_bound)<1e-12
 assert pair['ground_energy']>=pair_bound-1e-12
 B=selected('baryon',1,.5);M=selected('meson',1,0)
 baryon_gap=pair['ground_energy']+q['ground_energy']-B['ground_energy'];meson_gap=q['ground_energy']+a['ground_energy']-M['ground_energy']
 assert baryon_gap>=8*U/3-3*T and meson_gap>=8*U/3+3*J/2-2*T
 source={};receipts={}
 for package,pin_function in(('openwave_native_composites',composite_pins),('openwave_composite_scale',scale_pins)):
  path=ROOT/'model'/package/'results.json';data=json.loads(path.read_text());assert data['source_sha256']==pin_function()
  receipts[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest();source.update(data['source_sha256'])
  if package=='openwave_native_composites':
   for row in data['states']:
    r=row['complete_unconditional_records'];assert abs(sum(r['preparation_success_and_every_full_sector_occupation'])+r['preparation_failure']-1)<1e-9
  else:
   assert data['physical_particle_rest_mass_claim']is False
   assert all(float(x['ratio_to_comparison'])<.5 for x in data['physical_assignment_tests'])
 ownhash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 return {'verdict':'PASS_INDEPENDENT_NATIVE_COMPOSITE_SECTOR_SCALE_REVIEW',
  'source_sha256':{**source,**receipts,'model/openwave_composite_review/__main__.py':ownhash,
   'model/openwave_composite_review/__init__.py':hashlib.sha256((HERE/'__init__.py').read_bytes()).hexdigest()},
  'reviewed_receipt_sha256':receipts,'reviewed_source_sha256':source,
  'complete_fragment_spectra':[q,a,pair],'two_plus_one_finite_fragment_gap':baryon_gap,'one_plus_one_finite_fragment_gap':meson_gap,
  'fragment_meaning':'noninteracting copies of the same finite source in exact conserved native occupation sectors; no asymptotic scattering claim',
  'conditional_and_unconditional_records_separately_retained':True,
  'all_native_complement_sectors_controlled_by_positive_full_Fock_bound':True,
  'exact_weak_symmetry_and_family_number_offset_obstructions_retained':True,
  'physical_hadron_identification_claimed':False,
  'review_scope':'native Casimir normalization and intertwiners; singlet multiplicities120/64; invariance; positive full Fock floor; all fragment states; full preparation failures; common scale ancestry and rejected physical assignments',
  'review_code_sha256':ownhash}
def equivalent(a,b):
 if isinstance(b,dict):return isinstance(a,dict)and set(a)==set(b)and all(equivalent(a[k],v)for k,v in b.items())
 if isinstance(b,list):return isinstance(a,list)and len(a)==len(b)and all(equivalent(x,y)for x,y in zip(a,b))
 if isinstance(b,float):return isinstance(a,(float,int))and abs(a-b)<=max(2e-12,abs(b)*2e-10)
 return type(a)is type(b)and a==b
def main():
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=p.parse_args();r=certificate();path=HERE/'results.json'
 if args.write:path.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
 else:
  assert equivalent(json.loads(path.read_text()),r)
 print(r['verdict'])
if __name__=='__main__':main()
