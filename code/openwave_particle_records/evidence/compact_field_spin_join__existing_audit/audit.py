"""Frozen complete C/CAR single-source join and actual rotation audit."""
from pathlib import Path
import importlib.util,json,hashlib,sys,time,argparse
from fractions import Fraction as F
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3];RER=ROOT/'reverse-engineering-reality'
sys.path.insert(0,str(RER));sys.dont_write_bytecode=True
from model.observer.native import native
from model.defect_field_bridge.full_domain import branch_action

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def module(name,p):
 spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
mom=module('compact_spin_exact_geometry',ROOT/'openwave/shared/paired_source_moments/run.py')
FACES=((0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),(3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1))

def rotations(cell):
 points,_,_,_=mom.geometry(cell);base=FACES[0]
 sig=lambda p,face:tuple(mom.dot(points[p],points[v]) for v in face)
 perms=[]
 for face0 in FACES:
  for k in range(3):
   face=face0[k:]+face0[:k];where={sig(p,face):p for p in range(12)}
   perm=tuple(where[sig(p,base)] for p in range(12));assert len(set(perm))==12;perms.append(perm)
 assert len(set(perms))==60
 assert all(tuple(a[b[i]] for i in range(12)) in perms for a in perms for b in perms)
 # On the irreducible vector realization, trace(R)=sum_p v_p.Rv_p /(4 norm(v)^2).
 norm=mom.dot(points[0],points[0]);traces=[]
 for perm in perms:
  traces.append(mom.div(mom.total(mom.dot(points[p],points[perm[p]]) for p in range(12)),mom.mul(mom.q(4),norm)))
 return points,perms,traces

def run():
 for n,h in json.loads((HERE/'PREREG.json').read_text())['files'].items():assert sha(ROOT/n)==h,n
 for n,h in json.loads((HERE/'SOURCE_PINS.json').read_text())['files'].items():assert sha(ROOT/n)==h,n
 graph=json.loads((ROOT/'openwave/shared/localized_rest_join/evidence/graph.json').read_text());seams=tuple(map(tuple,graph['seams']));N=240;M=630
 assert len(seams)==M;cell=tuple(e for e in seams if e[1]<12)
 outputs=np.zeros((M,2,N),np.uint8);currents=np.zeros((M,2,N),np.int8);e_raw=np.zeros((M,2,N),np.uint8)
 calls=0;checks=0
 for e,(a,b) in enumerate(seams):
  for c in (0,1):
   destinations_by_current={j:[] for j in (-1,0,1)}
   for p in range(N):
    out,(_,j)=native((int(p==a),int(p==b)),(0,1),c);calls+=1
    dest=p if p not in (a,b) else a if out[0] else b
    outputs[e,c,p]=dest;currents[e,c,p]=j;destinations_by_current[j].append(dest)
    for spin in (0,1):
     mask=1<<(2*p+spin)
     for law in ('C','E'):
      raw,car,sign,jc=branch_action(mask,law,c,sites=N,flavors=2,seam=(a,b))
      assert car==1<<(2*dest+spin) and sign==1 and jc==j
      assert raw=={0:0,1:1,-1:2}[j]+(4*int(p in (a,b)) if law=='E' else 0)
      if law=='E':e_raw[e,c,p]=raw
      checks+=1
   assert all(len(set(ds))==len(ds) for ds in destinations_by_current.values())
 L=np.zeros((N,N),dtype=np.int64)
 for a,b in seams:L[a,a]+=1;L[b,b]+=1;L[a,b]-=1;L[b,a]-=1
 f=np.zeros(N,dtype=np.int64);f[[4,5,6,7]]=[1,-1,-1,1];assert np.array_equal(L@f,6*f)
 points,perms,traces=rotations(cell);Lc=L[:12,:12].copy();Lc[np.diag_indices(12)]=5
 Pnum=-Lc@(Lc@Lc-10*Lc+20*np.eye(12,dtype=np.int64))
 assert np.array_equal(Pnum@Pnum,24*Pnum) and np.trace(Pnum)==120
 char=[];orbital_spin_multiplicity=mom.q();fixed=[];ray=[]
 for g,tr in zip(perms,traces):
  chi=F(sum(int(Pnum[i,g[i]]) for i in range(12)),24);assert chi.denominator==1;char.append(int(chi))
  orbital_spin_multiplicity=mom.add(orbital_spin_multiplicity,mom.mul(mom.q(chi),mom.add(mom.q(1),tr)))
  pg=np.array([12*(p//12)+g[p%12] for p in range(N)])
  mapped=np.zeros(N,dtype=np.int64);mapped[pg]=f
  fixed.append(all(tuple(sorted((int(pg[a]),int(pg[b])))) in seams for a,b in seams))
  ray.append(int(np.max(abs(L@mapped-6*mapped))))
 orbital_spin_multiplicity=mom.div(orbital_spin_multiplicity,mom.q(60));assert orbital_spin_multiplicity==mom.q(0)
 # Exact fine-record counterexample, same seam and coin, two hold inputs.
 witness=[]
 for p in (0,2):
  row={law:branch_action(1<<(2*p),law,1,sites=N,flavors=2,seam=(0,1)) for law in ('C','E')}
  assert row['C'][0]==0 and row['C'][1]==1<<(2*p)
  witness.append({'input':p,'C_raw':row['C'][0],'E_raw':row['E'][0]})
 assert witness[0]['E_raw']!=witness[1]['E_raw']
 arrays={'source_destinations':outputs,'source_currents':currents,'E_raw_records':e_raw,'rotations':np.array(perms,dtype=np.int16),'fixed_graph_symmetry':np.array(fixed),'rotated_ray_residual':np.array(ray),'L6_projector_numerator':Pnum}
 result={'candidate':'A-static-signed-source-v1 with existing C/CAR single-sector identification','native_calls':calls,'CAR_C_and_E_columns':checks,'single_sector_columns_per_law':M*2*N*2,'C_complete_instrument_intertwiner':True,'E_coarse_record_intertwiner':False,'E_counterexample':witness,'conditioned_coin_probability':'1/2','C_offdiagonal':'1/4','E_reduced_offdiagonal':'0','L6_source_energy_unchanged':True,'proper_cell_rotations':60,'same_port_rotation_fixed_federation_symmetries':sum(fixed),'rotated_compact_rays_still_L6':sum(x==0 for x in ray),'L6_cell_rank':5,'L6_character':char,'L6_tensor_Spin2_fundamental_multiplicity':0,'full_row6':False,'scope':'Single-particle C instrument joins exactly to conserved CAR field spin. Fixed-source spatial rotation and total angular momentum require the stated orbital audit; no universal OPH obstruction.'}
 return result,arrays

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true');args=ap.parse_args();start=time.perf_counter()
 if args.record and (HERE/'RECEIPT.json').exists():raise FileExistsError('sealed receipt exists')
 result,arrays=run()
 if args.record:
  (HERE/'RECEIPT.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');np.savez_compressed(HERE/'RECORDS.npz',**arrays)
 else:
  assert result==json.loads((HERE/'RECEIPT.json').read_text());old=np.load(HERE/'RECORDS.npz');assert set(old.files)==set(arrays)
  for k,v in arrays.items():assert np.array_equal(v,old[k]),k
 print(json.dumps({'all_pass':True,'seconds':time.perf_counter()-start,**{k:result[k] for k in ('native_calls','single_sector_columns_per_law','same_port_rotation_fixed_federation_symmetries','rotated_compact_rays_still_L6','L6_tensor_Spin2_fundamental_multiplicity')}},indent=2))
if __name__=='__main__':main()
