"""Outgoing-volume refinement of the complete native meson current."""
from functools import lru_cache
from itertools import product
from fractions import Fraction as F
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply
import sympy as sy
from model.openwave_native_composites.spectrum import sector,one_body,norm
from model.native_flavour_transport.transport import source,word_unitary,numeric
from model.common_source import calibrated_matter as internal
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.openwave_weak_decay.decay import mass_data
from model.openwave_meson_weak_current.decay import BAR,ETA,G,_shifted_eigensystem

@lru_cache(None)
def hadron():
    native=sector('meson');subsets,_,_,_=representation();raising=numeric(internal.weak_raising_current(families=False))
    qi=[[subsets.index((c,w))for c in range(3)]for w in(3,4)]
    active=[i for i,(q,a)in enumerate(native['labels'])if a[2]==0]
    oldlabels=[native['labels'][i]for i in active];h=native['h'][np.ix_(active,active)]
    assert len(active)==32
    mlight,mheavy,_=mass_data();mass=np.array([1.,1/mlight,mheavy/mlight]);q=(mass[1]-2)/2-11*G/4
    ev,F=eigh(numeric(source()['y']));C=F.conj().T@numeric(word_unitary((2,3,10)))@F
    assert np.max(abs(C.conj().T@C-np.eye(3)))<2e-13
    hadron_labels=[(f,*qq[:2],qq[2],*aa[:2])for f in range(3)for qq,aa in oldlabels]
    Hh=sp.block_diag([sp.csr_matrix(h+(mf-1)*np.eye(32))for mf in mass],format='csr')
    full_diagonal=np.asarray(native['H'].diagonal()).real
    measured_hadron=np.asarray(abs(native['B'][:,active]).power(2).T@full_diagonal).ravel()
    measured_hadron=np.concatenate([measured_hadron+(mf-1)for mf in mass])
    qops=[]
    for x in range(2):
        row=[]
        for bar in BAR:
            aq=np.zeros((48,48),complex)
            for i,(species,cell,spin,w,c)in enumerate(native['modes']):
                if species or cell!=x:continue
                for j,(other,site,ss,ww,cc)in enumerate(native['modes']):
                    if other or site!=x:continue
                    aq[i,j]=bar[spin,ss]*raising[qi[w][c],qi[ww][cc]]
            full=one_body(aq,native['basis']);reduced=(native['B'].conj().T@full@native['B']).toarray()
            assert norm(full@native['B']-native['B']@sp.csr_matrix(reduced))<2e-13
            plus=reduced[np.ix_(active,active)]
            # The diagonal native current follows from the same full generator.
            a3=np.zeros((48,48),complex)
            T3=(raising@raising.conj().T-raising.conj().T@raising)/2
            for i,(species,cell,spin,w,c)in enumerate(native['modes']):
                if species or cell!=x:continue
                for j,(other,site,ss,ww,cc)in enumerate(native['modes']):
                    if other or site!=x:continue
                    a3[i,j]=bar[spin,ss]*T3[qi[w][c],qi[ww][cc]]
            zz=(native['B'].conj().T@one_body(a3,native['basis'])@native['B']).toarray()[np.ix_(active,active)]
            row.append((sp.csr_matrix(np.kron(C,plus)),sp.csr_matrix(np.kron(np.eye(3),zz))))
        qops.append(row)
    return dict(native=native,native_h=h,Hh=Hh,hadron_labels=hadron_labels,measured_hadron=measured_hadron,qops=qops,raising=raising,mass=mass,C=C,q_two=q)

def hadron_eigensystem(d,family,weak,*,driver=None):
 """Solve the complete native block and restore its known family energy.

 The full hadron matrix is checked against the original native block plus
 its scalar offset; no rounded subtraction or state projection is needed.
 """
 rows=[i for i,h in enumerate(d['hadron_labels'])if h[0]==family and h[3]==weak]
 native_rows=np.array(rows)%32
 H=d['Hh'].toarray()[np.ix_(rows,rows)];shift=d['mass'][family]-1
 e,U=_shifted_eigensystem(H,shift,native=d['native_h'][np.ix_(native_rows,native_rows)],driver=driver)
 return rows,e,U

def li(L,weak,x,spin,f):return weak*(6*L)+x*6+spin*3+f

def restricted(a,masks):
 lookup={mask:j for j,mask in enumerate(masks)};rows=[];cols=[];values=[]
 ai,aj=np.nonzero(a)
 for k,mask in enumerate(masks):
  for p,q in zip(ai,aj):
   p,q=int(p),int(q)
   if not(mask>>q)&1:continue
   before=mask^(1<<q)
   if(before>>p)&1:continue
   after=before|(1<<p);i=lookup.get(after)
   if i is None:continue
   sign=(-1)**((mask&((1<<q)-1)).bit_count()+(before&((1<<p)-1)).bit_count())
   rows.append(i);cols.append(k);values.append(sign*a[p,q])
 return sp.coo_matrix((values,(rows,cols)),shape=(len(masks),len(masks))).tocsr()

@lru_cache(maxsize=3)
def construction(L,recenter=False):
 assert L>=2 and L%2==0
 d=hadron();mass=d['mass'];q=d['q_two']if not recenter else (mass[1]-2)/2-G*(6-1/L)/2
 reference=(1<<(6*L))-1;masks=[reference];lepton=[None]
 for f,xe,se,xn,sn in product(range(3),range(L),range(2),range(L),range(2)):
  masks.append(reference^(1<<li(L,0,xn,sn,f))|(1<<li(L,1,xe,se,f)));lepton.append((f,xe,se,xn,sn))
 lap=2*np.eye(L)
 for x in range(L):lap[x,(x+1)%L]-=1;lap[x,(x-1)%L]-=1
 lep=np.zeros((12*L,12*L),complex)
 for weak,spin,f in product(range(2),range(2),range(3)):
  ids=[li(L,weak,x,spin,f)for x in range(L)];K=q*lap/4;lep[np.ix_(ids,ids)]=mass[f]*np.eye(L)+K if weak else -K
 nlep=len(masks);Hl=restricted(lep,masks)+sp.eye(nlep)*(3*q*L)
 assert abs(Hl[0,0])<1e-10
 j=d['raising'][np.ix_([10,11],[10,11])];t3=(j@j.conj().T-j.conj().T@j)/2
 lops=[]
 for x in range(2):
  one=[]
  for bar in BAR:
   pair=[]
   for weakmatrix in(j,t3):
    A=np.zeros((12*L,12*L),complex)
    for w,ww,s,ss,f in product(range(2),range(2),range(2),range(2),range(3)):
     A[li(L,w,x,s,f),li(L,ww,x,ss,f)]=bar[s,ss]*weakmatrix[w,ww]
    pair.append(restricted(A,masks))
   one.append(pair)
  lops.append(one)
 V=sp.csr_matrix((96*nlep,96*nlep),dtype=complex)
 for x,mu in product(range(2),range(4)):
  qp,qz=d['qops'][x][mu];lp,lz=lops[x][mu]
  V+=ETA[mu]*(sp.kron(qp,lp.conj().T)+sp.kron(qp.conj().T,lp)+2*sp.kron(qz,lz))
 H0=sp.kron(d['Hh'],sp.eye(nlep))+sp.kron(sp.eye(96),Hl);keep=[]
 for hi,(f,xq,sq,w,xa,sa)in enumerate(d['hadron_labels']):
  for k,ell in enumerate(lepton):
   Ne=ell is not None;sz=(1-2*sq)+(1-2*sa)
   if Ne:sz+=(1-2*ell[2])-(1-2*ell[4])
   if -w-Ne==-1 and sz==0:keep.append(hi*nlep+k)
 keep=np.array(keep);other=np.setdiff1d(np.arange(96*nlep),keep)
 assert len(keep)==24+216*L*L
 assert norm(H0.tocsr()[keep][:,other])<3e-12 and norm(V[keep][:,other])<3e-12
 H0=H0.tocsr()[keep][:,keep];V=V[keep][:,keep];H=H0+G*V
 where={int(v):i for i,v in enumerate(keep)};initial=np.zeros(len(keep),complex)
 parentrows=[i for i,h in enumerate(d['hadron_labels'])if h[0]==1 and h[3]==1]
 # Preserve the orbital seed rule independently of projector roundoff.
 _,seed_frame=eigh(d['Hh'].toarray()[np.ix_(parentrows,parentrows)])
 seed_hadron=parentrows[int(np.argmax(abs(seed_frame[:,0])**2))]
 parentrows,e,U=hadron_eigensystem(d,1,1)
 for k,z in zip(parentrows,U[:,0]):
  if k*nlep in where:initial[where[k*nlep]]=z
  else:assert abs(z)<1e-12
 finalrows,ef,Fh=hadron_eigensystem(d,0,0);fh=np.zeros(96,complex);fh[finalrows]=Fh[:,0];fl=np.zeros(nlep,complex)
 for k,label in enumerate(lepton):
  if label is None:continue
  f,xe,se,xn,sn=label
  if f==0 and se==sn:fl[k]=(-1)**(xe+xn+li(L,0,xn,sn,f)+6*L-1)/(L*np.sqrt(2))
 target=np.kron(fh,fl)[keep];assert abs(np.vdot(target,target)-1)<1e-12
 seed=where[seed_hadron*nlep];success=float(abs(initial[seed])**2)
 measured=np.repeat(d['measured_hadron'],nlep)+np.tile(Hl.diagonal().real,96)
 out=dict(d,L=L,recenter=recenter,q=q,nlep=nlep,keep=keep,lepton_labels=lepton,H=H,H0=H0,V=V,initial=initial,target=target,success=success,measured_energy=measured[keep]+G*V.diagonal().real,parent_energy=float(e[0]))
 return out

@lru_cache(None)
def singlets():
 spin=[sy.Matrix([[0,1],[1,0]])/2,sy.Matrix([[0,-sy.I],[sy.I,0]])/2,sy.diag(1,-1)/2];total=[]
 for axis in range(3):
  generators=[]
  for factor in range(4):
   factors=[sy.eye(2)for _ in range(4)];factors[factor]=-spin[axis].T if factor==3 else spin[axis];generators.append(sy.kronecker_product(*factors))
  total.append(sum(generators,sy.zeros(16)))
 S2=sum((a*a for a in total),sy.zeros(16));basis=sy.GramSchmidt(S2.nullspace(),orthonormal=True)
 return np.array([numeric(x).ravel()for x in basis])

def reduction(d):
 L=d['L'];nlep=d['nlep'];s=singlets();where={int(v):i for i,v in enumerate(d['keep'])};hd={x:i for i,x in enumerate(d['hadron_labels'])};ld={x:i for i,x in enumerate(d['lepton_labels'])};keys=[];columns=[]
 for f,xq,xa in product(range(3),range(2),range(2)):
  col={where[hd[(f,xq,sq,1,xa,sa)]*nlep]:(1 if sq==0 else -1)/np.sqrt(2)for sq,sa in((0,1),(1,0))};columns.append(col);keys.append(('parent',f,xq,xa))
 for f,ell,xq,xa,xe,xn,k in product(range(3),range(3),range(2),range(2),range(L),range(L),range(2)):
  col={}
  for sq,sa,se,sn in product(range(2),repeat=4):
   value=s[k,((2*sq+sa)*2+se)*2+sn]
   if value==0:continue
   h=hd[(f,xq,sq,0,xa,sa)];l=ld[(ell,xe,se,xn,sn)];col[where[h*nlep+l]]=(-1)**(li(L,0,xn,sn,ell)+6*L-1)*value
  columns.append(col);keys.append(('emission',f,ell,xq,xa,xe,xn,k))
 lookup={k:i for i,k in enumerate(keys)};rows=[];cols=[];data=[];count=0
 for i,key in enumerate(keys):
  other=('parent',key[1],1-key[2],1-key[3])if key[0]=='parent'else('emission',key[1],key[2],1-key[3],1-key[4],(1-key[5])%L,(1-key[6])%L,key[7])
  j=lookup[other]
  if i>=j:continue
  for q in(i,j):
   for row,v in columns[q].items():rows.append(row);cols.append(count);data.append(v/np.sqrt(2))
  count+=1
 B=sp.coo_matrix((data,(rows,cols)),shape=(len(d['keep']),count)).tocsr();h=(B.conj().T@d['H']@B).tocsr();error=norm(d['H']@B-B@h)
 assert count==6+36*L*L and norm(B.conj().T@B-sp.eye(count))<3e-12 and error<2e-9
 assert np.linalg.norm(d['initial']-B@(B.conj().T@d['initial']))<2e-12
 return B,h,error

def evolve(d,times):
 B,h,error=reduction(d);e,U=eigh(h.toarray());start=B.conj().T@d['initial'];frame=U.conj().T@start;rows=[]
 emitted=d['keep']%d['nlep']!=0
 for t in times:
  v=B@(U@(np.exp(-1j*t*e)*frame));p=abs(v)**2
  rows.append({'time':t,'conditional_pair_probability':float(p[emitted].sum()),'original_input_pair_probability':float(d['success']*p[emitted].sum()),'bound_pi_pi_mode_probability':float(abs(np.vdot(d['target'],v))**2),'full_energy':float(np.vdot(v,d['H']@v).real),'position_state':v if t==720 else None})
 return rows,(B,h,error)

def finite_volume_bound(L,time):
 R=(L-2)//2;x=F(103)*F(time);assert R>0 and x<R+1
 bound=2*x**R/F(__import__('math').factorial(R))/(1-x/F(R+1))
 return {'ring_sites':L,'source_time':str(time),'first_omitted_hop_order':R,'hopping_norm_upper':'103','prepared_domain_vector_error_upper_exact':str(bound),'prepared_domain_vector_error_upper':float(bound),'complete_probability_error_upper':float(2*bound),'floating_evolution_roundoff_included':False,'probability_transfer_domain':'common embedded/local instruments with the outside-line complement; no transfer from ring FFT to infinite-line momentum spectral instrument'}

def memory_kernel(d,B,h):
 # The six parent coordinates precede every emitted coordinate.
 A=h[:6,:6].toarray();C=h[6:,:6].toarray()/G;D=h[6:,6:].toarray();start=np.asarray(B.conj().T@d['initial']).ravel()[:6]
 coefficient=float(G*G*np.linalg.norm(C@start)**2)
 assert np.linalg.norm(C.conj().T@C-24*np.eye(6))<2e-11
 z=complex(d['parent_energy'],1);resolvent=np.linalg.inv(z*np.eye(len(h.toarray()))-h.toarray())[:6,:6]
 schur=np.linalg.inv(z*np.eye(6)-A-G*G*C.conj().T@np.linalg.solve(z*np.eye(len(D))-D,C))
 err=float(np.linalg.norm(resolvent-schur));assert err<2e-9
 return {'parent_dimension':6,'complete_emitted_dimension':len(D),'initial_quadratic_pair_probability_coefficient':coefficient,'memory_at_zero_trace':float(np.trace(C.conj().T@C).real),'exact_parent_Volterra_equation':'a_dot=-i A a-G^2 integral_0^t C^dagger exp(-iD(t-s)) C a(s) ds','reconstructed_emitted_state':'b(t)=-iG integral_0^t exp(-iD(t-s)) C a(s) ds','Schur_resolvent_numeric_residual':err,'D_includes_full_native_neutral_current':True,'irreversible_rate_claimed':False}

def fine_records(d,v):
 L=d['L'];nlep=d['nlep'];prob=abs(v)**2;success=d['success'];before=float(np.vdot(v,d['H']@v).real)
 # The native colour singlet splits each carrier occupation probability into
 # three equal matching-colour records. Other colour records have probability0.
 position_after=float(prob@d['measured_energy']);momentum=v.copy();energies=d['measured_energy'].copy();groups={}
 for index,original in enumerate(d['keep']):
  hi=int(original)//nlep;ell=d['lepton_labels'][int(original)%nlep]
  if ell is not None:
   f,xe,se,xn,sn=ell;groups.setdefault((hi,f,se,sn),[]).append((index,xe,xn))
 kinetic=d['q']/2*(1-np.cos(2*np.pi*np.arange(L)/L));vdiag=d['V'].diagonal().real
 for(hi,f,se,sn),indices in groups.items():
  assert len(indices)==L*L;matrix=np.zeros((L,L),complex)
  for index,xe,xn in indices:matrix[xe,xn]=v[index]*(-1)**(li(L,0,xn,sn,f)+6*L-1)
  transformed=np.fft.fft2(matrix,norm='ortho');neutral=float(np.mean(vdiag[[i for i,_,_ in indices]]))
  for index,ke,kn in indices:
   momentum[index]=transformed[ke,kn];energies[index]=d['measured_hadron'][hi]+d['mass'][f]+kinetic[ke]+kinetic[kn]+G*neutral
 mprob=abs(momentum)**2;assert abs(mprob.sum()-1)<1e-11
 momentum_after=float(mprob@energies)
 return {'record_domain':'one preparation failure quantum branch, then all native hadron occupations/three matching colours and all lepton occupation records; position and momentum are separate instruments','carrier_occupation_indices':d['keep'].tolist(),'conditional_position_probabilities':prob.tolist(),'conditional_momentum_probabilities':mprob.tolist(),'each_matching_colour_probability_factor':'1/3','unconditional_success_factor':success,'preparation_failure_probability':1-success,'all_other_source_records_probability':'0 on this preparation, retained by the orthogonal instrument complement','conditional_full_source_energy_before':before,'position_full_source_energy_after':position_after,'momentum_full_source_energy_after':momentum_after,'unconditional_position_readout_work':success*(position_after-before),'unconditional_momentum_readout_work':success*(momentum_after-before),'preparation_nonselective_mean_work':'0','no_postselection_renormalization':True,'both_complete_original_probability_sums':[float(1-success+success*prob.sum()),float(1-success+success*mprob.sum())]}

def certificate():
 rows=[]
 for L in(2,4,6,8):
  d=construction(L);series,(B,h,error)=evolve(d,(0,32,128,512,720,1024));state=next(x['position_state']for x in series if x['time']==720)
  for x in series:x.pop('position_state')
  coupling=complex(np.vdot(d['target'],d['V']@d['initial']));parent_v=float(np.vdot(d['initial'],d['V']@d['initial']).real);target_v=float(np.vdot(d['target'],d['V']@d['target']).real)
  assert abs(parent_v+3)<1e-10 and abs(target_v-(3-1/L))<1e-10
  expected=2*abs(d['C'][0,1])/L/np.sqrt(2);assert abs(abs(coupling)-expected)<2e-12
  initial_coefficient=float(np.linalg.norm((h[6:,:6]@(B.conj().T@d['initial'])[:6]))**2)
  rows.append({'ring_sites':L,'full_CAR_coordinates':180+12*(L-2),'charge_and_spin_sector_dimension':len(d['keep']),'exact_spin_zero_even_parity_dimension':h.shape[0],'source_box_q':d['q'],'prepared_parent_success_probability':d['success'],'full_source_reduction_numeric_residual':error,'native_bound_pi_pi_formfactor_absolute':abs(coupling),'exact_formfactor_formula':'2 abs(C_01)/(sqrt(2) L)','parent_contact_mean':parent_v,'bound_pi_pi_contact_mean':target_v,'fixed_box_contact_detuning':G*(.5-1/L),'declared_mean_recenter_box_control':(d['mass'][1]-2)/2-G*(6-1/L)/2,'initial_quadratic_pair_probability_coefficient':initial_coefficient,'complete_coarse_time_records':series,'complete_fine_position_and_momentum_records_at720':fine_records(d,state),'exact_memory_control':memory_kernel(d,B,h)if L in(2,4)else None})
 coeffs=[r['initial_quadratic_pair_probability_coefficient']for r in rows];assert max(coeffs)-min(coeffs)<3e-15
 # This sparse execution retains all36870 symmetry-reduced source states.
 d=construction(32);B,h,error=reduction(d);t=F(1,100);v=B@expm_multiply(-1j*float(t)*h,B.conj().T@d['initial']);emitted=d['keep']%d['nlep']!=0;p=float(np.linalg.norm(v[emitted])**2);assert abs(np.linalg.norm(v)-1)<2e-10
 large={'ring_sites':32,'full_CAR_coordinates':540,'full_charge_spin_sector_dimension':len(d['keep']),'exact_spin_zero_even_parity_dimension':h.shape[0],'source_time':str(t),'conditional_pair_probability_floating':p,'original_input_pair_probability_floating':d['success']*p,'complete_original_coarse_records':{'preparation_failure':1-d['success'],'prepared_no_pair':d['success']*(1-p),'prepared_pair':d['success']*p},'full_source_reduction_numeric_residual':error,'floating_propagation_norm_error':abs(float(np.linalg.norm(v))-1),'finite_volume_bound':finite_volume_bound(32,t),'floating_propagation_roundoff_certified':False}
 return {'schema':'oph.native_meson_outgoing_ring.v1','verdict':'PASS_COMPLETE_NATIVE_RING_RECORDS_AND_FINITE_VOLUME_BOUND','declared_source':'same two-site native bound meson, native semileptonic current and original180 CAR coordinates; add outgoing lepton ring coordinates, keeping every spectator and failed source preparation','coupling':str(F(1,64)),'fixed_source_box':'L2 resonance-selected box, held fixed through ring refinement','initial_quadratic_emission_coefficient_volume_independent':coeffs[0],'initial_quadratic_emission_coefficient_exact':'3/512','native_injection_Gram_exact':'C^dagger C=24 I_6','charge_reason_for_at_most_one_pair':'Q=-w-N_e=-1 and w in{0,1}; all charge-compatible emitted channels are retained','finite_rings':rows,'large_ring_short_time_execution':large,'finite_volume_theorem_examples':[finite_volume_bound(32,F(1,100)),finite_volume_bound(1024,F(1))],'exact_memory_kernel_parent_dimension':6,'infinite_volume_construction':'bounded same-source operator on finite hadron tensor zero/one outgoing e-hole pair over Z; local current at contacts0,1 and bounded nearest-neighbour hopping; filled neutral reference defines the hole Fock representation','observed_small_box_effect':'at source time720 the pair probability falls with L=2,4,6,8, so the L2 peak supplies no volume-stable irreversible rate','physical_strange_decay_lifetime_or_observed_meson_identity_claimed':False,'scope':'Exact native-current finite-volume refinement and complete observer records. The rigorous Dyson tail bounds source volume error, separately from floating propagation roundoff. Finite rings retain recurrences. No physical lifetime or infinite-time exponential is inferred.'}
