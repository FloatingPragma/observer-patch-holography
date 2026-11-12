"""Complete finite-block coloured fragment control, with an explicit copy convention."""
from itertools import combinations
from functools import lru_cache
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from model.openwave_calibrated_hadrons.source import native_colour48,angular4,integral
from model.openwave_native_composites.spectrum import one_body
from .source import sector,parameters,colour_global48,weights

@lru_cache(None)
def fragment_controls():
 p=parameters();modes=sector('baryon')['modes'];rows=[];mins={}
 for number in(1,2):
  basis=tuple(sum(1<<i for i in c)for c in combinations(range(24),number));n=len(basis);index={v:i for i,v in enumerate(modes)}
  k=np.zeros((24,24),dtype=int)
  for i,(species,x,s,f,c)in enumerate(modes):k[i,i]=1;k[i,index[species,1-x,s,f,c]]=-1
  n0=np.array([sum((b>>i)&1 for i,m in enumerate(modes)if m[1]==0)for b in basis]);nu=np.array([sum((b>>i)&1 for i,m in enumerate(modes)if m[3]==0)for b in basis])
  parts={'N':number*sp.eye(n,format='csr'),'D':sp.diags(number-2*nu,dtype=float),'K':integral(one_body(k,basis)),'C48':native_colour48(modes,basis),'Cg48':colour_global48(modes,basis),'S4':angular4(modes,basis,(0,),True)+angular4(modes,basis,(1,),True),'I4':angular4(modes,basis,(0,),False)+angular4(modes,basis,(1,),False),'Ig4':angular4(modes,basis,(0,1),False),'B9':sp.diags(n0*n0+(number-n0)**2,dtype=float)}
  for kappa in(0,1):
   H=sum(float(v)*parts[k]for k,v in weights({**p,'kappa_colour':kappa}).items()).toarray()
   for u in range(number+1):
    ids=np.flatnonzero(nu==u);e=eigh(H[np.ix_(ids,ids)],eigvals_only=True);mins[number,u,kappa]=float(e[0]);rows.append({'number':number,'up_number':u,'edge_colour_coefficient':kappa,'full_sector_dimension':len(ids),'all_eigenvalues':e.tolist()})
 baryon=sector('baryon');q=np.rint(np.diag(baryon['q']).real).astype(int);comparisons=[]
 for kappa in(0,1):
  for u in(1,2):
   e=float(eigh(baryon['h'][np.ix_(q==u-1,q==u-1)],eigvals_only=True)[0]);threshold=min(mins[1,a,kappa]+mins[2,u-a,kappa]for a in(0,1)if 0<=u-a<=2)
   comparisons.append({'baryon_up_number':u,'edge_colour_coefficient':kappa,'baryon_energy':e,'two_copy_fragment_threshold':threshold,'threshold_minus_baryon':threshold-e})
 return {'complete_fragment_spectra':rows,'copy_comparisons':comparisons,'arithmetic':'floating diagnostics; exact colour-Casimir lower shift per nontrivial 1Q/2Q copy is 4/3',
  'comparison_convention':'Each separated fragment occupies an independent complete two-cell copy, including that copy\'s edge-colour square. This is a comparison between disjoint copies, not deleting the edge of the original carrier.',
  'same_edge_separated_singlet_control':'The edge-colour square is zero on every global singlet in the original two-cell carrier, including a spatially separated quark and diquark. No two-cell spatial-binding repair follows.',
  'cut_convention':'Deleting the entire edge square removes its onsite and cross pieces together. Retaining onsite pieces while deleting its cross term defines a different cut action.'}
