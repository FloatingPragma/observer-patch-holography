"""One fixed spin-gradient strong action and its full finite contact source."""
from functools import lru_cache
from fractions import Fraction as F
import numpy as np
import scipy.sparse as sp
from model.openwave_spin_gradient_nuclei import source as gradient
from model.openwave_dressed_baryon_lines.source import common,parity_basis,exact_parity_matrix,G
from model.openwave_baryon_weak_formfactors.native import quartet
@lru_cache(None)
def construction(charge=0,coupling=G):
 provider=gradient;native=provider.sector('baryon');c=common();assert c['q']==(gradient.parameters()['delta']-1)/2;ne=c['ne'];qh=np.rint(np.diag(native['q']).real).astype(int);zh=np.rint(2*np.diag(native['z']).real).astype(int)
 charges=(qh[:,None]-ne).ravel();spins=(zh[:,None]+c['lz']).ravel();keep=np.flatnonzero((charges==charge)&(spins==1));outside=np.flatnonzero((charges!=charge)|(spins!=1))
 full0=sp.kron(sp.csr_matrix(native['h']),sp.eye(70))+sp.kron(sp.eye(120),c['Hl']);H0=full0.tocsr()[keep][:,keep];V=c['V'][keep][:,keep];H=H0+float(coupling)*V
 assert np.max(abs(full0.tocsr()[keep][:,outside].data),initial=0)<3e-10 and np.max(abs(c['V'][keep][:,outside].data),initial=0)<3e-12
 parity=c['parity'][keep][:,keep];assert np.max(abs((H@parity-parity@H).data),initial=0)<3e-10
 W,_=quartet(native['h']);input_full=np.kron(W[:,0 if charge==0 else 2],np.eye(70)[:,c['reference']]);initial=input_full[keep]
 assert np.linalg.norm(parity@initial-initial)<1e-10 and abs(np.linalg.norm(initial)-1)<1e-10
 return dict(c,member="spin_gradient",provider=provider,native=native,coupling=F(coupling),charge=charge,keep=keep,H0=H0,V=V,H=H,parity=parity,initial=initial,charge_dimension=int(np.sum(charges==charge)),spin_dimension=len(keep))
