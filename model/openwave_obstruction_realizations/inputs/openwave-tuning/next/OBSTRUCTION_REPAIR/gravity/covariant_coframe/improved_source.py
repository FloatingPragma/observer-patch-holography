"""Standard fourth-order differences; a new flat source, no entropy fit."""
import numpy as np
from audit import ALPHA,BETA,dagger

class ImprovedSource:
    def __init__(self,L,R):
        self.L=L
        n=np.arange(L)
        self.grid=np.stack(np.meshgrid(n,n,n,indexing='ij'),-1).reshape(-1,3)
        k=2*np.pi*(self.grid+.5)/L
        t=np.sin(k/2)**2
        d=np.sin(k)*(1+2*t/3)
        self.kin=sum(d[:,j,None,None]*ALPHA[j] for j in range(3))
        self.lap=np.sum(4*t*(1+t/3),axis=1)
        W=self.lap-np.sum(d*d,axis=1)
        exactW=np.sum((32/9)*t**3+(16/9)*t**4,axis=1)
        self.polynomial_identity_error=float(np.max(abs(W-exactW)))
        self.minimum_W=float(W.min())
        H=self.kin+R**3*self.lap[:,None,None]**2*BETA
        self.E,self.U=np.linalg.eigh(H)
        del H
        self.occ=self.E<0
        P=(self.U*self.occ[:,None,:])@dagger(self.U)
        self.kernel=np.fft.ifftn(P.reshape(L,L,L,4,4),axes=(0,1,2))
        del P

    def shifted(self,shift):
        idx=(self.grid+np.array(shift))%self.L
        return idx[:,0]*self.L*self.L+idx[:,1]*self.L+idx[:,2]
