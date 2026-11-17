"""Direct full625-register control of the H and H² block identities at Nq=1.

This control acts on all15,000 states; it is not the Nq=9 spectrum or a Gauss
projection. The generator identities used for the latter hold on this larger
unprojected control too. No matrix from the block formulas enters H_action.
"""
from pathlib import Path
import sys,json
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from local_auxiliary_interface import generators

def main():
    T,S=generators();I=np.eye(12);z=np.zeros((12,12));K=np.block([[z,-I/8],[-I/8,z]]);m=.000060300795404423236;M=np.diag([m]*12+[-m]*12);quarkG=[[],[]];quarkC=[[],[]]
    f=np.array([[[-2j*np.trace((T[a]@T[b]-T[b]@T[a])@T[c]) for c in range(8)] for b in range(8)] for a in range(8)]).real
    auxC=[]
    for a in range(8):
        c=np.zeros((25,25),complex);c[1:,1:]=np.kron(-1j*f[a],np.eye(3));auxC.append(c)
    star=[]
    for ai in range(24):
        x=np.zeros((25,25));x[0,ai+1]=x[ai+1,0]=1;star.append(x)
    for cell in range(2):
        for a in range(8):
            c=np.zeros((24,24),complex);c[12*cell:12*(cell+1),12*cell:12*(cell+1)]=np.kron(np.eye(4),T[a]);quarkC[cell].append(c)
            for i in range(3):
                g=np.zeros((24,24),complex);g[12*cell:12*(cell+1),12*cell:12*(cell+1)]=np.kron(np.kron(S[i],np.eye(2)),T[a]);quarkG[cell].append(g)
    number=np.diag([0]+[1]*24)
    def onq(h,v):return np.einsum('ab,xybk->xyak',h,v,optimize=True)
    def onl(h,v):return np.einsum('ab,byqk->ayqk',h,v,optimize=True)
    def onr(h,v):return np.einsum('ab,xbqk->xaqk',h,v,optimize=True)
    def H_action(v,kappa):
        out=onq(M+K,v)+3*kappa*(onl(number,v)+onr(number,v))
        for a in range(8):
            q=quarkC[0][a];c=auxC[a];out+=kappa*(onq(q@q,v)+2*onq(q,onl(c,v))+onl(c@c,v))
        for cell,fn in [(0,onl),(1,onr)]:
            for ai,g in enumerate(quarkG[cell]):out+=onq(g,fn(star[ai],v))/8
        return out
    P=np.zeros((25,25,24,24),complex);P[0,0]=np.eye(24);A=np.zeros_like(P)
    for ai in range(24):A[ai+1,0]=quarkG[0][ai];A[0,ai+1]=quarkG[1][ai]
    frame=np.concatenate([P,A],axis=3);flat=frame.reshape(15000,48);D=sum(g@g for cell in quarkG for g in cell);Dl=sum(g@g for g in quarkG[0]);Dr=sum(g@g for g in quarkG[1]);C=[K@g-g@K for cell in quarkG for g in cell];TC=sum(c.conj().T@c for c in C);results=[]
    for kappa in [.5,1,2]:
        H=M+K+kappa*sum(q@q for q in quarkC[0]);delta=3*kappa;Z=(H@D-D@H)/2-K;F=delta*D+(H@D+D@H)/2-K;J=H+delta*np.eye(24);W=J@D@J+J@Z+Z.conj().T@J+TC
        Hb=np.block([[H,D/8],[D/8,F]]);H2=np.block([[H@H+D/64,(H@D+F)/8],[(D@H+F)/8,D@D/64+W+Dl@Dr/16]])
        acted=H_action(frame,kappa).reshape(15000,48);resH=float(np.linalg.norm(flat.conj().T@acted-Hb));resH2=float(np.linalg.norm(acted.conj().T@acted-H2));assert max(resH,resH2)<1e-10
        results.append(dict(kappa=kappa,full_H_block_residual=resH,full_H_squared_block_residual=resH2))
        print('DIRECT FULL625 CONTROL',kappa,resH,resH2,flush=True)
    (HERE/'FULL_CARRIER_CONTROL.json').write_text(json.dumps(dict(schema='oph.auxiliary_current_block.full_carrier_control.v1',quark_number=1,full_carrier_dimension=15000,all625_auxiliary_states_retained=True,no_Gauss_postselection=True,baryon_spectrum_claim=False,results=results),indent=2)+'\n')
if __name__=='__main__':main()
