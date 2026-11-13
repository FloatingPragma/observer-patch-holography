"""Actual source-current transfer to the positive spin-gradient strong action."""
from fractions import Fraction as F
from pathlib import Path
import json
import numpy as np
from model.openwave_spin_gradient_nuclei.source import sector,parameters
from model.openwave_baryon_weak_formfactors.native import quartet,currents,PAULI
from .angular import source,records,explicit_spin_weight,native_spin_weight

def certificate():
    native=sector('baryon');h=native['h'];W,_=quartet(h);n=W[:,0];p=W[:,2:];d=source().copy()
    P={}
    for kind in ('plus','minus','third'):
        P[kind]=np.zeros((2,4,2),complex)
        for x in range(2):
            for mu in range(4):
                qp,qz=currents()['local'][x][mu];J=qp if kind=='plus' else qp.conj().T if kind=='minus' else qz
                P[kind][x,mu]=p.conj().T@J@n
    g=float(sum(P['plus'][x,3,0] for x in range(2)).real)
    enclosing=json.loads((Path(__file__).resolve().parents[1]/'openwave_spin_gradient_nuclei/results.json').read_text())['native_axial_ground']['gA_interval']
    lo,hi=map(F,enclosing);assert F('1.2752')<lo<hi<F('1.2754')
    for x in range(2):
        assert np.linalg.norm(P['plus'][x,0]-np.array([.5,0]))<2e-10
        for j in range(1,4):assert np.linalg.norm(P['plus'][x,j]-.5*g*PAULI[j][:,0])<2e-10
    En=float(np.vdot(n,h@n).real);Ep=float(np.vdot(p[:,0],h@p[:,0]).real)
    residual=max(np.linalg.norm(h@n-En*n),np.linalg.norm(h@p-Ep*p))
    assert residual<2e-8
    d.update(native=native,W=W,initial=n,proton=p,En=En,Ep=Ep,P=P,ga=g)
    a=np.array([.2,-.3,.4]);b=np.array([-.1,.7,.1]);err=abs(explicit_spin_weight(a,b,d)-native_spin_weight(a[None],b[None],g)[0,0]);assert err<2e-10
    den=1+3*g*g
    return {'source':'model/openwave_spin_gradient_nuclei/source.py','actual_strong_parameters':{k:str(v) for k,v in parameters().items()},'certified_actual_axial_interval':enclosing,'actual_native_axial':g,'actual_native_vector':float(sum(P['plus'][x,0,0] for x in range(2)).real),'native_matrix_polynomial_error':err,'actual_hadron_eigenstate_residual':float(residual),'coefficients':{'a':(1-g*g)/den,'A':2*g*(1-g)/den,'B':2*g*(1+g)/den},'finite_cube_records':[records(L,d=d) for L in (4,8)],'strong_free_channel_endpoint_GeV':(En-Ep-1)*.00051099895069,'scope':'The same full native current, actual source quartet, certified axial interval, preparation and Gaussian lepton source satisfy the angular proof. The full coloured complement and calibrated source conditions come from the bound owning source; no label-only transfer.'}
