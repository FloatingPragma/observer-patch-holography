"""Complete native beta instruments on a three-dimensional finite kinetic source."""
from fractions import Fraction as F
from collections import Counter
from itertools import product
import numpy as np
import sympy as sy
from model.openwave_beta_outgoing.response import source, G
from model.openwave_baryon_weak_formfactors.native import PAULI, currents


def spin_polynomial():
    g=sy.symbols('g',real=True);a=sy.symbols('ax ay az',real=True);b=sy.symbols('bx by bz',real=True)
    P=[sy.Matrix(p) for p in PAULI];P=[p.applyfunc(sy.nsimplify) for p in P]
    da=(P[0]-sum((a[j]*P[j+1] for j in range(3)),sy.zeros(2)))/2
    db=(P[0]-sum((b[j]*P[j+1] for j in range(3)),sy.zeros(2)))/2
    total=0
    for s in range(2):
        M=(P[0] if s==0 else sy.zeros(2))-g*sum((P[j][s,0]*P[j] for j in range(1,4)),sy.zeros(2))
        total+=sy.trace(da*M*db*M.conjugate().T)
    expected=(1+3*g*g+(1-g*g)*sum(x*y for x,y in zip(a,b))+2*g*(1-g)*a[2]+2*g*(1+g)*b[2])/2
    assert sy.expand(total-expected)==0
    return str(sy.expand(total))


def grid(L):
    assert L>=4 and L%2==0
    theta=np.array(list(product(2*np.pi*np.arange(L)/L,repeat=3)))
    p=2*np.sin(theta);p[np.abs(p)<1e-14]=0
    return theta,p



def exact_free_electron_energy_census(L):
    """Exact distinct free energies; all spatial momentum multiplicities retained."""
    # For these executed cubes, 4 sin²(2πn/L) has an exact integer cycle.
    if L==4:one=(0,4,0,4)
    elif L==8:one=(0,2,4,2,0,2,4,2)
    else:raise ValueError('Exact integer census is certified for sides4 and8 only')
    counts=Counter(1+sum(p2)for p2 in product(one,repeat=3))
    levels=sorted(counts)
    assert sum(counts.values())==L**3
    return {'squared_energy_source_units':levels,
            'spatial_momentum_multiplicities':[counts[x]for x in levels],
            'electron_energy_GeV':[float(np.sqrt(x)*.00051099895069)for x in levels],
            'scope':'Distinct free-energy levels; multiplicities count spatial momenta only. All spin, family, direction and probability records are unchanged.'}


def densities(p,m,sign):
    E=np.sqrt(np.sum(p*p,axis=1)+m*m)
    return (PAULI[0][None,:,:]+sign*np.einsum('ni,ijk->njk',p/E[:,None],np.asarray(PAULI[1:])))/2,E


def native_spin_weight(a,b,g):
    return (1+3*g*g+(1-g*g)*(a@b.T)+2*g*(1-g)*a[:,2,None]+2*g*(1+g)*b[None,:,2])/2


def explicit_spin_weight(a,b,d=None):
    """Independent native 120-state current contraction, without fitted spin weights."""
    d=source() if d is None else d;total=0j
    da=(PAULI[0]-sum(a[j]*PAULI[j+1] for j in range(3)))/2
    db=(PAULI[0]-sum(b[j]*PAULI[j+1] for j in range(3)))/2
    for sp in range(2):
        M=sum((1 if mu==0 else -1)*sum(d['P']['plus'][x,mu,sp] for x in range(2))*PAULI[mu] for mu in range(4))
        total+=np.trace(da@M@db@M.conj().T)
    assert abs(total.imag)<1e-12
    return float(total.real)


def wick_column(L,d=None):
    d=source() if d is None else d;theta,p=grid(L);N=len(p);rho=np.zeros((24,24),complex)
    def ix(x,s,w,f):return ((2*x+s)*2+w)*3+f
    for w in range(2):
        masses=d['neutral'] if w==0 else d['charged']
        bands=[densities(p,m,1)[0] for m in masses]
        R=d['R'] if w==0 else np.eye(3)
        for x,y,s,r in product(range(2),range(2),range(2),range(2)):
            values=[np.sum(np.exp(1j*(x-y)*theta[:,0])*b[:,s,r])/N for b in bands]
            block=R@np.diag(values)@R.conj().T
            for f,j in product(range(3),repeat=2):rho[ix(x,s,w,f),ix(y,r,w,j)]=block[f,j]
    assert np.linalg.norm(rho-rho.conj().T)<1e-13
    J=[];Q=[];fac=[];cur=currents()['local'];weak=[PAULI[j]/2 for j in range(1,4)]
    for x,mu,a in product(range(2),range(4),range(3)):
        matrix=np.zeros((24,24),complex)
        for s,r,w,v,f in product(range(2),range(2),range(2),range(2),range(3)):
            matrix[ix(x,s,w,f),ix(x,r,v,f)]=PAULI[mu][s,r]*weak[a][w,v]
        qp,qz=cur[x][mu];q=((qp+qp.conj().T)/2,(qp-qp.conj().T)/(2j),qz)[a]
        J.append(matrix);Q.append(q@d['initial']);fac.append(2*(1 if mu==0 else -1))
    means=np.array([np.trace(j@rho) for j in J]);assert np.max(abs(means))<3e-14
    value=sum(fac[a]*fac[b]*np.vdot(Q[a],Q[b])*np.trace(J[a]@(np.eye(24)-rho)@J[b]@rho) for a,b in product(range(24),repeat=2))
    assert abs(value.imag)<1e-10 and 0<value.real<54**2
    return float(value.real)


def records(L,t=F(1,10000),d=None):
    d=source() if d is None else d;theta,p=grid(L);N=len(p);Ee=np.sqrt(1+(p*p).sum(axis=1));ae=p/Ee[:,None]
    energy_census=exact_free_electron_energy_census(L)if L in(4,8)else None
    # Physical antineutral momentum is minus the canonical hole momentum.
    selected=np.zeros((N,N));neutral=[];g=d['ga']
    form=np.cos((theta[:,0,None]-theta[None,:,0])/2)**2
    minweight=1e9;maxdelta=0.
    for j,m in enumerate(d['neutral']):
        En=np.sqrt((p*p).sum(axis=1)+m*m);bn=-p/En[:,None]
        K=native_spin_weight(ae,bn,g);minweight=min(minweight,float(K.min()));assert K.min()>-2e-12
        mismatch=d['Ep']-d['En']+Ee[:,None]+En[None,:]
        maxdelta=max(maxdelta,float(abs(mismatch).max()))
        w=(float(G*t)/N)**2*abs(d['R'][0,j])**2*form*K*np.sinc(mismatch*float(t)/(2*np.pi))**2
        selected+=w;neutral.append(float(w.sum()))
    # Three independent 3D direction observables plus complete radial-energy bins.
    radius=np.linalg.norm(p,axis=1);angular=np.zeros((3,3));sign=np.sign(p[:,2]).astype(int)+1
    for a,b in product(range(3),repeat=2):angular[a,b]=selected[np.ix_(sign==a,sign==b)].sum()
    radial={str(v):float(selected[np.isclose(radius,v)].sum()) for v in np.unique(np.round(radius,12))}
    delta=float((54*G*t)**2/2);amp=float(G*t)*np.sqrt(2*(1+3*g*g));event_error=(2*amp+delta)*delta
    odd=float(angular[2].sum()-angular[0].sum());assert odd+event_error<0
    intervals={f'e_z_{a-1}_hole_z_{b-1}':[max(0,float(angular[a,b])-event_error)/120,min(1,float(angular[a,b])+event_error)/120] for a,b in product(range(3),repeat=2)}
    intervals['failed_preparation']=[119/120,119/120]
    intervals['prepared_full_quantum_complement']=[max(0,1-float(selected.sum())-event_error)/120,min(1,1-float(selected.sum())+event_error)/120]
    # A globally safe rational error used for the exact-source parity margin.
    assert amp<4*float(G*t)
    safe_error=float((8*G*t+(54*G*t)**2/2)*((54*G*t)**2/2))
    # Finite symmetric cube: reflecting only p_z leaves form/mismatch invariant;
    # sum over antineutral p_z cancels B and a cross-terms. Thus the odd
    # component has exact sign from 2g(1-g), with a simple selected-mode bound.
    # At least L/4 electron z modes have |p_z|>=sqrt(2), and Ee<=sqrt(13)<4.
    # Sum_h form=N/2; electron x,y are unrestricted. Spin factor yields
    # >= (Gt)^2 g(g-1)/(8 sqrt(2)); use 1/12 for a rational weaker bound.
    glo=F('1.2752');ghi=F('1.2754');lower=(G*t)**2*glo*(glo-1)/12
    # Sinc weighting differs from 1 by at most mismatch^2*t^2/12, all 24
    # current channels are in the full remainder. Selected mismatch<9.
    assert maxdelta<9
    weight_error=(G*t)**2*16*(9*t)**2/12
    # Finite neutral masses preserve p_z reflection exactly, so no mass limit.
    exact_margin=lower-weight_error-F(str(safe_error));assert exact_margin>0
    assert odd<-float(lower-weight_error)
    return {'cube_side':L,'spatial_sites':N,'native_CAR_coordinates':144+24*N,'source_time':str(t),'first_order_selected_probability':float(selected.sum()),'selected_probability_by_neutral_mass':neutral,'first_order_joint_z_sign_probabilities':angular.tolist(),'electron_energy_GeV':energy_census['electron_energy_GeV']if energy_census is not None else sorted(set((Ee*.00051099895069).tolist())),'exact_free_electron_energy_census':energy_census,'radial_momentum_source_weights':radial,'parity_odd_first_order_probability':odd,'exact_full_Fock_parity_upper_rational':str(-exact_margin),'exact_full_Fock_parity_upper':-float(exact_margin),'exact_full_Fock_probability_error_upper':safe_error,'complete_original_record_intervals':intervals,'full_current_Wick_column_squared_norm':wick_column(L,d),'momentum_velocity_distinction':'p_j=2 sin(k_j), v_j=4 sin(k_j) cos(k_j)/E; the instrument records momentum, not velocity','selected_hadron_state_count':2,'retained_hadron_state_count':120,'unresolved_quantum_complement':'all other hadrons, all charged/neutral families, all helicities, every multipair sector, and no emission'}


def certificate():
    d=source();g=d['ga'];den=1+3*g*g
    coeff={'a':(1-g*g)/den,'A':2*g*(1-g)/den,'B':2*g*(1+g)/den,'D':0.}
    errors=[]
    for a,b in [(np.array([.2,.3,.4]),np.array([-.1,.6,-.2])),(np.array([0.,0.,.8]),np.array([.6,0.,-.8])),(np.zeros(3),np.zeros(3))]:
        errors.append(abs(explicit_spin_weight(a,b)-native_spin_weight(a[None,:],b[None,:],g)[0,0]))
    assert max(errors)<2e-11
    from .transfer import certificate as transferred_certificate
    return {'positive_spin_gradient_same_current_transfer':transferred_certificate(),'schema':'oph.native_3d_beta_angular.v1','verdict':'PASS_COMPLETE_3D_NATIVE_PARITY_INSTRUMENT','native_spin_polynomial':spin_polynomial(),'native_axial':g,'allowed_point_contact_angular_coefficients':coeff,'native_matrix_polynomial_max_error':max(errors),'PDG_2025_zero_recoil_comparators':{'source':'https://pdg.lbl.gov/2025/listings/rpp2025-list-n.pdf','a0':{'value':-.1044,'uncertainty':.0007,'native_minus_central':coeff['a']+.1044},'A':{'value':-.11958,'uncertainty':.00021,'native_minus_central':coeff['A']+.11958},'B':{'value':.9807,'uncertainty':.003,'native_minus_central':coeff['B']-.9807},'role':'Comparison-only values; the calibrated gA average shares experimental ancestry with angular data, so these are conditional consistency comparisons, not independent validation. A and a0 use zero recoil; no experiment-specific recoil or radiative model is installed.'},'finite_cube_records':[records(L) for L in (4,8)],'calibration_inputs':['p,n,charged-pion masses and |gA| from openwave_axial_state','six neutral oscillation inputs and declared lightest mass from openwave_neutrino_calibration','electron energy anchor 0.00051099895069 GeV'],'source_geometry':'two native internal hadron cells at spatial sites (0,0,0),(1,0,0); 3D periodic Weyl kinetic p=2sin(k), every valley retained','native_current':'unchanged full 120-state baryon current, full lepton Fock space, G=1/64','strong_free_channel_endpoint_GeV':(d['En']-d['Ep']-1)*.00051099895069,'full_interacting_pole_endpoint_claimed':False,'angular_scope':'a,A,B are exact native point-contact spin polynomials. The executed finite two-contact cubes retain their anisotropic form factor; no physical continuum angular fit is asserted. gA is an empirical input.','work':{'nonselective_preparation_mean':'0','autonomous_mean_energy':'conserved','terminal_readout_absolute_work_source_units':'2G||V||<=27/16'},'physical_weak_clock_or_lifetime_claimed':False}
