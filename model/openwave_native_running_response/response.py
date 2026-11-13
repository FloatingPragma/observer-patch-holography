"""Full-current spectral sums, Ward contact term and a response coupling."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigh
from scipy.integrate import quad
from model.openwave_single_valley_leptons.spectrum import masses
from model.openwave_full_higgs_floor.spectrum import native_tensors
from model.openwave_lepton_mass_shell.spectrum import dirac_hamiltonian
from model.common_source.chiral_candidate import kinetic
ROOT=Path(__file__).resolve().parents[2]
INPUT=json.loads((ROOT/'model/openwave_lepton_widths/inputs.json').read_text())
ALPHA_INVERSE=F(INPUT['alpha_inverse'])
UNIT=F('0.00051099895069')


@lru_cache(None)
def gamma():
    beta=dirac_hamiltonian([0,0,0],1)
    alpha=[dirac_hamiltonian(np.eye(3)[j],1)-beta for j in range(3)]
    return alpha,beta


def native_factor():
    d=native_tensors();pairs=d['pairs'];weights=[F(x['left_charge'])**2 for x in pairs]
    assert sum(weights)==F(8,3)
    assert d['Q']*d['neutral']==d['neutral']*d['Q']==0*d['Q']
    return {'native_pairs':pairs,'charge_squared_weights_exact':list(map(str,weights)),
        'all_colour_and_lepton_weight_per_family_exact':str(sum(weights)),
        'neutral_current_exact_zero':True,'physical_CAR_per_site':90,'complete_Nambu_rows_per_momentum':180,
        'Nambu_duplication_factor_exact':'1/2'}


def block(angle,m,a):
    p=np.sin(angle)/a;M=m+np.sum(1-np.cos(angle))/a
    return dirac_hamiltonian(p,M),p,M


def full_native_control():
    d=native_tensors();C=np.asarray(d['D']+d['D'].T,float)
    Q=np.kron(np.eye(2),np.kron(np.asarray(d['Q'],float),np.eye(3)))
    Q=np.block([[Q,np.zeros_like(Q)],[np.zeros_like(Q),-Q.T]])
    A=kinetic.bdg_symbol([0,1,0],np.zeros((45,45)),speed=1.)
    B=kinetic.bdg_symbol([0,0,0],np.kron(C,np.eye(3)),speed=1.)
    ms=masses()/masses()[1];rows=[]
    for angle,a,nu in [(np.array([.3,.6,.4]),.2,1.),(np.array([.7,-.2,.5]),.4,3.)]:
        p=np.sin(angle)/a;w=sum(1-np.cos(angle))/a
        H=kinetic.bdg_symbol(p,np.kron(C,np.diag(ms+w)),speed=1.)
        J=Q@(np.cos(angle[1])*A+np.sin(angle[1])*B)
        e,U=eigh(H);minus=e<0;plus=e>0;assert minus.sum()==plus.sum()==90
        V=U[:,minus].conj().T@J@U[:,plus];gap=e[plus][None,:]-e[minus][:,None]
        actual=.5*np.sum(2*abs(V)**2*nu**2/(gap**3*(gap**2+nu**2)))
        expected=0.
        for m in ms:
            E=np.sqrt(p@p+(m+w)**2);dot=p[1]*np.cos(angle[1])+(m+w)*np.sin(angle[1])
            S=2*(1-dot*dot/(E*E));expected+=float(F(8,3))*S*nu**2/(4*E**3*(4*E**2+nu**2))
        assert abs(actual-expected)<2e-11*max(1,abs(expected))
        rows.append({'angle':angle.tolist(),'spacing':a,'frequency':nu,
            'complete180_spectral_response':float(actual),'native_pair_formula':float(expected),
            'absolute_difference':float(abs(actual-expected)),
            'all_negative_positive_transition_columns':8100})
    return rows


def ward_control(side=4,a=.5):
    alpha,beta=gamma();points=np.array(list(product(range(side),repeat=3)));V=len(points)
    angles=2*np.pi*points/side;Fmat=np.exp(1j*points@angles.T)/np.sqrt(V)
    hs=np.array([block(k,1,a)[0]for k in angles])
    H=np.einsum('xk,kab,yk->xayb',Fmat,hs,Fmat.conj()).reshape(4*V,4*V)
    theta=np.cos(2*np.pi*points[:,0]/side)+.3*np.sin(2*np.pi*points[:,1]/side)
    G=np.diag(np.repeat(theta,4));J=1j*(G@H-H@G);T=-(G@(G@H-H@G)-(G@H-H@G)@G)
    e,U=eigh(H);occupied=e<0;unoccupied=e>0
    Vj=U[:,occupied].conj().T@J@U[:,unoccupied];gap=e[unoccupied][None,:]-e[occupied][:,None]
    para=-2*np.sum(abs(Vj)**2/gap)
    dia=np.trace(U[:,occupied].conj().T@T@U[:,occupied]).real
    assert abs(para+dia)<3e-10
    assert abs(para)>1
    return {'side':side,'full_Dirac_matrix_dimension':len(H),'paramagnetic_second_derivative':float(para),
        'diamagnetic_contact_second_derivative':float(dia),'total_pure_gauge_energy_response':float(para+dia),
        'omitting_contact_fails':True,
        'operator_identity':'H(epsilon)=exp(i epsilon Q theta) H exp(-i epsilon Q theta), J=i[G,H], T=-[G,[G,H]]',
        'full_native_transfer':'Each exact charged pair and conjugate Nambu copy obeys the same identity; neutral current is zero.'}


def finite_response(side,a,frequencies,records=False):
    one=2*np.pi*np.arange(side)/side
    kx,ky,kz=np.meshgrid(one,one,one,indexing='ij')
    p2=(np.sin(kx)**2+np.sin(ky)**2+np.sin(kz)**2)/a**2
    w=(3-np.cos(kx)-np.cos(ky)-np.cos(kz))/a
    cy=np.cos(ky);sy=np.sin(ky);py=sy/a;V=side**3
    ms=masses()/masses()[1];sums=np.zeros(len(frequencies));family=[];norms=[];coarse=[]
    for f,m in enumerate(ms):
        E=np.sqrt(p2+(m+w)**2)
        # Positive algebraic form avoids cancellation near a current eigenvector.
        S=2*((np.sin(kx)**2+np.sin(kz)**2)/a**2+(py*sy-(m+w)*cy)**2)/E**2
        assert S.min()>=0
        mean=float(np.sum(-2*(py*cy+(m+w)*sy)/E))
        assert abs(mean)<1e-10*V
        values=[]
        for nu in frequencies:
            value=float(np.sum(S*nu**2/(4*E**3*(4*E**2+nu**2)))/(V*a**3))
            assert value>0
            values.append(value);sums[len(values)-1]+=float(F(8,3))*value
        family.append({'analyzer_character':f,'mass_electron_units':float(m),'unit_charge_response':values})
        norms.append(float(np.sum(S)))
        if records:coarse.append({'character':f,'all_momentum_spin_summed_pair_weights':S.ravel().tolist()})
    e2=4*np.pi/float(ALPHA_INVERSE);assert np.max(e2*sums)<1
    rows=[{'Euclidean_frequency_electron_units':float(nu),'Euclidean_frequency_GeV':float(nu)*float(UNIT),
        'native_subtracted_current_response':float(p),'response_alpha_inverse':float(ALPHA_INVERSE)-4*np.pi*float(p),
        'fractional_response_coupling_increase':float(1/(1-e2*p)-1)}for nu,p in zip(frequencies,sums)]
    for i in range(1,len(rows)):assert rows[i]['native_subtracted_current_response']>rows[i-1]['native_subtracted_current_response']
    r={'side':side,'spacing_electron_inverse_mass':a,'box_length_electron_inverse_mass':side*a,
        'all_native_Nambu_coordinates':180*V,'rows':rows,'all_three_families':family,
        'sampling_error':'none; complete finite Fourier sum','continuum_or_volume_error_certified':False}
    if records:
        scale=90*V;success=float(F(8,3))*sum(norms)/scale**2;assert 0<success<1
        r['complete_native_current_instrument']={'current_norm_bound_exact':str(scale),
            'Kraus_success':'J/(90 V), with exact positive square-root failure complement',
            'probability_rule':'pair_charge_squared times momentum_spin_weight divided by (90 V)^2; all seven native pairs and all three families',
            'all_momentum_addresses':[list(x)for x in product(range(side),repeat=3)],
            'families':coarse,'pair_charge_squared_exact':native_factor()['charge_squared_weights_exact'],
            'total_success_diagnostic':success,'failure_probability_diagnostic':1-success,
            'all_other_source_Fock_and_neutral_outputs_retained':True}
    return r


def continuum_reference(nu,m=1.):
    # This identity follows by radial integration of the same continuum
    # Dirac spectral sum; no beta-function value is supplied to the finite source.
    angular=lambda p:p*p*(1-p*p/(3*(p*p+m*m)))/((p*p+m*m)**1.5*(4*(p*p+m*m)+nu*nu))
    radial=nu*nu/(4*np.pi**2)*quad(angular,0,np.inf,epsabs=1e-12)[0]
    parameter=quad(lambda x:x*(1-x)*np.log1p(nu*nu*x*(1-x)/(m*m)),0,1,epsabs=1e-12)[0]/(2*np.pi**2)
    assert abs(radial-parameter)<3e-11
    return {'frequency':nu,'radial_current_integral':radial,'independently_integrated_parameter_identity':parameter}


def certificate():
    frequencies=[.05,.2,.5,1.,2.,4.,8.,16.]
    ms=masses()/masses()[1]
    return {'schema':'oph.native_Wilson_current_running_response.v1','verdict':'PASS_COMPLETE_NATIVE_CURRENT_RESPONSE_AND_WARD_CONTACT',
        'source':'complete 90-CAR native Higgs-Wilson source, filled negative bands and vacant neutral zeros; distinct from the192-CAR interacting transverse member',
        'input_ancestry':{'alpha_inverse_supplied_exact':str(ALPHA_INVERSE),'electron_anchor_GeV_exact':str(UNIT),
            'alpha_role':'zero-frequency normalization of an explicitly defined quadratic response kernel',
            'mass_source':'retrospective three-family analyzer and native common Higgs channels',
            'observed_quark_masses_fitted':False,'Standard_Model_beta_function_inserted':False},
        'exact_native_current_factor':native_factor(),'full180_native_spectral_controls':full_native_control(),
        'exact_gauge_Ward_and_diamagnetic_control':ward_control(),
        'response_definition':{'spectral_formula':'Pi_R(nu)=Vphys^-1 sum_occ,empty 2|J_ab|^2 nu^2/[Delta_ab^3(Delta_ab^2+nu^2)] with Nambu factor1/2',
            'coupling':'alpha_response(nu)^-1=alpha_IR^-1-4 pi Pi_R(nu)',
            'meaning':'exact second-order Gaussian-fermion current response and declared quadratic photon-probe kernel; not the all-orders interacting scattering coupling',
            'positivity_and_monotonicity':'Pi_R(0)=0 and derivative positive for nu>0 when pair current is nonzero',
            'finite_source_pair_threshold_electron_units':'2',
            'core_dictionary':'inverse pair threshold is a response time/length scale; no independently localized spatial quark core is identified'},
        'complete_current_records':finite_response(4,.5,frequencies,records=True),
        'fixed_box_cutoff_controls':[finite_response(n,8/n,frequencies)for n in (32,64,128)],
        'fixed_spacing_volume_controls':[finite_response(n,.25,frequencies)for n in (24,32,48)],
        'continuum_radial_identity_controls':[continuum_reference(n)for n in (.2,1.,4.,16.)],
        'derived_continuum_small_frequency_coefficient':{'Pi_R_over_nu_squared':float(F(8,3))*float(sum(1/(60*np.pi**2*m*m)for m in ms)),
            'formula':'sum_native_pairs,f q_pair^2/(60 pi^2 m_f^2)'},
        'scope':'A native quantum current response with supplied electromagnetic normalization, all retained massive charge branches and explicit finite-regulator effects. The equal up/down/lepton source masses preclude interpreting the complete curve as physical QED or QCD running. No continuum cutoff error, full interacting scattering theorem, particle-core radius, or Yang-Mills gap is asserted.'}
