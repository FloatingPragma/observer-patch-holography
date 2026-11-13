"""Native three-body widths with explicit units and complete family sums."""
import json,math
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh
from model.openwave_neutrino_calibration.calibration import calibration,sqrt_bounds
from model.openwave_inclusive_decay.decay import spin_certificate,exact_phase_space,phase_factor,source_data

HERE=Path(__file__).resolve().parent

def inputs():return json.loads((HERE/'inputs.json').read_text())

def units():
    x=inputs();hbar=float(F(x['h_J_s_exact'])/(F(x['elementary_charge_C_exact'])*10**9))/(2*np.pi)
    assert 6.582119568e-25<hbar<6.582119570e-25
    return {'hbar_GeV_s':hbar,'GF_GeV_minus2':float(x['GF_GeV_minus2']),'native_Weyl_G_GeV_minus2':np.sqrt(8)*float(x['GF_GeV_minus2'])}

def kallen_root(x,y,z):
    a=x-(np.sqrt(y)+np.sqrt(z))**2;b=x-(np.sqrt(y)-np.sqrt(z))**2
    return np.sqrt(max(0,a)*max(0,b))

def massive_integral(r,a,b):
    lo=(r+a)**2;hi=(1-b)**2
    if hi<=lo:return 0.
    def f(s):
        return (1+b*b-s)*(s-r*r-a*a)*kallen_root(1,s,b*b)*kallen_root(s,r*r,a*a)/s
    return quad(f,lo,hi,epsabs=2e-13,epsrel=2e-12,limit=200)[0]

def mass_error_certificate():
    # Uniform for every admitted daughter/parent charged ratio and neutrino pair.
    eps=F(1,20000);mu=F(1,2*10**9);rmax=F(1,16)
    floor=(1-8*rmax*rmax-rmax**8)/12;assert floor>F(2,25)
    assert 2*rmax*mu+mu*mu<eps and 2*mu<eps
    relative=(eps*eps/2+eps**3/3)/floor+9*mu*mu/(eps*eps)
    assert relative<F(17,10**9)
    return {'neutral_to_parent_mass_upper':str(mu),'daughter_to_parent_mass_upper':str(rmax),'central_invariant_cutoff':str(eps),'massless_integral_lower':str(floor),'relative_width_loss_upper':str(relative),'relative_width_loss_upper_decimal':float(relative),'direction':'0<=Gamma_massless-Gamma_massive<=bound*Gamma_massless; complete family weights sum to one'}

def spectral_transfer_certificate():
    # The exact coarea densities include every spin, family and angular fibre.
    # I0 >= D^5/30, I0-Iab <=21 mu^2 D^3 => squared Hellinger bound.
    coefficient=F(21**2*30*32,15*256*3**4)
    assert coefficient<4
    return {'massive_integral_loss_upper':'21*mu_max^2*(1-r^2)^3',
            'massless_integral_lower':'(1-r^2)^5/30',
            'weighted_spectral_Hellinger_norm_squared_coefficient_upper':str(coefficient),
            'inclusive_coupling_L2_difference_upper':'2*mu_max^2',
            'full_unitary_survival_transfer_error':'2*G*t*mu_max^2',
            'probability_transfer_error':'4*G*t*mu_max^2',
            'scope':'Exact angular coarea of the massive and massless positive-frequency sources, common daughter-energy carrier with all original dark complements retained and zero-amplitude auxiliary channels below the massive threshold. This is a survival-preserving dilation. Parent survival and complete inclusive energy-bin records transfer; individual angular records are not identified.'}

def kernel_certificate():
    # Uniform massless-neutral continuum comparison for all three charged channels.
    hi=F(1,16);lo=F(1,4000);assert sum(F(10)**j/math.factorial(j) for j in range(30))>2/lo
    S0=F(32,3);S1=32+64*hi**2+9*hi**5+528*hi**4
    S2=80+96*hi**2+251*hi**4;S3=160+96*hi**2+52*hi**3
    D=3072*3**4;A0=F(256,315)/D;A2=(16*S2+64*S1+48*S0)/D;A3=(16*S3+96*S2+144*S1+48*S0)/D
    B0=F(1,1000);B1=F(1,125)
    assert 4*A0*A2<B0*B0 and F(27,8)*A0*A3*A3<B1**3
    assert (1-8*hi*hi-hi**8)/(1536*F(22,7)**3)>F(1,50000)
    return {'mass_ratio_domain':[str(lo),str(hi)],'Fourier_bounds_A0_A2_A3':list(map(str,(A0,A2,A3))),'correlation_L1_upper':str(B0),'weighted_correlation_L1_upper':str(B1),'dimensionless_Gamma_over_G2_lower':'1/50000','comparison_source':'massless neutrals, native positive-frequency current, ultraviolet amplitude (2-E)^2 on daughter energy interval [r,2]','massive_neutral_long_time_transfer':'separate exact weighted spectral certificate supplies the transfer'}

def channel(parent,daughter):
    c=calibration();u=units();M=np.asarray(c['charged'])*float(c['middle_charged_unit_GeV_exact']);nu=c['physical_masses_eV']*1e-9
    mass=M[parent];r=M[daughter]/mass;g=u['native_Weyl_G_GeV_minus2'];weights=abs(c['R'][parent])[:,None]**2*abs(c['R'][daughter])[None,:]**2
    assert abs(weights.sum()-1)<2e-15
    all_rows=[];integral=0.
    for i,j in product(range(3),repeat=2):
        I=massive_integral(r,nu[i]/mass,nu[j]/mass);w=float(weights[i,j]);integral+=w*I
        all_rows.append({'neutrino_mass_index':i,'antineutrino_mass_index':j,'native_current_family_weight':w,'three_body_integral_diagnostic':I,'width_GeV':g*g*mass**5*w*I/(128*np.pi**3)})
    zero=g*g*mass**5*phase_factor(r)/(1536*np.pi**3);width=g*g*mass**5*integral/(128*np.pi**3);bound=mass_error_certificate()
    assert (1-float(bound['relative_width_loss_upper_decimal'])-2e-12)*zero<=width<=(1+2e-12)*zero
    G=g*mass*mass;kc=kernel_certificate();B0=float(F(kc['correlation_L1_upper']));B1=float(F(kc['weighted_correlation_L1_upper']))
    Gamma=zero/mass;window=[]
    exact_mass=c['charged_exact'][parent]*c['middle_charged_unit_GeV_exact'];mu2=max(c['mass_squared_eV2_exact'])/10**18/exact_mass**2
    G2=8*F(inputs()['GF_GeV_minus2'])**2*exact_mass**4;Glo,Ghi=sqrt_bounds(G2)
    rel_massive=F(630)*mu2/F(255,256)**2
    assert rel_massive<F(1,10**15)
    for k in (1,4):
        time=k/Gamma;error=2*G*G*B1*(1+G*G*B0*time);finite_bin=time*2**-110
        t_upper=F(k*50000)/G2;transfer=4*Ghi*t_upper*mu2
        window.append({'massive_exact_survival_additional_probability_error_upper_rational':str(transfer),'massive_exact_survival_additional_probability_error_upper':float(transfer),'massive_rate_vs_massless_relative_error_upper_rational':str(rel_massive),'massive_exact_survival_total_probability_error_upper':error+float(transfer)+finite_bin,'lifetimes':k,'source_time':time,'seconds':time*u['hbar_GeV_s']/mass,'massless_reference_survival':float(np.exp(-k)),'massless_exact_survival_error_upper':error,'finite_energy_bin_probability_error_upper':finite_bin})
    return {'parent_family':parent,'daughter_family':daughter,'parent_mass_GeV':float(mass),'daughter_mass_GeV':float(M[daughter]),'all9_massive_neutral_channels':all_rows,'massless_neutral_width_GeV':float(zero),'massive_neutral_width_GeV_diagnostic':float(width),'massive_width_GeV_certified_relative_to_massless':[float((1-bound['relative_width_loss_upper_decimal'])*zero),float(zero)],'hbar_over_width_seconds':u['hbar_GeV_s']/width,'physical_Weyl_contact_dimensionless_parent_units':G,'calibrated_massive_finite_lifetime_windows':window,'energy_release_GeV':float(mass-M[daughter]-sum(nu[:2]))}

def finite_records(parent,daughter,bins=64,time=10000.):
    """All bright bins and the survival complement of the actual weak-normalized source."""
    d=channel(parent,daughter);r=d['daughter_mass_GeV']/d['parent_mass_GeV'];G=d['physical_Weyl_contact_dimensionless_parent_units']
    edges=np.linspace(r,2,bins+1);energies=(edges[:-1]+edges[1:])/2
    rho=lambda E:E**5*phase_factor(r/E)*(2-E)**4/(3072*np.pi**4)
    w=np.array([quad(rho,a,b,epsabs=1e-15,epsrel=1e-11)[0] for a,b in zip(edges[:-1],edges[1:])]);H=np.diag(np.r_[1.,energies]);H[0,1:]=H[1:,0]=G*np.sqrt(w)
    E,U=eigh(H);psi=U@(np.exp(-1j*E*time)*U[0].conj());prob=abs(psi)**2
    assert abs(prob.sum()-1)<2e-13 and abs(np.vdot(psi,H@psi).real-1)<2e-13
    work=float(prob[1:]@energies+prob[0]-1)
    return {'parent_family':parent,'daughter_family':daughter,'bins':bins,'source_time':time,'all_survival_and_emission_probabilities':prob.tolist(),'bin_energies':energies.tolist(),'measurement_work_source_units':work,'energy_before':1.,'energy_after_free_readout':1+work,'scope':'physical weak-normalized massless-neutral comparison member; finite bin recurrences are retained'}

def certificate():
    c=calibration();x=inputs();u=units();rows=[channel(1,0),channel(2,0),channel(2,1)];ratio=rows[2]['massive_neutral_width_GeV_diagnostic']/rows[1]['massive_neutral_width_GeV_diagnostic'];obs=float(x['tau_mu_over_e_comparator']);sigma=float(x['tau_mu_over_e_uncertainty'])
    cGF=np.pi/float(x['alpha_inverse'])/(np.sqrt(2)*float(x['W_mass_GeV'])**2*(1-float(x['W_mass_GeV'])**2/float(x['Z_mass_GeV'])**2))
    mu_tau=rows[0]['hbar_over_width_seconds'];ew_tau=mu_tau*(u['GF_GeV_minus2']/cGF)**2
    return {'schema':'oph.native_leptonic_widths_GeV.v1','verdict':'PASS_NATIVE_MASSIVE_FAMILY_WIDTHS_AND_WEAK_UNIT_DICTIONARY','unit_dictionary':u,'input_ancestry':x,'native_current_coefficient_and_spin_checks':spin_certificate(),'phase_space_coefficient_reconstruction':exact_phase_space(),'massive_neutral_width_bound':mass_error_certificate(),'massless_continuum_kernel_certificate':kernel_certificate(),'massive_spectral_survival_transfer_certificate':spectral_transfer_certificate(),'all_leptonic_channels':rows,'held_out_tau_mu_over_e':{'native_ratio':ratio,'PDG2025_comparison':obs,'PDG2025_uncertainty':sigma,'difference':ratio-obs,'difference_in_reported_uncertainties_diagnostic':(ratio-obs)/sigma,'used_in_parameter_selection':False,'scope':'tree-level inclusive native ratio; electromagnetic/radiative corrections and observed beam identification are not supplied'},'muon_normalization_comparison':{'hbar_over_tree_width_seconds':mu_tau,'PDG2025_seconds':float(x['muon_lifetime_seconds_comparator']),'relative_difference':mu_tau/float(x['muon_lifetime_seconds_comparator'])-1,'independent_prediction':False,'reason':'GF has muon-lifetime ancestry; native source has no radiative electromagnetic correction'},'independent_tree_electroweak_control':{'GF_from_alpha0_W_Z_GeV_minus2':cGF,'muon_characteristic_seconds':ew_tau,'relative_to_observed':ew_tau/float(x['muon_lifetime_seconds_comparator'])-1,'matches_precision_muon_lifetime':False,'scope':'on-shell tree dictionary from declared alpha(0),W,Z values; electroweak loop matching is absent'},'finite_bright_sector_records':[finite_records(1,0),finite_records(2,0),finite_records(2,1)],'source_contract':'All momentum, spin, neutral-family and angular dark states of the declared positive-frequency current are retained; this is not an equivalence claim for the unfiltered native contact source. The weighted spectral transfer connects massive parent survival and inclusive energy-bin records to the Volterra comparison; no angular-record transfer is inferred.','tau_complete_leptonic_source':{'both_daughter_continua_retained':True,'all_massive_neutral_family_channels':18,'leptonic_partial_width_sum_GeV':rows[1]['massive_neutral_width_GeV_diagnostic']+rows[2]['massive_neutral_width_GeV_diagnostic'],'conditional_e_fraction_within_leptonic_outputs':1/(1+ratio),'conditional_mu_fraction_within_leptonic_outputs':ratio/(1+ratio),'total_tau_interpretation':'No hadronic tau vertices are supplied; this sum is a leptonic partial width, not the physical total width.'},'total_tau_width_or_lifetime_claimed':False,'neutrino_input_mass_GeV':(c['physical_masses_eV']*1e-9).tolist()}
