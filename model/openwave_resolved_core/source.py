"""Full-flux core and current certificates with an explicit regulator selection."""
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigh, solve
from model.openwave_dirac_maxwell_refinement.operators import (
    polynomial_member, charged_masses, root_tail, COORDINATES,
    ELECTRIC, MAGNETIC, SPEED_SQUARED)
from model.openwave_dressed_core_response.response import (
    exact_polynomial, finite_spectrum_bound, dyadic, norm_square, embed,
    exact_residual, sqrt_upper, inner, action)

CAP=24
DEGREE=6
CHARGE=-3
REGULATOR=F(1,2)
STATE_BUDGET=F(1,10**10)
RESPONSE_BUDGET=F(1,10**9)
WITNESS=Path(__file__).with_name('witnesses.json')
FREQUENCIES=tuple(map(F,('1/64','1/16','1/4','1')))
PROBE_STRENGTH=F(1,65536)
SPATIAL_FACTORS=tuple(map(F,('0','1','7/4','15/8','2')))


def generate_witnesses():
    """Generate candidates only; every consumer checks exact padded residuals."""
    result={}
    for name,mass in charged_masses():
        o=polynomial_member(CAP,CHARGE,mass,regulator=REGULATOR)
        n=len(o['basis'])//2;h=o['H_minus_mass'][:n,:n]
        es,u=eigh(h);v,vd=dyadic(u[:,0]);mu=F(round(float(es[0])*2**48),2**48)
        g=np.array([complex(a,b)/vd for a,b in v])
        x=np.array([COORDINATES[z[1]][0] for z in o['basis'][:n]])
        resolvents={}
        for nu in FREQUENCIES:
            z=solve(h-float(mu)*np.eye(n)+1j*float(nu)*np.eye(n),-x*g)
            resolvents[str(nu)]=dyadic(z)[0]
        spectra={}
        for spatial in SPATIAL_FACTORS:
            for sign in (-1,1):
                potential=sign*PROBE_STRENGTH*spatial
                hs=h+float(potential)*np.diag(x)
                spec,frame=eigh(hs);vec,_=dyadic(frame[:,0])
                spectra[f'{spatial}:{sign}']={'ground':vec,'centre':str(F(round(float(spec[0])*2**48),2**48))}
        result[name]={'ground':v,'centre':str(mu),'resolvents':resolvents,'scalar_probe_ground':spectra}
    return result


def member(name,mass):
    mass=F(mass);o=polynomial_member(CAP,CHARGE,mass,regulator=REGULATOR)
    n=len(o['basis'])//2;h=o['H_minus_mass'][:n,:n]
    assert np.array_equal(o['H_minus_mass'][:n,n:],np.zeros((n,n)))
    assert np.array_equal(h,o['H_minus_mass'][n:,n:].conj())
    es,u=eigh(h);basis,rows,den=exact_polynomial(CAP,CHARGE,mass,REGULATOR)
    reconstructed=np.zeros_like(h)
    for i,row in enumerate(rows):
        for j,(a,b) in row.items():reconstructed[i,j]=complex(float(F(a,den)),float(F(b,den)))
    assert np.max(abs(h-reconstructed))<3e-14
    second,finite=finite_spectrum_bound(basis,rows,den,(es,u))
    w=json.loads(WITNESS.read_text())[name];v=[tuple(map(int,z))for z in w['ground']]
    vd=1<<48;norm=norm_square(v,vd);mu=F(w['centre'])
    assert len(v)==len(basis) and abs(norm-1)<F(1,100)
    g=np.array([complex(a,b)/vd for a,b in v])/np.sqrt(float(norm))
    padded,prows,pden=exact_polynomial(CAP+2*DEGREE*abs(CHARGE),CHARGE,mass,REGULATOR)
    pv=embed(v,basis,padded);normlow=1-abs(norm-1)
    residual=exact_residual(prows,pden,pv,vd,mu)/normlow
    eps=root_tail(mass,DEGREE)
    outside=F((CAP+1-abs(CHARGE))**2,8)
    coupling=SPEED_SQUARED/(2*mass)+2*REGULATOR+MAGNETIC/2
    a1=second-eps;assert outside>a1
    infinite_second=a1-coupling*coupling/(outside-a1)
    gap=infinite_second-mu;assert gap>F(1,10)
    state_error=2*(residual+eps)/gap;assert state_error<=STATE_BUDGET
    x=np.array([COORDINATES[z[0]][0]for z in basis])
    p_exact=F(sum(a*a+b*b for (a,b),xx in zip(v,x)if xx),vd*vd)/norm
    p_interval=(p_exact-2*STATE_BUDGET,p_exact+2*STATE_BUDGET)
    assert F(1,4)<p_interval[0]<p_interval[1]<F(1,2)
    cert={'cap':CAP,'two_complete_spin_blocks_conjugate':True,
      'source_mass_electron_units_exact':str(mass),'Wilson_energy_coefficient_exact':str(REGULATOR),
      'complete_cap_spectrum':finite,'outside_electric_floor_exact':str(outside),
      'infinite_second_eigenvalue_lower_exact':str(infinite_second),
      'full_ground_energy_minus_bare_mass_interval_exact':[str(mu-residual-eps),str(mu+residual+eps)],
      'ground_vector_norm_error_upper_exact':str(STATE_BUDGET),
      'actual_state_error_upper_exact':str(state_error),
      'exact_padded_polynomial_residual_upper_exact':str(residual),
      'uniform_root_operator_error_exact':str(eps),
      'ground_radius_x_squared_cells_interval_exact':list(map(str,p_interval)),
      'variational_prepared_radius_x_squared_exact':str(p_exact)}
    return dict(name=name,mass=mass,regulator=REGULATOR,o=o,h=h,e=es,u=u,g=g,basis=basis,
      rows=rows,den=den,padded=padded,prows=prows,pden=pden,v=v,vd=vd,norm=norm,mu=mu,
      residual=residual,eps=eps,p_exact=p_exact,p_interval=p_interval,x=x,certificate=cert)


def response_records(d):
    x=d['x'];h=d['h'];u=d['u'];es=d['e'];g=d['g'];D=-x
    j=1j*(D[None,:]-D[:,None])*h;contact=-(D[None,:]-D[:,None])**2*h
    weights=abs(u.conj().T@(D*g))**2;jw=abs(u.conj().T@(j@g))**2;gaps=es-es[0]
    ward=float(np.vdot(g,contact@g).real-2*np.dot(gaps,weights));assert abs(ward)<1e-12
    v=d['v'];vd=d['vd'];norm=d['norm'];dv=[(-int(xx)*a,-int(xx)*b)for(a,b),xx in zip(v,x)]
    rhs=embed(dv,d['basis'],d['padded']);rhsnorm=sqrt_upper(norm_square(rhs,vd)/norm)
    w=json.loads(WITNESS.read_text())[d['name']];rows=[]
    for nu in FREQUENCIES:
        xv=[tuple(map(int,z))for z in w['resolvents'][str(nu)]];assert len(xv)==len(dv)
        px=embed(xv,d['basis'],d['padded'])
        residual=exact_residual(d['prows'],d['pden'],px,vd,d['mu'],rhs,nu)
        residual+=d['eps']*sqrt_upper(norm_square(px,vd))
        value=2*nu*nu*F(inner(dv,xv)[0],vd*vd)/norm
        error=2*nu*rhsnorm*residual/(1-abs(norm-1))+2*nu*STATE_BUDGET+2*(d['residual']+d['eps'])
        assert error<=RESPONSE_BUDGET
        diagnostic=float(2*float(nu)**2*np.sum(gaps*weights/(gaps*gaps+float(nu)**2)))
        assert abs(diagnostic-float(value))<1e-11
        rows.append({'imaginary_frequency_exact':str(nu),'response_interval_exact':[str(max(F(0),value-RESPONSE_BUDGET)),str(value+RESPONSE_BUDGET)],
          'actual_error_upper_exact':str(error),'response_diagnostic':diagnostic})
    assert all(F(a['response_interval_exact'][1])<F(b['response_interval_exact'][0])for a,b in zip(rows,rows[1:]))
    return {'static_Ward_residual_diagnostic':ward,'all_spectral_gaps_electron_units_diagnostic':gaps.tolist(),
       'all_dipole_weights_diagnostic':weights.tolist(),'all_current_weights_diagnostic':jw.tolist(),
       'individual_spectral_entries_are_diagnostics':True,'certified_infinite_responses':rows}


def probe_records(d):
    v=np.array([complex(a,b)/d['vd'] for a,b in d['v']]);v/=np.sqrt(float(d['norm']))
    h=d['h'];basis=d['basis'];x=d['x'];mass=float(d['mass'])
    initial=float(MAGNETIC);prepared=mass+float(np.vdot(v,h@v).real)
    preparation_work=.5*(prepared-initial)
    records=[]
    for k in (0.,np.pi/4,np.pi/2,np.pi):
        # The compensating reference remains at x=0: G=Q(cos(kx)-1).
        generator=1-np.cos(k*x);phase=np.exp(-1j*generator)
        spin_columns=[];spin_phase_energies=[];spin_final_energies=[]
        for state,hamiltonian in ((v,h),(v.conj(),h.conj())):
            columns=[]
            for axis in ('X','Y'):
                factor=1 if axis=='X' else -1j
                for sign in (-1,1):
                    columns.append((state+sign*factor*phase*state)/np.sqrt(8))
            nonselective=.5*np.outer(state,state.conj())+.5*np.outer(phase*state,(phase*state).conj())
            measured=sum((np.outer(w,w.conj())for w in columns),np.zeros_like(h))
            assert np.max(abs(nonselective-measured))<1e-15
            # The literal reported position/flux instrument dephases each
            # pointer branch. Its source-energy change is retained separately.
            dephased=np.diag(np.diag(measured))
            spin_phase_energies.append(mass+float(np.trace(measured@hamiltonian).real))
            spin_final_energies.append(mass+float(np.trace(dephased@hamiltonian).real))
            spin_columns.append(columns)
        probabilities=[sum(float(np.vdot(spin_columns[s][j],spin_columns[s][j]).real)for s in (0,1))/2 for j in range(4)]
        assert abs(sum(probabilities)-1)<1e-14
        phase_after=sum(spin_phase_energies)/2
        after=sum(spin_final_energies)/2
        phase_work=.5*(phase_after-prepared)
        position_work=.5*(after-phase_after)
        work=phase_work+position_work
        original_final=.5*initial+.5*after
        assert abs(original_final-initial-preparation_work-work)<2e-13
        records.append({'wavevector_radians_per_cell':k,
            'physical_probe_generator':'Q cos(k x)+R cos(0)=1-cos(k x); Q=-1,R=+1',
            'conditional_pointer_probabilities_Xminus_Xplus_Yminus_Yplus':probabilities,
            'complete_joint_position_flux_probabilities':[[ (.25*abs(spin_columns[s][j])**2).tolist() for s in (0,1)]for j in range(4)],
            'preparation_failure_probability_exact':'1/2',
            'spin_rule':'unpolarized two-spin mixture; each conjugate source branch is executed with the same probe phase; entries contain spin and preparation one-half factors',
            'original_input_phase_pointer_mean_work_electron_units':phase_work,
            'original_input_terminal_position_flux_mean_work_electron_units':position_work,
            'original_input_probe_mean_work_electron_units':work,
            'final_mean_source_energy_with_failure_electron_units':original_final,
            'actual_terminal_instrument':'phase pointer followed by all spin/position/flux projectors; the final source state is dephased in those reported outcomes'})
    p=d['p_interval'];rlo=sqrt_upper(p[0])-F(1,1<<63);rhi=sqrt_upper(p[1])
    assert F(0)<rlo<rhi<F(1)
    # x in {0,1}: exact all-wavevector formula; Taylor inequality is uniform.
    # Exact nested-radical endpoint brackets for a half-suppression crossing.
    from model.openwave_autonomous_clock_channel.clock import sqrt_bounds
    two_lo,two_hi=sqrt_bounds(F(2))
    cos_lo=sqrt_bounds(2+two_lo)[0]/2
    cos_hi=sqrt_bounds(2+two_hi)[1]/2
    lower_test=(p[0]*(1+two_lo/2),p[1]*(1+two_hi/2))
    upper_test=(p[0]*(1+cos_lo),p[1]*(1+cos_hi))
    assert lower_test[1]<F(1,2)<upper_test[0]
    from model.openwave_localized_rest_packet.packet import pi_interval
    pi_lo,pi_hi=pi_interval()
    onset_radius=(F(3,4)*pi_lo*rlo,F(7,8)*pi_hi*rhi)
    return {'position_flux_labels':[list(z)for z in basis],
        'preparation_input':'matter/reference vacuum and loop electric n=0; mean field energy1/8',
        'prepared_state':'normalized finite dyadic variational trial in each conjugate spin block, unpolarized; true infinite ground within certified norm error',
        'original_preparation_mean_work_electron_units':preparation_work,
        'input_failure_mean_field_energy_exact':'1/8',
        'preparation_external_above_true_field_vacuum':'1/8 minus the full rotor ground energy; separate from charged preparation work',
        'form_factor_exact':'F(k)=1-p_x1 (1-cos k), p_x1=<x^2>; charge response after retaining and separately reporting reference is Q F(k)',
        'reference_inclusive_response_exact':'Q[F(k)-1]; zero at k=0 because the complete state is neutral',
        'radius_x_cells_interval_exact':[str(max(F(0),rlo)),str(rhi)],
        'largest_first_Brillouin_zone_suppression_interval_exact':[str(2*p[0]),str(2*p[1])],
        'independent_radius_measurement':'complete native occupation/flux instrument, not inferred from a current gap',
        'small_wavevector_remainder':'|1-F(k)-k^2<x^2>/2| <= <x^2> k^4/24',
        'half_suppression_onset_resolved': bool(2*p[0]>F(1,2)),
        'half_suppression_wavevector_over_pi_interval_exact':['3/4','7/8'],
        'half_suppression_wavevector_times_measured_radius_interval_exact':list(map(str,onset_radius)),
        'suppression_at_three_quarter_pi_interval_exact':list(map(str,lower_test)),
        'suppression_at_seven_eighth_pi_interval_exact':list(map(str,upper_test)),
        'onset_test':'1-F(pi)=2<x^2>; compare its full ground-state interval with1/2',
        'records':records,'point_probabilities_and_work_are_numerical_diagnostics':True}



def with_scalar_potential(basis,rows,den,potential):
    potential=F(potential);scale=lcm(den,potential.denominator)
    result=[]
    for i,(coord,row) in enumerate(zip(basis,rows)):
        out={j:(a*(scale//den),b*(scale//den))for j,(a,b)in row.items()}
        a,b=out.get(i,(0,0));out[i]=(a+int(potential*scale)*COORDINATES[coord[0]][0],b)
        result.append(out)
    return result,scale


def scalar_probe_spectrum(d,spatial,sign):
    potential=sign*PROBE_STRENGTH*spatial
    w=json.loads(WITNESS.read_text())[d['name']]['scalar_probe_ground'][f'{spatial}:{sign}']
    v=[tuple(map(int,z))for z in w['ground']];vd=1<<48;mu=F(w['centre'])
    norm=norm_square(v,vd);assert abs(norm-1)<F(1,100)
    rows,den=with_scalar_potential(d['basis'],d['rows'],d['den'],potential)
    prows,pden=with_scalar_potential(d['padded'],d['prows'],d['pden'],potential)
    pv=embed(v,d['basis'],d['padded'])
    error=exact_residual(prows,pden,pv,vd,mu)/(1-abs(norm-1))+d['eps']
    second=F(d['certificate']['infinite_second_eigenvalue_lower_exact'])-abs(potential)
    assert mu+error<second
    assert error<F(1,10**10)
    prepared=F(inner(v,action(rows,v))[0],vd*vd*den)/norm+d['mass']
    terminal=sum(F((a*a+b*b)*rows[i][i][0],vd*vd*den)for i,(a,b)in enumerate(v))/norm+d['mass']
    preparation=(prepared-MAGNETIC)/2;readout=(terminal-prepared)/2
    final=(terminal+MAGNETIC)/2
    assert preparation+readout==final-MAGNETIC
    probabilities=[str(F(a*a+b*b,4*vd*vd)/norm)for a,b in v]
    assert 2*sum(map(F,probabilities))+F(1,2)==1
    return {'lambda_exact':str(sign*PROBE_STRENGTH),
      'energy_minus_bare_mass_interval_exact':[str(mu-error),str(mu+error)],
      'actual_energy_error_upper_exact':str(error),
      'both_spin_original_input_position_flux_probabilities_exact':[probabilities,probabilities],
      'preparation_failure_probability_exact':'1/2',
      'complete_native_and_flux_complement':'inherited full source; unprepared outcomes have explicit zero preparation probability',
      'preparation_work_interval_exact':[str(preparation-d['eps']/2),str(preparation+d['eps']/2)],
      'terminal_position_flux_work_interval_exact':[str(readout-d['eps']),str(readout+d['eps'])],
      'polynomial_work_balance_exact':'0',
      'full_source_total_work_minus_energy_change_exact':'0 by the source-energy definition of preparation and readout work',
      'intervention_work_on_input_vacuum_exact':'0; both particle and compensating-reference charge vanish on the failed vacuum input'}


def scalar_probe_coupling(d):
    base_upper=F(d['certificate']['full_ground_energy_minus_bare_mass_interval_exact'][1])
    gap=F(d['certificate']['infinite_second_eigenvalue_lower_exact'])-base_upper
    out=[]
    for spatial in SPATIAL_FACTORS:
        minus=scalar_probe_spectrum(d,spatial,-1);plus=scalar_probe_spectrum(d,spatial,1)
        em=list(map(F,minus['energy_minus_bare_mass_interval_exact']))
        ep=list(map(F,plus['energy_minus_bare_mass_interval_exact']))
        derivative=((ep[0]-em[1])/(2*PROBE_STRENGTH),(ep[1]-em[0])/(2*PROBE_STRENGTH))
        # Schur bound around the isolated ground, valid for both probe signs.
        assert gap>2*PROBE_STRENGTH*spatial
        curvature=PROBE_STRENGTH*spatial**2/(gap-2*PROBE_STRENGTH*spatial)
        hf=(d['p_interval'][0]*spatial,d['p_interval'][1]*spatial)
        assert derivative[0]-curvature<=hf[0]<=hf[1]<=derivative[1]+curvature
        coupling=(1-hf[1],1-hf[0])
        finite_coupling=(1-derivative[1]-curvature,1-derivative[0]+curvature)
        out.append({'one_minus_cos_wavevector_exact':str(spatial),
          'cos_wavevector_exact':str(1-spatial),'minus_intervention':minus,'plus_intervention':plus,
          'complete_source_energy_central_difference_interval_exact':list(map(str,derivative)),
          'finite_difference_to_derivative_error_upper_exact':str(curvature),
          'Hellmann_Feynman_total_charge_response_interval_exact':list(map(str,hf)),
          'particle_effective_charge_magnitude_interval_exact':list(map(str,coupling)),
          'independent_energy_difference_effective_charge_interval_exact':list(map(str,finite_coupling))})
    before=next(x for x in out if x['one_minus_cos_wavevector_exact']=='7/4')
    after=next(x for x in out if x['one_minus_cos_wavevector_exact']=='15/8')
    assert F(before['independent_energy_difference_effective_charge_interval_exact'][0])>F(1,2)
    assert F(after['independent_energy_difference_effective_charge_interval_exact'][1])<F(1,2)
    return {'same_source_static_native_charge_intervention':'H(lambda,k)=H+lambda[Q cos(kx)+R], Q=-1,R=1; on the charged sector this is H+lambda(1-cos k)x',
       'reference_term_reported_separately':True,
       'reference_subtracted_particle_derivative':'Q F(k); total neutral source derivative is Q F(k)+R',
       'small_k_onset':'1-F(k)=k^2*R_x^2/2+O(k^4*R_x^2), with the separately measured R_x^2=<x^2>',
       'scale_join_uses_no_frequency_to_momentum_dictionary':True,
       'complete_interventions':out,
       'static_effective_charge_is_universal_renormalized_charge':False}


def regulator_scan():
    rows=[]
    for regulator in map(F,('1/16','1/4','1/2','1','2')):
        o=polynomial_member(CAP,CHARGE,F(1),regulator=regulator)
        n=len(o['basis'])//2;h=o['H_minus_mass'][:n,:n];es,u=eigh(h)
        x=np.array([COORDINATES[z[1]][0]for z in o['basis'][:n]])
        p=float(abs(u[:,0])**2@x)
        rows.append({'Wilson_coefficient_exact':str(regulator),'finite_ground_energy_minus_mass_diagnostic':float(es[0]),
          'maximum_suppression_diagnostic':2*p,'radius_cells_diagnostic':p**.5})
    return {'all_tested_coefficients':rows,'selection':'First tested dyadic coefficient whose finite spatial suppression exceeds one half:1/2.',
      'selection_uses_resolution_outcome':True,'selection_uses_empirical_radius_or_running_data':False,
      'only_selected_member_has_infinite_ground_and_response_certificate':True}


def opposite_charge_control(mass):
    minus=polynomial_member(CAP,-3,mass,regulator=REGULATOR)
    plus=polynomial_member(CAP,3,mass,regulator=REGULATOR)
    index={z:i for i,z in enumerate(plus['basis'])}
    permutation=[index[(spin,x,-flux)]for spin,x,flux in minus['basis']]
    assert np.array_equal(plus['H_minus_mass'][permutation][:,permutation],minus['H_minus_mass'])
    return {'full_cap_dimension':len(permutation),'complete_permutation_checked':True,
      'infinite_intertwiner':'(spin,x,n,q) -> (spin,x,-n,-q); native opposite-charge copy retained',
      'electric_energy_preserved':'Every link electric flux changes sign; all squares are identical',
      'root_Wilson_and_magnetic_terms_preserved':True,
      'opposite_charge_has_identical_sector_energy_and_radius':True,
      'opposite_charge_scalar_probe':'Q and R reverse; lambda reverses, so both signs are in the complete intervention table',
      'whole_Standard_Model_charge_conjugation_symmetry_claimed':False}


def certificate():
    masses=charged_masses();rows=[]
    from model.openwave_common_mass_source.source import INPUTS
    unit=F(INPUTS['electron_unit_GeV'])
    from model.openwave_transverse_completion.source import field_floor
    vacuum=field_floor(CAP)
    vlo,vhi=map(F,vacuum['full_infinite_rotor_vacuum_interval_exact'])
    for name,mass in masses:
        d=member(name,mass);bounds=list(map(F,d['certificate']['full_ground_energy_minus_bare_mass_interval_exact']))
        rows.append({'family':name,'ground_certificate':d['certificate'],'opposite_charge_control':opposite_charge_control(mass),
          'absolute_charged_energy_GeV_interval_exact':[str(unit*(mass+x))for x in bounds],
          'reference_bound_energy_GeV_interval_exact':[str(unit*(mass+bounds[0]-vhi)),str(unit*(mass+bounds[1]-vlo))],
          'energy_reference':'true infinite compact-rotor vacuum; compensating quantum reference retained',
          'bare_mass_GeV_exact':str(unit*mass),'native_current_response':response_records(d),
          'complete_core_probe':probe_records(d),'same_source_scalar_probe_coupling':scalar_probe_coupling(d)})
    return {'schema':'oph.resolved_native_core.v1','verdict':'PASS_RESOLVED_INFINITE_FLUX_CORE_AND_COMPLETE_CURRENT_RECORDS',
      'architecture':{'parent':'openwave_dirac_maxwell_refinement','Maxwell_electric_coefficient_exact':str(ELECTRIC),
        'Maxwell_magnetic_coefficient_exact':str(MAGNETIC),'Wilson_coefficient_exact':str(REGULATOR),
        'all_integer_fluxes_both_spins_native_Fock_and_constraint_complements_retained':True,
        'positive_one_particle_reducing_sector_only':True,'electron_unit_GeV_exact':str(unit)},
      'full_rotor_vacuum_certificate':vacuum,'regulator_selection':regulator_scan(),'members':rows,
      'scope':{'half_suppression_core_onset_certified':True,'same_source_current_Ward_and_frequency_response':True,'same_source_static_charge_core_scale_join_certified':True,
        'physical_QED_or_QCD_beta_function_derived':False,'empirical_charge_radius_fit':False,
        'reference_bound_energies_are_observed_particle_poles':False,'unchanged_original_Wilson_member':False,
        'spatial_continuum_refinement_proved':False}}
