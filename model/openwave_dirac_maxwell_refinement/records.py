"""Full-flux source, native preparation/readout, and bounded finite executions."""
from fractions import Fraction as F
import numpy as np
from .operators import polynomial_member,charged_masses,ELECTRIC,MAGNETIC,SPEED_SQUARED
from .certificates import complete_eigensystem,actual_vector_certificate
from .bounds import infinite_flux_bound
from .free import refinement_records
from model.openwave_transverse_completion.source import flux,gauss_checks,field_floor
from model.openwave_common_mass_source.source import lepton_source,ix,INPUTS


def execute(cap,charge,mass,time_values=(2,4),family_index=0):
    member=polynomial_member(cap,charge,mass)
    n=len(member['basis'])//2
    H=member['H_minus_mass'];K=member['kinetic']
    assert np.count_nonzero(H[:n,n:])==np.count_nonzero(H[n:,:n])==0
    assert np.array_equal(H[n:,n:],H[:n,:n].conj())
    cert=complete_eigensystem(H[:n,:n],mass,6)
    Q=cert['Q'];lam=cert['eigenvalues'];index=member['index'][(0,0,0)]
    native=lepton_source();T=np.vstack((native['A'],native['B'].conj()))
    component=11 if charge==-3 else 9
    j=ix(0,0,component,family_index)
    native_prob=abs(T[:,j])**2
    assert abs(sum(native_prob)-1)<2e-15
    assert int(3*native['charge'][j])==charge
    electric=np.array([float(ELECTRIC*sum(v*v for v in flux(x,z,charge)))
                       for _,x,z in member['basis'][:n]])
    magnetic=H[:n,:n]-K[:n,:n]-np.diag(electric)
    initial=float(H[index,index].real)
    rows=[]
    for time in time_values:
        psi=Q@(np.exp(-1j*time*lam)*Q[index].conj())
        prob=abs(psi)**2
        assert abs(prob.sum()-1)<2e-13
        assert abs(np.vdot(psi,H[:n,:n]@psi).real-initial)<2e-12
        finite=actual_vector_certificate(cert,psi,index,time)
        infinite=infinite_flux_bound(cap,charge,mass,time)
        error=F(finite['conditional_vector_error_upper_exact'])+F(infinite['state_norm_error_before_numerical_diagonalization_exact'])
        rounding=sum((abs(F(float(p))-F(float(z.real))**2-F(float(z.imag))**2)
                      for p,z in zip(prob,psi,strict=True)),F(0))/4
        assert rounding<F(1,10**14)
        # Unconditional success probability1/2; the failure identity is retained.
        total_TV=error*(2+error)/4+F(1,10**14)
        ke=float(np.vdot(psi,K[:n,:n]@psi).real)
        me=float(np.vdot(psi,magnetic@psi).real)
        ee=float(np.dot(prob,electric))
        preparation=.5*(float(mass)+initial-float(MAGNETIC))
        terminal=.5*(-float(mass)-ke+float(MAGNETIC)-me)
        balance=preparation+terminal-.5*ee
        assert abs(balance)<3e-12
        rows.append({'cap':cap,'time_exact':str(time),'native_charge_in_thirds':charge,
                     'electron_charge_exact':str(F(charge,3)),
                     'initial_native_positive_mode':j,'initial_spin_exact':'up',
                     'complete_position_flux_labels':[list(z[1:]) for z in member['basis'][:n]],
                     'conditional_position_flux_probabilities':prob.tolist(),
                     'native_probe_probabilities_96':native_prob.tolist(),
                     'complete_joint_record_rule':'success probability1/2 times every position/flux probability times every native96-row weight; failure probability1/2; opposite spin and all native/Fock/Gauss complements retained with their exact zero or identity outcomes',
                     'preparation_failure_probability_exact':'1/2',
                     'finite_numerical_certificate':finite,'full_flux_truncation_certificate':infinite,
                     'original_input_position_flux_record_TV_upper_exact':str(total_TV),
                     'native_probe_numerical_scope':'fixed native analyzer weights are displayed floating diagnostics; certified TV concerns position/flux and preparation failure',
                     'conditional_source_energy_electron_units_diagnostic':float(mass)+initial,
                     'conditional_kinetic_energy_diagnostic':ke,'conditional_magnetic_energy_diagnostic':me,
                     'conditional_electric_energy_diagnostic':ee,
                     'preparation_mean_source_work_diagnostic':preparation,
                     'terminal_absorption_mean_source_work_diagnostic':terminal,
                     'initial_and_failure_field_energy_exact':str(MAGNETIC),
                     'final_mean_source_energy_diagnostic':float(MAGNETIC)+.5*ee,
                     'complete_work_balance_residual_diagnostic':balance})
    spectrum={k:v for k,v in cert.items() if k not in ('Q','eigenvalues','X','Y','E')}
    centre=F(int(cert['E'][0]),cert['energy_denominator'])
    err=F(cert['uniform_eigenvalue_error_upper_exact'])+F(member['uniform_root_operator_error_exact'])
    low,high=centre-err,centre+err
    outside=F((cap+1-abs(charge))**2,8)
    rho=SPEED_SQUARED/(mass*mass)
    boundary=mass*rho/(1-rho)+F(1,16)+MAGNETIC
    assert outside>high
    low-=boundary*boundary/(outside-high)
    vacuum=field_floor(cap)
    vlo,vhi=map(F,vacuum['full_infinite_rotor_vacuum_interval_exact'])
    gap=[mass+low-vhi,mass+high-vlo]
    unit=F(INPUTS['electron_unit_GeV'])
    return {'family_index':family_index,'native_charge_in_thirds':charge,'bare_mass_exact':str(mass),
            'cap':cap,'full_spin_dimension':2*n,'exact_spin_reduction':'H_down=conjugate(H_up), cross-spin entries exactlyzero',
            'complete_spin_block_spectrum_certificate':spectrum,
            'full_infinite_one_excitation_rest_minus_vacuum_interval_exact':list(map(str,gap)),
            'conditional_GeV_rest_interval_exact':[str(x*unit) for x in gap],
            'rest_energy_scope':'one-excitation ground relative to the full field vacuum with a compensating reference; no isolated physical pole or measured mass fit',
            'outside_electric_floor_exact':str(outside),'boundary_norm_upper_exact':str(boundary),
            'field_vacuum_certificate':vacuum,'records':rows}


def run():
    masses=charged_masses()
    cases=[]
    for f,(name,mass) in enumerate(masses):
        for charge in (-3,3):
            for cap in (24,48):
                value=execute(cap,charge,mass,family_index=f)
                value['family']=name;cases.append(value)
    return {'schema':'oph.native_Pauli_Maxwell_full_flux_refinement.v1',
            'verdict':'PASS_FULL_FLUX_PAULI_SOURCE_AND_COMMON_FREE_REFINEMENT',
            'construction':{'native_lepton_frame':'actual48-mode charged/neutral Bogoliubov frame;96-row probe',
                            'source':'positive number-preserving Pauli-root matter with positive Wilson energy plus compact Maxwell field',
                            'link_space':'L2(U1^4); every integer electric flux retained',
                            'one_particle_reduction':'fixed family/charge and compensating reference; two spin states and every cycle flux',
                            'kinetic':'sqrt(m²+(sigma_x Pi_x+sigma_y Pi_y)²)-m+(1/32)L_cov',
                            'Pauli_magnetic_commutator_retained':True,
                            'electric_coefficient_exact':str(ELECTRIC),'magnetic_coefficient_exact':str(MAGNETIC),
                            'common_free_speed_squared_exact':str(SPEED_SQUARED),
                            'pure_gradient_response':'unitary relative-dipole conjugation; first and second derivatives retained',
                            'threaded_probe':'external matter-link phase exp(iq theta); physical chargeQ=q/3',
                            'is_canonical_signed_Dirac_plus_Maxwell_unitary_transfer':False},
            'Gauss_exact_controls':gauss_checks(),'complete_executions':cases,
            'free_three_dimensional_refinement':refinement_records(masses),
            'scope_limits':{'vacuum_pair_creation_or_running_derived':False,
                            'interacting_Lorentz_covariant_limit_proved':False,
                            'physical_speed_or_length_calibrated':False,
                            'input_charged_mass_ratios_are_new_predictions':False},
            'input_ancestry':{'charged_masses':'frozen actual conditional charged analyzer/electron unit from openwave_common_mass_source',
                              'Wilson_coefficient':'declared1/16, zero-regulator doubler control retained'}}
