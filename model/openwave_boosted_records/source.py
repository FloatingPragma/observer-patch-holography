"""Positive square-root kinetic on literal native Dirac excitations."""
from fractions import Fraction as F
import numpy as np
from model.openwave_common_mass_source.source import lepton_source,ix,INPUTS


def source_packet():
    l=lepton_source();T=np.vstack((l['A'],l['B'].conj()))
    assert np.linalg.norm(T.conj().T@T-np.eye(48))<3e-15
    return {'schema':'oph.common_positive_kinetic_sector.v1',
        'matter_CAR_modes_per_source_block':192,'positive_Dirac_modes_per_block':48,
        'complete_native_particle_hole_probe_rows':96,
        'kinetic':'sqrt(m_j^2+4 sum_i sin(a p_i/2)^2/a^2)-m_j',
        'masses_electron_units':l['E'].tolist(),
        'native_positive_frame_real':T.real.tolist(),'native_positive_frame_imag':T.imag.tolist(),
        'charges':l['charge'].tolist(),
        'electron_unit_GeV_exact':INPUTS['electron_unit_GeV'],
        'cases':[{'name':name,'mode':ix(0,0,c,f),'momenta_in_grid_units':ks,
                  'amplitudes_real':weights,'time':2.0,'preparation_success_exact':'1/2'}
            for name,c,f,ks,weights in (
                ('rest_electron',11,0,[0],[1]),
                ('moving_electron',11,0,[1],[1]),
                ('electron_packet',11,0,[0,1,2],[1,2,1]),
                ('positron_packet',9,0,[0,1,2],[1,2,1]),
                ('muon_packet',11,1,[0,1,2],[1,2,1]),
                ('tau_packet',11,2,[0,1,2],[1,2,1]))],
        'spatial_execution':'one periodic direction; other two momentum coordinates zero',
        'circumference_exact':'8*pi/3','sides':[16,32,64],
        'source_domain':'Product of complete longitudinal192-mode blocks, plus positive gauge-dressed Dirac hopping. Actual one-lepton/quark-vacuum reducing sector; not a numerical simulation of the interacting full Fock space.'}


def dispersion(m,p,a):
    p=np.asarray(p);return np.sqrt(m*m+np.sum((2*np.sin(a*p/2)/a)**2,axis=-1))


def boost_checks():
    rows=[];unit=F(INPUTS['electron_unit_GeV'])
    for n in (16,32,64,128):
        a=8*np.pi/(3*n);p=.75;m=1.;E=float(dispersion(m,[p],a));v=np.sin(a*p)/(a*E)
        delta=a*a*p**4/(24*m);gamma=F(5,4);beta=F(3,5)
        defect=abs(1-float(gamma)*(E-float(beta)*p))
        assert defect<=float(gamma)*delta
        velocity_bound=a*a*(p**3/(6*m)+p**5/(24*m**3))
        assert abs(v-.6)<=velocity_bound
        grid=np.linspace(-np.pi/a,np.pi/a,1001)
        velocities=np.sin(a*grid)/(a*np.sqrt(m*m+(2*np.sin(a*grid/2)/a)**2))
        assert np.max(abs(velocities))<1
        rows.append({'side':n,'spacing_electron_inverse_mass_diagnostic':a,
            'moving_energy_GeV_diagnostic':E*float(unit),
            'continuum_boost_target_energy_GeV_exact':str(F(5,4)*unit),
            'continuum_beta_exact':str(beta),'continuum_gamma_exact':str(gamma),
            'lattice_group_velocity_diagnostic':float(v),
            'group_velocity_error_upper_diagnostic':velocity_bound,
            'boost_energy_shell_defect_electron_units_diagnostic':defect,
            'boost_energy_shell_defect_upper_diagnostic':float(gamma)*delta,
            'full_regulator_frame_retained':True})
    return rows


def theorem_metadata():
    return {'positive_kinetic_full_Fock_extension':True,
        'dressed_hopping_changes_all_affected_flux_and_reference_labels':True,
        'rest_eigenvalues_and_exact_full_charge_threshold_preserved':True,
        'nonlocal_square_root_kinetic_declared':True,
        'spatial_hopping_satisfies_strict_microcausality':False,
        'photon_coupling_or_full_interacting_Lorentz_covariance':False,
        'three_dimensional_velocity_norm_bound':'|grad_p E_a|<=1 for every retained momentum',
        'all_momentum_minima':'One per massive branch in the Brillouin torus for this positive lattice-Laplacian regulator',
        'energy_error_on_momentum_ball':'0<=sqrt(m^2+|p|^2)-E_a(p)<=a^2 P^4/(24m)',
        'same_input_complete_record_TV_bound':'min(1,s |t| a^2 P^4/(24m)), with preparation success s and retained failure',
        'boost_energy_shell_defect_bound':'|E_a(p_boost)-gamma(E_a(p)-beta.p)|<=delta(p_boost)+gamma delta(p)',
        'boost_bound_scope':'Continuum boosted on-shell momentum p_boost and both momenta inside the admitted Brillouin domain; not a unitary finite-lattice Lorentz representation',
        'continuum_record_comparator':'Same finite Fourier amplitudes and native complete probe on the same periodic band subspace, replacing E_a by sqrt(m^2+p^2). No independently prepared reference transfer.'}
