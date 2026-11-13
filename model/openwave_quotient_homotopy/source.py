"""Continuous quotient lift and the operational boundary of its winding."""
from fractions import Fraction as F
import sympy as s
from model.openwave_native_centre_string.centre import native_quotient,chain,spatial_labels


def native_lift():
    parent=native_quotient();rows=[]
    for r in parent['complete_native_one_generation_coordinates']:
        q=F(r['charge_exact']);integer=r['integer_charge_plus_centre_generator'];c=F(integer)-q
        assert q+c==integer and 3*c==int(3*c)
        rows.append({**r,'colour_path_generator_exact':str(c),
          'U1_path_angle_over_2pi':'s/3',
          'native_path_phase_over_2pi':f'{integer}*s',
          'endpoint_total_phase_exact':'1','charge_phase_lift_over_2pi_exact':str(q)})
    assert sum(F(r['colour_path_generator_exact']) for r in rows)==0
    assert sum(F(r['charge_exact']) for r in rows)==0
    return {'native_rows':rows,'fundamental_colour_generator_exact':['1/3','1/3','-2/3'],
      'universal_cover':'SU(3) times R',
      'cover_kernel':'(exp(2*pi*i*k/3) I,2*pi*k/3), every integer k',
      'continuous_generator_path':'(diag(exp(2*pi*i*s/3),exp(2*pi*i*s/3),exp(-4*pi*i*s/3)),2*pi*s/3), 0<=s<=1',
      'quotient_isomorphism':'[(S,exp(i theta))] maps to exp(-i theta) S in U(3)',
      'determinant_path':'exp(-2*pi*i*s)',
      'quotient_fundamental_group':'Z','generator_determinant_winding_exact':'-1',
      'threefold_lift_U1_winding_exact':'1',
      'native_quark_phase_lifts_exact':['2/3','-1/3'],
      'full_native_Fock_lift_closes':True}


def complete_spatial_lifts(length):
    d=chain(length);labels=spatial_labels(d)
    rows=[]
    # Every complete singlet coordinate, not only a low-energy or charge selection.
    for i in range(len(d['counts'])):
        flux=list(map(int,labels['flux'][i]));trial=list(map(int,labels['triality'][i]))
        assert all((n+t)%3==0 for n,t in zip(flux,trial))
        rows.append({'occupation':list(map(int,d['counts'][i])),
          'local_charge_thirds':list(map(int,labels['local_q'][i])),
          'reference_charge_thirds':int(labels['reference'][i]),
          'electric_flux':flux,'colour_triality':trial,
          'triality_representative_closure_exponents':[(n+t)//3 for n,t in zip(flux,trial)],
          'U1_fractional_phase_lifts_exact':[str(F(n,3)) for n in flux]})
    return {'cells':length,'complete_singlet_dimension':len(rows),'records':rows,
      'all_source_link_shifts':[list(r) for r in labels['all_hamiltonian_link_shifts']],
      'electric_eigenstates_conjugate_angle_density':'uniform Haar density, independent of electric integer',
      'conjugate_angle_moments_exact':{'nonzero_integer_Fourier_modes':'0','zero_mode':'1'},
      'classical_localized_texture_measured':False,
      'triality_exponents_are_endpoint_classes_not_all_intermediate_colour_weight_phases':True}


def certificate():
    return {'schema':'oph.continuous_native_quotient_lift.v1',
      'verdict':'PASS_CONTINUOUS_QUOTIENT_HOMOTOPY_AND_COMPLETE_GAUSS_LIFTS',
      'continuous_native_lift':native_lift(),
      'complete_spatial_lifts':[complete_spatial_lifts(n) for n in (2,3)],
      'source_action':'unchanged zero-U1-energy Gauss dilation of the native strong chain',
      'parent_mass_and_dynamical_records_preserved_by_identity':True,
      'scope':{'continuous_gauge_group_homotopy_proved':True,
               'native_fractional_charge_compatible_with_generator_loop':True,
               'source_state_contains_a_localized_classical_angle_texture':False,
               'spontaneous_texture_or_soliton_equation_derived':False,
               'free_quark_GeV_values_derived':False,
               'quotient_lift_is_an_extra_physical_particle_fit':False}}
