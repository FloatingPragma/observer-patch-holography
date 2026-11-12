"""Prepared packet bounds, unresolved dispersion and explicit GeV calibration."""
import sympy as s
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import outward,sqrt_interval,divide,mul
from model.openwave_particle_readouts.readouts import ANCHOR_GEV
from .isolation import data
from .filter import design


def calibration():
    d=data();x=s.Symbol('x')
    ec=s.Poly(x**4-2*x**3-212*x**2+213*x+11114,x)
    root=tuple(map(s.Rational,ec.intervals(eps=s.Rational(1,10**110))[0][0]))
    anchor=(root[0]-d['e0'][1],root[1]-d['e0'][0])
    scale=divide((ANCHOR_GEV,ANCHOR_GEV),anchor)
    energy=mul(d['rest'],scale)
    return {'empirical_electron_anchor_GeV':str(ANCHOR_GEV),
        'native_ec_local_energy_gap_interval':outward(anchor,160),
        'GeV_per_native_energy_interval':outward(scale,160),
        'neutral_true_ground_relative_zero_momentum_band_GeV_interval':outward(energy,180),
        'neutral_GeV_decimal_interval':[float(z) for z in energy],
        'neutral_keV_decimal_interval':[float(10**6*z) for z in energy],
        'calibration_is_supplied_empirical_energy_mapping':True,
        'neutral_energy_or_observed_neutrino_mass_used_to_select_parameters':False,
        'physical_relativistic_rest_mass_identification':False}


def certificate():
    d=data();f=design();eta2=d['eta']**2
    initial_amplitude=sqrt_interval((eta2/(1-eta2),)*2,180)[1]
    filter_trace=4*f['epsilon']/(1-eta2)
    rows=[]
    for time in (0,1,10**6,10**9,10**12):
        escape=min(s.Integer(1),(initial_amplitude+d['width']*time)**2+filter_trace/2)
        rows.append({'source_time':str(time),
            'same_atomic_position_record_lower':str(1-escape),
            'all_other_positions_and_full_profile_complement_upper':str(escape),
            'escape_decimal_upper':float(escape),
            'positive_motion_lower_bound_certified':False})
    c2=2*d['rest'][1]*d['width']+d['width']**2
    c=sqrt_interval((c2,c2),160)[1]
    # Exact eight-site measurement incompatibility, in a translated native
    # orbit basis; all native CAR phases are retained by translations.py.
    p=s.zeros(8);p[0,0]=1;k=s.ones(8)/8
    comm=p*k-k*p
    assert s.trace(comm.H*comm)==s.Rational(7,32)
    assert p*k*p==p/8
    from itertools import product
    characters=[s.Matrix([(-1)**sum(a*b for a,b in zip(momentum,point)) for point in product((0,1),repeat=3)])/s.sqrt(8)
                for momentum in product((0,1),repeat=3)]
    assert sum((v*v.T for v in characters),s.zeros(8))==s.eye(8)
    return {'localized_preparation':'known P_y/128 followed by the complete full-H band filter, with success and failure retained',
        'atomic_localized_preparation_is_actual_interacting_ground':False,
        'ideal_band_herald_success_lower':str(1-eta2),
        'finite_filter_conditional_trace_error_upper':str(filter_trace),
        'complete_momentum_readout':'all eight characters of the three periodic translations',
        'each_momentum_probability_after_localized_band_herald_interval':outward(((1-eta2)/8,1/(8*(1-eta2))),160),
        'position_momentum_noncommuting_control':{'atomic_P_y_P_k_P_y':'P_y/8',
            'commutator_Hilbert_Schmidt_squared':'7/32','simultaneous_sharp_position_and_momentum_claimed':False},
        'source_spatial_spacing':'pi','momentum_phase_coordinates':'0 or pi per translation; model wave numbers 0 or 1',
        'every_momentum_band_width_upper':str(d['width']),
        'difference_of_any_two_momentum_band_energies_upper':str(d['width']),
        'individual_dispersion_eigenvalues_resolved':False,
        'flat_atomic_dispersion_control':'at lambda=0 the complete band has exactly one energy for all eight momenta',
        'fixed_position_observable_velocity_bound':'|d<X_j>/dt| <= (pi/2) band_width for X_j=sum_y centered_coordinate_j(y) P_y, with the complement separately recorded',
        'conditional_relativistic_speed_squared_upper_if_same_band_has_mass_shell':str(c2),
        'conditional_relativistic_speed_decimal_upper':float(c),
        'unit_model_speed_relativistic_dispersion_compatible':bool(c2>=1),
        'group_velocity_or_physical_light_speed_measured':False,
        'localization_rows':rows,'calibration':calibration(),
        'scope':'The same source supplies stationary zero-momentum energy and controlled narrow-band packets. Native motion, a nonflat dispersion, relativistic spin, a physical neutrino/dark-matter dictionary and an independent physical time/length scale are not derived.'}
