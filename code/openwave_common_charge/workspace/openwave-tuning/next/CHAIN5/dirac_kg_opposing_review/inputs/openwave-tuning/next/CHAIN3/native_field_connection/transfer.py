"""State-dependent tensor-box transfer for coherent and free harmonic modes."""
from fractions import Fraction as F
import math
from providers import bounds, total_cap
from geometry import data, D, A, B


def up(x):
    return bounds.sqrt_upper(F(x), bits=100)


def outward(x, bits=160):
    return total_cap.outward_dyadic(F(x),bits)


def parameters(macros, particles=1, initial_cells=8):
    inherited = total_cap.parameters(macros,particles,initial_cells)
    spec = data()['spectral']
    eps = F(spec['frame_Gram_error_upper'])
    # Rational Q has norm <=sqrt(1+eps). Its coefficient current defines the
    # diagonal comparison source. Every signed interval contributes.
    length = up(1+eps)*F(inherited['cumulative_Weyl_length_upper_exact'])
    current_length = macros*sum((up(x) for x in inherited['axis_current_norm_squared_max']),F())
    return inherited, length, current_length


def coexact_box(macros, cap, particles=1, initial_cells=8, extra_length=F(0)):
    if type(cap) is not int or cap < 0:
        raise ValueError('nonnegative tensor cap required')
    p,L,S = parameters(macros,particles,initial_cells)
    L+=F(extra_length)
    lam=L*L; branches=p['all_branches']*(2 if extra_length else 1)
    if not lam:
        point=tail=error2=F()
    elif cap+2 <= lam:
        point=tail=F(1); error2=F(4)
    else:
        point=lam**(cap+1)/math.factorial(cap+1)
        tail=point*F(cap+2)/(cap+2-lam)
        # Box leakage is contained in total N>R. Each of D distinct lowering
        # faces has norm sqrt(R+1), and sum|gamma_j|<=sqrt(D)||gamma||.
        error2=min(F(4),2*branches*(tail+D*L*L*(cap+1)*point))
    return {'cap':cap,'modes':D,'length_upper_exact':str(L),
            'occupation_upper_exact':str(lam), 'prefix_branches':branches,
            'point_upper_exact':str(outward(point)), 'tail_upper_exact':str(outward(tail)),
            'state_error_squared_upper_exact':str(outward(error2)),
            'state_error_upper_exact':str(up(outward(error2))),
            'physical_dimension_expression':f'({cap+1})^190',
            'binary_wires':D*cap.bit_length(),
            'flux_formula':'sqrt(B)*sqrt(190)*L*sqrt(R+1)*sqrt(point)',
            'domain':'vacuum coexact input, all routed branches, arbitrary spectator/reference'}


def harmonic_box(macros, cap, order=6, particles=1, initial_cells=8, extra_current_length=F(0)):
    if type(cap) is not int or cap < 2 or type(order) is not int or order < 2:
        raise ValueError('harmonic cap>=2 and integer moment order>=2 required')
    p,L,S=parameters(macros,particles,initial_cells)
    S+=F(extra_current_length)
    T=F(macros); w=up(p['all_branches']*(2 if extra_current_length else 1))
    # Under free shear and any prefix of the actual kicks, |u|+|v|<=1+aT
    # and |alpha|<=(1+aT)S for each auxiliary Hermite oscillator.
    linear=1+A*T; displacement=linear*S
    moment=((linear*up(2*order)+displacement)**2+1)**order
    tail=moment/F(cap+2)**order
    flux_kick=up(3)*S*up(F(cap+1,2))*tail
    flux_free=3*T*A*up((cap+1)*(cap+2))*tail/4
    error=min(F(2),w*(up(3)*tail+flux_kick+flux_free))
    return {'cap':cap,'modes':3,'moment_order':order,
            'N_plus_one_power_norm_upper_exact':str(moment),
            'per_mode_tail_norm_upper_exact':str(tail),
            'terminal_tail_upper_exact':str(w*up(3)*tail),
            'kick_boundary_flux_upper_exact':str(w*flux_kick),
            'free_boundary_flux_upper_exact':str(w*flux_free),
            'state_error_upper_exact':str(outward(error)),
            'physical_dimension_expression':f'({cap+1})^3','binary_wires':3*cap.bit_length(),
            'Hamiltonian':'P_R (a p^2/2) P_R, not (a/2)(P_R p P_R)^2',
            'initial_packet':'three Hermite vacua; Var(q)=Var(p)=1/2',
            'domain':'all finite-time sheared/displaced Gaussian branches and arbitrary reference'}


def geometry_transfer(macros, particles=1, initial_cells=8):
    """Actual K vacuum/dynamics versus frozen rational Q,Lambda oscillator."""
    p,L,S=parameters(macros,particles,initial_cells)
    spec=data()['spectral']; eta=F(spec['actual_curvature_operator_error_upper'])
    low=F(spec['minimum_eigenvalue'])-eta
    omega= bounds.sqrt_lower(A*B*low,bits=100)
    v=A/(2*omega); W=up(p['all_branches']); lam=L*L
    delta=F(spec['frame_residual_upper'])/F(spec['minimum_eigenvalue'])+F(spec['polar_Gram_error_upper'])
    vacuum=up(D*A*B*eta/(omega*omega))
    quartic=up((4*lam+D)**2+16*lam+2*D)
    dynamics=macros*W*B*eta*v*quartic/2
    currents=W*delta*S*up(v*(4*lam+1))
    return {'polar_frame_distance_upper_exact':str(delta),
            'curvature_error_upper_exact':str(eta),
            'initial_vacuum_norm_error_upper_exact':str(vacuum),
            'all_free_interval_Duhamel_upper_exact':str(dynamics),
            'all_current_coordinate_errors_upper_exact':str(currents),
            'state_error_upper_exact':str(vacuum+dynamics+currents),
            'harmonic_frame_rotation':'exact orthonormal change, isotropic packet and p^2 unchanged',
            'global_phase':'all oscillator half-frequency phases retained, no scalar subtraction'}


def choose_caps(macros, tolerance=F(1,10**6), particles=1, initial_cells=8):
    p,L,S=parameters(macros,particles,initial_cells)
    co=4*math.ceil(L*L)+128
    while F(coexact_box(macros,co,particles,initial_cells)['state_error_upper_exact'])>tolerance:
        co*=2
    harmonic=32
    while F(harmonic_box(macros,harmonic,6,particles,initial_cells)['state_error_upper_exact'])>tolerance:
        harmonic*=2
    return co,harmonic


def edge_transfer(macros,coexact_cap,harmonic_cap,time=F(1),particles=1,initial_cells=8):
    p,L,S=parameters(macros,particles,initial_cells)
    spec=data()['spectral']; eps=F(spec['frame_Gram_error_upper'])
    omega_max=up(A*B*(F(spec['maximum_eigenvalue'])+F(spec['actual_curvature_operator_error_upper'])))
    vp=omega_max/(2*A)
    extra=abs(F(time))*up((1+eps)*vp)
    co=coexact_box(macros,coexact_cap,particles,initial_cells,extra)
    harmonic=harmonic_box(macros,harmonic_cap,6,particles,initial_cells,abs(F(time)))
    geometry=geometry_transfer(macros,particles,initial_cells)
    frame=abs(F(time))*F(geometry['polar_frame_distance_upper_exact'])*up(2*p['all_branches'])*up(vp*(4*(L+extra)**2+1))
    total=F(geometry['state_error_upper_exact'])+frame+F(co['state_error_upper_exact'])+F(harmonic['state_error_upper_exact'])
    return {'edge':'uniform over all288 actual edges','Fourier_parameter_exact':str(time),
            'coexact_including_probe':co,'harmonic_including_probe':harmonic,
            'measurement_frame_error_upper_exact':str(frame),
            'source_state_and_measurement_state_error_upper_exact':str(total),
            'observable':'E_edge=P_q(edge)+sum U(edge,j)*sqrt(omega_j/(2a))*i(a†_j-a_j)+sum H(edge,k)*p_h,k',
            'readout':'Ramsey binary law after controlled exp(i t E_edge); all C/E fine records and failures retained',
            'scope':'one edge characteristic, not an inferred joint sharp-edge PVM'}


def static_harmonic_box(time,cap,order=6):
    time=F(time)
    if time<0 or cap<2 or order<2:
        raise ValueError('nonnegative time, cap>=2, order>=2 required')
    moment=(((1+A*time)*up(2*order))**2+1)**order
    tail=moment/F(cap+2)**order
    flux=3*time*A*up((cap+1)*(cap+2))*tail/4
    return {'cap':cap,'time_exact':str(time),'moment_order':order,
            'state_error_upper_exact':str(outward(min(F(2),up(3)*tail+flux))),
            'vacuum_packet':'three harmonic Hermite vacua, no current kicks',
            'retained_free_boundary_flux_exact':str(flux)}
