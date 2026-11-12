from fractions import Fraction as F
import numpy as np
from scipy.linalg import expm
from .packet import one_axis_moments,moment_upper_bounds,certify,finite_execution,pi_interval,clock_controls
from model.openwave_lepton_mass_shell.spectrum import dirac_hamiltonian


def test_exact_boundary_moments_on_multiple_finite_volumes():
    for length in (4,5,6,9,17):
        for side in (length+5,2*(length+5)):
            g=np.zeros(side);j=np.arange(1,length+1);g[j]=np.sqrt(8/(3*(length+1)))*np.sin(np.pi*j/(length+1))**2
            D=(np.roll(np.eye(side),1,axis=0)-np.roll(np.eye(side),-1,axis=0))/(2j)
            exact=one_axis_moments(length);direct=(np.vdot(D@g,D@g).real,np.vdot(D@D@g,D@D@g).real)
            assert np.max(abs(np.array(exact)-direct))<2e-15
            bounds=moment_upper_bounds(length)
            assert exact[0]<=float(bounds['one_axis_second'])and exact[1]<=float(bounds['one_axis_fourth'])
            # Removing outside collars strictly loses the complete derivative norm.
            assert np.linalg.norm((D@g)[1:length+1])**2<direct[0]


def test_complete_native_small_box_energy_position_and_clock():
    for time in (0.,1.3,4.7):
        result=finite_execution(length=4,side=10,time=time)
        positions=np.array(result['complete_position_spinor_probabilities'])
        weights=np.array(result['complete_energy_momentum_probabilities']);energies=np.array(result['energy_momentum_values'])
        failure=result['original_preparation_failure'];assert abs(positions.sum()+failure-1)<1e-14
        assert abs(weights.sum()+failure-1)<1e-14
        assert abs(failure-result['positive_projection_failure_finite_identity'])<1e-14
        C=np.sum(weights*np.exp(-1j*time*energies))/weights.sum()
        assert abs(C-complex(*result['conditional_complex_return']))<1e-14
        assert abs(sum(r['original_probability']for r in result['complete_clock_records'])-1)<1e-14
        assert min(r['original_probability']for r in result['complete_clock_records'])>-1e-14


def test_independent_dense_source_matches_fft_evolution():
    length=4;side=10;mass=.4;time=.7
    result=finite_execution(length,side,time,mass)
    g=np.zeros(side);j=np.arange(1,length+1);g[j]=np.sqrt(8/(3*(length+1)))*np.sin(np.pi*j/(length+1))**2
    scalar=g[:,None,None]*g[None,:,None]*g[None,None,:];fhat=np.fft.fftn(scalar,norm='ortho')
    f=np.array([1,0,1,0])/np.sqrt(2);momentum=np.zeros((side,side,side,4),complex)
    for address in np.ndindex((side,side,side)):
        p=np.sin(2*np.pi*np.array(address)/side);H=dirac_hamiltonian(p,mass);E=np.sqrt(p@p+mass*mass)
        momentum[address]=expm(-1j*time*H)@(np.eye(4)+H/E)@f*fhat[address]/2
    probs=abs(np.fft.ifftn(momentum,axes=(0,1,2),norm='ortho'))**2
    assert np.max(abs(probs.reshape((-1,4))-result['complete_position_spinor_probabilities']))<1e-15


def test_large_exact_certificate_and_same_width_volume_control():
    a=certify();b=certify(side=8*(2**30+1))
    assert a['moment_upper_bounds']==b['moment_upper_bounds']
    assert F(a['mean_relative_rest_energy_excess_upper'])<F(3014,10**12)
    assert F(a['conditional_outside_cube_uniform_time_upper'])<F(1103,10**10)
    assert F(a['autonomous_clock']['visibility_uniform_lower'])>F('0.9997728')
    assert F(a['finite_full_Fock_compiler']['full_source_unhalved_trace_bound'])<=F(1,10**12)
    assert a['finite_full_Fock_compiler']['cycles']>10**52
    compiler=a['finite_full_Fock_compiler']
    assert F(a['mean_source_energy_excess_lower'])>F(compiler['conditional_on_initial_success_physical_energy_mean_transfer_error_upper'])
    assert F(compiler['conditional_on_initial_success_physical_energy_mean_transfer_error_upper'])<=F(a['mean_source_energy_excess_upper'])/10**6
    assert F(compiler['full_source_unhalved_trace_bound'])*F(compiler['physical_excitation_energy_operator_norm_upper'])==F(compiler['conditional_on_initial_success_physical_energy_mean_transfer_error_upper'])


def test_clock_drive_reference_inverse_and_coherent_herald():
    import sympy as s
    from model.openwave_kinetic_phase_cycle.compiler import annihilator
    # Literal two-mode source/reference CAR product; full source complement.
    c,d=[annihilator(2,j)for j in range(2)];A=c.H*d.H;w=s.Rational(3,5)
    K=s.diag(s.sqrt(w)*s.eye(4),A);G=K.H*K
    failure=s.diag(*[s.sqrt(1-G[j,j])for j in range(8)])
    assert G+failure.H*failure==s.eye(8)
    assert A*A==s.zeros(4)
    c0,c1,ref=[annihilator(3,j)for j in range(3)]
    occupied=s.zeros(8,1);occupied[1]=s.sqrt(1-w);occupied[2]=s.sqrt(w)
    A3=c0.H*ref.H;assert (occupied.H*A3.H*A3*occupied)[0]==w
    K3=s.diag(s.sqrt(w)*s.eye(8),A3);input_state=s.Matrix.vstack(occupied,occupied)/s.sqrt(2)
    assert (input_state.H*K3.H*K3*input_state)[0]==w
    rows=clock_controls()['independent_finite_two_line_controls'];assert rows[-1]['dephased_preparation_XY_contrast']==[0.,0.]
    assert abs(complex(*rows[-1]['same_source_then_actual_inverse_return'])-1)<1e-15
    C=complex(*rows[-1]['free_return']);driven=complex(*rows[-1]['differential_native_charge_drive_return'])
    assert abs(driven/C-np.exp(3j))<1e-14


def test_moment_bounds_imply_persistence_on_an_independent_packet():
    length=12;side=18;mass=3.;time=2.
    result=finite_execution(length,side,time,mass);mom=moment_upper_bounds(length)
    eta=float(mom['kinetic_second'])/(4*mass*mass);w=1-eta
    b=np.sqrt(float(mom['kinetic_fourth'])/(4*mass*mass*w))
    conditional=result['original_probability_outside_cube']/(1-result['original_preparation_failure'])
    assert conditional<=(np.sqrt(eta/w)+time*b)**2
    energy=result['conditional_mean_energy']
    assert 0<=energy-mass<=float(mom['kinetic_second'])/(2*mass*w)


def test_same_packet_phase_and_frequency_bounds():
    m=3.;time=2.;result=finite_execution(12,18,time,m)
    weights=np.array(result['complete_energy_momentum_probabilities']);weights/=weights.sum()
    E=np.array(result['energy_momentum_values']);C=np.sum(weights*np.exp(-1j*time*E))
    d=float(moment_upper_bounds(12)['kinetic_second'])/(2*m*(1-float(moment_upper_bounds(12)['kinetic_second'])/(4*m*m)))
    derivative=np.sum(-1j*E*weights*np.exp(-1j*time*E))
    frequency=-(derivative/C).imag
    assert abs(C-np.exp(-1j*m*time))<=time*d
    assert abs(frequency-m)<=d/(1-time*d)
