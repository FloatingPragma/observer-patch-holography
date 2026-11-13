"""A complete signed-band rotor diagnostic of the dynamical field connection.

This is an explicit analytical diagnostic member, not a source-equivalence
claim for the full native Pauli construction or a physical particle fit.
"""
from fractions import Fraction as F
from math import sqrt, sin
import numpy as np
from model.openwave_dirac_maxwell_refinement.operators import charged_masses,ELECTRIC
from model.openwave_capacity_stress.source import sqrt_interval,enclosure

HOP=F(1,16)

def spectral_intervals(n,q,m,hop=HOP,alpha=ELECTRIC):
    m,hop,alpha=map(F,(m,hop,alpha))
    if type(n)is not int or type(q)is not int or m<=0 or alpha<=0:raise ValueError('invalid rotor domain')
    eps=sqrt_interval(m*m+hop*hop,24)
    centre=alpha*F(n*n+(n+q)**2,2)
    delta=m+alpha*F(n*n-(n+q)**2,2)
    rad=sqrt_interval(delta*delta+hop*hop,24)
    return ((centre+eps[0]-rad[1],centre+eps[1]-rad[0]),
            (centre+eps[0]+rad[0],centre+eps[1]+rad[1]))

def block(n,q,m,hop=HOP,alpha=ELECTRIC):
    m,hop,alpha=map(float,(m,hop,alpha));eps=sqrt(m*m+hop*hop)
    c=sqrt((1+m/eps)/2);s=hop/(2*eps*c)
    rotation=np.array([[c,-s],[s,c]])
    bare=np.diag([alpha*n*n+m,alpha*(n+q)**2-m])+np.array([[0,hop],[hop,0]])+eps*np.eye(2)
    # Exact all-Fourier-block form of U*(E+A)^2 U, not a clipped flux shift.
    correct=rotation.T@bare@rotation
    naive=np.diag([alpha*n*n+2*eps,alpha*(n+q)**2])
    off=alpha*c*s*((n+q)**2-n*n)
    assert np.allclose(correct[0,1],off,rtol=3e-12,atol=1e-12)
    return bare,rotation,correct,naive

def records(n,q,m,time):
    bare,u,h,naive=block(n,q,m)
    w,v=np.linalg.eigh(h)
    initial=np.array([1.,0.],complex)
    evolved=v@(np.exp(-1j*time*w)*(v.T@initial))
    probabilities=abs(evolved)**2
    before=float(np.vdot(evolved,h@evolved).real)
    after=float(np.dot(probabilities,np.diag(h)))
    initial_energy=float(h[0,0])
    assert abs(sum(probabilities)-1)<1e-12 and abs(before-initial_energy)<2e-11
    m=F(m);eps=sqrt_interval(m*m+HOP*HOP,24)
    cterm=ELECTRIC*F(n*n-(n+q)**2,2)*m
    mean=ELECTRIC*F(n*n+(n+q)**2,2)
    initial_interval=(2*eps[0]+mean+cterm/(eps[1] if cterm>=0 else eps[0]),
                      2*eps[1]+mean+cterm/(eps[0] if cterm>=0 else eps[1]))
    return {'flux_block':[n,n+q],'time':time,
            'complete_band_records':[{'band':b,'probability':float(p)} for b,p in zip(('+','-'),probabilities)],
            'all_other_integer_flux_blocks_probability_exact':'0',
            'initial_source_energy':initial_energy,'before_readout_energy':before,
            'initial_source_energy_interval_exact':list(map(str,initial_interval)),
            'floating_probabilities_and_readout_work_are_diagnostics':True,
            'after_readout_energy':after,'terminal_readout_work':after-before,
            'free_band_dephasing_is_energy_free':False,
            'missing_field_off_diagonal':float(h[0,1]),
            'naive_field_spectrum_error':float(np.max(abs(np.linalg.eigvalsh(h)-np.linalg.eigvalsh(naive))))}

def run():
    cases=[]
    for name,m in charged_masses():
        for q in (-3,3):
            window=range(-16,17)
            levels=[spectral_intervals(n,q,m)[0] for n in window]
            # All omitted n have min(|n|,|n+q|)>=14; D+eps is positive.
            outside=ELECTRIC*14**2
            ground=(min(x[0] for x in levels),min(x[1] for x in levels))
            assert ground[1]<outside
            rows=[records(n,q,m,t) for n in (-1,0,1) for t in (.5,2.)]
            for r in rows:
                r['ground_reference_interval']=enclosure(ground)
                ei=tuple(map(F,r['initial_source_energy_interval_exact']))
                r['preparation_work_interval']=enclosure((ei[0]-ground[1],ei[1]-ground[0]))
            cases.append({'family':name,'inherited_mass_parameter_exact':str(m),'charge_integer':q,
                          'full_rotor_ground_interval':enclosure(ground),
                          'outside_flux_floor_exact':str(outside),'window_block_labels':list(window),
                          'records':rows})
    return {'schema':'oph.dynamical_field_connection.v1','verdict':'PASS_FULL_ROTOR_BAND_CONNECTION',
            'source':{'carrier':'L2(U1) tensor C2, every integer flux; signed internal bands',
                      'Hamiltonian':'alpha E^2 + m sigma_z + t(cos(q theta) sigma_x+sin(q theta) sigma_y) + sqrt(m^2+t^2) I',
                      'alpha_exact':str(ELECTRIC),'hopping_exact':str(HOP),
                      'declared_diagnostic_member':True,'native_Fock_source_equivalence_proved':False},
            'exact_connection':'A=-i U* derivative_theta U=q [[s^2,cs exp(-iq theta)],[cs exp(iq theta),-s^2]]',
            'correct_transformed_field':'alpha(E+A)^2, including derivative A and A^2',
            'all_flux_blocks_retained_analytically':True,'cases':cases,
            'scope':{'new_physical_mass_prediction':False,'interacting_Lorentz_covariance':False,
                     'fermionic_pair_creation_dictionary':False,'field_connection_omission_refuted':True}}
