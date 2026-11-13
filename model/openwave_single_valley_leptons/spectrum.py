"""Native Higgs-Wilson spectrum and complete charged-lepton records."""
from fractions import Fraction
from functools import lru_cache
from itertools import product
import numpy as np
import sympy as sy
from model.openwave_full_higgs_floor.spectrum import native_tensors, exact_symbol
from model.openwave_lepton_mass_shell.spectrum import native_frames, dirac_hamiltonian
from model.common_source.chiral_candidate import kinetic
from model.lepton_mass_compiler.analyzer import rational_effect_intervals
from model.lepton_mass_compiler.compiler import electron_calibration


def masses():
    return np.array([float(sum(map(Fraction, row), Fraction(0))/2)
                     for row in rational_effect_intervals()['effect_intervals_k0_k1_k2']])


def symbol(angle, spacing=1., wilson=1.):
    angle=np.asarray(angle,float)
    if angle.shape!=(3,) or not np.isfinite(angle).all() or spacing<=0 or not np.isfinite(spacing) or wilson<0:
        raise ValueError('finite three-momentum, positive spacing and nonnegative Wilson coefficient required')
    p=np.sin(angle)/spacing
    w=wilson*np.sum(1-np.cos(angle))/spacing
    d=native_tensors()['D'];C=np.asarray(d+d.T,float)
    K=np.kron(C,np.diag(masses()+w))
    return kinetic.bdg_symbol(p,K,speed=1.),p,w


@lru_cache(None)
def minimum_theorem():
    # Independent symbolic identity proves monotonicity on each coordinate
    # between 0 and pi; reflection supplies the entire Brillouin torus.
    x=sy.symbols('x0:3',real=True);a,m=sy.symbols('a m',positive=True)
    w=sum(1-sy.cos(t) for t in x)
    energy2=(sum(sy.sin(t)**2 for t in x)+(a*m+w)**2)/a**2
    for j in range(3):
        expected=2*sy.sin(x[j])*(a*m+1+sum(1-sy.cos(x[k]) for k in range(3)if k!=j))/a**2
        assert sy.simplify(sy.diff(energy2,x[j])-expected)==0
    native=exact_symbol()
    rep=native_frames()['representation'];ids=list(rep['selected'])
    color=rep['operators']['color_C2'].extract(ids,ids);Q=native_tensors()['Q']
    colourless_charged=[i for i in range(15)if color[i,i]==0 and Q[i,i]!=0]
    assert len(colourless_charged)==2
    lepton=next(p for p in native_tensors()['pairs']if p['channel']=='charged_lepton')
    assert set(colourless_charged)=={lepton['left'],lepton['right']}
    return {
        'energy_squared':'sum sin(a*p_j)^2/a^2 + [m + sum(1-cos(a*p_j))/a]^2',
        'coordinate_derivative':'2 sin(x_j) [a*m+1+sum_(k!=j)(1-cos(x_k))]/a^2',
        'strictly_positive_derivative_for_0_to_pi':True,
        'unique_charged_branch_minimum':'p=0 on the Brillouin torus',
        'rest_Hessian':'(1/m+a) delta_ij',
        'other_corners':'at least one negative Hessian direction; no other charged local minimum',
        'exact_colourless_charged_internal_coordinates':colourless_charged,
        'charged_colourless_Dirac_families':3,
        'massive_charged_Dirac_channels_all_families':21,
        'neutral_internal_coefficients_retained':3,
        'full_native_symbol':native,
        'full_Fock_charged_floor':'H-E_reference >= min(m_k)*abs(Q)',
        'finite_range':'onsite plus nearest-neighbour kinetic and Higgs pairing',
        'old_mean_action_compiler_transferred':False,
        'chiral_gauge_continuum_or_anomaly_transfer_claimed':False,
    }


def complete_census(side=4,spacing=1.,wilson=1.):
    if side!=4:raise ValueError('complete literal execution uses the four-site periodic cube')
    values=masses();frames=[np.asarray(v,complex)for v in native_frames()['frames']]
    rows=[];ground=0.;charged_minima=[[] for _ in range(3)];maxerr=0.;intertwiner=0.
    # Exact trigonometric values avoid turning a mathematical zero into a
    # vendor-dependent sine(pi) residual. No eigenstate or momentum is removed.
    sins=(0,1,0,-1);coss=(1,0,-1,0)
    C=np.asarray(native_tensors()['D']+native_tensors()['D'].T,float)
    Q=np.kron(np.eye(2),np.kron(np.asarray(native_tensors()['Q'],float),np.eye(3)))
    charge=np.block([[Q,np.zeros_like(Q)],[np.zeros_like(Q),-Q.T]])
    for address in product(range(side),repeat=3):
        p=np.array([sins[i]/spacing for i in address]);w=wilson*sum(1-coss[i]for i in address)/spacing
        H=kinetic.bdg_symbol(p,np.kron(C,np.diag(values+w)),speed=1.)
        eig=np.linalg.eigvalsh(H);positive=np.r_[np.repeat(np.sqrt(p@p+(values+w)**2),28),np.repeat(np.sqrt(p@p),6)]
        expected=np.sort(np.r_[-positive,positive]);err=float(np.max(abs(eig-expected)))
        assert err<2e-12 and np.linalg.norm(H@charge-charge@H)<2e-12
        for k,U in enumerate(frames):
            diff=float(np.linalg.norm(H@U-U@dirac_hamiltonian(p,values[k]+w)))
            assert diff<2e-12;intertwiner=max(intertwiner,diff)
            if p@p==0 and w==0:charged_minima[k].append(list(address))
        maxerr=max(maxerr,err);ground-=float(sum(positive)/2)
        rows.append({'momentum_address':list(address),'sin_momentum':p.tolist(),'Wilson_mass_shift':w,
                     'all_180_Nambu_energies':eig.tolist(),
                     'charged_lepton_positive_energies':np.sqrt(p@p+(values+w)**2).tolist()})
    return {'side':side,'spacing':spacing,'wilson_coefficient':wilson,'physical_CAR_modes':90*side**3,
            'Nambu_rows':180*side**3,'complete_momentum_records':rows,
            'charged_minimum_addresses_by_family':charged_minima,
            'minimum_count_per_charged_family':[len(x)for x in charged_minima],
            'unshifted_ground_energy':ground,'full_eigenvalue_max_error':maxerr,
            'full_native_lepton_intertwiner_max_error':intertwiner,
            'neutral_zero_CAR_modes':48,'reference_degeneracy':'2^48 from retained neutral zeros only'}


def complete_arrival_records(character,time,side=4,spacing=1.):
    m=masses()[character];points=list(product(range(side),repeat=3));n=len(points)
    f=np.array([1.,0,0,0],complex);psis=[];angles=[]
    for address in points:
        angle=2*np.pi*np.array(address)/side;angles.append(angle)
        sins=(0,1,0,-1);coss=(1,0,-1,0)
        p=np.array([sins[i]/spacing for i in address]);w=sum(1-coss[i]for i in address)/spacing
        h=dirac_hamiltonian(p,m+w);E=np.sqrt(p@p+(m+w)**2)
        initial=(np.eye(4)+h/E)@f/(2*np.sqrt(n))
        psis.append((np.cos(time*E)*np.eye(4)-1j*np.sin(time*E)*h/E)@initial)
    real=np.array([sum((np.exp(1j*np.dot(a,x))*v for a,v in zip(angles,psis)),np.zeros(4,complex))/np.sqrt(n)for x in points])
    probs=np.abs(real)**2;success=float(probs.sum());assert abs(success-.5)<2e-12
    return {'character':character,'source_time':time,'complete_256_arrival_probabilities':probs.tolist(),
            'creation_failure_probability':.5,'original_input_total':float(probs.sum()+.5),
            'source_readout':'complete positive-band CAR addition instrument and its reservoir failure',
            'reference':'filled negative bands, vacant neutral zeros, compensating idle charge reference',
            'all_other_native_modes_retained':True}


def continuum_bound(mass,spacing,momentum_radius,time):
    # ||alpha.(sin(ap)/a-p)+beta W|| <= a^2 P^3/6 + a P^2/2.
    if min(mass,spacing)<=0 or min(momentum_radius,time)<0:raise ValueError('invalid positive parameters')
    d=spacing**2*momentum_radius**3/6+spacing*momentum_radius**2/2
    return {'mass':mass,'spacing':spacing,'momentum_radius':momentum_radius,'source_time':time,
            'operator_error_upper':d,'complete_record_TV_upper':min(1.,time*d),
            'domain':'same initial reference/instrument and momenta with Euclidean norm at most P',
            'claim':'uniform finite-time continuum Dirac comparison on this restricted momentum domain'}


def certificate():
    theorem=minimum_theorem();census=complete_census();control=complete_census(wilson=0.)
    assert census['minimum_count_per_charged_family']==[1]*3
    assert control['minimum_count_per_charged_family']==[8]*3
    cal=electron_calibration();values=masses()
    boxes=[tuple(map(Fraction,b))for b in rational_effect_intervals()['effect_intervals_k0_k1_k2']]
    mids=[sum(b,Fraction(0))/2 for b in boxes]
    anchor=Fraction('0.00051099895069');unit_exact=anchor/mids[1];unit=float(unit_exact)
    rows=[]
    for r in cal['rows']:
        k=r['character'];low,high=map(Fraction,rational_effect_intervals()['effect_intervals_k0_k1_k2'][k])
        rows.append({'species':r['species'],'character':k,'source_rest_energy':values[k],
            'GeV':float(mids[k]*unit_exact),'GeV_exact':str(mids[k]*unit_exact),
            'source_effect_interval':[str(low),str(high)],
            'analyzer_enclosure_GeV':([str(anchor),str(anchor)]if k==1 else [str(anchor*low/boxes[1][1]),str(anchor*high/boxes[1][0])]),
            'rest_energy_unchanged_by_Wilson_term':True,'particle_and_antiparticle_spin_degeneracy':4,
            'unique_colourless_charged_minimum':True,
            'comparison_role':'electron anchor'if k==1 else 'retrospective analyzer comparison'})
    corners=[]
    for count in range(1,4):
        corners.append({'pi_coordinates':count,'Wilson_shift_source':2*count,
                        'lowest_charged_corner_energy_GeV':float((values[1]+2*count)*unit),
                        'local_minimum':False})
    return {'schema':'oph.single_valley_native_leptons.v1','verdict':'PASS_COMPLETE_NATIVE_SINGLE_VALLEY_CHARGED_LEPTONS',
        'source':'all90 native CAR coordinates per cell; three Higgs channels, declared r=1 nearest-neighbour Higgs-Wilson action',
        'architecture':{'phase':'2/9, retrospective native analyzer choice','wilson_coefficient':1,
                        'Higgs_background_norm':1,'same_up_down_lepton_family_coefficients':True,
                        'unique_architecture_or_physical_species_identification_claimed':False},
        'minimum_theorem':theorem,'complete_source_spectrum':census,
        'retained_zero_Wilson_control':{'coefficient':0,'charged_minima_per_family':control['minimum_count_per_charged_family'],
                'full_spectrum_error':control['full_eigenvalue_max_error']},
        'charged_lepton_GeV':rows,'electron_anchor_GeV':float(cal['rows'][0]['electron_anchored_GeV']),
        'source_energy_unit_GeV':unit,'source_energy_unit_GeV_exact':str(unit_exact),'former_corner_modes_retained':corners,
        'complete_arrival_records':[complete_arrival_records(k,t)for k in (1,2,0)for t in (0.,.5,1.)],
        'continuum_comparisons':[continuum_bound(values[1],a,.1,1.)for a in (1.,.5,.25,.125)],
        'scope':{'physical_interacting_lepton_identification':False,'neutral_valleys_removed':False,
                 'physical_quark_mass_fit':False,'old_compiler_transferred':False,
                 'neutral_spectators_and_pairing_reference_retained':True,
                 'exact_Lorentz_covariance_at_finite_spacing':False}}
