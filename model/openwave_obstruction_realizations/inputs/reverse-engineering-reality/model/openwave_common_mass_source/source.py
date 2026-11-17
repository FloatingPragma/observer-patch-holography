"""Native192-CAR strong/Dirac source and dimension-free moment reductions."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import numpy as np
import sympy as sy
from model.openwave_axial_state.source import parameters,sector,INPUTS
from model.openwave_baryon_weak_formfactors.native import currents,quartet,PAULI
from model.openwave_neutrino_calibration.calibration import calibration
from model.openwave_native_composites.spectrum import one_body
from model.common_source import calibrated_matter as internal

G=F(1,64)
COMPONENTS=(9,10,11,15)


def ix(x,s,c,f):
    return x*24+s*12+COMPONENTS.index(c)*3+f


@lru_cache(None)
def lepton_source():
    data=calibration()
    charged=np.array([float(x/data['charged_exact'][0]) for x in data['charged_exact']])
    neutral=data['physical_masses_eV']*1e-9/float(data['empirical_electron_anchor_GeV'])
    R=data['R'];A=np.zeros((48,48),complex);B=np.zeros_like(A);E=np.zeros(48)
    for x,s,f in product(range(2),range(2),range(3)):
        sign=1-2*s
        for c,r,mass in ((11,9,charged[f]),(10,15,neutral[f])):
            a=ix(x,s,c,f);b=ix(x,1-s,r,f);E[a]=E[b]=mass
            A[b,b]=1/np.sqrt(2);B[b,a]=sign/np.sqrt(2)
            if c==11:A[a,a]=1/np.sqrt(2);B[a,b]=-sign/np.sqrt(2)
        for alpha in range(3):
            a=ix(x,s,10,alpha);b=ix(x,1-s,15,f)
            A[a,ix(x,s,10,f)]=R[alpha,f]/np.sqrt(2)
            B[a,b]=-sign*R[alpha,f]/np.sqrt(2)
    T=np.block([[A,B],[B.conj(),A.conj()]])
    assert np.linalg.norm(T@T.conj().T-np.eye(96))<4e-15
    H=T@np.diag(np.r_[E,-E])@T.conj().T
    left=[ix(x,s,c,f) for x,s,c,f in product(range(2),range(2),(10,11),range(3))]
    rho=B@B.conj().T;pair=A@B.T
    assert np.linalg.norm(rho[np.ix_(left,left)]-.5*np.eye(24))<1e-15
    assert np.linalg.norm(pair[np.ix_(left,left)])==0
    charge=np.array([1 if c==9 else -1 if c==11 else 0
                     for x,s,c,f in product(range(2),range(2),COMPONENTS,range(3))])
    number=np.array([1 if c in (10,11) else -1
                     for x,s,c,f in product(range(2),range(2),COMPONENTS,range(3))])
    for q in (charge,number):
        Q=np.diag(np.r_[q,-q])
        assert np.linalg.norm(H@Q-Q@H)<1e-12
    native=np.array(internal.weak_raising_current(families=False)).astype(complex)
    raising=native[np.ix_([10,11],[10,11])]
    weak=((raising+raising.conj().T)/2,(raising-raising.conj().T)/(2j),
          (raising@raising.conj().T-raising.conj().T@raising)/2)
    JL=[]
    for x,mu,a in product(range(2),range(4),range(3)):
        J=np.zeros((48,48),complex)
        for s,z,w,v,f in product(range(2),range(2),range(2),range(2),range(3)):
            J[ix(x,s,(10,11)[w],f),ix(x,z,(10,11)[v],f)]=PAULI[mu][s,z]*weak[a][w,v]
        JL.append(J)
    means=np.array([np.trace(J@rho) for J in JL])
    gram=np.array([[np.trace(J@K)/4 for K in JL]for J in JL])
    assert np.linalg.norm(means)<1e-15
    assert np.linalg.norm(gram-.75*np.eye(24))<1e-15
    return dict(A=A,B=B,E=E,H=H,rho=rho,pair=pair,left=left,JL=JL,
                charged=charged,neutral=neutral,charge=charge,number=number,
                current_means=means,current_covariance=gram,calibration=data)


@lru_cache(None)
def local_norm_certificate():
    """Exact one-quark spin/weak state times the full four-mode lepton Fock space."""
    pauli=[sy.eye(2),sy.Matrix([[0,1],[1,0]]),sy.Matrix([[0,-sy.I],[sy.I,0]]),sy.diag(1,-1)]
    creators=[]
    for mode in range(4):
        c=sy.zeros(16)
        for mask in range(16):
            if not(mask>>mode)&1:c[mask|(1<<mode),mask]=(-1)**((mask&((1<<mode)-1)).bit_count())
        creators.append(c)
    V=sy.zeros(64)
    for mu,a in product(range(4),range(1,4)):
        J=sy.kronecker_product(pauli[mu],pauli[a]/2)
        L=sy.zeros(16)
        for (i,j),v in J.todok().items():L+=v*creators[i]*creators[j].H
        V+=2*(1 if mu==0 else -1)*sy.kronecker_product(J,L)
    assert V==V.H
    x=sy.Symbol('x')
    polynomial=sy.factor(V.charpoly(x).as_expr())
    expected=x**28*(x-1)**17*(x+1)**3*(x-2)**3*(x-3)*(x+3)**3*(x+6)*(x*x+2*x-9)**4
    assert sy.expand(polynomial-expected)==0
    return {'dimension':64,'characteristic_polynomial':str(polynomial),
            'exact_operator_norm':'6','three_family_norm_per_quark':'18',
            'three_quark_all_family_norm':'54'}


def baryon_moments():
    d=sector('baryon');Q=currents();W,_=quartet(d['h'])
    operators=[]
    for x,mu in product(range(2),range(4)):
        raising,third=Q['local'][x][mu]
        operators.extend(((raising+raising.conj().T)/2,(raising-raising.conj().T)/(2j),third))
    variance=3*sum(J@J for J in operators)
    assert np.linalg.norm(variance-variance.conj().T)<1e-13
    rows=[]
    unit=F(INPUTS['electron_unit_GeV'])
    for i,name in enumerate(('neutron_up','neutron_down','proton_up','proton_down')):
        v=W[:,i];energy=float(np.vdot(v,d['h']@v).real)
        residual=float(np.linalg.norm(d['h']@v-energy*v))
        assert residual<3e-9
        kappa=float(np.vdot(v,variance@v).real)
        seed_index=int(np.argmax(abs(v)**2));seed=np.eye(120)[:,seed_index]
        success=float(abs(v[seed_index])**2);failure=seed-v*np.vdot(v,seed)
        prep_work=success*energy+float(np.vdot(failure,d['h']@failure).real)-float(d['h'][seed_index,seed_index].real)
        windows=[]
        for width in (F(1,2),F(1),F(2)):
            # Numerical trial residual is retained separately from the exact eigenstate theorem.
            upper=float(G)**2*kappa/float(width)**2
            windows.append({'half_width_electron_units':str(width),
                            'half_width_GeV':str(width*unit),
                            'exact_eigenstate_outside_probability_bound_diagnostic':upper,
                            'original_inside_probability_lower_diagnostic':success*max(0,1-upper)})
        rows.append({'state':name,'strong_energy_electron_units':energy,
                     'full_source_energy_mean_GeV_diagnostic':energy*float(unit),
                     'trial_strong_eigenstate_residual':residual,
                     'native_contact_mean_exact':'0','contact_second_moment_diagnostic':kappa,
                     'contact_variance_electron_units_squared_diagnostic':float(G)**2*kappa,
                     'source_preparation_success':success,'preparation_failure':1-success,
                     'preparation_mean_work_electron_units_diagnostic':prep_work,
                     'seed_colour_label':seed_index,'seed_filter_success_vector_real':v.real.tolist(),
                     'seed_filter_success_vector_imag':v.imag.tolist(),
                     'failure_vector_real':failure.real.tolist(),'failure_vector_imag':failure.imag.tolist(),
                     'complete_energy_window_records':windows})
    return rows,variance


def exact_norm_squared(C,v):
    C=C.tocsr();total=0
    for i in range(C.shape[0]):
        real=imag=0
        for k in range(C.indptr[i],C.indptr[i+1]):
            z=C.data[k];a=round(float(z.real));b=round(float(z.imag))
            assert abs(z-a-1j*b)<1e-12
            real+=a*v[C.indices[k]];imag+=b*v[C.indices[k]]
        total+=real*real+imag*imag
    return total


def exact_baryon_window():
    from model.openwave_axial_state.certificate import ground_certificate,exact_action
    from model.openwave_calibrated_hadrons.source import weights,integral
    certificate=ground_certificate();v=list(map(int,certificate['integer_native_occupation_trial']))
    norm=sum(x*x for x in v);c=currents();d=sector('baryon');p=parameters()
    numerator=0
    for x,mu in product(range(2),range(4)):
        raising,third=c['full'][x][mu]
        for J in ((raising+raising.conj().T)/2,(raising-raising.conj().T)/(2j),third):
            numerator+=exact_norm_squared(2*J,v)
    kappa=F(3*numerator,4*norm);assert 36<kappa<37
    variance=F(certificate['trial_residual_squared_rational'])+G*G*kappa
    mean=F(certificate['trial_energy_rational'])
    maximal=max(range(len(v)),key=lambda i:abs(v[i]));success=F(v[maximal]**2,norm)
    hv=[F(0) for _ in v]
    for name,coefficient in weights(p).items():
        for i,value in enumerate(exact_action(d['parts'][name],v)):hv[i]+=coefficient*value
    if 'Cg48' in d['parts']:
        assert not any(exact_action(d['parts']['Cg48'],v))
    work=2*F(v[maximal],norm)*(v[maximal]*mean-hv[maximal])
    # Literal native weak raising generates the proton trial without changing its norm.
    raising=integral(sum(c['full'][x][0][0] for x in range(2)))
    proton=exact_action(raising,v)
    assert sum(x*x for x in proton)==norm
    proton_max=max(range(len(proton)),key=lambda i:abs(proton[i]))
    assert proton[proton_max]**2==v[maximal]**2
    proton_hv=[F(0) for _ in v]
    for name,coefficient in weights(p).items():
        for i,value in enumerate(exact_action(d['parts'][name],proton)):proton_hv[i]+=coefficient*value
    assert sum(x*y for x,y in zip(proton,proton_hv))/norm==mean-p['delta']
    assert sum((a-(mean-p['delta'])*x)**2 for a,x in zip(proton_hv,proton))/norm==F(certificate['trial_residual_squared_rational'])
    proton_numerator=0
    for x,mu in product(range(2),range(4)):
        raising,third=c['full'][x][mu]
        for J in ((raising+raising.conj().T)/2,(raising-raising.conj().T)/(2j),third):
            proton_numerator+=exact_norm_squared(2*J,proton)
    assert proton_numerator==numerator
    proton_work=2*F(proton[proton_max],norm)*(proton[proton_max]*(mean-p['delta'])-proton_hv[proton_max])
    assert proton_work==work
    unit=F(INPUTS['electron_unit_GeV'])
    outcomes=[]
    for name,energy in (('neutron',mean),('proton',mean-p['delta'])):
        for width in (F(1,2),F(1),F(2)):
            outside=variance/(width*width);assert 0<outside<1
            outcomes.append({'assigned_channel':name,'mean_energy_electron_units':str(energy),
                'mean_energy_GeV':str(energy*unit),'half_width_electron_units':str(width),
                'half_width_GeV':str(width*unit),
                'conditional_inside_probability_lower':str(1-outside),
                'original_inside_probability_lower':str(success*(1-outside)),
                'original_outside_probability_upper':str(success*outside),
                'preparation_failure_probability':str(1-success),
                'total_probability_exact':'1',
                'record_semantics':'Three records: preparation failure, energy inside, energy outside. Inside/outside probabilities are rigorous bounds, not separately calculated exact probabilities.'})
    return {'integer_trial_norm_squared':str(norm),'neutron_integer_trial':list(map(str,v)),
        'proton_integer_trial':list(map(str,proton)),
        'native_occupation_seeds':{
            'neutron':{'index':maximal,'mask':str(d['basis'][maximal])},
            'proton':{'index':proton_max,'mask':str(d['basis'][proton_max])}},
        'exact_preparation_success':str(success),'exact_preparation_mean_work_electron_units':str(work),
        'contact_mean_exact':'0','contact_second_moment_exact':str(kappa),
        'literal_proton_current_partner_mean_residual_and_second_moment_equalities_checked':True,
        'full_source_energy_variance_electron_units_squared_exact':str(variance),
        'trial_strong_variance_electron_units_squared_exact':certificate['trial_residual_squared_rational'],
        'complete_energy_window_instrument_bounds':outcomes,
        'energy_window_instrument_mean_work_exact':'0: its effects are spectral projectors of the full interacting Hamiltonian',
        'actual_trial_preparation_not_assumed_exact_stationary_ground':True}


def full_index(x,s,c,f):
    return x*96+s*48+c*3+f


def sparse_snapshot(A):
    return [[int(i),int(j),float(A[i,j].real),float(A[i,j].imag)]
            for i,j in zip(*np.nonzero(A))]


@lru_cache(None)
def full_source_terms():
    """Executable operator polynomial on the entire192-mode CAR algebra."""
    from model.coherent_defect_wave.physical_attachment.chiral import representation
    from model.common_source.chiral_candidate.action import ordering
    p=parameters();rows=ordering()['all_left_components']
    colourful=[i for i,row in enumerate(rows)if row['internal_type']in('Q','u^c','d^c')]
    Qids=[i for i,row in enumerate(rows)if row['internal_type']=='Q']
    projector=np.diag([int(i in colourful)for i in range(16)])
    qprojector=np.diag([int(i in Qids)for i in range(16)])
    def extend(A):
        out=np.zeros((16,16),complex);out[:15,:15]=np.array(A).astype(complex);return out
    generators=representation()[2]
    colours=[extend(a)for a in generators[:7]]+[extend(generators[6]+2*generators[7])]
    cm=[F(1)]*7+[F(1,3)]
    raising=extend(internal.weak_raising_current(families=False))
    weak=[(raising+raising.conj().T)/2,(raising-raising.conj().T)/(2j),
          (raising@raising.conj().T-raising.conj().T@raising)/2]
    baryon=np.diag([F(1,3)if i in Qids else -F(1,3)if i in colourful else F(0)for i in range(16)]).astype(complex)
    def embed(A,spin=np.eye(2),cells=(0,1)):
        M=np.zeros((192,192),complex)
        block=np.kron(spin,np.kron(A,np.eye(3)))
        for x in cells:M[x*96:(x+1)*96,x*96:(x+1)*96]=block
        return M
    masses=np.zeros((192,192),complex)
    for x,s,c,f in product(range(2),range(2),colourful,range(3)):
        row=rows[c];down=row['internal_type']=='d^c' or (row['internal_type']=='Q'and F(row['T3'])<0)
        masses[full_index(x,s,c,f),full_index(x,s,c,f)]=float(p['mbar']+(p['delta']/2 if down else -p['delta']/2))
        j=full_index(1-x,s,c,f);i=full_index(x,s,c,f)
        masses[i,i]+=float(p['t']);masses[i,j]-=float(p['t'])
    terms=[{'kind':'one_body','name':'colourful_number_and_hopping','weight':'1','matrix':masses}]
    def square(name,A,weight):
        terms.append({'kind':'one_body_square','name':name,'weight':str(weight),'matrix':A})
    for x in range(2):
        for a,C in enumerate(colours):square(f'local_colour_{x}_{a}',embed(C,cells=(x,)),p['U']*cm[a])
        for a in range(3):
            square(f'local_colourful_spin_{x}_{a}',embed(projector,PAULI[a+1]/2,(x,)),p['Js'])
            square(f'local_quark_weak_{x}_{a}',embed(qprojector@weak[a]@qprojector,cells=(x,)),p['Jw'])
        square(f'local_baryon_density_{x}',embed(baryon,cells=(x,)),p['A'])
    for a,C in enumerate(colours):square(f'edge_colour_{a}',embed(C),p.get('kappa_colour',F(1))*cm[a])
    for a in range(3):square(f'edge_quark_weak_{a}',embed(qprojector@weak[a]@qprojector),p['eta'])
    left=np.diag([int(i in (10,11))for i in range(16)])
    for x,mu,a in product(range(2),range(4),range(3)):
        terms.append({'kind':'one_body_product','name':f'semileptonic_{x}_{mu}_{a}',
            'weight':str(2*G*(1 if mu==0 else -1)),
            'left_matrix':embed(qprojector@weak[a]@qprojector,PAULI[mu],(x,)),
            'right_matrix':embed(left@weak[a]@left,PAULI[mu],(x,))})
    l=lepton_source()
    ids=[full_index(x,s,c,f) for x,s,c,f in product(range(2),range(2),COMPONENTS,range(3))]
    for j,E in enumerate(l['E']):
        terms.append({'kind':'positive_linear_fermion_square','name':f'Dirac_mode_{j}',
                      'energy_electron_units_diagnostic':float(E),
                      'annihilation_coefficients':[[ids[i],float(l['A'][i,j].conjugate().real),float(l['A'][i,j].conjugate().imag)]
                          for i in np.flatnonzero(l['A'][:,j])],
                      'creation_coefficients':[[ids[i],float(l['B'][i,j].real),float(l['B'][i,j].imag)]
                          for i in np.flatnonzero(l['B'][:,j])]})
    # Rebuild the interacting native baryon restriction from this full polynomial.
    d=sector('baryon')
    lookup={tuple(row['exterior_subset']):i for i,row in enumerate(rows)}
    q_modes=[full_index(x,s,lookup[(c,w+3)],0)for species,x,s,w,c in d['modes']]
    assert len(q_modes)==24
    total=None
    for term in terms:
        if term['kind']not in('one_body','one_body_square'):continue
        A=one_body(term['matrix'][np.ix_(q_modes,q_modes)],d['basis'])
        if term['kind']=='one_body_square':A=A@A
        A=float(F(term['weight']))*A
        total=A if total is None else total+A
    error=float(np.max(abs(total@d['B']-d['B']@d['h'])))
    assert error<2e-9
    snapshots=[]
    for term in terms:
        item=dict(term)
        for key in ('matrix','left_matrix','right_matrix'):
            if key in item:
                assert np.linalg.norm(item[key]-item[key].conj().T)<1e-12
                item[key]=sparse_snapshot(item[key])
        snapshots.append(item)
    return {'native_CAR_modes':192,'colourful_modes':144,'Dirac_lepton_modes':48,
            'exact_operator_definition':'Native current polynomial with stated exact parameters; floating sparse snapshots below are reproducibility diagnostics.',
            'terms':snapshots,'complete_first_family_baryon_intertwiner_error':error,
            'purely_leptonic_contact_or_transverse_photons_installed':False,
            'three_quark_spectator_family_vacuum_reducing':True}
