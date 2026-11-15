"""Gauge-invariant quantum Yukawa memory and its declared background response."""
from functools import lru_cache
from collections import Counter
import sympy as s
from ..physical_attachment.chiral import representation

EPS=s.ImmutableMatrix([[0,1],[-1,0]])
PAULI=(s.ImmutableMatrix([[0,1],[1,0]]),s.ImmutableMatrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))


def clean(m):
    return s.ImmutableSparseMatrix(m.applyfunc(s.simplify))


@lru_cache(maxsize=1)
def higgs_memory():
    """Total scalar occupation <=1: vacuum, two particles, two antiparticles."""
    fields=[]
    for k in range(2):
        phi=s.zeros(5);phi[0,1+k]=1;phi[3+k,0]=1
        fields.append(s.ImmutableSparseMatrix(phi))
    _,g5,_,_=representation()
    currents=tuple(s.diag(s.zeros(1),t[3:,3:],-t[3:,3:].T) for t in g5)
    n=s.diag(0,1,1,1,1)
    checks=0
    for q,t in zip(currents,g5,strict=True):
        t=t[3:,3:]
        for k,phi in enumerate(fields):
            assert clean(q*phi-phi*q+sum((t[k,j]*fields[j] for j in range(2)),s.zeros(5)))==s.zeros(5)
            checks+=1
    return tuple(fields),currents,n,checks


def pair_matrix(channel,component):
    vertices=representation()[3]
    return clean(s.kronecker_product(EPS,vertices[channel]['components'][component]))


def neutral_matrix(yu,yd,ye,v):
    values=tuple(map(s.sympify,(yu,yd,ye,v)))
    if any(x.has(s.Float) for x in values):raise ValueError('exact coefficients required')
    yu,yd,ye,v=values;vertices=representation()[3]
    return clean(v*(yu*vertices['up']['components'][1]+yd*vertices['down']['components'][1]
                    +ye*vertices['charged_lepton']['components'][1]))


def weyl_symbol(momentum):
    p=tuple(map(s.sympify,momentum))
    if len(p)!=3 or any(x.has(s.Float) for x in p):raise ValueError('three exact momentum components required')
    return clean(-sum((x*t for x,t in zip(p,PAULI,strict=True)),s.zeros(2)))


def bdg_symbol(momentum,yu,yd,ye,v):
    """Nambu (c(p),c^dagger(-p)); duplication is an operator identity."""
    p=tuple(map(s.sympify,momentum));m=neutral_matrix(yu,yd,ye,v)
    delta=s.kronecker_product(EPS,m)
    h=s.kronecker_product(weyl_symbol(p),s.eye(15))
    hm=s.kronecker_product(weyl_symbol(tuple(-x for x in p)),s.eye(15))
    return clean(h.row_join(-delta.conjugate()).col_join(delta.row_join(-hm.T)))


@lru_cache(maxsize=1)
def tensor_checks():
    _,g5,currents,vertices=representation();fields,higgs_q,_,field_checks=higgs_memory()
    ward=0;antisym=0
    for name,vertex in vertices.items():
        ds=tuple(pair_matrix(name,k) for k in range(2))
        for d in ds:
            assert d.T==-d
            antisym+=1
        for t,g in zip(currents,g5,strict=True):
            matter=s.kronecker_product(s.eye(2),t)
            scalar=g[3:,3:]
            if vertex['conjugate_Higgs']:scalar=-scalar.T
            for k,d in enumerate(ds):
                residual=matter.T*d+d*matter+sum((scalar[j,k]*ds[j] for j in range(2)),s.zeros(30))
                assert clean(residual)==s.zeros(30)
                ward+=1
    spin_checks=0
    # SL(2,C) singlet checks include nonunitary boosts, not only rotations.
    for b in [s.diag(2,s.Rational(1,2)),s.Matrix([[s.Rational(5,3),s.Rational(4,3)],[s.Rational(4,3),s.Rational(5,3)]]),
              s.Matrix([[s.Rational(3,5),s.Rational(4,5)],[-s.Rational(4,5),s.Rational(3,5)]])]:
        assert b.det()==1 and clean(b.T*EPS*b)==EPS
        spin_checks+=1
    # A gauge-invariant finite reference has zero charged Higgs expectation.
    hyper=higgs_q[-1]
    for phi in fields:assert clean(hyper*phi-phi*hyper+s.Rational(1,2)*phi)==s.zeros(5)
    return {'antisymmetric_Weyl_pair_tensors':antisym,'full_twelve_current_pair_Ward_checks':ward,
            'quantum_scalar_commutator_checks':field_checks,'SL2C_spin_singlet_checks':spin_checks,
            'scalar_memory_dimension':5,'scalar_memory_total_number_cutoff':1,
            'canonical_scalar_CCR_at_cutoff':False,'gauge_invariant_reference_Higgs_expectation':'0'}


@lru_cache(maxsize=6)
def annihilator(index):
    entries={}
    for mask in range(64):
        if mask & (1<<index):entries[mask^(1<<index),mask]=(-1)**((mask & ((1<<index)-1)).bit_count())
    return s.ImmutableSparseMatrix(64,64,entries)


@lru_cache(maxsize=1)
def quantum_lepton_block():
    """Complete 6-mode lepton CAR space times the quantum doublet memory."""
    vertices=representation()[3];indices=[9,10,11]
    cs=[annihilator(i) for i in range(6)];fields,hq,number,_=higgs_memory()
    aa=[]
    for component in vertices['charged_lepton']['components']:
        c=component.extract(indices,indices)
        pair=sum((c[i,j]*cs[i]*cs[3+j] for i in range(3) for j in range(3)),s.zeros(64))
        aa.append(clean(pair))
    # Both Higgs components and the full particle/antiparticle memory are retained.
    interaction=sum((s.kronecker_product(aa[k],fields[k].H)+s.kronecker_product(aa[k].H,fields[k])
                     for k in range(2)),s.zeros(320))
    h=clean(s.kronecker_product(s.eye(64),number)+interaction)
    charges=[]
    for t,q in zip(representation()[2],hq,strict=True):
        tl=t.extract(indices,indices)
        matter=sum((tl[i,j]*cs[spin*3+i].H*cs[spin*3+j] for spin in range(2)
                    for i in range(3) for j in range(3)),s.zeros(64))
        charges.append(clean(s.kronecker_product(matter,s.eye(5))+s.kronecker_product(s.eye(64),q)))
    assert h==h.H
    ward=0
    for q in charges:
        assert clean(q*h-h*q)==s.zeros(320);ward+=1
    # Graph closure is checked in the full 320-state matrix, not assumed.
    adjacency={i:set() for i in range(320)}
    for (i,j),value in h.todok().items():
        if value:adjacency[j].add(i)
    seen={0};front={0}
    while front:
        more=set().union(*(adjacency[j] for j in front))-seen
        seen|=more;front=more
    support=sorted(seen);block=h.extract(support,support)
    assert all(i in seen for j in seen for i in adjacency[j])
    return h,charges,support,clean(block),ward


def quantum_response():
    h,charges,support,block,ward=quantum_lepton_block()
    return {'full_dimension':320,'complete_native_lepton_CAR_modes':6,
            'spectator_quark_CAR_modes':24,'full_R15_and_scalar_dimension':5*2**30,
            'every_active_energy_multiplicity_spectator_factor':2**24,
            'all_gauge_current_commutators_checked':ward,'cyclic_block_dimension':len(support),
            'cyclic_source_basis':support,'cyclic_H':[[str(x) for x in row] for row in block.tolist()],
            'cyclic_characteristic_polynomial':str(block.charpoly().as_expr()),
            'source_vacuum_is_stationary':False}


def background_checks():
    yu,yd,ye,v=s.symbols('y_u y_d y_e v',positive=True)
    m=neutral_matrix(yu,yd,ye,v);square=clean(m.H*m)
    assert square==s.diag(*square.diagonal())
    spectrum=dict(Counter(map(str,square.diagonal())))
    p=s.symbols('p_x p_y p_z',real=True)
    b=bdg_symbol(p,yu,yd,ye,v)
    assert b==b.H
    target=s.kronecker_product(s.eye(4),square)+(sum(x*x for x in p))*s.eye(60)
    assert clean(b*b-target)==s.zeros(60)
    # Nambu redundancy; this does not introduce 30 independent extra Weyl modes.
    tau=s.zeros(60)
    for i in range(30):tau[i,30+i]=tau[30+i,i]=1
    assert clean(tau*b.conjugate()*tau+bdg_symbol(tuple(-x for x in p),yu,yd,ye,v))==s.zeros(60)
    values=bdg_symbol((0,0,0),s.Integer(1),s.Integer(2),s.Integer(3),s.Integer(1)).eigenvals()
    return {'canonical_Weyl_modes_per_site':30,'Nambu_matrix_dimension':60,
            'Nambu_particle_hole_redundancy_checked':True,
            'symbol_square_identity':'B(p)^2=p^2 I+I4 tensor M^dagger M (real coefficients)',
            'M_dagger_M_spectrum':spectrum,'independent_complex_Yukawa_channels':3,
            'neutral_background_selected_by_source':False,
            'declared_test_coefficients':{'y_u':'1','y_d':'2','y_e':'3','v':'1'},
            'analytic_zero_momentum_symbol_eigenvalues':{str(k):int(n) for k,n in values.items()},
            'massless_minimal_neutrino_Weyl_species':1,'massive_Dirac_color_components':{'up':3,'down':3,'charged_lepton':1},
            'source_selected_measured_masses':False}


@lru_cache(maxsize=1)
def full_quantum_spectrum():
    """Every 320-state energy, via the exact 64-by-64 singular-value problem."""
    h,qs,support,block,ward=quantum_lepton_block()
    vacuum=list(range(0,320,5));excited=[j for j in range(320) if j%5]
    w=h.extract(excited,vacuum);b=clean(w.H*w)
    assert h.extract(vacuum,vacuum)==s.zeros(64)
    assert h.extract(excited,excited)==s.eye(256)
    projectors=[]
    for r in range(7):
        q=s.eye(64)
        for k in range(7):
            if k!=r:q=q*(b-k*s.eye(64))/(r-k)
        q=clean(q);assert clean(q*q-q)==s.zeros(64)
        assert clean(b*q-r*q)==s.zeros(64)
        projectors.append(q)
    assert sum(projectors,s.zeros(64))==s.eye(64)
    multiplicities=[int(s.trace(q)) for q in projectors]
    assert multiplicities==[6,20,22,8,2,4,2]
    f=s.zeros(320,64)
    for i,x in enumerate(vacuum):f[x,i]=1
    for (i,j),value in w.todok().items():f[excited[i],j]=-value/3
    pg=clean(s.Rational(3,5)*f*projectors[6]*f.H)
    assert pg.H==pg and clean(pg*pg-pg)==s.zeros(320) and s.trace(pg)==2
    assert clean(h*pg+2*pg)==s.zeros(320)
    assert all(clean(q*pg)==s.zeros(320) for q in qs)
    parity=s.diag(*[(-1)**(j//5).bit_count() for j in range(320)])
    assert clean(parity*pg+pg)==s.zeros(320)
    rho=pg/2
    fields,_,_,_=higgs_memory()
    assert all(s.trace(s.kronecker_product(s.eye(64),phi)*rho)==0 for phi in fields)
    return h,w,tuple(projectors),rho,vacuum,excited


def spectral_weights(density):
    _,w,qs,_,vacuum,excited=full_quantum_spectrum()
    vv=density.extract(vacuum,vacuum);vs=density.extract(vacuum,excited)
    sv=density.extract(excited,vacuum);ss=density.extract(excited,excited)
    wssw=clean(w.H*ss*w);cross=clean(vs*w+w.H*sv)
    result=[(s.Integer(0),s.simplify(s.trace(qs[0]*vv)),6)]
    plus_inverse=sum((qs[r]/r for r in range(1,7)),s.zeros(64))
    result.append((s.Integer(1),s.simplify(s.trace(ss)-s.trace(plus_inverse*wssw)),198))
    for r in range(1,7):
        for sign in (-1,1):
            energy=(1+sign*s.sqrt(1+4*r))/2;d=energy-1
            value=s.simplify((s.trace(qs[r]*vv)+s.trace(qs[r]*cross)/d+s.trace(qs[r]*wssw)/d**2)/(1+r/d**2))
            result.append((energy,value,int(s.trace(qs[r]))))
    assert s.simplify(sum(x[1] for x in result)-s.trace(density))==0
    assert sum(x[2] for x in result)==320
    return result


def quantum_ground_and_readout():
    h,w,qs,rho,_,_=full_quantum_spectrum()
    source=s.kronecker_product(annihilator(0),s.eye(5));n=source.H*source
    assert s.simplify(s.trace(n*rho))==s.Rational(1,2)
    alternatives=[('add',0,source.H),('add',1,-n),('remove',0,s.eye(320)-n),('remove',1,source)]
    assert sum((op.H*op/2 for _,_,op in alternatives),s.zeros(320))==s.eye(320)
    records=[]
    for setting,flag,op in alternatives:
        output=clean(op*rho*op.H/2)
        for energy,prob,rank in spectral_weights(output):
            assert prob.is_nonnegative is True
            records.append({'setting':setting,'reservoir_out':flag,'source_energy':str(energy),
                            'excitation_energy_above_full_ground':str(s.simplify(energy+2)),
                            'energy_multiplicity':rank,'probability':str(prob)})
    assert s.simplify(sum(s.sympify(r['probability']) for r in records))==1
    # Gauge-singlet, zero lepton-charge sector is a distinct reference contract.
    lepton=[sum((-1 if i%3==0 else 1) for i in range(6) if mask>>i&1) for mask in range(64)]
    zero=[i for i in range(320) if lepton[i//5]==0]
    _,currents,_,_,_=quantum_lepton_block()
    constraints=s.Matrix.vstack(*(currents[j].extract(zero,zero) for j in (8,9,10,11)))
    basis=s.Matrix.hstack(*constraints.nullspace());gram=basis.H*basis
    restricted=basis.H*h.extract(zero,zero)*basis
    e=s.Symbol('e');characteristic=s.factor((restricted-e*gram).det()/gram.det())
    assert basis.shape[1]==5 and characteristic==-(e-1)**3*(e**2-e-4)
    # Vacuum cyclic execution: t=pi/sqrt17 gives a complete 320-outcome PVM law.
    _,_,support,block,_=quantum_lepton_block();rad=s.sqrt(17)
    initial=s.zeros(5,1);initial[0]=1
    centered=block-s.eye(5)/2
    assert clean(centered*centered*initial-s.Rational(17,4)*initial)==s.zeros(5,1)
    # On the vacuum cyclic two-dimensional space, exp(-it H) at t=pi/sqrt17
    # equals a global phase times -2i(H-I/2)/sqrt17.
    final=clean(-2*s.I*centered*initial/rad)
    assert final[0]==s.I/rad and all(final[k]==2*s.I/rad for k in range(1,5))
    assert s.simplify((final.H*final)[0])==1
    probabilities=[s.Rational(1,17) if j==0 else s.Rational(4,17) if j in support else s.Integer(0) for j in range(320)]
    assert sum(probabilities)==1
    phi=higgs_memory()[0][0];chi=s.Matrix([1,1,0,0,0])/s.sqrt(2)
    mean=(chi.H*phi*chi)[0];variance=s.simplify((chi.H*phi.H*phi*chi)[0]-s.conjugate(mean)*mean)
    assert mean==s.Rational(1,2) and variance==s.Rational(3,4)
    return {'full_energy_multiplicities':{str(e):rank for e,_,rank in spectral_weights(rho)},
            'full_ground_energy':'-2','full_ground_rank':2,'full_ground_reference':'normalized exact ground projector',
            'full_ground_is_gauge_singlet':True,'full_ground_fermion_parity':'odd',
            'full_ground_excitation_gap':'(5-sqrt(21))/2','full_reference_is_physical_Lorentz_vacuum':False,
            'complete_driven_reservoir_energy_records':records,'number_of_complete_records':len(records),
            'reference_preparation_service':'complete energy PVM on I320/320; no dynamical vacuum selection',
            'ground_preparation_success_probability':'1/160',
            'complete_reference_preparation_records':[{'source_energy':str(e),'probability':str(s.Rational(rank,320)),
                    'ground_success':e==-2} for e,_,rank in spectral_weights(rho)],
            'driven_apparatus_work_computed':False,'reservoir_carried_hypercharge':'1',
            'neutral_gauge_singlet_reference':{'lepton_charge':0,'sector_dimension':5,
                'ground_energy':'(1-sqrt(17))/2','ground_rank':1,'ground_is_full_H_ground':False,
                'characteristic_polynomial':str(characteristic)},
            'vacuum_input_interaction':{'time':'pi/sqrt(17)','complete_basis_probabilities':list(map(str,probabilities)),
                'failure_or_no_conversion_probability':'1/17','pair_plus_scalar_probability':'16/17',
                'ordinary_zero_fermion_state_stationary':False},
            'nonzero_mean_not_classical_field_control':{'scalar_state':'(vacuum+particle0)/sqrt2',
                'field_mean':str(mean),'field_variance':str(variance)},
            'onsite_gap_is_relativistic_particle_mass':False}


def certificate():
    return {'schema':'oph.coherent_defect_wave.chiral_higgs.v1',
            'source':'native-chiral-Weyl-R15-quantum-Higgs-v1',
            'source_spin_order':'spin*15+internal','all_left_Weyl_kinetic_symbol':'-sigma.p',
            'tensor_and_memory':tensor_checks(),'neutral_background':background_checks(),
            'quantum_interaction':quantum_response(),'quantum_spectrum':quantum_ground_and_readout(),
            'scope':{'exact_finite_quantum_Yukawa_interaction':True,
                     'actual_R15_source_current_tensors':True,'full_local_lepton_CAR_and_scalar_domain':True,
                     'conditional_broken_phase_BdG_mass_relation':True,
                     'background_Higgs_state_selected':False,'Yukawa_coefficients_selected':False,
                     'quantum_Higgs_continuum_constructed':False,'interacting_chiral_continuum':False,
                     'measured_species_identified':False,'GeV_masses':False}}
