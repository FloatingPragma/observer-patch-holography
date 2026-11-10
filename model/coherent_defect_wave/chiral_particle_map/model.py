"""Exact broken-phase field dictionary for one declared native Weyl family."""
from functools import lru_cache
import sympy as s

from ..physical_attachment.chiral import representation


def clean(matrix):
    return s.ImmutableSparseMatrix(matrix.applyfunc(s.expand))


def pauli():
    return (s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))


@lru_cache(maxsize=1)
def dictionary():
    subsets,_,currents,vertices=representation()
    y=currents[-1];t3=currents[-2];charge=y+t3
    rows=[]
    for channel in ('up','down','charged_lepton'):
        c=vertices[channel]['components'][1]
        for left in range(15):
            if t3[left,left]==0:continue
            for right in range(15):
                if c[left,right]==0:continue
                assert t3[right,right]==0 and abs(c[left,right])==1
                assert charge[left,left]==-charge[right,right]
                rows.append({'channel':channel,'left':left,'conjugate_right':right,
                             'sign':c[left,right], 'charge':charge[left,left],
                             'left_subset':subsets[left],'right_subset':subsets[right]})
    paired={i for row in rows for i in (row['left'],row['conjugate_right'])}
    neutrino=next(i for i in range(15) if i not in paired)
    assert len(rows)==7 and len(paired)==14 and charge[neutrino,neutrino]==0
    assert y[neutrino,neutrino]==-s.Rational(1,2)
    epsilon=s.Matrix([[0,1],[-1,0]])
    # psi_phys=A c+B c^dagger. Output: seven four-component Dirac fields,
    # then the two left-handed neutral components. This is a CAR automorphism.
    A=s.zeros(30);B=s.zeros(30)
    for row_index,row in enumerate(rows):
        for spin in range(2):
            A[4*row_index+spin,15*spin+row['left']]=1
            for other in range(2):
                B[4*row_index+2+spin,15*other+row['conjugate_right']]=-row['sign']*epsilon[spin,other]
    for spin in range(2):A[28+spin,15*spin+neutrino]=1
    assert A*A.H+B*B.H==s.eye(30)
    assert A*B.T+B*A.T==s.zeros(30)
    assert A.H*A+B.T*B.conjugate()==s.eye(30)
    S=A.row_join(B)
    assert S*S.H==s.eye(30)
    return rows,neutrino,clean(A),clean(B),clean(S)


def physical_currents():
    _,_,currents,_=representation();rows,neutrino,A,B,S=dictionary()
    output=[]
    for t in currents:
        G=s.kronecker_product(s.eye(2),t)
        # The source annihilator/current convention fixes the minus transpose
        # on charge-conjugate right fields; constant terms are retained below.
        out=clean(A*G*A.H-B*G.T*B.H)
        assert clean(A*G*B.T-B*G.T*A.T)==s.zeros(30)
        output.append(out)
    q=output[-1]+output[-2]
    expected=s.diag(*[row['charge'] for row in rows for _ in range(4)],0,0)
    assert q==expected
    left=s.diag(*[v for _ in rows for v in (1,1,0,0)],1,1)
    weak=output[8:11]
    assert all(t==left*t*left for t in weak)
    # Color acts equally on the physical L/R halves after conjugation.
    for t in output[:8]:
        for i in range(7):
            for j in range(7):
                assert t.extract([4*i,4*i+1],[4*j,4*j+1])==t.extract([4*i+2,4*i+3],[4*j+2,4*j+3])
    constants=[]
    for t in currents:
        G=s.kronecker_product(s.eye(2),t)
        constants.append(s.trace(G*B.T*B.conjugate()))
    assert constants==[0]*12
    # Charged weak current makes three up/down and one neutral/lepton pair.
    raising=clean(output[8]+s.I*output[9])
    entries=[(i,j,raising[i,j]) for i in range(30) for j in range(30) if raising[i,j]]
    assert len(entries)==8 and all(abs(value)==1 for _,_,value in entries)
    return output,{'physical_electric_charge_diagonal':[str(x) for x in q.diagonal()],
                   'weak_generators_supported_only_on_left_components':True,
                   'color_generators_equal_on_physical_left_and_right':True,
                   'charged_current_nonzero_spin_resolved_entries':[[i,j,str(v)] for i,j,v in entries],
                   'current_reordering_constants':[str(x) for x in constants]}


def spectrum_identity():
    rows,nu,A,B,S=dictionary();_,_,_,vertices=representation()
    px,py,pz=s.symbols('p_x p_y p_z',real=True)
    yu,yd,ye,v=s.symbols('y_u y_d y_e v',positive=True)
    parameters={'up':yu,'down':yd,'charged_lepton':ye}
    M=sum((v*parameters[name]*vertices[name]['components'][1] for name in parameters),s.zeros(15))
    epsilon=s.Matrix([[0,1],[-1,0]]);delta=s.kronecker_product(epsilon,M)
    h=-s.kronecker_product(sum((p*t for p,t in zip((px,py,pz),pauli())),s.zeros(2)),s.eye(15))
    bdg=h.row_join(-delta.conjugate()).col_join(delta.row_join(h.T))
    assert bdg==bdg.H
    physical=clean(S*bdg*S.H)
    assert clean(S*bdg-physical*S)==s.zeros(30,60)
    expected=s.zeros(30)
    two=-sum((p*t for p,t in zip((px,py,pz),pauli())),s.zeros(2))
    for i,row in enumerate(rows):
        m=v*parameters[row['channel']]
        block=two.row_join(m*s.eye(2)).col_join((m*s.eye(2)).row_join(-two))
        expected[4*i:4*i+4,4*i:4*i+4]=block
    expected[28:30,28:30]=two
    assert physical==expected
    momentum=px**2+py**2+pz**2
    assert clean(physical*physical)==s.diag(*[momentum+(v*parameters[row['channel']])**2 for row in rows for _ in range(4)],momentum,momentum)
    return {'Dirac_fields_including_color':7,'massless_left_Weyl_fields':1,
            'physical_Cauchy_rows':30,'redundant_Nambu_rows':60,
            'Nambu_rows_are_independent_mirror_species':False,
            'exact_source_current_and_BdG_intertwiner':True,
            'background_mass_energy':{'up':'abs(v*y_u)','down':'abs(v*y_d)','charged_lepton':'abs(v*y_e)','neutral':'0'},
            'squared_branch_energy':'|p|^2+abs(v*y_channel)^2',
            'background_is_a_quantum_Higgs_vacuum_theorem':False}


def fock_pair_check():
    # Exhaustive four-mode realization, independent of the Nambu calculation.
    def annihilator(mode):
        c=s.zeros(16)
        for mask in range(16):
            if mask&(1<<mode):c[mask^(1<<mode),mask]=(-1)**((mask&((1<<mode)-1)).bit_count())
        return c
    c=[annihilator(i) for i in range(4)]
    m=s.Symbol('m',real=True)
    pair=c[0]*c[3]-c[1]*c[2]
    checks=0
    for sign in (-1,1):
        b=[c[0],c[1],-sign*c[3].H,sign*c[2].H]
        native=sign*m*(pair+pair.H)
        physical=m*(b[0].H*b[2]+b[1].H*b[3]+b[2].H*b[0]+b[3].H*b[1])
        assert clean(native-physical)==s.zeros(16)
        for i in range(4):
            for j in range(4):
                assert b[i]*b[j].H+b[j].H*b[i]==(s.eye(16) if i==j else s.zeros(16))
                checks+=1
    return {'full_Fock_dimension':16,'CAR_operator_checks':checks,'both_source_pair_signs_match_Dirac_mass_operator':True}


def shared_operator_checks():
    from ..chiral_source.model import currents, native_basis
    from ..chiral_regulator.model import symbol, spectral_momenta
    from ..chiral_higgs.model import bdg_symbol, neutral_matrix
    subsets,_,native,_=representation()
    assert tuple(native_basis())==tuple(subsets)
    assert all(a==s.kronecker_product(s.eye(2),b) for a,b in zip(currents(),native,strict=True))
    rows,nu,A,B,S=dictionary()
    yu,yd,ye,v=map(s.Integer,(2,3,5,1))
    points=[(s.Rational(1,2),)*3,(s.Rational(-3,2),s.Rational(1,2),s.Rational(-1,2))]
    menu=spectral_momenta(4)
    for momentum in points:
        assert all(p in menu and -p in menu for p in momentum)
        actual=bdg_symbol(momentum,yu,yd,ye,v)
        h=clean(S*actual*S.H)
        assert clean(S*actual-h*S)==s.zeros(30,60)
        params={'up':yu,'down':yd,'charged_lepton':ye}
        for i,row in enumerate(rows):
            two=symbol(momentum);mass=params[row['channel']]*v
            block=two.row_join(mass*s.eye(2)).col_join((mass*s.eye(2)).row_join(-two))
            assert h[4*i:4*i+4,4*i:4*i+4]==block
        assert h[28:30,28:30]==symbol(momentum)
    return {'same_native_basis_and_all_twelve_source_currents':True,
            'actual_regulator_momenta_and_Higgs_BdG_points_checked':len(points),
            'quantum_gauge_links_replaced_by_identity':False,
            'domain':'specified neutral-background free kinetic benchmark; quantum link/scalar reference is separate'}


def certificate():
    rows,nu,*_=dictionary();_,current_result=physical_currents()
    labels={'up':'up-type quark field','down':'down-type quark field','charged_lepton':'charged-lepton field'}
    mapped=[]
    for i,row in enumerate(rows):
        mapped.append({'id':f'native-chiral-pair-{i}',
             'source_left_subset':list(row['left_subset']),
             'source_conjugate_right_subset':list(row['right_subset']),
             'source_Yukawa_sign':str(row['sign']),
             'conditional_role':labels[row['channel']],
             'spin':'1/2','electric_charge':str(row['charge']),
             'antiparticle_charge':str(-row['charge']),
             'color':'triplet' if row['channel']!='charged_lepton' else 'singlet',
             'mass_energy_source_units':{'up':'abs(v*y_u)','down':'abs(v*y_d)','charged_lepton':'abs(v*y_e)'}[row['channel']],
             'empirically_identified_species':None})
    subsets,_,_,_=representation()
    mapped.append({'id':'native-chiral-neutral','source_left_subset':list(subsets[nu]),
        'source_conjugate_right_subset':None,'conditional_role':'left neutral-lepton field',
        'spin':'1/2','electric_charge':'0','antiparticle_charge':'0','color':'singlet',
        'mass_energy_source_units':'0','empirically_identified_species':None})
    return {'schema':'oph.coherent_defect_wave.chiral_particle_map.v1',
            'shared_operators':shared_operator_checks(),'conditional_field_rows':mapped,'physical_currents':current_result,
            'background_spectrum':spectrum_identity(),'full_Fock_pair_control':fock_pair_check(),
            'parameter_controls':{'v1_yu1_yd1_ye1':{'up':'1','down':'1','charged_lepton':'1'},
                                  'v1_yu2_yd3_ye5':{'up':'2','down':'3','charged_lepton':'5'},
                                  'both_satisfy_same_charge_and_chirality_constraints':True},
            'identification':{'family_count':1,'generation_selected':False,
               'massive_neutrino_sector_present':False,'measured_GeV_predictions':{},
               'calibration_data_used':False,'identity_probability_assigned':None,
               'reason_no_probability':'No likelihood, prior or experimental discrimination is supplied.',
               'source_parameters_selected_by_native_idle':False,
               'interacting_relativistic_particle_identification':False}}
