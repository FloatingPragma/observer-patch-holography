"""Finite affine-Gauss source on the committed cone, with an explicit cutoff boundary."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import itertools
import math
import sys

import sympy as s

HERE=Path(__file__).resolve().parent
RER=HERE.parents[1]
sys.path.insert(0,str(RER/'code/electromagnetism'))
import cone_whitney_bridge as cone
import discrete_coulomb_green as green

X=s.Symbol('x')
ROOT5=s.sqrt(5)
ROOTS=(s.Integer(0),4-ROOT5,s.Integer(3),6-ROOT5,s.Integer(4),s.Integer(6),4+ROOT5,s.Integer(7),6+ROOT5)
MULTIPLICITIES=(12,3,5,3,4,4,3,5,3)
TOTAL_CUTOFF=3
DEGREE=TOTAL_CUTOFF+1
REFERENCE_VERTEX=1


def simple(value):
    return s.simplify(s.expand(value))


def serialize(value):
    if isinstance(value,dict):return {str(k):serialize(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [serialize(v) for v in value]
    if isinstance(value,s.Basic):return str(value)
    return value


@lru_cache(None)
def geometry():
    # Reuse the committed orientation and cone extension, not a new mesh.
    edges,faces,d,c,vertices=cone.load_geometry()
    maps=cone.cone_maps(d,c)
    D,C=maps['D'],maps['C'];L=D.T*D;K=C.T*C
    assert D.shape==(42,13) and C.shape==(50,42)
    assert C*D==s.zeros(50,13) and D.rank()==12 and C.rank()==30
    # The inherited inverse routine is dimension independent. Its twelve-port
    # green_matrix wrapper is not used on the thirteen-vertex cone.
    inverse=green.invert([[F(v+1) for v in row] for row in L.tolist()])
    G=s.Matrix(inverse)-s.ones(13)/169
    assert L*G==s.eye(13)-s.ones(13)/13 and G==G.T and G*s.ones(13,1)==s.zeros(13,1)
    PT=s.eye(42)-D*G*D.T
    assert PT==PT.T and PT*PT==PT and D.T*PT==s.zeros(13,42)
    assert PT*K==K and K*PT==K
    m=s.Poly(X*(X-3)*(X-4)*(X-6)*(X-7)*(X*X-12*X+31)*(X*X-8*X+11),X)
    powers=[s.eye(42)]
    for _ in range(m.degree()):powers.append(powers[-1]*K)
    assert sum((coef*powers[i] for (i,),coef in m.terms()),s.zeros(42))==s.zeros(42)
    expected=s.prod((X-lam)**mult for lam,mult in zip(ROOTS,MULTIPLICITIES))
    assert s.Poly(s.expand(expected),X)==K.charpoly(X).as_poly()
    polynomials=[]
    for lam in ROOTS:
        numerator=s.div(m.as_expr(),X-lam,X,extension=ROOT5)[0]
        denominator=s.diff(m.as_expr(),X).subs(X,lam)
        polynomial=s.Poly(numerator/denominator,X,extension=ROOT5)
        assert s.rem(polynomial*polynomial-polynomial,m,extension=ROOT5).is_zero
        assert s.rem((s.Poly(X,X)-lam)*polynomial,m,extension=ROOT5).is_zero
        polynomials.append(polynomial)
    assert sum(polynomials,s.Poly(0,X)).as_expr()==1
    p0=sum((coef*powers[i] for (i,),coef in polynomials[0].terms()),s.zeros(42))
    assert p0==s.eye(42)-PT
    # Symmetric orthogonal projector polynomials certify nonnegative weights.
    return edges,faces,vertices,D,C,G,K,PT,powers,polynomials


def projected_weight(vector,polynomial,powers):
    values=[(vector.T*power*vector)[0] for power in powers]
    return simple(sum(coef*values[i] for (i,),coef in polynomial.terms()))


@lru_cache(None)
def charge_geometry():
    edges,faces,vertices,D,C,G,K,PT,powers,polynomials=geometry()
    Q=s.eye(13)-s.eye(13)[:,REFERENCE_VERTEX]*s.ones(1,13)
    electric=D*G*Q
    V=Q.T*G*Q
    assert D.T*electric==Q and C*electric==s.zeros(50,13)
    assert electric.T*electric==V
    assert PT*electric==s.zeros(42,13)
    distances=green.bfs_distances([u for u,v in edges],[v for u,v in edges])
    points=[s.zeros(3,1)]+vertices
    rows=[]
    for y in range(13):
        rho=Q[:,y];field=electric[:,y]
        ec=V[y,y]/2
        assert (field.T*field)[0]/2==ec
        assert ec>=0
        if y!=REFERENCE_VERTEX:assert ec>0
        rows.append({'source_vertex':y,'reference_vertex':REFERENCE_VERTEX,
            'boundary_surface_hop_distance':None if y==0 else distances[REFERENCE_VERTEX-1][y-1],
            'cone_graph_hop_distance':0 if y==REFERENCE_VERTEX else (1 if y==0 or distances[REFERENCE_VERTEX-1][y-1]==1 else 2),
            'declared_coordinate_distance_squared':simple((points[y]-points[REFERENCE_VERTEX]).dot(points[y]-points[REFERENCE_VERTEX])),
            'charge_vector':list(rho),'longitudinal_electric_field':list(field),
            'static_electric_energy':ec,'matter_energy':s.Integer(1),'reference_energy':s.Integer(1),
            'total_energy_above_actual_neutral_ground':2+ec,
            'source_energy_excluding_reference':1+ec,
            'complete_pump_probability':s.Rational(1,4),'pump_failure_probability':s.Rational(3,4),
            'conditional_energy_line_probability':s.Integer(1),'conditional_all_other_energy_probability':s.Integer(0)})
    return rows,V


def cutoff_control(d=DEGREE):
    if d<2:raise ValueError('at least two oscillator levels are required')
    a=s.zeros(d)
    for n in range(1,d):a[n-1,n]=s.sqrt(n)
    number=s.diag(*range(d));top=s.zeros(d);top[d-1,d-1]=1
    # Frequency-one calculation scales exactly to every positive omega.
    q=(a+a.H)/s.sqrt(2);p=s.I*(a.H-a)/s.sqrt(2)
    bare=(p*p+q*q)/2
    correction=s.Rational(d,2)*top
    H=bare+correction
    assert a*a.H-a.H*a==s.eye(d)-d*top
    assert H==number+s.eye(d)/2
    assert H*q-q*H==-s.I*p and H*p-p*H==s.I*q
    vacuum=s.eye(d)[:,0]
    assert (q*vacuum).norm()**2==s.Rational(1,2)
    assert bare!=H and s.trace(q*p-p*q)==0
    if d==2:assert bare==s.eye(2)/2
    # Every vacuum word through degree six agrees with the unbounded
    # oscillator at d=4: a return path cannot rise above three quanta.
    moments=[]
    for k in range(1,min(d,4)):
        moment=simple((vacuum.T*(q**(2*k))*vacuum)[0])
        expected=s.Rational(math.prod(range(1,2*k,2)),2**k)
        assert moment==expected
        moments.append({'order':2*k,'vacuum_q_moment_at_frequency_one':moment})
    return {'levels_per_mode':d,'finite_CCR':'[a,a^dag]=I-d|d-1><d-1|',
        'boundary_energy_per_mode':'(d omega/2)|d-1><d-1|',
        'completed_H_per_mode':'omega (n+1/2)',
        'bare_quadrature_diagonal_at_frequency_one':list(bare.diagonal()),
        'completed_diagonal_at_frequency_one':list(H.diagonal()),
        'exact_Heisenberg_equations_for_finite_q_p':True,
        'untruncated_CCR_claim':False,'checked_vacuum_moments':moments}


@lru_cache(None)
def field_records():
    edges,faces,vertices,D,C,G,K,PT,powers,polynomials=geometry()
    # omega_max<3, omega_min>1; a row of C has squared norm three.
    # Cauchy-Schwarz over 30 modes gives norm^2 <= 2(d-1)*30*3=540.
    norm_upper=s.Integer(24)
    assert 2*(DEGREE-1)*30*3<norm_upper**2
    channels=[]
    for kind,indices in (('transverse_electric_edge',(0,12)),('magnetic_face',(0,20))):
        for index in indices:
            vector=s.eye(42)[:,index] if kind=='transverse_electric_edge' else C[index,:].T
            lines=[];success=0
            for lam,mult,polynomial in zip(ROOTS[1:],MULTIPLICITIES[1:],polynomials[1:]):
                weight=projected_weight(vector,polynomial,powers)
                assert weight>=0
                residue=simple(s.sqrt(lam)*weight/2 if kind=='transverse_electric_edge' else weight/(2*s.sqrt(lam)))
                probability=residue/norm_upper**2
                lines.append({'squared_frequency':lam,'frequency':s.sqrt(lam),'degeneracy':mult,
                    'spatial_projector_weight':weight,'vacuum_residue':residue,
                    'unconditional_pump_and_energy_line_probability':probability})
                success+=probability
            success=simple(success)
            assert success>0 and success<1
            channels.append({'kind':kind,'index':index,'full_finite_operator_norm_upper':norm_upper,
                'energy_lines':lines,'pump_success_probability':success,
                'pump_failure_probability':1-success,'success_other_energy_probability':s.Integer(0),
                'probabilities_sum_exactly_one':simple(1-success+sum(row['unconditional_pump_and_energy_line_probability'] for row in lines))==1,
                'time_correlation':'sum_line residue * exp(-i frequency time)',
                'time_scope':'all real times on the same completed finite H'})
    return channels


def total_cutoff_control():
    count=0;coefficient_checks=0
    for total in range(TOTAL_CUTOFF+1):
        for word in itertools.combinations_with_replacement(range(30),total):
            n=[word.count(j) for j in range(30)]
            for value in n:
                a_dag_a=value
                a_a_dag=value+1 if total<TOTAL_CUTOFF else 0
                correction=s.Rational(value+1,2) if total==TOTAL_CUTOFF else 0
                assert s.Rational(a_dag_a+a_a_dag,2)+correction==value+s.Rational(1,2)
                coefficient_checks+=1
            count+=1
    assert count==math.comb(30+TOTAL_CUTOFF,TOTAL_CUTOFF)==5456
    minimum=4-ROOT5
    vacuum_energy=sum(mult*s.sqrt(lam)/2 for lam,mult in zip(ROOTS[1:],MULTIPLICITIES[1:]))
    wrong_top_energy=TOTAL_CUTOFF*s.sqrt(minimum)/2
    assert vacuum_energy>wrong_top_energy
    return {'total_photon_cutoff':TOTAL_CUTOFF,'photon_dimension':count,
        'mode_boundary_coefficient_checks':coefficient_checks,
        'source_representation':'direct sum Sym^n(C^30), n=0,...,3',
        'basis_independent_in_degenerate_frequency_spaces':True,
        'boundary_energy':'(1/2) P_(n_total=3) [dGamma(omega)+Tr(omega)]',
        'completed_H':'Tr(omega)/2 + dGamma(omega)',
        'bare_top_state_energy':wrong_top_energy,'actual_completed_vacuum_energy':vacuum_energy,
        'bare_cutoff_top_state_strictly_below_bare_vacuum':True,
        'all_vacuum_field_words_through_degree_six_exact':True,
        'error_bound_on_unrestricted_Weyl_or_large_occupation_states_claimed':False}


def joint_charge_field_records():
    charges,V=charge_geometry()
    # Every longitudinal bit contributes a dipole with Euclidean norm <1.
    assert all(V[j,j]<1 for j in range(13))
    # ||E_L,e||<12; ||E_T,e||<24; choose a single full-carrier bound48.
    assert 12+24<48
    by_edge={row['index']:row for row in field_records() if row['kind']=='transverse_electric_edge'}
    rows=[]
    for charge in charges:
        for edge,field in by_edge.items():
            static=charge['longitudinal_electric_field'][edge]
            elastic=static**2/s.Integer(48**2)
            inelastic=[{'frequency':line['frequency'],
                'total_energy_above_neutral_ground':charge['total_energy_above_actual_neutral_ground']+line['frequency'],
                'conditional_field_pump_and_energy_probability':line['vacuum_residue']/48**2,
                'unconditional_probability':line['vacuum_residue']/(4*48**2)} for line in field['energy_lines']]
            field_success=simple(elastic+sum(line['conditional_field_pump_and_energy_probability'] for line in inelastic))
            charge_failure=s.Rational(3,4);field_failure=(1-field_success)/4
            total=charge_failure+field_failure+elastic/4+sum(line['unconditional_probability'] for line in inelastic)
            assert simple(total)==1 and field_success>0 and field_success<1
            rows.append({'source_vertex':charge['source_vertex'],'reference_vertex':REFERENCE_VERTEX,
                'electric_edge':edge,'full_electric_operator_norm_upper':48,
                'source_state_after_charge_success':'same-H static charge ground and transverse vacuum',
                'charge_pump_failure_probability':charge_failure,
                'charge_success_field_pump_failure_probability':field_failure,
                'elastic_total_energy_above_neutral_ground':charge['total_energy_above_actual_neutral_ground'],
                'elastic_unconditional_probability':elastic/4,
                'inelastic_lines':inelastic,'other_energy_unconditional_probability':s.Integer(0),
                'sum_of_all_twelve_records':'1','retained_pointer_dimension':12,
                'electric_mean_in_charged_ground':static,
                'charged_ground_full_electric_correlation':'electric_mean^2 + sum_line vacuum_residue * exp(-i frequency time)'})
    return rows


def charge_car_check():
    # V embeds every one of 2^13 source occupations with R=N and twelve
    # matched longitudinal bits. Verify the exact dressed ladder action.
    checked=0
    for mask in range(1<<13):
        occupied=[i for i in range(13) if mask>>i&1]
        count=len(occupied)
        charge=[int(mask>>i&1)-count*int(i==REFERENCE_VERTEX) for i in range(13)]
        longitudinal=[int(mask>>i&1) if i!=REFERENCE_VERTEX else -sum(int(mask>>j&1) for j in range(13) if j!=REFERENCE_VERTEX) for i in range(13)]
        assert longitudinal==charge
        for y in range(13):
            if mask>>y&1:continue
            after=mask|(1<<y)
            delta=[int(i==y)-int(i==REFERENCE_VERTEX) for i in range(13)]
            target=[int(after>>i&1)-(count+1)*int(i==REFERENCE_VERTEX) for i in range(13)]
            assert [a+b for a,b in zip(charge,delta)]==target
            assert count+1<=13
            sign=(-1)**sum(mask>>i&1 for i in range(y))
            assert sign in (-1,1)
            checked+=1
    # Independent Jordan-Wigner anticommutation by applying both orders,
    # including occupied-state zeros, on every thirteen-mode basis vector.
    def create(mask,y):
        return None if mask>>y&1 else (mask|(1<<y),(-1)**((mask&((1<<y)-1)).bit_count()))
    pairs=0
    for mask in range(1<<13):
        for x,y in ((0,1),(1,2),(0,12),(5,9)):
            a=create(mask,x);b=create(mask,y)
            ax=None if a is None else create(a[0],y)
            by=None if b is None else create(b[0],x)
            assert (ax is None)==(by is None)
            if ax is not None:assert ax[0]==by[0] and a[1]*ax[1]+b[1]*by[1]==0
            pairs+=1
    return {'complete_source_CAR_modes':13,'all_source_occupation_states':1<<13,
        'Gauss_preserving_creation_columns':checked,'independent_CAR_basis_pair_controls':pairs,
        'reference_charge_states':14,'reference_charge':'-R at vertex 1',
        'reference_energy':'mu R, mu=1','reference_parity':'(-1)^R',
        'longitudinal_bits':12,'physical_Gauss_constraint':'ell_y=n_y (y!=1), R=N',
        'dressed_creation':'V c_y^dag V^dag, full extension c_y^dag S_R^+ S_ell,y^+',
        'joint_dressed_creation_is_even':True,
        'matter_creation_is_odd':True,'physical_spin_identification':False}


def exact_ground_gaps():
    _,V=charge_geometry()
    sites=[j for j in range(13) if j!=REFERENCE_VERTEX]
    diagonal=[int(V[j,j]*217/2) for j in sites]
    pairs={(i,j):int(V[sites[i],sites[j]]*217) for i in range(12) for j in range(i+1,12)}
    assert all(s.Integer(diagonal[i])==V[j,j]*217/2 for i,j in enumerate(sites))
    assert all(s.Integer(value)==V[sites[i],sites[j]]*217 for (i,j),value in pairs.items())
    energies=[]
    for mask in range(1<<12):
        numerator=sum(diagonal[i] for i in range(12) if mask>>i&1)
        numerator+=sum(value for (i,j),value in pairs.items() if mask>>i&1 and mask>>j&1)
        energies.append(s.Rational(numerator,217))
    assert energies[0]==0 and min(energies[1:])==s.Rational(26,217)
    gap=s.sqrt(4-ROOT5)
    assert s.Rational(26,217)<1<gap<2
    a=s.zeros(4)
    for j in range(1,4):a[j-1,j]=s.sqrt(j)
    q=(a+a.H)/s.sqrt(2)
    finite=(q**8)[0,0]
    assert finite==s.Rational(81,16)!=s.Rational(105,16)
    return {'complete_longitudinal_patterns_checked':len(energies),
        'minimum_nonzero_longitudinal_energy':s.Rational(26,217),
        'actual_full_kinematic_ground_gap':s.Rational(26,217),
        'actual_physical_ground_gap':gap,
        'frequency_one_single_mode_vacuum_eighth_moment':finite,
        'unrestricted_Gaussian_eighth_moment':s.Rational(105,16),
        'higher_Weyl_or_all_moment_transfer_rejected':True}


def pins():
    names=('code/electromagnetism/cone_whitney_bridge.py',
        'code/electromagnetism/discrete_coulomb_green.py',
        'code/electromagnetism/runtime/discrete_coulomb_green_receipt.json',
        'Lean/Screen/SeamCurrentCarrierQuotient.lean',
        'Lean/ObserverPatchHolography/CoreAxioms.lean',
        'Lean/Screen/SeamCurrentEdge30Moment.lean',
        'Lean/Screen/DiscreteCoulombGreen.lean','Lean/Screen/ConeCochainBridge.lean',
        'Lean/Screen/WhitneyConeModes.lean','Lean/Screen/WhitneyNormalModeConstruction.lean',
        'Lean/Screen/WhitneyQuantumBridge.lean','paper/tex_fragments/WHITNEY_MAXWELL_DYNAMICS.tex')
    return {name:hashlib.sha256((RER/name).read_bytes()).hexdigest() for name in names}


def calculate():
    inherited=green.verify_committed_receipt()
    edges,faces,vertices,D,C,G,K,PT,powers,polynomials=geometry()
    charges,V=charge_geometry()
    result={'schema':'oph.quantum_worker.finite_coulomb_source.v1',
        'architecture':{'source':'finite affine-Gauss abelian source on the committed contractible cone',
            'relation_to_PW119':'different declared source; no vacuum or current transfer',
            'spatial_geometry':'inherited 13 vertices, 42 edges, 50 faces, 20 tetrahedra',
            'electric_and_magnetic_metrics':'identity on edges and faces; declared counting metrics, not Whitney mass forms',
            'charge_unit':'1','matter_energy_per_fermion':'1','reference_energy_per_unit':'1',
            'reference_vertex':REFERENCE_VERTEX,'transverse_modes':30,'total_photon_number_cutoff':TOTAL_CUTOFF,
            'clock':'exp(-i t H), hbar=1',
            'Hamiltonian':'N+R+||E||^2/2+||B||^2/2+(1/2)P_top[dGamma(omega)+Tr(omega)]',
            'ground':'all matter, reference and longitudinal registers zero; every mode vacuum',
            'physical_ground_is_unique':True,'preparation_success_exact':'1',
            'finite_physical_dimension_expression':'2^13 * 5456 = 44695552',
            'finite_kinematic_dimension_expression':'2^13 * 14 * 2^12 * 5456',
            'physical_source_dimension':2**13*5456,
            'kinematic_source_reference_dimension':2**13*14*2**12*5456,
            'retained_complete_joint_detector_pointer_dimension':12,
            'kinematic_source_reference_and_pointer_dimension':12*2**13*14*2**12*5456,
            'other_controller_storage_included':False,
            'compact_local_link_realization':False,'mobile_matter_hopping':False,
            'thermodynamic_Coulomb_phase':False,'measured_particle_assignment':False},
        'geometry':{'gradient_rank':12,'curl_rank':30,'constant_longitudinal_zero_mode_removed':True,
            'transverse_harmonic_zero_modes':0,'neutral_charge_constraint_exact':True,
            'green_matrix':G.tolist(),'scalar_product':'Euclidean counting metrics',
            'electric_transverse_projector_diagonal':list(PT.diagonal()),
            'positive_squared_frequencies':[{'value':lam,'multiplicity':m} for lam,m in zip(ROOTS[1:],MULTIPLICITIES[1:])],
            'minimum_polynomial':str(s.prod(X-lam for lam in ROOTS).expand()),
            'full_stiffness_characteristic_polynomial':str(s.factor(K.charpoly(X).as_expr()))},
        'charge_reference_and_Gauss':charge_car_check(),
        'exact_ground_gaps_and_moment_boundary':exact_ground_gaps(),
        'cutoff_boundary':total_cutoff_control(),'single_mode_cutoff_corner_check':cutoff_control(),
        'bare_two_level_obstruction':cutoff_control(2),
        'complete_static_charge_records':charges,'complete_field_records':field_records(),
        'complete_same_source_charge_then_full_electric_records':joint_charge_field_records(),
        'inherited_twelve_port_Green_receipt_freshly_replayed':inherited['receipt_sha256'],
        'input_sha256':pins(),
        'owned_source_sha256':{path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(HERE.iterdir()) if path.suffix in ('.py','.md')}}
    return serialize(result)
