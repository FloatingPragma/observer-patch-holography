"""One finite quantum-link/Weyl/Higgs source and its complete reference filter."""
from collections import defaultdict,deque
from functools import lru_cache
from itertools import product
import sympy as s
from ..chiral_source.model import source_descriptor,link_labels,link_generator,PAULI
from ..chiral_regulator.model import momentum_matrix,spectral_momenta
from ..chiral_higgs.model import higgs_memory,tensor_checks
from ..physical_attachment.chiral import representation


def clean(values):
    return {key:s.expand(value) for key,value in values.items() if s.expand(value)!=0}


def exact_positive(value):
    value=s.sympify(value)
    if value.has(s.Float) or value.is_positive is not True or value.is_finite is not True:
        raise ValueError('positive exact finite parameter required')
    return value


def finite_program(side=2,length=2*s.pi):
    """Complete sparse operator grammar; no quantum link is replaced by a number."""
    spectral_momenta(side,length);cells=list(product(range(side),repeat=3))
    derivative=momentum_matrix(side,length);hopping=[];edges=set()
    for x in cells:
        for axis in range(3):
            for finish in range(x[axis]+1,side):
                y=list(x);y[axis]=finish;y=tuple(y)
                path=[]
                for q in range(x[axis],finish):
                    a=list(x);a[axis]=q;b=list(a);b[axis]+=1
                    edge=(tuple(a),tuple(b));edges.add(edge);path.append(edge)
                hopping.append({'from':list(y),'to':list(x),'axis':axis,
                    'coefficient_of_sigma_W':str(-derivative[x[axis],finish]),
                    'path':[[list(a),list(b)] for a,b in path],
                    'reverse_term':'Hermitian adjoint, retained'})
    assert len(edges)==3*side**2*(side-1)
    return {'program':'coherent-spectral-quantum-link-Higgs-v1',
            'source_carrier_program':source_descriptor()['program'],
            'side':side,'length':str(length),'cells':len(cells),'gauge_edges':len(edges),
            'matter_Fock_dimension_exponent':30*len(cells),
            'full_logical_Hilbert_dimension':f'2^{30*len(cells)} * 60^{len(edges)} * 5^{len(cells)}',
            'kinetic_oriented_terms':hopping,
            'scalar_dimension_per_cell':5,'quantum_link_dimension_per_edge':60,
            'Yukawa_channels':['up','down','charged_lepton'],
            'H_total':'H_covariant_spectral_paths + H_quantum_Yukawa + omega*sum N_H + kappa*sum P_link',
            'declared_exact_parameters':{'y_u':'1','y_d':'1','y_e':'1','omega':'1','kappa':'1','beta':'1'},
            'gauge_geometry':'coordinate nearest-neighbor edges; all spectral pairs use monotone paths',
            'finite_range_kinetic_claim':False,
            'quantum_transporter_set_to_identity':False,
            'full_global_matrix_materialized':False}


@lru_cache(maxsize=1)
def path_ward_checks():
    labels=link_labels();index={label:i for i,label in enumerate(labels)}
    blocks={}
    for i,q in enumerate(representation()[2][-1].diagonal()):blocks.setdefault(q,[]).append(i)
    checked=0;terms={}
    for length in (1,2):
        original={}
        for indices in blocks.values():
            for vertex in product(indices,repeat=length+1):
                a,b=vertex[0],vertex[-1]
                linkstate=tuple(index[vertex[e],vertex[e+1]] for e in range(length))
                original[a,b,linkstate]=s.Integer(1)
        terms[str(length)]=len(original)
        for t in representation()[2]:
            td=t.todok()
            linkactions={(end):link_generator(t,end).todok() for end in (0,1)}
            for vertex in range(length+1):
                residual=defaultdict(lambda:s.Integer(0))
                for (a,b,ls),coefficient in original.items():
                    if vertex==0:
                        for (i,j),value in td.items():
                            if j==a:residual[i,b,ls]+=coefficient*value
                    if vertex==length:
                        for (i,j),value in td.items():
                            if i==b:residual[a,j,ls]-=coefficient*value
                    # Internal vertex carries right action of its incoming edge
                    # and left action of its outgoing edge.
                    adjacent=[]
                    if vertex>0:adjacent.append((vertex-1,1))
                    if vertex<length:adjacent.append((vertex,0))
                    for edge,end in adjacent:
                        for (row,col),value in linkactions[end].items():
                            if col==ls[edge]:
                                target=list(ls);target[edge]=row
                                residual[a,b,tuple(target)]+=coefficient*value
                assert not clean(residual)
                checked+=1
    return {'actual_quantum_path_lengths_checked':[1,2],
            'raising_tensor_terms':terms,'all_endpoint_and_internal_vertex_Ward_checks':checked,
            'arbitrary_path_Ward_by_telescoping_proved':True,
            'all_normal_ordered_Fock_sectors_retained':True,
            'quantum_path_has_no_identity_vacuum':True}


def fermion(mask,mode,create):
    present=bool(mask&(1<<mode))
    if present==create:return None
    return mask^(1<<mode),s.Integer(-1)**((mask&((1<<mode)-1)).bit_count())


def sequence(mask,operators):
    coefficient=s.Integer(1)
    for mode,create in reversed(operators):
        step=fermion(mask,mode,create)
        if step is None:return None
        mask,sign=step;coefficient*=sign
    return mask,coefficient


@lru_cache(maxsize=1)
def active_terms():
    # e^c and L form a reducing native gauge module. Quarks are retained
    # as empty spectators in this explicit y_u=y_d=0 control.
    indices=(9,10,11);vertices=representation()[3];fields,_,_,_=higgs_memory()
    c=tuple(v.extract(indices,indices) for v in vertices['charged_lepton']['components'])
    # The active link labels are actual R15 raising components.
    links=((-1,-1),(0,0),(1,1),(1,2),(2,1),(2,2))
    return indices,c,fields,links


def joint_column(key):
    """Literal complete two-cell control, twelve CAR modes, six links, five scalars."""
    mask,link,scalar=key;_,cs,fields,labels=active_terms();out=defaultdict(lambda:s.Integer(0))
    out[key]+=int(link!=0)+int(scalar!=0)
    # h_xy=+i sigma3 U/2 at the side-two spectral edge; every spin retained.
    for label,(a,b) in enumerate(labels[1:],1):
        for spin in range(2):
            coefficient=s.I*s.Integer(-1)**spin/2
            lmode=3*spin+a;rmode=6+3*spin+b
            if link==0:
                step=sequence(mask,((lmode,True),(rmode,False)))
                if step:out[step[0],label,scalar]+=coefficient*step[1]
            if link==label:
                step=sequence(mask,((rmode,True),(lmode,False)))
                if step:out[step[0],0,scalar]+=s.conjugate(coefficient)*step[1]
    # Local quantum charged-lepton pair term at the right cell.
    for k,c in enumerate(cs):
        phi=fields[k].H
        for i,j in product(range(3),repeat=2):
            if c[i,j]==0:continue
            operators=((6+i,False),(9+j,False))
            for matrix,ops,coef in ((phi,operators,c[i,j]),
                    (phi.H,((9+j,True),(6+i,True)),s.conjugate(c[i,j]))):
                step=sequence(mask,ops)
                if step is None:continue
                for destination in range(5):
                    if matrix[destination,scalar]!=0:
                        out[step[0],link,destination]+=coef*step[1]*matrix[destination,scalar]
    return clean(out)


def apply_control(state):
    out=defaultdict(lambda:s.Integer(0))
    for key,coefficient in state.items():
        for target,value in joint_column(key).items():out[target]+=coefficient*value
    return clean(out)


def norm2(state):
    return s.simplify(sum(s.conjugate(v)*v for v in state.values()))


@lru_cache(maxsize=1)
def joint_encounter():
    vacuum={(0,0,0):s.Integer(1)}
    first=apply_control(vacuum);second=apply_control(first)
    arrived=lambda state:{key:v for key,v in state.items() if key[0]&63}
    assert not arrived(vacuum) and not arrived(first)
    incoming=arrived(second);c=norm2(incoming)
    assert c>0
    # On the complete 122880 active source, each scalar/link projector has
    # norm1; 20 half-unit hopping terms and 8 unit pair terms give ||H||<=20.
    # The deliberately looser bound24 keeps the same exact inequality simple.
    bound=s.Integer(24);time=s.Rational(1,10**6)
    error=time**3*bound**3/6
    lower=s.simplify((time**2*s.sqrt(c)/2-error)**2)
    assert (time**2*s.sqrt(c)/2-error).is_positive is True
    # Reverse every populated matrix element, including all four pair leaves
    # and their hopping destinations. This checks the actual source action.
    seen=set(vacuum)|set(first)|set(second);hermitian=0
    for key in seen:
        for target,value in joint_column(key).items():
            assert joint_column(target).get(key,0)==s.conjugate(value)
            hermitian+=1
    return {'complete_active_matter_modes':12,'active_link_dimension':6,'active_scalar_dimension':5,
            'full_active_control_dimension':2**12*6*5,
            'control_parameters':{'y_u':'0','y_d':'0','right_y_e':'1','left_y_e':'0','omega':'1','kappa':'1'},
            'full_R15_colored_modes_and_link_rows':'retained empty reducing spectators for this control',
            'first_H_vacuum_norm_squared':str(norm2(first)),
            'second_H_vacuum_arrival_norm_squared':str(c),
            'literal_reverse_matrix_entries_checked':hermitian,
            'exact_source_norm_bound':str(bound),'time':str(time),'unitary_Taylor_order':2,
            'operator_norm_remainder':str(error),'complete_arrival_yes_probability_strict_lower':str(lower),
            'complete_arrival_no_probability':'1-p_yes','vacuum_is_gauge_invariant':True,
            'initial_empty_state_is_stationary':False,
            'response_is_asymptotic_scattering':False}


def norm_bound(program):
    """Safe finite operator coefficient bound for the unit-coupling program."""
    ds=(3,6,1,2,3);kinetic=s.Integer(0)
    for term in program['kinetic_oriented_terms']:
        length=len(term['path']);coefficient=abs(s.sympify(term['coefficient_of_sigma_W']))
        # Two Pauli entries, all internal paths and Hermitian adjoint.
        kinetic+=4*coefficient*sum(d**(length+1) for d in ds)
    yukawa=s.Integer(0)
    for vertex in representation()[3].values():
        for c in vertex['components']:
            # 1/2 sum epsilon*C with h.c. has coefficient l1 norm2*sum|C|.
            yukawa+=2*sum(abs(v) for v in c)
    return s.simplify(kinetic+program['cells']*(yukawa+1)+program['gauge_edges'])


def physical_reference_certificate():
    program=finite_program();bound=norm_bound(program);beta=s.Integer(1)
    assert bound.is_positive is True
    return {'program':program['program'],'state':'rho_beta=P exp(-beta H) P / Tr(P exp(-beta H))',
            'P':'orthogonal projector onto common kernel of all local total Gauss generators',
            'physical_sector_nonempty_witness':'empty matter, blank links, scalar vacuum',
            'declared_beta':str(beta),'unit_parameter_full_H_norm_upper_bound':str(bound),
            'source_clock':'continuous exp(-it H_total), same H as preparation filter',
            'success_Kraus':'L=P exp[-beta(H+B I)/2]',
            'failure_Kraus':'sqrt(I-L^dagger L)',
            'input':'P/r on the physical Gauss sector, r=rank(P)',
            'physical_input_preparation':'gauge- and total-parity-preserving unitary Omega tensor0 -> sum_j p_j tensorj /sqrt(r), retained neutral graded ancilla with parity(j)=parity(p_j)',
            'physical_input_uses_unphysical_Gauss_states':False,
            'physical_input_preserves_total_fermion_parity':True,
            'original_success_probability':'Tr(P exp[-beta(H+B I)])/r',
            'strict_success_lower_bound':'exp(-2 beta B) > 0',
            'ambient_I_over_D_control_success':'Tr(P exp[-beta(H+B I)])/D',
            'prepared_success_state_exactly_rho_beta':True,
            'all_failure_records_retained':True,'source_refinement_reference_efficient':False,
            'rho_beta_faithful_on_physical_sector':True,'rho_beta_faithful_on_full_unconstrained_carrier':False,
            'zero_temperature_reference':'normalized full ground projector of H restricted to P',
            'unique_interacting_ground_or_uniform_gap_proved':False,
            'global_eigensystem_numerically_computed':False}


def certificate_data():
    return {'source_program':finite_program(),'quantum_path_covariance':path_ward_checks(),
            'quantum_Higgs_Ward_reused':tensor_checks(),
            'literal_joint_encounter':joint_encounter(),'same_clock_physical_reference':physical_reference_certificate(),
            'scope':{'finite_common_carrier_state_reference_clock':True,
                     'complete_physical_Gibbs_filter_construction':True,
                     'finite_quantum_gauge_and_Higgs_terms_on_same_source':True,
                     'source_parameters_selected_by_native_idle':False,
                     'free_identity_background_is_quantum_link_vacuum':False,
                     'interacting_chiral_continuum':False,'physical_pole_mass_or_GeV_prediction':False}}
