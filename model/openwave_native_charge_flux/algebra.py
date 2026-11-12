"""Complete endpoint algebra and exact six-route charge/flux attachment."""
from collections import Counter
from functools import lru_cache
import sympy as s
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory, tensor_checks
from model.coherent_defect_wave.chiral_joint.oriented_links.model import endpoint_current, labels
from model.coherent_defect_wave.chiral_joint.peter_weyl_links.checks import matrix_checks
from model.coherent_defect_wave.chiral_joint.peter_weyl_links import model as pw
from model.coherent_defect_wave.chiral_joint.interacting_charge_localization.geometry import graph
from model.coherent_defect_wave.chiral_joint.interacting_charged_transport.algebra import routes
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import FAMILIES, operator_certificate
from model.coherent_defect_wave.chiral_joint.periodic_family.model import Configuration, kinetic_terms, cells


def current(current_index, endpoint):
    return endpoint_current(current_index, endpoint)


@lru_cache(maxsize=2)
def charge_weights(endpoint):
    y=current(11,endpoint)
    assert y==s.diag(*y.diagonal())
    return tuple(6*x for x in y.diagonal())


@lru_cache(maxsize=256)
def cut_endpoints(mask):
    if type(mask) is not int or mask not in range(256):
        raise ValueError('one subset of the eight actual cells is required')
    return tuple((edge,0 if mask>>left&1 else 1)
                 for edge,left,right in graph()[0]
                 if ((mask>>left)^(mask>>right))&1)


def cut_flux_six(mask, link_labels):
    if len(link_labels)!=24 or any(type(i) is not int or i not in range(119) for i in link_labels):
        raise ValueError('all twenty-four complete PW119 link coordinates are required')
    return -sum(charge_weights(endpoint)[link_labels[edge]] for edge,endpoint in cut_endpoints(mask))


@lru_cache(maxsize=2)
def constrained_arrival_link(origin_endpoint):
    """Kernel of all origin Gauss currents with the actual ec reference."""
    if origin_endpoint not in (0,1):raise ValueError('endpoint zero or one required')
    ec=FAMILIES['ec'][0]; gauges=representation()[2]
    y=current(11,origin_endpoint)
    candidates=[i for i in range(119) if y[i,i]==-gauges[11][ec,ec]]
    if len(candidates)!=1:raise ArithmeticError('the complete origin hypercharge constraint is not rank one')
    index=candidates[0];vector=s.eye(119)[:,index]
    for a,t in enumerate(gauges):
        if current(a,origin_endpoint)*vector+t[ec,ec]*vector!=s.zeros(119,1):
            raise ArithmeticError('full origin Gauss current failed')
        if current(a,1-origin_endpoint)*vector-t[ec,ec]*vector!=s.zeros(119,1):
            raise ArithmeticError('arrival endpoint does not carry the retained native charge')
    assert gauges[10][ec,ec]==0 and gauges[11][ec,ec]==1
    assert (current(10,1-origin_endpoint)+current(11,1-origin_endpoint))*vector==vector
    return index


@lru_cache(maxsize=25)
def boundary_spectrum(edge_count):
    """Exact full tensor-carrier spectral multiplicities, without its matrix."""
    if type(edge_count) is not int or edge_count not in range(25):
        raise ValueError('boundary edge count must lie between zero and twenty-four')
    single=Counter(map(int,charge_weights(0)))
    assert single==Counter(map(int,charge_weights(1)))
    counts=Counter({0:1})
    for _ in range(edge_count):
        target=Counter()
        for first,m in counts.items():
            for second,n in single.items():target[first-second]+=m*n
        counts=target
    assert sum(counts.values())==119**edge_count
    return counts


@lru_cache(maxsize=1)
def gauss_cut_census():
    """Independent actual Gauss compiler on a complete additive factor basis."""
    y=representation()[2][11];ys=higgs_memory()[1][11]
    edges={edge:(left,right) for edge,left,right in graph()[0]}
    blank=(0,(0,)*24,(0,)*8)
    cases=[(blank,[s.Integer(0)]*8)]
    for cell in range(8):
        for mode in range(30):
            expected=[s.Integer(0)]*8;expected[cell]=y[mode%15,mode%15]
            cases.append(((1<<(30*cell+mode),blank[1],blank[2]),expected))
        for scalar in range(1,5):
            values=list(blank[2]);values[cell]=scalar
            expected=[s.Integer(0)]*8;expected[cell]=ys[scalar,scalar]
            cases.append(((0,blank[1],tuple(values)),expected))
    for edge,(left,right) in edges.items():
        for label in range(1,119):
            links=list(blank[1]);links[edge]=label
            expected=[s.Integer(0)]*8
            expected[left]=current(11,0)[label,label]
            expected[right]=current(11,1)[label,label]
            cases.append(((0,tuple(links),blank[2]),expected))
    for key,expected in cases:
        for cell in range(8):
            actual=pw.gauss_column(key,cell,11)
            assert actual==({key:expected[cell]} if expected[cell] else {})
    assert len(cases)==1+240+32+24*118
    return {'independent_actual_Gauss_column_factor_basis_inputs':len(cases),
        'independent_cell_Gauss_comparisons':8*len(cases),
        'all_matter_modes':240,'all_scalar_nonvacuum_factor_states':32,
        'all_link_nonvacuum_factor_states':24*118,
        'extension':'The actual Y Gauss compiler is the sum of these diagonal factor operators; arbitrary occupations and simultaneous nonblank factors follow by additivity.',
        'reference_extension':'Add the same native one-particle reference Y at cell zero; the actual full-Fock pump commutes with all twelve joint currents.'}


@lru_cache(maxsize=1)
def endpoint_certificate():
    gauges=representation()[2]; scalar_gauges=higgs_memory()[1]
    y=gauges[11];ys=scalar_gauges[11]
    # Central hypercharge commutes with every native gauge generator, including
    # the full scalar and link memories. Its cut measurement is gauge invariant.
    for a in range(12):
        assert y*gauges[a]==gauges[a]*y
        assert ys*scalar_gauges[a]==scalar_gauges[a]*ys
        for end in (0,1):
            for other in (0,1):
                assert current(11,end)*current(a,other)==current(a,other)*current(11,end)
    assert current(11,0)+current(11,1)==s.zeros(119)
    electric_total=current(10,0)+current(11,0)+current(10,1)+current(11,1)
    witness=next(i for i in range(119) if electric_total[i,i]!=0)
    assert electric_total==s.diag(*electric_total.diagonal())
    ward=matrix_checks();scalar=tensor_checks();pump=operator_certificate()
    assert ward['complete_endpoint_Ward_identities']==1416
    return {'native_current':'Y, native current index 11; Q=T3+Y on the attached ec arrival sector',
        'all_native_currents':12,'link_carrier_dimension':119,
        'full_endpoint_Ward_identities':ward['complete_endpoint_Ward_identities'],
        'complete_scalar_Yukawa_checks':scalar,
        'complete_static_reference_pump_checks':pump,
        'matter_one_spin_six_Y':list(map(str,y.diagonal()*6)),
        'scalar_six_Y':list(map(str,ys.diagonal()*6)),
        'link_endpoint_six_Y':[list(map(str,charge_weights(end))) for end in (0,1)],
        'hypercharge_is_central_on_matter_scalar_and_both_link_endpoints':True,
        'complete_119_state_interior_hypercharge_endpoint_sum':'0',
        'noncentral_Q_cancellation_control':{'link_basis_index':witness,
            'link_label':list(labels()[witness]),'Q_endpoint_sum':str(electric_total[witness,witness]),
            'discarding_internal_weak_gauge_charge_is_valid':False},
        'flux_convention':'F_S=-sum of native Y endpoint currents on links crossing the cut',
        'full_Gauss_cut_identity':'sum_(x in S) G_x^Y = rho_source,S^Y + 1_(0 in S)Y_ref - F_S',
        'physical_joint_record_identity':'F_S=rho_source,S^Y+1_(0 in S)Y_ref',
        'matter_current_only_inferred_by_omitting_scalar_memory':False}


@lru_cache(maxsize=1)
def route_certificate():
    terms={term.path[0]:term for term in kinetic_terms(Configuration())}
    edge_data={edge:(left,right) for edge,left,right in graph()[0]}
    half_cuts=[sum(1<<i for i,xyz in enumerate(cells(2)) if xyz[axis]==1) for axis in range(3)]
    records=[];cut_rows=[]
    for route in routes():
        edge=route['edge'];left,right=edge_data[edge]
        origin_end=0 if left==0 else 1
        label=constrained_arrival_link(origin_end)
        # Independent actual full-source compiler column fixes the orientation
        # and scalar factors of the refilling branch, rather than copying a label.
        ec=FAMILIES['ec'][0]
        key=(1<<(30*route['arrival']+ec),(0,)*24,(0,)*8)
        image=pw.hopping_column(key,(terms[edge],))
        assert image and sum(s.conjugate(v)*v for v in image.values())==s.Rational(1,16)
        for (mask,links,scalars),amplitude in image.items():
            assert links[edge]==label and sum(i!=0 for i in links)==1
            assert mask.bit_count()==1 and (mask.bit_length()-1)//30==0
            assert scalars==(0,)*8 and amplitude!=0
        link_labels=[0]*24;link_labels[edge]=label
        all_cuts=[]
        for mask in range(256):
            actual=cut_flux_six(mask,link_labels)
            expected=6*((mask&1)-((mask>>route['arrival'])&1))
            assert actual==expected
            all_cuts.append(str(s.Rational(actual,6)))
        records.append({**route,'origin_endpoint':origin_end,'forced_full_link_basis_index':label,
            'forced_link_label':list(labels()[label]),'full_joint_origin_Gauss_kernel_dimension':1,
            'all_twelve_current_constraints_checked_at_both_endpoints':True,
            'origin_reference_Y_and_Q':'1','origin_source_matter_plus_scalar_Y_and_Q':'0',
            'arrival_source_matter_plus_scalar_Y_and_Q':'-1',
            'raw_origin_link_endpoint_Y_and_Q':'-1','raw_arrival_link_endpoint_Y_and_Q':'1',
            'arrival_enclosing_outward_Gauss_flux':'-1',
            'three_half_cube_fluxes':[str(s.Rational(cut_flux_six(mask,link_labels),6)) for mask in half_cuts],
            'all_256_cut_fluxes':all_cuts})
    for mask in range(256):
        boundary=cut_endpoints(mask);counts=boundary_spectrum(len(boundary))
        cut_rows.append({'cell_mask':mask,'boundary':list(map(list,boundary)),
            'boundary_edges':len(boundary),'spectrum_key':str(len(boundary)),
            'spectral_multiplicity_sum':str(sum(counts.values()))})
    assert len(records)==6
    return {'actual_routes':records,'all_route_cut_identities':6*256,
        'half_cube_cell_masks':half_cuts,'all_256_cuts':cut_rows,
        'full_boundary_flux_spectral_multiplicities_by_edge_count':{
            str(n):{str(k):str(v) for k,v in sorted(boundary_spectrum(n).items())}
            for n in sorted({row['boundary_edges'] for row in cut_rows})},
        'spectral_outcome_convention':'integer k records F_S=k/6; disconnected boundary components and both periodic routes are retained',
        'boundary_spectra_are_ambient_admissible_outcomes_not_physical_probability_predictions':True}
