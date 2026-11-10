"""Native hopping routes and exact atomic-spin arrival coefficients."""
from functools import lru_cache
import sympy as s
from ..periodic_family.model import Configuration,kinetic_terms,cells
from ..peter_weyl_links.clebsch import link
from ..atomic_reference_control.ground_exchange import actual_spin_coefficients
from ..atomic_spin_geometry.cube import spectrum
from ..family_cell.quotient import add,mul,outward
from ..interacting_charged_response.model import FAMILIES
from ...physical_attachment.chiral import representation
from ...chiral_source.model import PAULI


def routes():
    rows=[]
    for term in kinetic_terms(Configuration()):
        if 0 not in (term.left_cell,term.right_cell):continue
        assert len(term.path)==1 and abs(term.coefficient)==s.Rational(1,4)
        target=term.right_cell if term.left_cell==0 else term.left_cell
        assert PAULI[term.axis].H==PAULI[term.axis]
        rows.append({'edge':term.path[0],'origin':0,'arrival':target,
            'arrival_cell':list(cells(2)[target]),'axis':term.axis,
            'actual_term_coefficient':str(term.coefficient),
            'refilling_branch':'forward' if term.left_cell==0 else 'adjoint',
            'selected_edge_full_norm_upper':'15'})
    assert len(rows)==6 and sorted(row['arrival'] for row in rows)==[1,1,2,2,4,4]
    return rows


@lru_cache(maxsize=1)
def spin_effects():
    b=s.Symbol('b',real=True)
    def one_body(i,j):
        return s.eye(2)*s.Rational(int(i==j),2)+b*sum((a[j,i]*a/2 for a in PAULI),s.zeros(2))
    effects=[]
    for axis,sigma in enumerate(PAULI):
        effect=s.zeros(4)
        for spin in range(2):
            for alpha,beta in sigma.todok():
                for alphap,betap in sigma.todok():
                    effect+=s.conjugate(sigma[alphap,betap])*sigma[alpha,beta]*s.kronecker_product(
                        one_body(alphap,spin).H*one_body(alpha,spin),one_body(betap,beta))/32
        exchange=sum(((1 if j==axis else -1)*s.kronecker_product(a/2,a/2)
            for j,a in enumerate(PAULI)),s.zeros(4))
        expected=(1+3*b*b)*s.eye(4)/128+b*b*(1-b)*exchange/16
        assert (effect-expected).applyfunc(s.expand)==s.zeros(4)
        effects.append(effect.applyfunc(s.expand))
    return b,tuple(effects)


def coefficient(family):
    if family not in ('ec','L'):raise ValueError('the proved ec/L charged source channels are required')
    actual=actual_spin_coefficients();beta=tuple(map(s.Rational,actual['beta_ec_interval']))
    b=beta if family=='ec' else ((1-beta[1])/2,(1-beta[0])/2)
    x=tuple(4*s.Rational(value) for value in spectrum()['ground_energy_interval'])
    b2=mul(b,b)
    result=add(mul((s.Rational(1,128),)*2,add((1,1),mul((3,3),b2))),
        mul((s.Rational(1,768),)*2,mul(mul(b2,(1-b[1],1-b[0])),x)))
    assert result[0]>0
    return {'b_interval':b,'cube_4H_ground_interval':x,'moment_interval':result}


def certificate():
    spin_effects();records=[];hyper=representation()[2][-1]
    # Every native pair has nonzero total hypercharge, excluding a second
    # annihilation at the origin from its neutral local ground projector.
    assert all(hyper[a,a]+hyper[b,b]!=0 for a in range(15) for b in range(15))
    for family in ('ec','L'):
        indices=FAMILIES[family];d=len(indices);columns=[]
        for a in indices:
            for b in indices:columns.extend((link(a,b)[:,0],link(a,b).H[:,0]))
        matrix=s.ImmutableSparseMatrix.hstack(*columns)
        assert matrix.H*matrix==s.eye(2*d*d)/d
        data=coefficient(family)
        charge=int(6*hyper[indices[0],indices[0]])
        records.append({'family':family,'native_internal_dimension':d,
            'vacuum_link_Gram_dimension':2*d*d,'exact_vacuum_link_Gram':'I/d_r',
            'origin_reference_six_Y':charge,'arrival_matter_six_Y':-charge,
            'flux_origin_endpoint_six_Y':-charge,'flux_arrival_endpoint_six_Y':charge,
            'source_raw_L_change':1 if family=='ec' else -1,
            'origin_source_raw_L_on_arrival':1,'arrival_source_raw_L':2 if family=='ec' else 0,
            'all_other_source_cells_raw_L':1,'reference_raw_L':-1 if family=='ec' else 1,
            'reference_location':'origin throughout the clock evolution',
            'source_arrival_gauge_representation':'conjugate of the native removal field irrep',
            'flux_origin_representation':'conjugate of the retained reference irrep',
            'flux_arrival_representation':'native removal field irrep',
            'b_definition':'beta_ec' if family=='ec' else '(1-beta_ec)/2',
            'b_interval':outward(data['b_interval'],60),
            'exact_two_spin_arrival_effect':'(1+3b^2)I/128 + b^2(1-b) S_origin dot R_axis S_arrival/16',
            'exact_cube_ground_arrival_moment':'(1+3b^2)/128 + b^2(1-b)x/768, x=4 epsilon_cube',
            'arrival_moment_interval':outward(data['moment_interval'],64)})
    return {'actual_kinetic_routes':routes(),'complete_three_axis_spin_matrix_identities':3,
        'native_hypercharge_excludes_all_origin_pair_annihilation_compressions':True,
        'only_same_native_family_can_refill_the_origin_ground':True,
        'reference_moved_to_arrival_cell':False,
        'same_source_flux_route_coefficients':records}
