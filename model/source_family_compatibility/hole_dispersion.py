"""Exact aligned-hole observables of the declared signed Hubbard source.

The sector and centered clock are specified inputs. The native repair is
expressed in the same signed-Hubbard occupation chart. No physical particle
assignment or full-ground selection is made.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, prod

import numpy as np
from scipy import sparse

from model.field_symmetry import symmetry
from model.full_fock_field.field import FieldGeometry
from model.strong_coupling_refinement import algebra
from model.defect_field_bridge.full_domain import branch_action, RAW_SLOTS


def box(sides):
    if (type(sides) is not tuple or len(sides) != 3
            or any(type(n) is not int or n < 1 for n in sides)):
        raise ValueError('three positive integer side lengths required')
    points = tuple(product(*(range(n) for n in sides)))
    if len(points) > 4096:
        raise ValueError('explicit one-body allocation exceeds 4096 sites')
    lookup = {x: i for i, x in enumerate(points)}
    edges = []
    for i, x in enumerate(points):
        for j in range(3):
            y = list(x); y[j] += 1
            if tuple(y) in lookup:
                edges.append((i, lookup[tuple(y)], symmetry.link_sign(x, j), j))
    return points, tuple(edges)


def directional_hopping(sides):
    """Integer matrices A_j with h=(c/2) sum_j A_j, from actual link_sign."""
    points, edges = box(sides)
    result = []
    for axis in range(3):
        rows, cols, entries = [], [], []
        for x, y, sign, j in edges:
            if j == axis:
                rows.extend((x, y)); cols.extend((y, x)); entries.extend((sign, sign))
        result.append(sparse.csr_matrix((np.array(entries, dtype=np.int64), (rows, cols)),
                                       shape=(len(points), len(points))))
    return tuple(result)


def _zero(matrix, message):
    matrix = matrix.tocsr(); matrix.eliminate_zeros()
    if matrix.nnz:
        raise ArithmeticError(message)


def box_identity(sides):
    points, edges = box(sides)
    directions = directional_hopping(sides)
    for a in range(3):
        for b in range(a):
            _zero(directions[a] @ directions[b]+directions[b] @ directions[a],
                  'actual signed directional hoppings do not anticommute')
    h2 = sum(directions, sparse.csr_matrix((len(points),)*2, dtype=np.int64))
    square = sum((a @ a for a in directions), sparse.csr_matrix(h2.shape, dtype=np.int64))
    _zero(h2 @ h2-square, 'the source square differs from its separated sum')
    parity = sparse.diags([(-1)**sum(x) for x in points], dtype=np.int64)
    _zero(parity @ h2+h2 @ parity, 'the source bipartite spectral symmetry failed')
    plaquettes = 0
    table = {(x,y): s for x,y,s,_ in edges}
    lookup = {x:i for i,x in enumerate(points)}
    for x in points:
        for a,b in combinations(range(3),2):
            xa=list(x); xa[a]+=1; xb=list(x); xb[b]+=1
            xab=xa.copy(); xab[b]+=1
            if all(tuple(y) in lookup for y in (xa,xb,xab)):
                i,j,k,l=(lookup[tuple(y)] for y in (x,xa,xb,xab))
                if table[i,j]*table[j,l]*table[i,k]*table[k,l] != -1:
                    raise ArithmeticError('an actual source plaquette lost its pi flux')
                plaquettes += 1
    return {'sides':list(sides),'sites':len(points),'edges':len(edges),
            'negative_flux_plaquettes':plaquettes,'directional_anticommutators_zero':True,
            'h_squared_separates_exactly':True,'bipartite_spectral_pairing':True}


def _hole_column(sites, holes):
    """f_i1* ... f_iq* |filled>, f_i*=d_(i,up), with increasing labels."""
    mask = sum(1 << (2*x) for x in range(sites)); sign = 1
    for x in reversed(holes):
        mode = 2*x
        sign *= (-1)**(mask & ((1 << mode)-1)).bit_count()
        mask ^= 1 << mode
    return mask, sign


def hole_lift(sides, holes):
    """Exact all-column reduction of the original full fixed-number source."""
    sites, edges = algebra.geometry(sides)
    if type(holes) is not int or not 0 <= holes <= sites:
        raise ValueError('a hole count in the complete aligned sector required')
    masks, T2, doublons = algebra.operators(sites, edges, sites-holes)
    index={mask:i for i,mask in enumerate(masks)}
    labels=tuple(combinations(range(sites),holes))
    rows, signs=zip(*(_hole_column(sites, label) for label in labels),strict=True)
    inclusion=sparse.csr_matrix((signs,([index[mask] for mask in rows],range(len(labels)))),
                               shape=(len(masks),len(labels)),dtype=np.int64)
    one=sum(directional_hopping(sides),sparse.csr_matrix((sites,)*2,dtype=np.int64))
    expected=sparse.lil_matrix((len(labels),)*2,dtype=np.int64)
    hole_masks=tuple(sum(1 << x for x in label)for label in labels)
    inverse={mask:i for i,mask in enumerate(hole_masks)}
    one=one.tocoo()
    for j,mask in enumerate(hole_masks):
        for destination,origin,sign in zip(one.row,one.col,one.data,strict=True):
            image=algebra.hop(mask,int(origin),int(destination))
            if image:
                output,car=image
                expected[inverse[output],j] -= int(sign)*car
    expected=expected.tocsr()
    _zero(T2@inclusion-inclusion@expected,'full Hubbard hopping differs from hole exterior power')
    _zero(sparse.diags(doublons,dtype=np.int64)@inclusion,'contact does not vanish on the aligned sector')
    _zero(inclusion.T@inclusion-sparse.eye(len(labels),dtype=np.int64),
          'hole CAR inclusion is not an isometry')
    return {'sides':list(sides),'holes':holes,'source_number':sites-holes,
            'whole_number_dimension':len(masks),'aligned_dimension':len(labels),
            'full_source_invariant':True,'contact_restriction_zero':True,
            'relative_centered_energy_shift':'q U/2',
            'relative_hopping_is_minus_h_transpose_exterior_lift':True}


def frame_audit():
    """The actual Weyl frame transports the Hamiltonian and must transport repair."""
    g=FieldGeometry((2,1,1))
    source=symmetry.operators(g)
    if not (source.frame@source.H@source.frame.adjoint()-source.hubbard_H).is_zero():
        raise ArithmeticError('actual original Weyl action is not the signed Hubbard action')
    # d-up at site 1 is original c-down. A bare c-flavour hop reaches c-down
    # at site 0, which is d-down, outside the selected d-up sector.
    frame=source.frame.matrix().toarray()
    aligned=[1,4]
    original=frame.conj().T[:,4]
    mask=int(np.flatnonzero(original)[0])
    raw,out,sign,_=branch_action(mask,'E',1,sites=2,flavors=2,seam=(0,1))
    transported=sign*original[mask]*frame[:,out]
    if raw!=6 or np.linalg.norm(transported[aligned])!=0:
        raise ArithmeticError('bare original-flavour repair unexpectedly preserves hidden alignment')
    return {'actual_original_Weyl_frame_identity':True,
            'original_input_mask':mask,'bare_original_repair_output_mask':out,'raw':raw,
            'hidden_aligned_output_weight':0,
            'repair_chart_attachment':'signed Hubbard and full_domain in the same declared occupation chart',
            'untransformed_original_Weyl_repair_claimed':False}


def gamma_matrices():
    """Actual eight parity-cell Clifford matrices, with every cell component."""
    parities=tuple(product(range(2),repeat=3)); lookup={x:i for i,x in enumerate(parities)}
    result=[]
    for j in range(3):
        g=np.zeros((8,8),dtype=np.int64)
        for i,x in enumerate(parities):
            y=list(x);y[j]^=1
            g[i,lookup[tuple(y)]]=symmetry.link_sign(x,j)
        result.append(g)
    for i,g in enumerate(result):
        if not np.array_equal(g@g,np.eye(8,dtype=int)):
            raise ArithmeticError('source parity-cell Clifford square failed')
        for h in result[:i]:
            if np.any(g@h+h@g):raise ArithmeticError('source Clifford cross relation failed')
    return tuple(result)


def bloch_h(k,speed=1.):
    k=np.asarray(k,dtype=float)
    if k.shape!=(3,)or not np.isfinite(k).all() or not np.isfinite(speed)or speed<0:
        raise ValueError('three finite momenta and nonnegative speed required')
    return speed*sum(np.cos(k[j])*g for j,g in enumerate(gamma_matrices()))


def band_identification(k):
    """Canonical Clifford rotor from the k=0 positive band to the k band."""
    k=np.asarray(k,dtype=float)
    if k.shape!=(3,) or not np.isfinite(k).all() or np.linalg.norm(k)>.5:
        raise ValueError('finite momentum in the radius-one-half chart required')
    gamma=gamma_matrices();v=np.cos(k);n=v/np.linalg.norm(v);n0=np.ones(3)/np.sqrt(3)
    a=sum(n[j]*g for j,g in enumerate(gamma));a0=sum(n0[j]*g for j,g in enumerate(gamma))
    rotor=(np.eye(8)+a@a0)/np.sqrt(2*(1+np.dot(n,n0)))
    return rotor,(np.eye(8)+a)/2,(np.eye(8)+a0)/2


def finite_gap_interval(sides,speed=F(1,8),interaction=F(1)):
    """Rational enclosure of the infinite gap and its finite-box excess."""
    if type(sides)is not tuple or len(sides)!=3 or any(type(n)is not int or n<1 for n in sides):
        raise ValueError('three positive side lengths required')
    speed,interaction=F(speed),F(interaction)
    if speed<0 or interaction<0:raise ValueError('nonnegative source coefficients required')
    # Squaring the endpoints proves this outward rational sqrt(3) enclosure.
    lo,hi=F(1732050807568877,10**15),F(1732050807568878,10**15)
    if not lo*lo<3<hi*hi:raise ArithmeticError('sqrt(3) enclosure failed')
    excess=speed/lo*sum((F(22,7*(n+1)))**2 for n in sides)
    return {'bulk_gap_lower':interaction/2-speed*hi,
            'bulk_gap_upper':interaction/2-speed*lo,
            'finite_box_gap_minus_bulk_upper':excess,
            'exact_finite_gap':'U/2-c sqrt(sum_j cos(pi/(L_j+1))^2)',
            'physical_mass':None,'reference':'filled hidden-up aligned state; K=H-U N/2'}


def continuum_budget(a,radius,duration,speed=F(1,8),holes=1):
    a,radius,duration,speed=map(F,(a,radius,duration,speed))
    if (a<=0 or radius<0 or a*radius>F(1,2)or speed<0
            or type(holes)is not int or holes<0):
        raise ValueError('small positive spacing, bounded momentum support and hole count required')
    return {'one_hole_dynamics_operator_error':min(F(2),abs(duration)*speed*a*a*radius**4/3),
            'q_hole_dynamics_operator_error':min(F(2),holes*abs(duration)*speed*a*a*radius**4/3),
            'one_hole_position_derivative_error':2*a*a*radius,
            'kinetic_coefficient':'c/(2 sqrt(3))',
            'clock':'source model duration T/a^2; fixed-q scalar gap phase removed',
            'relativistic_or_interacting_limit':False}


def translation_identity(side=6):
    """Exact interior action of actual magnetic translations on an open box.

    Boundary shifts are partial isometries. Only columns padded by two
    layers enter the identities; no finite periodic source is substituted.
    """
    if type(side)is not int or not 5<=side<=12:
        raise ValueError('an explicit padded open cube of side five through twelve required')
    points,_=box((side,)*3);lookup={x:i for i,x in enumerate(points)};v=len(points)
    shifts=[]
    for j in range(3):
        rows,cols,values=[],[],[]
        for i,x in enumerate(points):
            y=list(x);y[j]+=1
            if tuple(y)in lookup:
                rows.append(lookup[tuple(y)]);cols.append(i);values.append((-1)**sum(x[:j]))
        shifts.append(sparse.csr_matrix((values,(rows,cols)),shape=(v,v),dtype=np.int64))
    inside=[i for i,x in enumerate(points)if all(1<=a<=side-3 for a in x)]
    embed=sparse.csr_matrix((np.ones(len(inside),dtype=np.int64),(inside,range(len(inside)))),
                           shape=(v,len(inside)))
    h=sum(directional_hopping((side,)*3),sparse.csr_matrix((v,v),dtype=np.int64))
    positions=[sparse.diags([x[j]for x in points],dtype=np.int64)for j in range(3)]
    for j,u in enumerate(shifts):
        _zero((h@u-u@h)@embed,'magnetic translation does not preserve the actual interior action')
        for k,x in enumerate(positions):
            _zero((u.T@x@u-x-(sparse.eye(v,dtype=np.int64)if j==k else 0))@embed,
                  'operational site position fails its interior translation covariance')
        for w in shifts[:j]:
            _zero((u@w+w@u)@embed,'one-step magnetic translation flux sign failed')
            _zero((u@u@w@w-w@w@u@u)@embed,'two-step translations do not commute')
    return {'open_cube_side':side,'interior_columns':len(inside),
            'magnetic_translation_preserves_source_on_interior':True,
            'position_covariance':'U_j* X_i U_j = X_i + delta_ij',
            'one_step_magnetic_translations_anticommute':True,
            'two_step_translations_commute':True,'finite_boundary_discarded':False}


def native_dimer_analyzer():
    """A supplied coherent energy analyzer, complete E records and exports.

    It acts on all sixteen full source input columns. Numerical native
    execution is an admission diagnostic, separate from the exact algebra.
    """
    from model.common_source.chiral_candidate import gaussian_source as gaussian
    from model.observer.source_scalar import compiler
    from model.field_thermal.thermal import execute_with_exports
    from model.field_symmetry.readout import pulse_budget
    swap=gaussian.local_word('fswap')
    hop=gaussian.local_word('hopping',1j,np.pi/4)
    word=[]
    for local,left in((swap,1),(hop,0),(swap,1)):
        word += [{**event,'targets':[left+x for x in event['targets']]}for event in local]
    actual,failures=execute_with_exports(word,4,np.eye(16,dtype=complex))
    if (type(actual)is not np.ndarray or actual.shape!=(16,16)
            or not np.isfinite(actual).all()or type(failures)is not tuple or len(failures)!=len(word)):
        raise ArithmeticError('native analyzer omitted complete output columns or events')
    for event,outputs in zip(word,failures,strict=True):
        dimension=1<<len(event['targets'])
        if (type(outputs)is not tuple or len(outputs)!=6-dimension
                or any(type(x)is not np.ndarray or x.shape!=(16//dimension,16)
                       or not np.isfinite(x).all()for x in outputs)):
            raise ArithmeticError('native analyzer omitted an export complement')
    # Independent exact CAR exponential: the generator satisfies H^3=H.
    from model.full_fock_field.field import FullFock
    fock=FullFock(4);a,b=fock.annihilator(0),fock.annihilator(2)
    generator=(1j*(a.conj().T@b)-1j*(b.conj().T@a)).toarray()
    target=np.eye(16)-1j*np.sin(np.pi/4)*generator+(np.cos(np.pi/4)-1)*(generator@generator)
    complements=[block for event in failures for block in event]
    difference=np.concatenate([actual-target,*complements],axis=0)
    error=float(np.linalg.norm(difference,2));budget,_=pulse_budget(word)
    if error>float(2*budget):raise ArithmeticError('complete native analyzer fails its admitted dilation budget')
    actual_isometry=np.concatenate([actual,*complements],axis=0)
    mass_error=float(np.linalg.norm(actual_isometry.conj().T@actual_isometry-np.eye(16),2))
    if mass_error>2*error+error*error+1e-14:
        raise ArithmeticError('complete native analyzer probability mass is inconsistent')
    hole=np.zeros((16,2),complex)
    for x in range(2):
        mask,sign=_hole_column(2,(x,));hole[mask,x]=sign
    eigenvectors=np.array([[1,1],[1,-1]],complex)/np.sqrt(2)
    evolved=target@hole@eigenvectors
    weights=np.zeros((2,2))
    rows=configuration_scan(2,1,((0,1),))
    for word_record,row in rows.items():
        initial=next(iter(row));hole_position=0 if initial==4 else 1
        weights[hole_position] += np.abs(evolved[initial])**2/2
    if not np.allclose(weights,np.eye(2),atol=2e-15,rtol=0):
        raise ArithmeticError('coherent source analyzer does not read the signed hole energy modes')
    counts={key:value for key,value in compiler.resources(word).items()if type(value)in(int,bool)}
    return {'native_pulses':len(word),'full_source_input_columns':16,
            'all_first_failure_slots':len(complements),'terminal_E_raw_slots':2*RAW_SLOTS,
            'full_dilation_operator_error_interval':['0',str(2*budget)],
            'complete_probability_operator_error_interval':['0',str(4*budget+4*budget*budget)],
            'numerical_admission_passed':True,
            'lower_and_upper_one_hole_mode_position_probabilities':[[1,0],[0,1]],
            'resources':counts,'coherent_analyzer_service_supplied':True,
            'floating_native_execution_is_outward_certificate':False}


def configuration_scan(sites,holes,edges):
    """All nonzero full-domain E words; absent 32^m entries are zero maps."""
    if type(sites)is not int or not 2<=sites<=8 or type(holes)is not int or not 0<=holes<=sites:
        raise ValueError('bounded complete aligned sector required')
    edges=tuple(edges)
    known={0}
    for a,b in edges:
        if not 0<=a<b<sites or (a in known)==(b in known):
            raise ValueError('a rooted spanning-tree traversal required')
        known.update((a,b))
    if len(known)!=sites:raise ValueError('all sites must be covered by the tree')
    masks=tuple(_hole_column(sites,label)[0]for label in combinations(range(sites),holes))
    rows={}
    for initial in masks:
        for coins in product(range(2),repeat=len(edges)):
            mask,sgn,word=initial,1,[]
            for edge,coin in zip(edges,coins,strict=True):
                raw,mask,s,_=branch_action(mask,'E',coin,sites=sites,flavors=2,seam=edge)
                sgn*=s;word.append((coin,raw))
            row=rows.setdefault(tuple(word),{})
            row[initial]=(mask,sgn)
    for word,row in rows.items():
        if len(row)!=1:raise ArithmeticError('complete tree did not identify initial aligned configuration')
        initial=next(iter(row)); decoded=decode_configuration(sites,holes,edges,word)
        if initial!=decoded:raise ArithmeticError('record-only position decoder disagrees with source')
    return rows


def decode_configuration(sites,holes,edges,word):
    """Read total occupation from the declared sector, then solve original records."""
    if len(edges)!=len(word):raise ValueError('one complete record per tree seam required')
    shift=[0]*sites;offset={0:0}
    for (a,b),(coin,raw)in zip(edges,word,strict=True):
        if coin not in (0,1)or raw not in (0,4,5,6):raise ValueError('a supported E record required')
        current={0:0,4:0,5:1,6:-1}[raw]
        difference=0 if raw==0 else 2*(current+coin)-1
        delta=difference-shift[a]+shift[b]
        if a in offset:offset[b]=offset[a]-delta
        elif b in offset:offset[a]=offset[b]+delta
        else:raise ValueError('tree traversal lost its known endpoint')
        shift[a]-=current;shift[b]+=current
    base=F(sites-holes-sum(offset.values()),sites)
    occupancy=[base+offset[x]for x in range(sites)]
    if any(x not in(0,1)for x in occupancy):raise ValueError('record has no input in this sector')
    return sum(int(n)<<(2*x)for x,n in enumerate(occupancy))


def many_hole_dilation(sites,holes,edges,e_weight=F(1,2)):
    """Exact complete (V-V0)* (V-V0) on every q-hole aligned input column."""
    e_weight=F(e_weight);edges=tuple(edges)
    if not 0<=e_weight<=1 or not edges:raise ValueError('valid law weight and nonempty seams required')
    labels=tuple(combinations(range(sites),holes))
    masks=tuple(_hole_column(sites,label)[0]for label in labels)
    gram={}
    for edge in edges:
        for law,weight in(('C',1-e_weight),('E',e_weight)):
            for coin in range(2):
                deltas={}
                for mask in masks:
                    raw,out,sign,_=branch_action(mask,law,coin,sites=sites,flavors=2,seam=edge)
                    deltas.setdefault((raw,out),{})[mask]=sign
                    row=deltas.setdefault((0,mask),{});row[mask]=row.get(mask,0)-1
                for row in deltas.values():
                    for a,sa in row.items():
                        for b,sb in row.items():
                            gram[a,b]=gram.get((a,b),F())+weight*F(sa*sb,2*len(edges))
    gram={key:v for key,v in gram.items()if v}
    expected={}
    for mask,label in zip(masks,labels,strict=True):
        crossing=sum((a in label)!=(b in label)for a,b in edges)
        if crossing:expected[mask,mask]=(1+e_weight)*F(crossing,len(edges))
    if gram!=expected:raise ArithmeticError('many-hole complete source dilation differs from cut count')
    return {'holes':holes,'dimension':comb(sites,holes),'maximum_squared_step_difference':max(gram.values(),default=F()),
            'exact_diagonal_cut_formula':True}


def certificate():
    defaults=FieldGeometry()
    if defaults.speed!=.125 or defaults.interaction!=1.:
        raise ArithmeticError('readable source defaults changed')
    scans=[]
    for sites in(2,3,4):
        edges=tuple((x,x+1)for x in range(sites-1))
        for q in range(sites+1):
            rows=configuration_scan(sites,q,edges)
            scans.append({'sites':sites,'holes':q,'input_dimension':comb(sites,q),
                          'allocated_records':(2*RAW_SLOTS)**len(edges),'positive_records':len(rows),
                          'each_initial_configuration_total_weight':1,
                          **many_hole_dilation(sites,q,edges)})
    return {'source_chart':frame_audit(),
            'pi_flux_box_identities':[box_identity(s)for s in((2,1,1),(2,2,1),(2,2,2),(3,4,5))],
            'full_source_exterior_reductions':[hole_lift(s,q)for s in((2,1,1),(2,2,1))for q in range(prod(s)+1)],
            'eight_cell_clifford_dimension':len(gamma_matrices()[0]),
            'operational_position_translation':translation_identity(),
            'native_complete_energy_analyzer':native_dimer_analyzer(),
            'chosen_reference_gap':finite_gap_interval((17,17,17)),
            'complete_position_scans':scans,
            'continuum_budget':continuum_budget(F(1,16),F(2),F(3)),
            'physical_particle_labels':[],'physical_masses':[],
            'full_interacting_ground_reference_selected':False,
            'source_capacity_or_time_trajectory_selected':False}
