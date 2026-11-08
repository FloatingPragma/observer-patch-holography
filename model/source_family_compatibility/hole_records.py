"""Complete native repair on a mobile hole and its internal spin register."""
from fractions import Fraction as F
from itertools import product

from model.defect_field_bridge.full_domain import branch_action, RAW_SLOTS
from model.strong_coupling_refinement.algebra import operators


def hole_basis(sites):
    if type(sites) is not int or not 2 <= sites <= 8:
        raise ValueError('two through eight explicit sites required')
    rows = []
    for hole in range(sites):
        occupied = [x for x in range(sites) if x != hole]
        for spins in product(range(2), repeat=sites-1):
            mask = sum(1 << (2*x+s) for x, s in zip(occupied, spins))
            rows.append((mask, hole))
    return tuple(rows)


def validate_edges(sites, edges):
    edges = tuple(edges)
    if not edges or len(set(edges)) != len(edges):
        raise ValueError('a nonempty set of undirected seams required')
    for edge in edges:
        if (type(edge) is not tuple or len(edge) != 2
                or any(type(x) is not int for x in edge)
                or not 0 <= edge[0] < edge[1] < sites):
            raise ValueError('ordered source endpoints required')
    return edges


def branches(sites, edge, law):
    validate_edges(sites, (edge,))
    rows = {(coin, raw): {} for coin in range(2) for raw in range(RAW_SLOTS)}
    basis = dict(hole_basis(sites))
    for mask, hole in basis.items():
        for coin in range(2):
            raw, out, sign, _ = branch_action(mask, law, coin, sites=sites,
                                             flavors=2, seam=edge)
            if out not in basis or sign not in (-1, 1):
                raise ArithmeticError('the full source left the one-hole sector')
            rows[coin, raw][mask] = (out, sign)
    for row in rows.values():
        if len({out for out, _ in row.values()}) != len(row):
            raise ArithmeticError('a complete record branch is not injective')
    return rows


def cross_factor(x, y, edge, law):
    if x == y or law not in ('C', 'E'):
        raise ValueError('two different hole positions and a native law required')
    a, b = int(x in edge), int(y in edge)
    return 1-F(a+b, 2) if law == 'C' else F((1-a)*(1-b))


def cross_matrix_unit(sites, edge, law, left, right):
    """Unconditioned output of the actual complete instrument, exactly."""
    out = {}
    for row in branches(sites, edge, law).values():
        if left in row and right in row:
            a, sa = row[left]; b, sb = row[right]
            out[a, b] = out.get((a, b), F())+F(sa*sb, 2)
    return {key: value for key, value in out.items() if value}


def coherence_rate(sites, edges, x, y, gamma=F(1), e_weight=F(1)):
    """Rate in the explicitly supplied independent-seam Poisson clock."""
    edges = validate_edges(sites, edges)
    gamma, e_weight = F(gamma), F(e_weight)
    if (x == y or x not in range(sites) or y not in range(sites)
            or gamma < 0 or not 0 <= e_weight <= 1):
        raise ValueError('valid sites, nonnegative rate and law probability required')
    return gamma*sum((1-e_weight)*(1-cross_factor(x,y,e,'C'))
                     +e_weight*(1-cross_factor(x,y,e,'E')) for e in edges)


def counting_word_multiplier(sites,edges,x,y,attempts):
    """The conditional full-alphabet counting schedule, in native attempts."""
    edges=validate_edges(sites,edges)
    if type(attempts)is not int or attempts<0:
        raise ValueError('a nonnegative literal attempt count required')
    decrement=coherence_rate(sites,edges,x,y,F(1),F(1,2))/len(edges)
    return {'one_attempt_decrement':decrement,
            'attempts':attempts,'complete_coherence_multiplier':(1-decrement)**attempts}


def identity_dilation_gram(sites,edges,e_weight=F(1,2)):
    """Exact (V-V0)* (V-V0), including the seam/law/coin/raw registers."""
    edges=validate_edges(sites,edges);e_weight=F(e_weight)
    if not 0<=e_weight<=1:raise ValueError('a source-law probability required')
    basis=dict(hole_basis(sites));gram={}
    for edge in edges:
        for law,law_weight in (('C',1-e_weight),('E',e_weight)):
            for (coin,raw),row in branches(sites,edge,law).items():
                delta={}
                for initial in basis:
                    if initial in row:
                        output,sign=row[initial]
                        delta.setdefault(output,{})[initial]=sign
                    if raw==0:
                        entries=delta.setdefault(initial,{})
                        entries[initial]=entries.get(initial,0)-1
                for entries in delta.values():
                    for left,sl in entries.items():
                        for right,sr in entries.items():
                            key=(left,right)
                            gram[key]=gram.get(key,F())+law_weight*F(sl*sr,2*len(edges))
    gram={key:value for key,value in gram.items()if value}
    expected={(mask,mask):(1+e_weight)*F(sum(hole in e for e in edges),len(edges))
              for mask,hole in basis.items() if any(hole in e for e in edges)}
    if gram!=expected:
        raise ArithmeticError('complete source dilation differs from its degree formula')
    return gram


def global_attempt_bound(n):
    """Count-only open-cube family; no large source allocation or clock selection."""
    if type(n)is not int or n<2:raise ValueError('integer n at least two required')
    side=n*n+1;edges=3*(side-1)*side*side;attempts=n*n
    squared_step_bound=F(12,edges)
    if squared_step_bound>F(4,n**6):raise ArithmeticError('cube count bound failed')
    return {'n':n,'spatial_scaling':str(F(1,n*n)),'cube_side':side,
            'cube_edges':edges,'maximum_degree':6,'global_attempts':attempts,
            'one_step_dilation_squared_bound':str(squared_step_bound),
            'complete_source_archive_trace_bound':str(min(F(2),F(4,n))),
            'physical_clock_or_capacity_trajectory_selected':False}


def actual_aligned_hubbard_sector():
    rows=[]
    for sites in(2,3,4):
        edges=tuple((x,x+1,1)for x in range(sites-1))
        masks,two_hopping,doubles=operators(sites,edges,sites-1)
        matrix=two_hopping.toarray()
        aligned={sum(1<<(2*x)for x in range(sites)if x!=hole)for hole in range(sites)}
        selected=[i for i,mask in enumerate(masks)if mask in aligned]
        outside=[i for i,mask in enumerate(masks)if mask not in aligned]
        if matrix[outside][:,selected].any() or any(doubles[i]for i in selected):
            raise ArithmeticError('the actual Hubbard source left the aligned hole sector')
        rows.append({'sites':sites,'whole_fixed_occupation_dimension':len(masks),
                     'invariant_aligned_hole_dimension':len(selected),
                     'interaction_D_restriction_zero':True})
    masks,two_hopping,_=operators(3,((0,1,1),(1,2,1)),2)
    entry=int(two_hopping[masks.index(12),masks.index(9)])
    if not entry:raise ArithmeticError('the source doublon leakage control vanished')
    return {'aligned_cases':rows,'general_one_hole_leakage':{
        'input_mask':9,'doublon_output_mask':12,'T_over_c_half_matrix_element':entry},
        'aligned_sector_identified_as_interacting_ground':False}


def scan(sites, edges):
    """Every fine coin/raw record and signed spin transport of an E scan."""
    edges = validate_edges(sites, edges)
    if set(sum(edges, ())) != set(range(sites)):
        raise ValueError('the prescribed seam sequence must cover every site')
    rows = {}
    basis = dict(hole_basis(sites))
    for initial in basis:
        for coins in product(range(2), repeat=len(edges)):
            mask, sign, history = initial, 1, []
            for seam, coin in zip(edges, coins):
                raw, mask, local_sign, _ = branch_action(
                    mask, 'E', coin, sites=sites, flavors=2, seam=seam)
                sign *= local_sign; history.append((coin, raw))
            row = rows.setdefault(tuple(history), {})
            row[initial] = (mask, sign)
    for row in rows.values():
        if len({basis[x] for x in row}) != 1:
            raise ArithmeticError('one complete history failed to identify the initial hole')
        if len({out for out, _ in row.values()}) != len(row):
            raise ArithmeticError('spin transport is not a signed partial permutation')
    return rows


def corrected_coefficient(sites, edges, left, right):
    """Coefficient after inverse signed transport, with every record retained."""
    # The inverse corrects each row's actual CAR sign as well as its permutation.
    return sum((F(1, 2**len(edges)) for row in scan(sites, edges).values()
                if left in row and right in row), F())


def scan_branch(sites,edges,history):
    """Zero-default access to every allocated coin/raw history slot."""
    edges=validate_edges(sites,edges)
    if (len(history)!=len(edges) or any(type(pair)is not tuple or len(pair)!=2
            or type(pair[0])is not int or pair[0]not in (0,1)
            or type(pair[1])is not int or pair[1]not in range(RAW_SLOTS)
            for pair in history)):
        raise ValueError('one exact coin/raw pair per attempted seam required')
    return scan(sites,edges).get(tuple(history),{})


def certificate():
    cases = []
    for sites, edges in ((2,((0,1),)), (3,((0,1),(1,2))),
                         (4,((0,1),(1,2),(2,3)))):
        basis = dict(hole_basis(sites)); checked = 0
        for edge in edges:
            for law in ('C','E'):
                for a,x in basis.items():
                    for b,y in basis.items():
                        if x == y: continue
                        factor = cross_factor(x,y,edge,law)
                        expected = {(a,b):factor} if factor else {}
                        if cross_matrix_unit(sites,edge,law,a,b) != expected:
                            raise ArithmeticError('actual repair violates the coherence formula')
                        checked += 1
        rows = scan(sites,edges)
        # No matrix element or independent reference coefficient is discarded.
        for a,x in basis.items():
            for b,y in basis.items():
                if corrected_coefficient(sites,edges,a,b) != int(x==y):
                    raise ArithmeticError('complete corrected detector differs from hole pinching')
        cases.append({'sites':sites,'source_one_hole_dimension':len(basis),
            'all_input_cross_hole_matrix_units_checked':checked,
            'scan_seams':[list(e) for e in edges],
            'all_allocated_scan_record_slots':(2*RAW_SLOTS)**len(edges),
            'nonzero_scan_records':len(rows),
            'probability_per_coin_word':str(F(1,2**len(edges))),
            'full_internal_spin_dimension_per_hole':2**(sites-1),
            'complete_uniform_repair_identity_dilation_gram':[
                {'mask':a,'squared_difference':str(value)}for(a,b),value in
                sorted(identity_dilation_gram(sites,edges).items())]})
    return {'source':'model.defect_field_bridge.full_domain.branch_action',
        'cases':cases,
        'E_scan_corrected_complete_channel':'sum_x P_x rho P_x',
        'record_erasure_restores_hole_coherence':False,
        'arbitrary_reference_preserved_in_each_internal_spin_block':True,
        'correction_native_pulse_tape_executed':False,
        'zero_history_representation':'absent sparse entries are zero maps; scan_branch covers every allocated slot',
        'conditional_counting_schedule_controls':[
            {k:str(v)if type(v)is F else v for k,v in counting_word_multiplier(
                4,((0,1),(1,2),(2,3),(0,3)),0,y,7).items()}for y in (1,2)],
        'complete_word_simplex_and_counting_reference_premises_discharged':False,
        'global_fair_attempt_complete_record_workaround':[
            global_attempt_bound(n)for n in(2,4,8,16,64)],
        'actual_full_Hubbard_propagation_scope':actual_aligned_hubbard_sector(),
        'Poisson_rate_supplied_not_selected_by_repair':True,
        'one_hole_ground_or_relativistic_particle_assigned':False}
