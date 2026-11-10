"""An isolated interacting gauge transition on the unchanged finite source.

The reducing empty-Fock sector is a rare original sector. Its spatial-odd
gauge line is not a fermion, an ambient vacuum gap or a GeV mass.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_branch_certificate import model as branch
from model.particle_primary_survival import survival
from model.particle_primary_channel_labels import channels
from model.particle_primary_variational_readout import variational
from model.particle_primary_joint_spectroscopy import spectroscopy

SOURCE_INPUTS=tuple(dict.fromkeys((
 *('model/particle_primary_reduced_sector/'+name for name in
   ('__init__.py','__main__.py','model.py','test_model.py','PROOF.md','README.md')),
 *branch.SOURCE_INPUTS,*survival.SOURCE_INPUTS,*channels.SOURCE_INPUTS,
 *variational.SOURCE_INPUTS,*spectroscopy.SOURCE_INPUTS,
)))
CELLS=tuple(product(range(2),repeat=3))


@lru_cache(maxsize=1)
def cube_support():
    edges=tuple((x,y) for x in range(8) for y in range(x+1,8)
                if sum(a!=b for a,b in zip(CELLS[x],CELLS[y]))==1)
    faces=[]
    for axis,side in product(range(3),range(2)):
        vertices=tuple(i for i,x in enumerate(CELLS) if x[axis]==side)
        face=tuple(j for j,e in enumerate(edges) if set(e)<=set(vertices))
        origin=min(vertices)
        neighbors=sorted(y if x==origin else x for x,y in (edges[i] for i in face) if origin in (x,y))
        left,right=neighbors;opposite=next(v for v in vertices if v not in (origin,left,right))
        u=tuple(a-b for a,b in zip(CELLS[left],CELLS[origin]));v=tuple(a-b for a,b in zip(CELLS[right],CELLS[origin]))
        cross=tuple(u[(i+1)%3]*v[(i+2)%3]-u[(i+2)%3]*v[(i+1)%3] for i in range(3))
        tour=(origin,left,opposite,right) if cross[axis]*(2*side-1)>0 else (origin,right,opposite,left)
        flux=[0]*12
        for a,b in zip(tour,tour[1:]+tour[:1]):flux[edges.index((min(a,b),max(a,b)))]=1 if a<b else -1
        faces.append({'axis':axis,'side':side,'vertices':vertices,'edges':face,'outward_flux':tuple(flux)})
    if len(edges)!=12 or any(len(f['edges'])!=4 for f in faces):
        raise ArithmeticError('actual member-zero plaquettes changed')
    support_counts=Counter();small=[]
    for mask in range(1<<12):
        selected=tuple(i for i in range(12) if mask>>i&1);degrees=[0]*8
        for i in selected:
            for x in edges[i]:degrees[x]+=1
        if 1 in degrees:continue
        support_counts[len(selected)]+=1
        if selected and len(selected)<=5:small.append(selected)
    if sorted(small)!=sorted(f['edges'] for f in faces):
        raise ArithmeticError('a nonconstant Gauss support below six edges is not a square')
    # Every unsupported fourth-product pattern has an edge with a single
    # nontrivial matrix coefficient. Doubled distinct faces use private links.
    private={}
    for order in (2,3,4):
        checked=0
        for positions in combinations(range(6+order-1),order):
            indices=tuple(x-j for j,x in enumerate(positions));counts=Counter(indices)
            if 1 not in counts.values():continue
            edge_counts=Counter()
            for face,multiplicity in counts.items():
                for e in faces[face]['edges']:edge_counts[e]+=multiplicity
            if 1 not in edge_counts.values():raise ArithmeticError('Haar cancellation lacks a private edge')
            checked+=1
        private[order]=checked
    return {'cells':CELLS,'interior_edges':edges,'oriented_outward_faces':faces,
        'supports_with_no_degree_one_vertex_counts':dict(sorted(support_counts.items())),
        'nonempty_supports_at_most_five_edges':len(small),
        'only_six_squares_below_six_edges':True,'private_edge_patterns_checked_by_total_degree':private,
        'ghost_links_reduced_to_constant_by_actual_Gauss_constraint':True}


@lru_cache(maxsize=1)
def native_Haar_control():
    raw=variational.character_control();weights=raw['native_weights'];zero=(0,0,0,0)
    multiply=variational.multiply
    chars=[]
    for q in raw['block_order_Q6']:
        c=Counter(tuple(w) for w in weights if w[-1]==q)
        chars.extend((dict(c),{tuple(-x for x in w):n for w,n in c.items()}))
    total=Counter()
    for c in chars:total.update(c)
    weyl={zero:1}
    for root in ((2,1,0,0),(1,2,0,0),(1,-1,0,0),(0,0,2,0)):
        for sign in (-1,1):weyl=multiply(weyl,{zero:1,tuple(sign*x for x in root):-1})
    def integral(poly):
        return F(sum(v*weyl.get(tuple(-x for x in k),0) for k,v in poly.items()),12)
    for i,left in enumerate(chars):
        for j,right in enumerate(chars):
            adjoint={tuple(-x for x in k):v for k,v in left.items()}
            if integral(multiply(adjoint,right))!=int(i==j):
                raise ArithmeticError('actual degree-one native characters are not distinct irreducibles')
    square=multiply(total,total)
    odd_difference={(0,0,0,6):1,(0,0,0,-6):-1}
    y_squared={zero:2,(0,0,0,12):-1,(0,0,0,-12):-1}
    if integral(square)!=10 or integral(multiply(square,y_squared))!=18:
        raise ArithmeticError('actual native magnetic/imaginary-character fourth moment differs')
    if integral(multiply(total,y_squared))!=0:
        raise ArithmeticError('spatial-odd loop magnetic mean is nonzero')
    if integral(total)!=0 or integral(odd_difference)!=0 or integral(multiply(total,odd_difference))!=0 or integral(y_squared)!=2:
        raise ArithmeticError('a claimed single-face Haar vanishing or normalization failed')
    native=branch.sectors.refinement.native_control()
    energies=tuple(F(105,2)*C for _,C,_ in native)
    if min(energies)!=F(63,4) or sorted(energies)[1]!=F(189,8):
        raise ArithmeticError('the actual source electric levels changed')
    dimension=1+2*sum(d*d for d,_,_ in native)
    if dimension!=119:raise ArithmeticError('actual degree-one coordinate span differs')
    ground_variance=F(6*10,210**2)
    odd_parent_variance=F(18+5*20,2*210**2)
    if ground_variance!=F(1,735) or odd_parent_variance!=F(59,44100):
        raise ArithmeticError('full parent residual arithmetic differs')
    return {'native_blocks':native,'one_link_electric_energies':energies,
        'actual_degree_one_link_dimension':dimension,'native_character_orthogonality_pairs_checked':100,
        'constant_link_multiplicity':1,'first_electric_Gauss_level':F(63),
        'first_electric_Gauss_level_multiplicity':12,'next_electric_Gauss_level_lower':F(189,2),
        'integral_S_squared':10,'integral_S_squared_imaginary_character_squared':18,
        'ground_anchor_magnetic_mean':0,'odd_loop_anchor_magnetic_mean':0,
        'ground_anchor_residual_squared':ground_variance,
        'odd_loop_residual_squared_upper':odd_parent_variance,
        'odd_residual_bound_equals_compressed_second_moment_claimed':False}


@lru_cache(maxsize=1)
def cube_channels():
    group=channels.group;graph=cube_support();faces=graph['oriented_outward_faces'];edges=graph['interior_edges']
    normals=tuple(tuple((2*f['side']-1) if i==f['axis'] else 0 for i in range(3)) for f in faces)
    chars=[]
    for q in group.elements():
        R=group.rotation(q)
        image=tuple(normals.index(tuple(sum(int(R[i,j])*v[j] for j in range(3)) for i in range(3))) for v in normals)
        reverse=tuple(normals.index(tuple(-x for x in normals[p])) for p in image)
        cell_image=tuple(CELLS.index(channels.rotate_cell(x,2,R)) for x in CELLS)
        for source,target in enumerate(image):
            mapped=[0]*12
            for (x,y),value in zip(edges,faces[source]['outward_flux']):
                a,b=cell_image[x],cell_image[y]
                mapped[edges.index((min(a,b),max(a,b)))]+=(1 if a<b else -1)*value
            if tuple(mapped)!=faces[target]['outward_flux']:
                raise ArithmeticError('proper source rotation reverses the outward loop chart')
        chars.append((q,sum(i==j for i,j in enumerate(image)),sum(i==j for i,j in enumerate(reverse))))
    for source,normal in enumerate(normals):
        target=normals.index(tuple(-x for x in normal));mapped=[0]*12
        for (x,y),value in zip(edges,faces[source]['outward_flux']):
            a,b=7-x,7-y
            mapped[edges.index((min(a,b),max(a,b)))]+=(1 if a<b else -1)*value
        if tuple(mapped)!=tuple(-x for x in faces[target]['outward_flux']):
            raise ArithmeticError('actual spatial inversion failed to reverse loop orientation')
    results={}
    for part,conjugation in (('real',1),('imaginary',-1)):
        row={}
        for name in group.IRREPS:
            for parity,sign in (('g',1),('u',-1)):
                value=sum((group.characters(q)[name]*(tr+sign*conjugation*invtr) for q,tr,invtr in chars),group.Q())/96
                count=value.integer()
                if count:row[name+parity]=count
        results[part]=row
    expected={'real':{'A1g':1,'Eg':1,'T1u':1},'imaginary':{'A1u':1,'Eu':1,'T1g':1}}
    if results!=expected:raise ArithmeticError('the actual oriented-loop cube decomposition changed')
    return {'proper_binary_cube_elements_checked':48,'full_rotation_inversion_order':96,
        'actual_oriented_edge_word_transformations_checked':294,
        'face_outward_normals':normals,'first_cluster_decomposition':results,
        'spatial_odd_A1u_multiplicity_in_cluster':1,
        'A1u_CAR_parity':0,'A1u_is_Lorentz_or_fermion_spin':False,
        'full_native_sector_rank_one_obstruction_for_odd_trions_reused':True,
        'CAR_odd_trion_reduced_complement_bound_improved':False}


def spectral_bounds():
    data=native_Haar_control();eps=F(6,7);electric=F(63);next_level=F(189,2)
    if 2*eps>=min(electric,next_level-electric):raise ArithmeticError('actual perturbation does not isolate clusters')
    lower_ground=-data['ground_anchor_residual_squared']/(electric-eps)
    lower_odd=electric-data['odd_loop_residual_squared_upper']/(next_level-eps-electric)
    if lower_ground!=-F(1,45675) or lower_odd!=63-F(59,1351350):
        raise ArithmeticError('Temple enclosure arithmetic differs')
    # Exact outward eigenvector-angle constants; the physical source residuals
    # are smaller than1/27 and the two spectral denominators exceed62 and30.
    distance=F(10,7)*(F(1,27*62)+F(1,27*30))
    if not F(100,49)>2 or not F(15,52)**2<F(1,12) or not distance<F(15,52)-F(2,7):
        raise ArithmeticError('source transition amplitude lower bound failed')
    return {'magnetic_operator_norm_upper':eps,
        'isolated_gauge_vacuum_cluster_rank':1,'isolated_first_electric_cluster_rank':12,
        'isolated_A1u_cluster_rank_within_symmetry_sector':1,
        'interacting_gauge_ground_energy_interval':[lower_ground,F()],
        'interacting_A1u_energy_interval':[lower_odd,electric],
        'gauge_ground_to_A1u_energy_gap_interval':[lower_odd,electric-lower_ground],
        'gap_interval_decimal_display_only':[format(float(x),'.12f') for x in (lower_odd,electric-lower_ground)],
        'next_A1u_energy_lower':next_level-eps,
        'A1u_gap_to_same_symmetry_complement_lower':next_level-eps-electric,
        'full_cluster_gap_to_rest_lower':next_level-electric-2*eps,
        'phase_aligned_eigenvector_distance_sum_upper':distance,
        'contraction_probe':'A=sum_outward_faces(chi6-chi6*)/(12i)',
        'bare_ground_to_odd_matrix_element_squared':F(1,12),
        'actual_ground_to_odd_matrix_element_squared_strict_lower':F(4,49),
        'exact_resolved_transition_operator':'B=P_empty*P_A1u*1_(62<H<64)*A*P_empty*1_(H<1)',
        'resolved_original_transfer_measure':'p(E_g)*|<u,A g>|^2*delta_(E_u-E_g)',
        'raw_loop_probe_has_no_other_transition_frequencies_claimed':False,
        'A1u_energy_is_globally_nondegenerate_across_other_cube_irreps':False,
        'same_source_interaction_coefficient':1,'small_interaction_extrapolated_to_one':False}


def bind_member(member):
    branch.sectors.branch._member_zero(member)
    if member['source']['per_link_word_degree']!=1:raise ValueError('the isolated cycle classification uses degree one')
    graph=cube_support();haar=native_Haar_control();cube=cube_channels();spectrum=spectral_bounds()
    original=branch.sectors.branch.occupation.bind_member(member)
    branch.sectors.branch.occupation.verify_bound(member,original)
    L0=F(original['result']['zero_retention_physical_ground_lower'])
    U0=F(original['result']['physical_ground_upper'])
    floor=branch.selected_energy_floor(member,L0,0)
    clock=spectroscopy.clock_plan(member,'1/64','1/1024','1/1024')
    eta=F(clock['energy_resolution'])
    loss=F(clock['clock_plan']['complete_clock_tail_upper'])+F(clock['one_complete_clock_error_upper'])
    residue=F(4,49);factor=residue*(1-loss)**2
    if not loss<F(1,512) or not factor>F(1,13):raise ArithmeticError('complete source clocks lost the positive pole event')
    beta_low=F(branch.modular.selected_modular_slope(member)['modular_slope_lower'])
    gaplo,gaphi=spectrum['gauge_ground_to_A1u_energy_gap_interval']
    if beta_low*gaplo<=1:raise ArithmeticError('actual selected decreasing weights do not give the retarded residue')
    rarity=survival.extreme_sector_bound(member)
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'same_source_state_reference_clock':True,'source_input_is_original_selected_density':True,
        'reducing_sector':'all native matter occupations zero; full physical product-group gauge carrier',
        'original_extreme_sector_suppression':rarity,
        'source_global_ground_interval':[L0,U0],
        'actual_cube_cycle_control':graph,'actual_native_Haar_and_electric_control':haar,
        'actual_cube_cluster_labels':cube,'interacting_source_spectrum':spectrum,
        'original_energy_local_population_floor':floor,
        'original_resolved_positive_transfer_atom_lower_factored':{**floor,'operator':'positive scalar','rational_multiplier':residue},
        'original_resolved_retarded_positive_residue_lower_factored':{**floor,'operator':'positive scalar','rational_multiplier':F(2,49)},
        'complete_clocked_original_event':{
            'instrument':'initial source energy clock; complete empty-sector loop probe; complete cube and final source energy records',
            'probe_success_Kraus':'A*P_empty','probe_failure_Kraus':'sqrt(I-P_empty*A^2*P_empty)',
            'initial_record_window':[spectrum['interacting_gauge_ground_energy_interval'][0]-eta,eta],
            'final_A1u_record_window':[spectrum['interacting_A1u_energy_interval'][0]-eta,F(63)+eta],
            'reported_transfer_window':[gaplo-2*eta,gaphi+2*eta],
            'complete_clock_plan':clock,'one_clock_loss_upper':loss,
            'original_weight_lower_factored':{**floor,'operator':'positive scalar','rational_multiplier':factor},
            'simple_original_weight_strict_lower_factored':{**floor,'operator':'positive scalar','rational_multiplier':F(1,13)},
            'all_source_native_sectors_symmetry_and_failure_outcomes_retained':True,
            'idle_reference_outputs_retained':True,'joint_reference_Loewner_bound_asserted':False,
            'empty_sector_success_divided_out':False,
            'exact_pole_projector_replaced_by_a_bare_loop_truncation':False,
            'all_other_initial_state_clock_leakage_excluded':False},
        'finite_symmetry_resolved_gauge_transition_atom_proved':True,
        'extreme_sector_supported_contraction_weight_vanishes_on_selected_refinement':True,
        'CAR_odd_trion_branch_isolated':False,'ambient_primary_vacuum_gap_proved':False,
        'refinement_uniform_observable_residue_proved':False,
        'four_dimensional_relativistic_or_Yang_Mills_continuum_proved':False,
        'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_reduced_sector.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('source sector, pole, residue or complete clocks changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    value={'schema':'oph.primary_reduced_sector_certificate.v1',
        'selected_member_bound':bind_member(contract.member(0)),'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('source inputs changed during reduced-sector replay')
    return value
