"""Exact admitted two-link virtual packets improve a full odd-sector upper.

Only the specified real-coefficient trial is evaluated. No invariant
compression, interacting pole or prepared-state substitution is inferred.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_composite_transport import model as transport
from model.particle_primary_sector_threshold import model as sectors
from model.particle_primary_odd_resolvent import model as resolvent
from model.particle_spatial_modular_response import regional
from model.particle_primary_gravity_premises import model as modular
from model.particle_primary_composite_spin import model as composite
from model.particle_primary_charge_conjugation import model as conjugation

SOURCE_INPUTS=tuple(dict.fromkeys((
 *('model/particle_primary_branch_certificate/'+name for name in
   ('__init__.py','__main__.py','model.py','test_model.py','PROOF.md','README.md')),
 *transport.SOURCE_INPUTS,*sectors.SOURCE_INPUTS,*resolvent.SOURCE_INPUTS,*regional.SOURCE_INPUTS,
 *modular.SOURCE_INPUTS,*composite.SOURCE_INPUTS,*conjugation.SOURCE_INPUTS,
)))
DELTA=(0,2,0,0,1)


def _move(mask,i,j):
    if not(mask>>j&1) or mask>>i&1:return None
    intermediate=mask^(1<<j)
    parity=(mask&((1<<j)-1)).bit_count()+(intermediate&((1<<i)-1)).bit_count()
    return intermediate|(1<<i),(-1 if parity%2 else 1)


def _gaussian(value):
    """Read bounded Gaussian integers, with no floating tolerance."""
    if not value.real.is_integer() or not value.imag.is_integer() or abs(value.real)+abs(value.imag)>=2**20:
        raise ArithmeticError('exact small Gaussian-integer arithmetic domain exceeded')
    return int(value.real),int(value.imag)


@lru_cache(maxsize=1)
def source_packets():
    cells=tuple(product(range(2),repeat=3))
    edges=tuple((x,y,next(a for a in range(3) if cells[x][a]!=cells[y][a]))
                for x in range(8) for y in range(x+1,8) if sum(a!=b for a,b in zip(cells[x],cells[y]))==1)
    blocks=sectors.native_blocks();native=sectors.refinement.native_control()
    dimensions=tuple(len(indices) for _,indices in blocks);casimirs=tuple(row[1] for row in native)
    signed_charges={sign*q for q,_ in blocks for sign in (-1,1)}
    if len(edges)!=12 or len(signed_charges)!=10:raise ArithmeticError('actual native Haar channels or cube changed')
    omega=sum(1<<(60*x+j) for x in range(8) for j in range(30,60))
    local=sectors.local_singlet(DELTA)
    if len(local['state'])!=1 or local['norm_squared']!=1:raise ArithmeticError('unit physical Slater witness required')
    phi=next(iter(local['state']))|(omega&~((1<<60)-1))
    positive=sum(1<<(60*x+j) for x in range(8) for j in range(30))
    def defects(mask):return (mask&positive).bit_count()+(~mask&omega).bit_count()
    terms=[]
    for x,y,axis in edges:
        row=[];hop=transport.spin_blocks()[axis]*2
        for r,(_,indices) in enumerate(blocks):
            for i,j in product(indices,repeat=2):
                for a,b in product(range(4),repeat=2):
                    value=hop[a,b]
                    if value:
                        z=complex(value)
                        if z not in (1,-1,1j,-1j):raise ArithmeticError('literal doubled Wilson entry is not a Gaussian unit')
                        row.append(((r,1,i,j),60*x+15*a+i,60*y+15*b+j,z))
                        row.append(((r,-1,i,j),60*y+15*b+j,60*x+15*a+i,z.conjugate()))
        terms.append(tuple(row))
    one=[];first_norms=defaultdict(F)
    for row in terms:
        out=defaultdict(complex)
        for label,i,j,z in row:
            moved=_move(phi,i,j)
            if moved:out[label,moved[0]]+=moved[1]*z
        out={key:value for key,value in out.items() if value};one.append(out)
        for (label,mask),value in out.items():
            real,imag=_gaussian(value);band=(label[0],defects(mask)-3)
            first_norms[band]+=F(real*real+imag*imag,4*dimensions[label[0]])
    inherited=sectors.hopping_packets(DELTA)
    for r,row in enumerate(inherited):
        for k,key in ((0,'unchanged_defect_packet_norm_squared'),(2,'two_added_defect_packet_norm_squared')):
            if first_norms[r,k]!=row[key]:raise ArithmeticError('independent full480-mode first-action norm differs')
    norms=defaultdict(F);couplings=defaultdict(F);imaginary=defaultdict(F);pair_rows=[]
    for b in range(12):
        for c in range(b+1,12):
            total=defaultdict(complex);partial=defaultdict(lambda:defaultdict(complex))
            for first,second in ((b,c),(c,b)):
                for (label,mask),coefficient in one[first].items():
                    band=(label[0],defects(mask)-3)
                    for other,i,j,z in terms[second]:
                        moved=_move(mask,i,j)
                        if moved:
                            lb,lc=(label,other) if first==b else (other,label)
                            key=(lb,lc,moved[0]);value=coefficient*z*moved[1]
                            total[key]+=value;partial[band][key]+=value
            total={key:value for key,value in total.items() if value};pair_norm=F();pair_cost=F()
            for (lb,lc,mask),value in total.items():
                real,imag=_gaussian(value);group=(*sorted((lb[0],lc[0])),defects(mask)-3)
                if group[2] not in (0,2,4):raise ArithmeticError('two-link action has an unaccounted defect count')
                weight=F(real*real+imag*imag,16*dimensions[lb[0]]*dimensions[lc[0]])
                norms[group]+=weight;pair_norm+=weight
                pair_cost+=weight*(3*group[2]+F(105,2)*(casimirs[group[0]]+casimirs[group[1]]))
                if mask==phi:raise ArithmeticError('distinct-link second action returned the anchor matter profile')
            for band,vector in partial.items():
                for (lb,lc,mask),value in vector.items():
                    real,imag=_gaussian(value.conjugate()*total.get((lb,lc,mask),0))
                    group=(*sorted((lb[0],lc[0])),defects(mask)-3)
                    denominator=16*dimensions[lb[0]]*dimensions[lc[0]]
                    couplings[group,band]+=F(real,denominator);imaginary[group,band]+=F(imag,denominator)
            adjacent=bool(set(edges[b][:2])&set(edges[c][:2]))
            if pair_norm!=(F(1785,2) if adjacent else F(900)):
                raise ArithmeticError('source two-link pair norm differs from its independently grouped geometry')
            pair_rows.append({'edges':[b,c],'share_vertex':adjacent,'norm_squared':pair_norm,
                              'relative_onsite_electric_numerator':pair_cost,'Haar_CAR_coordinates':len(total)})
    total_norm=sum(norms.values(),F());total_cost=sum((x['relative_onsite_electric_numerator'] for x in pair_rows),F())
    if total_norm!=59220 or total_cost!=F(11350341,2) or sum(couplings.values(),F())!=total_norm:
        raise ArithmeticError('full two-link norm, H0 mean or adjoint-coupling sum differs')
    if any(imaginary.values()):
        raise ArithmeticError('an aggregate actual one-link/two-link hopping entry is not real')
    groups=[]
    for group,norm in sorted(norms.items()):
        if not norm:continue
        cost=F(105,2)*(casimirs[group[0]]+casimirs[group[1]])+3*group[2]
        groups.append({'group':group,'norm_squared':norm,'relative_H0_energy':cost,
          'one_packet_couplings':[{'packet':band,'real':couplings[group,band],'imaginary':imaginary[group,band]}
                                 for band in sorted(first_norms) if first_norms[band]]})
    if len(groups)!=30:raise ArithmeticError('admitted two-link packet count differs')
    return {'native_dimension_order':dimensions,'native_Casimir_order':casimirs,
        'actual_source_modes':480,'actual_native_signed_Haar_channels':10,
        'source_first_packet_norms':[{'packet':band,'norm_squared':n} for band,n in sorted(first_norms.items()) if n],
        'two_distinct_link_pairs':pair_rows,'second_action_groups':groups,
        'total_second_action_norm_squared':total_norm,'total_second_action_relative_H0_numerator':total_cost,
        'all_first_and_second_actions_within_actual_product_cutoff':True,
        'zero_net_matter_profile_in_second_action':False,
        'only_real_trial_coefficients_used':True,
        'uncomputed_full_source_histogram_replaced_by_control':False}


def _solve(matrix,rhs):
    A=[list(map(F,row))+[F(b)] for row,b in zip(matrix,rhs)];n=len(A)
    for j in range(n):
        pivot=next((i for i in range(j,n) if A[i][j]),None)
        if pivot is None:raise ArithmeticError('singular exact Ritz system')
        A[j],A[pivot]=A[pivot],A[j];scale=A[j][j];A[j]=[x/scale for x in A[j]]
        for i in range(n):
            if i!=j:
                factor=A[i][j]
                A[i]=[x-factor*y for x,y in zip(A[i],A[j])]
    return tuple(row[-1] for row in A)


def _positive_pivots(matrix):
    A=[list(map(F,row)) for row in matrix];pivots=[]
    for j in range(len(A)):
        d=A[j][j]
        if d<=0:raise ArithmeticError('the exact excited variational block is not positive')
        pivots.append(d)
        for i in range(j+1,len(A)):
            for k in range(j+1,len(A)):
                A[i][k]-=A[i][j]*A[j][k]/d
    return tuple(pivots)


@lru_cache(maxsize=1)
def trial_certificate():
    controls=source_packets();first=controls['source_first_packet_norms'];groups=controls['second_action_groups']
    native=sectors.refinement.native_control()
    n=[r['norm_squared'] for r in first]
    delta=[F(105,2)*native[r['packet'][0]][1]+3*r['packet'][1] for r in first]
    coupling=[[next(x['real'] for x in group['one_packet_couplings'] if x['packet']==row['packet']) for group in groups] for row in first]
    cn=[group['norm_squared'] for group in groups]
    cd=[group['relative_H0_energy']+F(6,7) for group in groups]
    char=sectors.branch.variational.character_control()['cubic_loop_singlet_matrix']
    loop=[[F(210)*native[i][1]*(i==j)-F(char[i][j],420) for j in range(5)] for i in range(5)]
    g=F(6,22050)
    def coefficients(z):
        z=F(z)
        schur=[[n[i]*(delta[i]+z)*(i==j)-sum((coupling[i][k]*coupling[j][k]/(cn[k]*(cd[k]+z)) for k in range(len(groups))),F())
                for j in range(len(n))] for i in range(len(n))]
        h=_solve(schur,n)
        second=tuple(-sum((coupling[i][k]*h[i] for i in range(len(n))),F())/(cn[k]*(cd[k]+z)) for k in range(len(groups)))
        ell=_solve([[loop[i][j]+z*(i==j) for j in range(5)] for i in range(5)],[1]*5)
        q=sum((n[i]*h[i] for i in range(len(n))),F())+g*sum(ell)
        return h,second,ell,q,schur
    zero=coefficients(0);pivots=_positive_pivots(zero[-1]);_positive_pivots(loop)
    lo,hi=F(),F(8)
    if coefficients(hi)[3]>=hi:raise ArithmeticError('variational root lacks its upper bracket')
    while hi-lo>F(1,2**24):
        mid=(lo+hi)/2
        if coefficients(mid)[3]>=mid:lo=mid
        else:hi=mid
    h,second,ell,q,_=coefficients(lo)
    square=sum((n[i]*h[i]**2 for i in range(len(n))),F())+sum((cn[i]*second[i]**2 for i in range(len(cn))),F())+g*sum(x*x for x in ell)
    shift=(-q-lo*square)/(1+square)
    if not -hi<=shift<=-lo:raise ArithmeticError('rational Rayleigh upper left its exact scalar bracket')
    return {'real_trial_dimension':1+len(n)+len(cn)+5,'excited_block_positive_pivots':pivots,
        'magnetic_upper_on_second_action':F(6,7),'relative_variational_root_bracket':[-hi,-lo],
        'one_packet_real_coefficients':h,'second_action_real_coefficients':second,'loop_real_coefficients':ell,
        'trial_norm_squared':1+square,'source_Rayleigh_quotient_upper':-711+shift,
        'source_Rayleigh_shift_upper':shift,'trial_is_an_eigenvector':False,
        'full_source_omitted_space_assumed_invariant':False,'source_action_or_input_density_changed':False}


@lru_cache(maxsize=1)
def symmetry_reduction():
    """Literal source G1/H matrix coefficients and the complete odd-sector menu."""
    group=composite.group;elements=group.elements();clean=composite.clean
    local=composite.local_multiplet(DELTA);cube=composite.cube_restriction(DELTA)
    if cube['proper_binary_cube_irreps']!={'G1':1,'H':1}:
        raise ArithmeticError('the actual physical local odd multiplet changed')
    central={name:(group.characters(group.negative(group.ONE))[name]/dimension).integer()
             for name,dimension in zip(group.IRREPS,group.DIMENSIONS)}
    odd={name:dimension for name,dimension in zip(group.IRREPS,group.DIMENSIONS) if central[name]==-1}
    if odd!={'G1':2,'G2':2,'H':4}:raise ArithmeticError('the full source odd cube representations changed')
    bases={};representations={}
    for name,value in (('G1',s.Rational(3,4)),('H',s.Rational(15,4))):
        vectors=local['casimir_projectors'][value].columnspace()
        B=s.Matrix.hstack(*s.GramSchmidt(vectors,orthonormal=True));bases[name]=B
        representations[name]=tuple(clean(B.H*D*B) for D in cube['matrices'])
    representations['G2']=tuple(clean(composite.channels._sym(group.characters(q)['A2'])*D)
                                for q,D in zip(elements,representations['G1']))
    # Schur orthogonality is checked on every coefficient of all three
    # spinorial representations. Their source extensions are group sums.
    coefficients=[]
    for name,d in odd.items():
        for a,b in product(range(d),repeat=2):
            coefficients.append((name,a,b,tuple(D[a,b] for D in representations[name])))
    for i,(name,a,b,left) in enumerate(coefficients):
        for j,(other,c,d,right) in enumerate(coefficients):
            inner=s.simplify(s.expand(sum(s.conjugate(x)*y for x,y in zip(left,right))))
            if inner!=(F(48,odd[name]) if i==j else 0):
                raise ArithmeticError('exact binary-cube matrix coefficient orthogonality failed')
    matrices={};diagonal_sum=s.zeros(6)
    for name,d in odd.items():
        for a,b in product(range(d),repeat=2):
            E=clean(sum((s.conjugate(D[a,b])*U for D,U in zip(representations[name],cube['matrices'])),s.zeros(6))*s.Rational(d,48))
            expected=s.zeros(6) if name=='G2' else bases[name][:,a]*bases[name][:,b].H
            if clean(E-expected)!=s.zeros(6):raise ArithmeticError('group sum is not the actual source matrix unit')
            matrices[name,a,b]=E
            if a==b:diagonal_sum+=E
    if clean(diagonal_sum)!=s.eye(6):raise ArithmeticError('the diagonal source matrix-unit menu is incomplete')
    for (name,a,b),E in matrices.items():
        if clean(E.H-matrices[name,b,a])!=s.zeros(6):raise ArithmeticError('matrix-unit adjoint law failed')
    return {'actual_source_proper_group_order':48,'source_anchor_CAR_parity':0,
        'sector_increment_CAR_parity':1,'full_odd_sector_irrep_dimensions':odd,
        'every_full_sector_energy_multiplicity_lower':2,
        'full_sector_rank_one_complement_bottom_equals_sector_bottom':True,
        'actual_local_gauge_kernel_dimension':6,'actual_local_irrep_multiplicities':cube['proper_binary_cube_irreps'],
        'literal_matrix_coefficient_count':len(coefficients),
        'matrix_coefficient_orthogonality_pairs_checked':len(coefficients)**2,
        'source_matrix_unit_formula':'E^lambda_ab=(d_lambda/48) sum_g conjugate(D^lambda_ab(g)) U_g',
        'source_matrix_unit_multiplication':'E^lambda_ab E^mu_cd=delta_lambda_mu delta_bc E^lambda_ad',
        'complete_odd_sector_diagonal_projector_outcomes':sum(odd.values()),
        'actual_local_diagonal_projector_ranks':{name:([1]*d if name!='G2' else [0]*d) for name,d in odd.items()},
        'complete_source_menu_requires_all_native_sectors_and_no_discarded_complement':True,
        'symmetry_matrix_units_commute_with_actual_H':True,
        'compressed_multiplicity_operator':'H|isotypic(lambda)=I_dlambda tensor h_lambda',
        'multiplicity_Feshbach_operator':'PHP-zP-PHQ*(QHQ-zQ)^(-1)*QHP, Q=E^lambda_aa-P',
        'rank_one_obstruction_applies_after_diagonal_symmetry_reduction':False,
        'an_isolated_multiplet_or_scalar_cyclic_pole_excluded':False,
        'reduced_complement_gap_or_pole_computed':False,
        'representations':{name:[[[str(x) for x in D.row(i)] for i in range(D.rows)] for D in representations[name]] for name in odd}}


def selected_energy_floor(member,ground_lower,energy_upper):
    """Actual selected polynomial density, on its lower-energy spectral subspace."""
    slope=modular.selected_modular_slope(member);beta=F(slope['modular_slope_upper'])
    low,high=F(ground_lower),F(energy_upper)
    if high<low:raise ValueError('energy-local floor requires an ordered energy interval')
    bits=member['preparation_plan']['source']['dimension_bits_upper']
    exponent=bits+(2*beta*(high-low)).__ceil__()
    global_floor=regional.common_contract_payload(member,1)['physical_density_floor_factored']
    if exponent>=global_floor['negative_exponent']:raise ArithmeticError('energy-local selected floor failed to sharpen the full-space floor')
    return {'base':2,'negative_exponent':exponent,'operator':'P_G 1_(H<=energy_upper)',
        'energy_upper':high,'global_ground_lower':low,'ambient_dimension_bits_upper':bits,
        'selected_modular_slope_upper':beta,'global_physical_density_floor_negative_exponent':global_floor['negative_exponent'],
        'exact_Gibbs_substituted_for_selected_polynomial':False,
        'joint_source_reference_Loewner_inequality_asserted':False}


def bind_member(member):
    sectors.branch._member_zero(member)
    inherited=sectors.bind_member(member);sectors.verify_bound(member,inherited)
    C=conjugation.bind_member(member);conjugation.verify_bound(member,C)
    controls=source_packets();trial=trial_certificate();symmetry=symmetry_reduction()
    row=inherited['result']['sector_results'][0]['signed_sectors'][0]
    if tuple(row['delta'])!=DELTA:raise ArithmeticError('the actual selected odd occupation sector changed')
    lower=F(row['physical_sector_energy_lower']);upper=trial['source_Rayleigh_quotient_upper']
    old=F(row['trial']['sector_ground_energy_upper'])
    L0=F(inherited['result']['global_ground_lower']);U0=F(inherited['result']['global_ground_upper'])
    if not U0<lower<=upper<old:raise ArithmeticError('the enlarged admitted trial failed the full-source comparison')
    floor=selected_energy_floor(member,L0,upper)
    clock=sectors.spectroscopy.clock_plan(member,'1/64','1/1024','1/1024')
    eta=F(clock['energy_resolution'])
    loss=F(clock['clock_plan']['complete_clock_tail_upper'])+F(clock['one_complete_clock_error_upper'])
    if not loss<F(1,512):raise ArithmeticError('the actual complete one-clock implementation lost its budget')
    ceilings=resolvent.retained_method_ceiling()
    deficit=upper-max(x['comparison_floor_strict_upper'] for x in ceilings['cases'])
    if deficit<=0:raise ArithmeticError('retained scalar complement method comparison changed')
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'same_source_state_reference_clock':True,'source_sector_receipt_sha256':contract.digest(inherited),
        'source_charge_conjugation_receipt_sha256':contract.digest(C),
        'physical_native_increment':DELTA,'opposite_increment':tuple(-n for n in DELTA),
        'full_source_ground_interval':[L0,U0],'full_odd_sector_bottom_interval':[lower,upper],
        'odd_sector_bottom_excitation_interval':[lower-U0,upper-L0],
        'thirteen_coordinate_Rayleigh_upper':old,'Rayleigh_upper_improvement':old-upper,
        'actual_admitted_two_link_controls':controls,'forty_three_coordinate_real_trial':trial,
        'source_binary_cube_reduction':symmetry,
        'retained_scalar_complement_method':ceilings,'Ritz_upper_minus_largest_method_ceiling':deficit,
        'original_selected_energy_local_density_floor':floor,
        'complete_original_final_energy_event':{
            'input_density':'original selected polynomial density',
            'instrument':'all native occupation projectors and one complete final energy clock',
            'exact_source_energy_window':[lower,upper],'reported_source_energy_window':[lower-eta,upper+eta],
            'reported_excitation_window':[lower-eta-U0,upper+eta-L0],
            'complete_clock_plan':clock,'one_complete_clock_loss_upper':loss,
            'per_signed_sector_weight_lower_factored':{**floor,'rational_multiplier':2*(1-loss),'operator':'positive scalar'},
            'either_signed_sector_weight_lower_factored':{**floor,'rational_multiplier':4*(1-loss),'operator':'positive scalar'},
            'all_native_sector_and_energy_records_retained':True,
            'original_quantum_reference_outputs_retained':True,
            'source_marginal_only_probability_bound':True,
            'initial_energy_measurement_performed':False,
            'event_is_an_energy_transfer_or_particle_creation_record':False,
            'probability_divided_by_occupation_success':False},
        'analysis_trial_used_as_original_input':False,
        'finite_selected_sector_final_energy_event_proved':True,
        'isolated_branch_or_resolved_mass_proved':False,
        'uniform_physical_refinement_residue_proved':False,
        'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_branch_certificate.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('actual source trial, symmetry, selected floor or complete event changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    value={'schema':'oph.primary_branch_certificate_receipt.v1',
        'selected_member_bound':bind_member(contract.member(0)),'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('source changed during branch certificate replay')
    return value
