"""Exact open-edge obstructions and a declared two-link periodic completion.

The periodic source is a distinct full carrier. Its two oriented links retain
all six actual quantum states each; neither is replaced by an identity value.
"""
from collections import defaultdict
from functools import lru_cache
import sympy as s
from ..model import active_terms, sequence, clean
from . import hamiltonian as opened
from .basis import current_data

PERIODIC_DIMENSION = 2**12 * 6**2 * 5**2


def apply_column(state, action):
    output = defaultdict(lambda: s.Integer(0))
    for key, weight in state.items():
        for target, value in action(key).items():
            output[target] += weight * value
    return clean(output)


def subtract(left, right):
    output = defaultdict(lambda: s.Integer(0), left)
    for key, value in right.items():
        output[key] -= value
    return clean(output)


def translate_mask(mask):
    """Gamma(T): c_left -> c_right, c_right -> -c_left, including CAR signs."""
    if type(mask) is not int or not 0 <= mask < 4096:
        raise ValueError('twelve-mode occupation required')
    left, right = mask & 63, mask >> 6
    phase = (-1)**(right.bit_count() + left.bit_count()*right.bit_count())
    return right | (left << 6), phase


def open_translation(key):
    """A tested candidate only: transpose the single link; it is not a symmetry."""
    mask, link, hl, hr = opened.validate_key(key)
    labels = active_terms()[3]
    a, b = labels[link]
    newmask, phase = translate_mask(mask)
    return {(newmask, labels.index((b, a)), hr, hl): s.Integer(phase)}


def open_edge_obstruction():
    operators = []
    for label in range(1, 6):
        u = s.zeros(6); u[label, 0] = 1
        operators.append(u)
    gram_in = sum((u.H*u for u in operators), s.zeros(6))
    gram_out = sum((u*u.H for u in operators), s.zeros(6))
    assert gram_in == s.diag(5, 0, 0, 0, 0, 0)
    assert gram_out == s.diag(0, 1, 1, 1, 1, 1)
    seed = (1 << 6, 0, 0, 0)
    difference = subtract(apply_column(open_translation(seed), opened.column),
                          apply_column(opened.column(seed), open_translation))
    assert difference[(1 << 6, 1, 0, 0)] == -s.I/2
    return {'link_common_domain_rank':1, 'link_total_range_rank':5,
            'input_Gram_eigenvalues':['5','0','0','0','0','0'],
            'output_Gram_eigenvalues':['0','1','1','1','1','1'],
            'single_link_unitary_or_antiunitary_orientation_reversal':False,
            'literal_AP_exchange_counterexample_seed':list(seed),
            'counterexample_H_T_minus_T_H_coefficient':'-I/2',
            'scope':'translations implementing the stated local CAR shift and scalar swap; arbitrary re-encoded observables are not excluded'}


def spin_raising_column(key):
    mask, link, hl, hr = opened.validate_key(key)
    out = defaultdict(lambda:s.Integer(0))
    for site in range(2):
        for internal in range(3):
            result=sequence(mask,((6*site+internal,True),(6*site+3+internal,False)))
            if result: out[result[0],link,hl,hr] += result[1]
    return clean(out)


def spin_obstruction():
    seed=(1 << 9,0,0,0)
    difference=subtract(apply_column(opened.column(seed),spin_raising_column),
                        apply_column(spin_raising_column(seed),opened.column))
    assert difference[(1,1,0,0)] == -s.I
    return {'conserved_labels':['L','2Sz'],
            'full_internal_spin_SU2_commutes_H':False,
            'Jplus_H_minus_H_Jplus_seed':list(seed),
            'nonzero_output':[1,1,0,0], 'coefficient':'-I',
            'physical_little_group_spin_identified':False}


def validate_periodic_key(key):
    if len(key)!=5 or any(type(v) is not int for v in key):
        raise ValueError('periodic key is (mask,link01,link10,hleft,hright)')
    mask,l01,l10,hl,hr=key
    if not (0<=mask<4096 and 0<=l01<6 and 0<=l10<6 and 0<=hl<5 and 0<=hr<5):
        raise ValueError('key outside full periodic carrier')
    return key


@lru_cache(maxsize=49152)
def oriented_hopping(mask, link, site):
    """+i sigma3/4 for the 01 edge, -i sigma3/4 for the 10 edge."""
    labels=active_terms()[3];out=defaultdict(lambda:s.Integer(0))
    for label,(a,b) in enumerate(labels[1:],1):
        for spin in range(2):
            coefficient=s.I*(-1)**(spin+site)/4
            target_mode=6*site+3*spin+a
            origin_mode=6*(1-site)+3*spin+b
            if link==0:
                result=sequence(mask,((target_mode,True),(origin_mode,False)))
                if result:out[result[0],label]+=coefficient*result[1]
            if link==label:
                result=sequence(mask,((origin_mode,True),(target_mode,False)))
                if result:out[result[0],0]+=s.conjugate(coefficient)*result[1]
    return clean(out)


def periodic_column(key):
    mask,l01,l10,hl,hr=validate_periodic_key(key)
    out=defaultdict(lambda:s.Integer(0))
    out[key]+=int(l01!=0)+int(l10!=0)+int(hl!=0)+int(hr!=0)
    for (target,link),value in oriented_hopping(mask,l01,0).items():
        out[target,link,l10,hl,hr]+=value
    for (target,link),value in oriented_hopping(mask,l10,1).items():
        out[target,l01,link,hl,hr]+=value
    lm,rm=mask&63,mask>>6
    for (target,h),value in opened.local_yukawa(lm,hl).items():
        out[target+(rm<<6),l01,l10,h,hr]+=value
    for (target,h),value in opened.local_yukawa(rm,hr).items():
        out[lm+(target<<6),l01,l10,hl,h]+=value
    return clean(out)


def periodic_translation(key):
    mask,l01,l10,hl,hr=validate_periodic_key(key)
    target,phase=translate_mask(mask)
    return {(target,l10,l01,hr,hl):s.Integer(phase)}


def periodic_current_column(key,site,current,boundary=(0,0)):
    mask,l01,l10,hl,hr=validate_periodic_key(key)
    if site not in (0,1) or current not in range(12):raise ValueError('current index')
    matter,scalar,edge=current_data();out=defaultdict(lambda:s.Integer(0))
    for (a,b),value in matter[current].todok().items():
        for spin in range(2):
            result=sequence(mask,((6*site+3*spin+a,True),(6*site+3*spin+b,False)))
            if result:out[result[0],l01,l10,hl,hr]+=value*result[1]
    for (a,b),value in scalar[current].todok().items():
        if b==(hl if site==0 else hr):
            out[mask,l01,l10,a if site==0 else hl,hr if site==0 else a]+=value
    for link_number,link in enumerate((l01,l10)):
        endpoint=site if link_number==0 else 1-site
        for (a,b),value in edge[endpoint][current].todok().items():
            if b==link:out[mask,a if link_number==0 else l01,l10 if link_number==0 else a,hl,hr]+=value
    if current==11 and boundary[site]:out[key]+=s.sympify(boundary[site])
    return clean(out)


def momentum_projector(state, label):
    """Four complete phase effects for T^4=I; eigenphase k=pi*label/2."""
    if label not in range(4):raise ValueError('momentum phase label 0..3')
    total=defaultdict(lambda:s.Integer(0));power=dict(state)
    for m in range(4):
        coefficient=s.I**(-label*m)/4
        for key,value in power.items():total[key]+=coefficient*value
        power=apply_column(power,periodic_translation)
    return clean(total)


def norm2(state):
    return s.simplify(sum(s.conjugate(v)*v for v in state.values()))


@lru_cache(maxsize=1)
def hopping_grammar_checks():
    """Every actual link matrix unit and CAR bilinear, with spectators retained."""
    labels=active_terms()[3];terms=defaultdict(lambda:s.Integer(0))
    for site in range(2):
        for label,(a,b) in enumerate(labels[1:],1):
            for spin in range(2):
                i,j=6*site+3*spin+a,6*(1-site)+3*spin+b
                value=s.I*(-1)**(spin+site)/4
                for spectator in range(6):
                    row=(label,spectator) if site==0 else (spectator,label)
                    col=(0,spectator) if site==0 else (spectator,0)
                    terms[i,j,row,col]+=value
                    terms[j,i,col,row]+=s.conjugate(value)
    terms=clean(terms);transformed=defaultdict(lambda:s.Integer(0))
    for (i,j,row,col),value in terms.items():
        phase=(-1)**(int(i>=6)+int(j>=6))
        transformed[(i+6)%12,(j+6)%12,row[::-1],col[::-1]]+=phase*value
    assert clean(transformed)==terms
    return {'complete_link_matrix_unit_bilinear_terms':len(terms),
            'all_term_translation_intertwiners_exact':True}


def periodic_checks():
    grammar=hopping_grammar_checks()
    for mask in range(4096):
        m,p=translate_mask(mask);mm,pp=translate_mask(m)
        assert mm==mask and p*pp==(-1)**mask.bit_count()
    # Complete local Yukawa occupation/scalar domain under the full common H,
    # plus occupied spectators and every two-link configuration.
    fixtures={(mask,0,0,h,0) for mask in range(64) for h in range(5)}
    fixtures.update((mask,a,b,(mask+a)%5,(mask+b+1)%5)
                    for mask in (0,1,8,65,521,2047,4095) for a in range(6) for b in range(6))
    for key in fixtures:
        assert apply_column(periodic_translation(key),periodic_column)==apply_column(periodic_column(key),periodic_translation)
        for target,value in periodic_column(key).items():
            assert periodic_column(target).get(key,0)==s.conjugate(value)
    gauge_checks=0
    for key in sorted(fixtures)[::max(1,len(fixtures)//12)]:
        for site in range(2):
            for current in range(12):
                left=apply_column(periodic_current_column(key,site,current),periodic_translation)
                right=apply_column(periodic_translation(key),lambda k:periodic_current_column(k,1-site,current))
                assert left==right;gauge_checks+=1
    seed=(1,0,0,0,0)
    fixed_left=apply_column(periodic_current_column(seed,0,11,(0,-1)),periodic_translation)
    fixed_right=apply_column(periodic_translation(seed),lambda k:periodic_current_column(k,1,11,(0,-1)))
    assert fixed_left!=fixed_right
    moving_right=apply_column(periodic_translation(seed),lambda k:periodic_current_column(k,1,11,(-1,0)))
    assert fixed_left==moving_right
    records=[]
    for label in range(4):
        output=momentum_projector({seed:s.Integer(1)},label)
        assert momentum_projector(output,label)==output
        assert all(not momentum_projector(output,k) for k in range(4) if k!=label)
        assert apply_column(output,periodic_column)==momentum_projector(periodic_column(seed),label)
        records.append(str(norm2(output)))
    assert records==['0','1/2','0','1/2']
    assert sum(s.sympify(p) for p in records)==1
    return {'declared_periodic_full_dimension':PERIODIC_DIMENSION,
            'two_actual_oriented_link_dimensions':[6,6],**grammar,
            'T_squared_equals_matter_parity_masks_checked':4096,
            'complete_local_and_occupied_column_controls':len(fixtures),
            'exact_local_Gauss_translation_controls':gauge_checks,
            'complete_four_phase_readout_on_odd_one_particle_input':records,
            'fixed_charged_boundary_preserved_by_translation':False,
            'translated_boundary_restores_Gauss_intertwiner':True,
            'full_spatial_rotations_or_relativistic_mass_shell_proved':False,
            'momentum_phase_readout_is_finite_and_complete':True,
            'periodic_source_equals_single_open_edge_source':False}



@lru_cache(maxsize=4096)
def conjugating_mask(mask):
    """Unitary part of A=U K, with A c_x A^-1 = i c_(1-x)^dagger."""
    if type(mask) is not int or not 0<=mask<4096:raise ValueError('twelve-mode mask')
    occupied=[j for j in range(12) if mask>>j&1]
    result=sequence(4095,tuple(((j+6)%12,False) for j in occupied))
    assert result is not None
    return result[0],s.expand((-s.I)**len(occupied)*result[1])


def antiunitary_column(key):
    mask,link,hl,hr=opened.validate_key(key)
    labels=active_terms()[3];a,b=labels[link]
    scalar_conjugate=(0,3,4,1,2)
    target,phase=conjugating_mask(mask)
    return {(target,labels.index((b,a)),scalar_conjugate[hr],scalar_conjugate[hl]):phase}


def apply_antiunitary(state):
    output=defaultdict(lambda:s.Integer(0))
    for key,value in state.items():
        for target,phase in antiunitary_column(key).items():
            output[target]+=s.conjugate(value)*phase
    return clean(output)


def source_probe(key,mode,create):
    mask,link,hl,hr=opened.validate_key(key)
    result=sequence(mask,((mode,create),))
    return {} if result is None else {(result[0],link,hl,hr):result[1]}


def anti_hopping_grammar():
    labels=active_terms()[3];terms=defaultdict(lambda:s.Integer(0))
    for label,(a,b) in enumerate(labels[1:],1):
        for spin in range(2):
            i,j=3*spin+a,6+3*spin+b
            value=s.I*(-1)**spin/2
            terms[i,j,label,0]+=value
            terms[j,i,0,label]+=s.conjugate(value)
    transformed=defaultdict(lambda:s.Integer(0))
    transpose=lambda label:labels.index(labels[label][::-1])
    for (i,j,row,col),value in terms.items():
        # The reordered c_(other i) c_(other j)^dagger contributes the minus.
        transformed[(j+6)%12,(i+6)%12,transpose(row),transpose(col)]-=s.conjugate(value)
    assert clean(transformed)==clean(terms)
    return len(terms)


def antiunitary_checks():
    from .basis import local_current_column
    checked=0
    for mask in range(4096):
        target,phase=conjugating_mask(mask)
        back,second=conjugating_mask(target)
        assert back==mask and s.expand(s.conjugate(phase)*second)==(-1)**mask.bit_count()
        assert s.conjugate(phase)*phase==1
        # Exhaust every CAR annihilator column; adjoints follow by isometry.
        for j in range(12):
            source=sequence(mask,((j,False),))
            image=sequence(target,(((j+6)%12,True),))
            assert (source is None)==(image is None)
            if source is not None:
                actual,p=conjugating_mask(source[0])
                assert actual==image[0]
                assert s.expand(source[1]*p-s.I*phase*image[1])==0
            checked+=1
    fixtures={(mask,0,h,0) for mask in range(64) for h in range(5)}
    fixtures.update((mask,link,hl,hr) for mask in (0,1,8,65,521,2047,4095)
                    for link in range(6) for hl,hr in ((0,0),(1,2),(3,4),(2,3)))
    for key in fixtures:
        assert apply_antiunitary(opened.column(key))==apply_column(antiunitary_column(key),opened.column)
        ell,spin=opened.conserved_labels(key)
        target=next(iter(antiunitary_column(key)))
        assert opened.conserved_labels(target)==(4-ell,-spin)
        probe=apply_antiunitary(source_probe(key,6,True))
        reverse=apply_column(antiunitary_column(key),lambda k:source_probe(k,0,False))
        assert probe==clean({k:-s.I*v for k,v in reverse.items()})
    gauge_checks=0
    for key in sorted(fixtures)[::max(1,len(fixtures)//10)]:
        for boundary in ((0,0),(0,-1),(0,1),(-1,0),(1,0)):
            transformed_boundary=(-boundary[1],-boundary[0])
            for site in range(2):
                for current in range(12):
                    left=apply_antiunitary(local_current_column(key,site,current,boundary))
                    right=apply_column(antiunitary_column(key),lambda k:
                        local_current_column(k,1-site,current,transformed_boundary))
                    assert left==clean({k:-v for k,v in right.items()})
                    gauge_checks+=1
    return {'same_open_two_cell_H_antiunitary_symmetry':True,
            'annihilator_map':'A c_(x,a) A^-1 = i c_(1-x,a)^dagger',
            'scalar_map':'particle/antiparticle exchange and site exchange',
            'link_map':'|a,b><0| -> |b,a><0|; raising remains raising',
            'complete_CAR_basis_columns_checked':checked,
            'A_squared_is_matter_parity_masks_checked':4096,
            'all_open_hopping_grammar_terms_checked':anti_hopping_grammar(),
            'full_H_and_probe_fixture_columns_checked':len(fixtures),
            'actual_local_Gauss_with_boundary_checks':gauge_checks,
            'sector_map':{'L':'4-L','2Sz':'-2Sz','boundary':'(q_left,q_right) -> (-q_right,-q_left)'},
            'probe_map':'right addition -> left removal, up to phase -I',
            'same_side_addition_removal_equality_implied':False,
            'matched_transformed_reference_and_complete_records_identical':True,
            'physical_CPT_or_relativistic_particle_identification_claim':False}



@lru_cache(maxsize=4096)
def spin_reversal_mask(mask):
    """Unitary spin part of Theta: up->down, down->-up, all occupied signs."""
    if type(mask) is not int or not 0<=mask<4096:raise ValueError('twelve-mode mask')
    occupied=[j for j in range(12) if mask>>j&1]
    image=[6*(j//6)+(j%6+3)%6 for j in occupied]
    down=sum(j%6>=3 for j in occupied)
    inversions=sum(image[i]>image[j] for i in range(len(image)) for j in range(i+1,len(image)))
    return sum(1<<j for j in image),s.Integer((-1)**(down+inversions))


def spin_reversal_column(key):
    mask,link,hl,hr=opened.validate_key(key)
    target,phase=spin_reversal_mask(mask)
    return {(target,link,hl,hr):phase}


def apply_spin_reversal(state):
    output=defaultdict(lambda:s.Integer(0))
    for key,value in state.items():
        for target,phase in spin_reversal_column(key).items():
            output[target]+=s.conjugate(value)*phase
    return clean(output)


def spin_reversal_checks():
    from .basis import local_current_column
    checked=0
    for mask in range(4096):
        target,phase=spin_reversal_mask(mask);back,second=spin_reversal_mask(target)
        assert back==mask and phase*second==(-1)**mask.bit_count()
        ell=opened.conserved_labels((mask,0,0,0))[0]
        assert (-1)**(ell%2)==(-1)**mask.bit_count()
        for j in range(12):
            image=6*(j//6)+(j%6+3)%6
            sign=1 if j%6<3 else -1
            left=sequence(mask,((j,False),));right=sequence(target,((image,False),))
            assert (left is None)==(right is None)
            if left is not None:
                result,p=spin_reversal_mask(left[0])
                assert result==right[0] and left[1]*p==sign*phase*right[1]
            checked+=1
    fixtures={(mask,0,h,0) for mask in range(64) for h in range(5)}
    fixtures.update((mask,link,hl,hr) for mask in (0,1,8,65,521,2047,4095)
                    for link in range(6) for hl,hr in ((0,0),(1,2),(3,4),(2,3)))
    for key in fixtures:
        assert apply_spin_reversal(opened.column(key))==apply_column(spin_reversal_column(key),opened.column)
        ell,sz=opened.conserved_labels(key)
        target=next(iter(spin_reversal_column(key)))
        assert opened.conserved_labels(target)==(ell,-sz)
    matter,scalar,edges=current_data()
    for current in range(12):
        sign=-1 if current==9 else 1
        for matrix in (matter[current],scalar[current],edges[0][current],edges[1][current]):
            assert matrix.conjugate()==sign*matrix
    gauge_checks=0
    for key in sorted(fixtures)[::max(1,len(fixtures)//10)]:
        for boundary in ((0,0),(0,-1),(0,1),(-1,0),(1,0)):
            for site in range(2):
                for current in range(12):
                    left=apply_spin_reversal(local_current_column(key,site,current,boundary))
                    right=apply_column(spin_reversal_column(key),lambda k:
                        local_current_column(k,site,current,boundary))
                    sign=-1 if current==9 else 1
                    assert left==clean({k:sign*v for k,v in right.items()})
                    gauge_checks+=1
    return {'same_boundary_uniform_H_antiunitary_symmetry':True,
            'source_action':'up -> down; down -> -up; coefficients conjugated',
            'complete_CAR_columns_checked':checked,
            'Theta_squared_equals_fermion_parity_masks_checked':4096,
            'full_H_fixture_columns_checked':len(fixtures),
            'actual_fixed_boundary_Gauss_checks':gauge_checks,
            'current_conjugation_signs':[(-1 if j==9 else 1) for j in range(12)],
            'sector_map':{'L':'L','2Sz':'-2Sz','boundary':'unchanged'},
            'odd_L_energy_degeneracy_at_least_two_exact':True,
            'finite_odd_sector_Kramers_pairs_are_physical_Lorentz_spin_proof':False,
            'real_Yukawa_parameters_required_for_this_unmodified_operator':True}


def physical_periodic_momentum_control():
    """Complete phase readout on a gauge-singlet local pair/scalar preparation."""
    state={(mask,0,0,scalar,0):value/2
           for (mask,scalar),value in opened.local_yukawa(0,0).items()}
    assert norm2(state)==1
    for site in range(2):
        for current in range(12):
            assert not apply_column(state,lambda k:periodic_current_column(k,site,current))
    values=[str(norm2(momentum_projector(state,k))) for k in range(4)]
    assert values==['1/2','0','1/2','0']
    return {'preparation':'normalized local quantum Yukawa action on empty physical input',
            'source_state_is_ground_or_stationary':False,
            'all_24_local_Gauss_actions_annihilate_state':True,
            'complete_phase_probabilities':values,
            'phase_labels':['0','pi/2','pi','3pi/2'],
            'particle_rest_mass_inferred':False}


def certificate():
    return {'open_edge':open_edge_obstruction(),'spin':spin_obstruction(),
            'periodic_completion':periodic_checks(),
            'physical_periodic_readout':physical_periodic_momentum_control(),
            'boundary_charge_conjugation':antiunitary_checks(),
            'fixed_boundary_spin_reversal':spin_reversal_checks()}
