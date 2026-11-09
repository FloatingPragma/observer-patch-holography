"""Exact primary pair multiplets, complete channel records and original bands."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_channel_labels import channels
from model.particle_primary_response_moments import model as response
from model.particle_primary_occupation_bound import model as occupation
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_primary_mass_readout import readout


SOURCE_INPUTS=(
    'model/particle_primary_pair_channels/__init__.py',
    'model/particle_primary_pair_channels/__main__.py',
    'model/particle_primary_pair_channels/model.py',
    'model/particle_primary_pair_channels/test_model.py',
    'model/particle_primary_pair_channels/PROOF.md',
    'model/particle_primary_response_moments/PROOF.md',
    'model/particle_primary_occupation_bound/PROOF.md',
    'model/particle_primary_channel_labels/PROOF.md',
    'model/particle_primary_momentum/PROOF.md',
    'model/particle_primary_joint_spectroscopy/PROOF.md',
    'model/particle_primary_mass_readout/PROOF.md',
    'model/particle_gauss_ground_state/PROOF.md',
    'Lean/QFT/PeriodicMomentum.lean',
    'Lean/QFT/GroundRecordTail.lean',
)


@lru_cache(maxsize=1)
def pair_basis():
    """Actual Clifford chart; each two-spinor coefficient is unitary/sqrt2."""
    pairs=response.beta_pairs()
    positive=s.Matrix.hstack(*(p for p,_ in pairs))
    negative=s.Matrix.hstack(*(m for _,m in pairs))
    gamma=response.density.gamma_data()
    sigma=(s.eye(2),)+tuple(response.apparatus.clean(-s.I*positive.H*gamma[(j+1)%3]*gamma[(j+2)%3]*positive)
                            for j in range(3))
    matrices=tuple(response.apparatus.clean(positive*t*negative.H/s.sqrt(2)) for t in sigma)
    for t in sigma:
        if t.H*t!=s.eye(2):raise ArithmeticError('the spin coefficients are not unitary')
    for a in range(4):
        for b in range(4):
            if s.simplify(s.trace(matrices[a].H*matrices[b]))!=int(a==b):
                raise ArithmeticError('the pair coefficients are not orthonormal')
    return sigma,matrices


def car_bilinear(modes,i,j):
    """Exact small Fock operator, retaining every number sector."""
    if type(modes) is not int or not 1<=modes<=8:raise ValueError('small exact CAR control required')
    if not (0<=i<modes and 0<=j<modes):raise ValueError('CAR index outside the control')
    out=s.zeros(1<<modes)
    for mask in range(1<<modes):
        if not mask&(1<<j):continue
        sign=(-1)**((mask&((1<<j)-1)).bit_count());reduced=mask^(1<<j)
        if reduced&(1<<i):continue
        sign*=(-1)**((reduced&((1<<i)-1)).bit_count())
        out[reduced|(1<<i),mask]=sign
    return out


@lru_cache(maxsize=1)
def local_controls():
    sigma,_=pair_basis();identity=s.eye(16)
    number=[car_bilinear(4,j,j) for j in range(4)]
    defect=number[0]+number[1]+2*identity-number[2]-number[3]
    probes=tuple(response.apparatus.clean(sum((t[a,b]*car_bilinear(4,a,2+b)/s.sqrt(2)
                           for a in range(2) for b in range(2)),s.zeros(16))) for t in sigma)
    for b in probes:
        effect=response.apparatus.clean(b.H*b)
        if any(value<0 for value in effect.eigenvals()) or any(value<0 for value in (identity-effect).eigenvals()):
            raise ArithmeticError('a full-Fock pair is not a contraction')
        if any(value<0 for value in (effect-identity+defect/2).eigenvals()):
            raise ArithmeticError('the pair defect inequality failed')
        if response.apparatus.clean(b.H*b-b*b.H)!=identity-defect/2:
            raise ArithmeticError('the full-Fock pair commutator differs')
    return probes,defect


@lru_cache(maxsize=1)
def coefficient_representation():
    """Physical proper cube action and beta inversion on all32 pair slots."""
    cells=tuple(product(range(2),repeat=3));locate={x:j for j,x in enumerate(cells)}
    _,pairs=pair_basis();gamma=response.density.gamma_data();beta=gamma[3]
    rotations=[];seen=set()
    for q in channels.group.elements():
        raw=channels.group.rotation(q);key=tuple(int(x) for x in raw.flat)
        if key in seen:continue
        seen.add(key);R=s.Matrix(raw);spin=channels.spin_lift(q)
        internal=s.diag(1,R)
        for a in range(4):
            actual=response.apparatus.clean(spin*pairs[a]*spin.H)
            expected=sum((internal[b,a]*pairs[b] for b in range(4)),s.zeros(4))
            if actual!=expected:raise ArithmeticError('actual Clifford pair covariance differs')
        representation=s.zeros(32)
        for x,cell in enumerate(cells):
            y=locate[channels.rotate_cell(cell,2,raw)]
            representation[4*y:4*y+4,4*x:4*x+4]=internal
        rotations.append((q,representation))
    inversion=s.zeros(32)
    for x,cell in enumerate(cells):
        y=locate[tuple(1-z for z in cell)]
        inversion[4*y:4*y+4,4*x:4*x+4]=-s.eye(4)
    if any(response.apparatus.clean(beta*b*beta)!=-b for b in pairs):
        raise ArithmeticError('pair inversion parity is not odd')
    for g in gamma[:3]:
        hop=(-beta+beta*g)/2
        if response.apparatus.clean(beta*hop*beta-hop.H)!=s.zeros(4):
            raise ArithmeticError('inversion does not map a hopping to its reversed edge')
    if len(rotations)!=24 or inversion**2!=s.eye(32):raise ArithmeticError('cube action coverage differs')
    return rotations,inversion


@lru_cache(maxsize=1)
def cubic_projectors():
    rotations,inversion=coefficient_representation();sigma,_=pair_basis()
    projectors={}
    for parity,sign in (('g',1),('u',-1)):
        for name,dimension in zip(channels.group.IRREPS[:5],channels.group.DIMENSIONS[:5]):
            proper=sum((s.Integer(channels.group.characters(q)[name].integer())*R
                        for q,R in rotations),s.zeros(32))*s.Rational(dimension,24)
            p=(s.eye(32)+sign*inversion)*proper/2
            rank=int(s.trace(p))
            if p.T!=p or p*p!=p or rank not in (1,2,6):
                raise ArithmeticError('the exact isotypic projector changed')
            # This diagonal contraction is the full aggregate CAR commutator.
            for x in range(8):
                local=sum((p[4*x+a,4*x+b]*sigma[a].H*sigma[b]/2
                           for a in range(4) for b in range(4)),s.zeros(2))
                if response.apparatus.clean(local)!=s.Rational(rank,16)*s.eye(2):
                    raise ArithmeticError('the isotypic defect coefficient is not uniform')
            projectors[name+parity]=p
    if sum(projectors.values(),s.zeros(32))!=s.eye(32):raise ArithmeticError('cubic channels omit a slot')
    for name,p in projectors.items():
        for other,q in projectors.items():
            if name!=other and p*q!=s.zeros(32):raise ArithmeticError('cubic channels overlap')
    zero_scalar=s.kronecker_product(s.ones(8)/8,s.diag(1,0,0,0))
    if projectors['A1u']!=zero_scalar:raise ArithmeticError('the zero-wavevector scalar channel differs')
    return projectors


@lru_cache(maxsize=2)
def coefficient_transform(kind='cubic'):
    """Exact32-slot unitary with every success and failure label retained."""
    if kind=='cubic':
        rows=[];labels=[]
        for name,p in cubic_projectors().items():
            basis=s.GramSchmidt(p.columnspace(),orthonormal=True)
            rows.extend(v.T for v in basis)
            labels.extend((name,j) for j in range(len(basis)))
        transform=s.Matrix.vstack(*rows)
    elif kind=='walsh':
        cells=tuple(product(range(2),repeat=3));transform=s.zeros(32);labels=[]
        for k,frequency in enumerate(cells):
            for a in range(4):
                labels.append((frequency,a))
                for x,position in enumerate(cells):
                    transform[4*k+a,4*x+a]=(-1)**sum(u*v for u,v in zip(frequency,position))/s.sqrt(8)
    else:raise ValueError('cubic or walsh coefficient transform required')
    if response.apparatus.clean(transform*transform.H)!=s.eye(32):
        raise ArithmeticError('complete channel coefficient transform is not unitary')
    return transform,tuple(labels)


def symmetry_summary():
    rows=[]
    for name,p in cubic_projectors().items():
        irrep=name[:-1];dimension=channels.group.DIMENSIONS[channels.group.IRREPS.index(irrep)]
        rank=int(s.trace(p))
        rows.append({'channel':name,'irrep_dimension':dimension,'isotypic_rank':rank,
            'multiplicity':rank//dimension,'diagonal_defect_coefficient':F(rank,16),
            'projector_entries':[[str(x) for x in row] for row in p.tolist()]})
    return {'proper_rotation_order':24,'rotation_and_inversion_order':48,
        'local_pair_internal_representation':'A1u + T1u',
        'complete_coefficient_dimension':32,'isotypic_channels':rows,
        'zero_Walsh_scalar_channel':'A1u','zero_Walsh_vector_channel':'T1u subspace',
        'coefficient_to_operator_map_injective_on_physical_anchor':True,
        'physical_anchor_is_selected_state':False,
        'operator_irrep_is_initial_or_final_state_irrep':False,
        'continuum_spin_or_species_assigned':False}


def passive_window(commutator_lower,beta,center,moment,lower,upper,mass_upper=1):
    """Original spectral band from the actual polynomial log-slope premise.

    Each input bound needs an authenticated source consumer. This arithmetic
    does not establish those premises for an unrelated source or operator.
    """
    values=(commutator_lower,beta,center,moment,lower,upper,mass_upper)
    if any(type(x) not in (int,str,F) for x in values):raise ValueError('exact rational spectral bounds required')
    chi,b,c,m,lo,hi,cap=map(F,values)
    if min(chi,b,lo,cap)<=0 or m<0 or chi>cap or hi<=max(lo,c):
        raise ValueError('positive commutator/slope/window and valid moment/mass bounds required')
    low=b*lo*cap;high=m/(hi-c)**2
    return {'original_transfer_interval':[lo,hi],'commutator_lower':chi,
            'polynomial_log_slope_upper':b,'original_total_mass_upper':cap,
            'centered_second_moment_center':c,'centered_second_moment_upper':m,
            'low_positive_frequency_allowance_upper':low,
            'high_positive_frequency_allowance_upper':high,
            'original_transfer_weight_lower':max(F(),chi-low-high),
            'success_probability_normalized':False}


def selected_channels(member):
    contract.verify_member(member)
    if member['source']['level']!=0:raise ValueError('the occupation-supported channel bands are member-zero results')
    inherited=occupation.bind_member(member);occupation.verify_bound(member,inherited)
    response.source_probe_control();local_controls();coefficient_transform('cubic');coefficient_transform('walsh')
    w=F(inherited['result']['average_matched_q6_pair_weight_lower'])
    beta=F(member['state']['beta']);epsilon=F(member['state']['target_ground_excess'])
    if w<F(11,16) or beta<=0:raise ArithmeticError('actual selected-state weight or cooling changed')
    delta=F(11,2048)/beta;upper=F(40);moment=F(6)
    low=beta*delta;high=moment/(upper-6)**2
    if high>F(11,2048):raise ArithmeticError('the high-frequency allowance failed')
    clock=spectroscopy.clock_plan(member,delta/4,'1/4096','1/4096')
    clockloss=2*F(clock['clock_plan']['complete_clock_tail_upper'])+F(clock['two_complete_clock_error_upper'])
    if clockloss>F(1,1024):raise ArithmeticError('the channel clock budget failed')
    rows=[]
    for name,p in cubic_projectors().items():
        rank=int(s.trace(p));comm=rank*w/32
        window=passive_window(comm,beta,6,moment,delta,upper)
        band=window['original_transfer_weight_lower']
        simple=F(11*rank,1024)
        if band<simple:raise ArithmeticError('the original isotypic band does not pass its rational floor')
        final=readout.excitation_window((delta,upper),band,epsilon,8)
        rows.append({'channel':name,'isotypic_rank':rank,
            'passive_window':window,
            'aggregate_original_commutator_lower':comm,
            'original_success_weight_lower':comm,
            'original_positive_transfer_interval':[delta,upper],
            'original_transfer_weight_lower':band,'original_transfer_weight_simple_lower':simple,
            'original_final_excitation':final,
            'original_final_excitation_weight_simple_lower':F(11*rank-8,1024),
            'finite_reported_transfer_interval':[delta/2,upper+delta/2],
            'finite_reported_original_weight_lower':band-clockloss,
            'finite_reported_original_weight_simple_lower':F(11*rank-1,1024),
            'weight_is_isotypic_aggregate':True,'multiplicity_resolved_weight_claimed':rank in (1,2)})
    # Evaluate the existing unit adapter at a unit tick only to obtain the
    # coefficient multiplying 1/tau_seconds; no absolute tick is selected.
    unit=readout.declared_SI_comparison((delta,upper),seconds_per_source_time=(1,1))
    return contract.serial({'schema':'oph.primary_pair_channels.v1',
        'member_id':member['member_id'],**member['ids'],
        'occupation_receipt_sha256':contract.digest(inherited),
        'selected_density_log_slope_upper':beta,'log_slope_formula':'beta*T_(m-1)(x)/T_m(x), x>=0',
        'exact_Gibbs_state_substituted':False,'polynomial_order':member['state']['Taylor_order'],
        'symmetry':symmetry_summary(),'per_local_Pauli_probe_norm_upper':1,
        'spatially_averaged_each_Pauli_success_lower':w,
        'spatially_averaged_each_Pauli_centered_second_moment_upper':moment,
        'spatially_averaged_each_Pauli_transfer_interval':[2,10],
        'spatially_averaged_each_Pauli_transfer_weight_lower':w-F(3,8),
        'spatially_averaged_each_Pauli_transfer_weight_simple_lower':F(5,16),
        'each_Walsh_pair_record_original_success_lower':w/32,
        'each_Walsh_pair_record_positive_frequency_weight_lower':w/32,
        'each_Walsh_pair_record_passive_window':passive_window(w/32,beta,6,moment,delta,upper),
        'each_Walsh_pair_record_positive_band_weight_simple_lower':F(11,1024),
        'coherent_instrument_success_records':32,'coherent_instrument_failure_records':32,
        'fine_channel_Kraus_formula':'K_(b,j)=sum_i W_(j,i) C_(b,i)/sqrt(32); C_(0,i)=C_i, C_(1,i)=sqrt(I-C_i^dagger C_i)',
        'cubic_parity_channels':rows,'channel_low_frequency_allowance_upper':low,
        'channel_high_frequency_allowance_upper':high,
        'finite_clock_plan':clock,'complete_energy_channel_history_outcomes':64*clock['clock_plan']['clock_size']**2,
        'energy_interval_GeV_times_seconds_per_source_tick_coefficient':unit['energy_GeV_interval'],
        'converted_quantity_is_transfer_energy_not_rest_mass':True,
        'absolute_seconds_per_source_tick':None,'unit_conversion_owner':'particle_primary_mass_readout.readout.declared_SI_comparison',
        'reference_idle_and_retained':True,'initial_energy_readout_preserves_source_marginal':True,
        'initial_energy_readout_preserves_joint_purification':False,
        'same_source_action_and_energy_zero':True,'original_weights_conditioned_on_success':False,
        'full_primary_spectral_matrix_or_histogram_executed':False,
        'Fourier_label_is_conserved_particle_momentum':False,
        'isolated_pole_threshold_or_dispersion_proved':False,
        'uniform_refinement_weight_or_fixed_lower_endpoint_proved':False,
        'physical_mass_table_GeV':[]})


def bind_member(member):
    value=selected_channels(member)
    return {'result':value,'result_sha256':contract.digest(value)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('the actual pair channels or original bounds changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*occupation.SOURCE_INPUTS,*response.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    value={'schema':'oph.primary_pair_channel_certificate.v1','selected_member':bind_member(contract.member(0)),
           'source_sha256':pins,'physical_mass_table_GeV':[]}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('pair-channel inputs changed during replay')
    return value
