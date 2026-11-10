"""Source energy-current response to a complete local lapse-control pulse."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import numpy as np

from model.particle_common_contract import contract, provenance
from model.particle_primary_covariant_stress import model as covariant
from model.particle_primary_gravity_premises import model as modular
from model.particle_primary_sector_threshold import model as sectors
from model.particle_primary_composite_transport import model as transport
from model.particle_spatial_modular_response import regional

SOURCE_INPUTS=tuple('model/particle_primary_energy_response/'+n for n in
    ('__init__.py','__main__.py','model.py','test_model.py','PROOF.md','README.md'))+(
    'model/particle_primary_covariant_stress/PROOF.md',
    'model/particle_primary_gravity_premises/PROOF.md',
    'model/particle_primary_locality/PROOF.md',
    'Lean/QFT/CoarseEnergyWard.lean',
    'Lean/QFT/PathTimeSliceInterface.lean')


def move(mask,i,j):
    if not(mask>>j&1) or mask>>i&1:return None
    left=mask^(1<<j)
    parity=(mask&((1<<j)-1)).bit_count()+(left&((1<<i)-1)).bit_count()
    return left|(1<<i),(-1 if parity%2 else 1)


@lru_cache(maxsize=1)
def current_witness():
    """Literal120-mode one-word projection of [H,h'_0] on a physical flip."""
    blocks=sectors.native_blocks();q6=next(indices[0] for charge,indices in blocks if charge==6)
    anchor=sum(1<<(60*x+j) for x in range(2) for j in range(30,60))
    phi=(anchor|(1<<q6))^(1<<(30+q6))
    rows=[]
    for axis in range(3):
        values=defaultdict(complex);hop=2*transport.spin_blocks()[axis]
        for r,(charge,indices) in enumerate(blocks):
            for i,j in product(indices,repeat=2):
                for a,b in product(range(4),repeat=2):
                    z=complex(hop[a,b])
                    if not z:continue
                    if z not in (1,-1,1j,-1j):raise ArithmeticError('source doubled hopping entry changed')
                    for target,origin,coefficient,orientation in (
                        (15*a+i,60+15*b+j,z,1),
                        (60+15*b+j,15*a+i,z.conjugate(),-1)):
                        result=move(phi,target,origin)
                        if result is None:continue
                        state,sign=result
                        # onsite eigenvalues are 3 beta per particle. The
                        # half-link electric contribution cancels exactly.
                        change0=3*((1 if target//15<2 else -1) if target<60 else 0)
                        change0-=3*((1 if origin//15<2 else -1) if origin<60 else 0)
                        change1=3*((1 if (target-60)//15<2 else -1) if target>=60 else 0)
                        change1-=3*((1 if (origin-60)//15<2 else -1) if origin>=60 else 0)
                        factor=F(change1-change0,2)
                        values[(r,orientation,i,j,state)]+=complex(factor)*coefficient*sign
        values={k:v for k,v in values.items() if v}
        norm=F()
        for (r,orientation,i,j,state),z in values.items():
            if not z.real.is_integer() or not z.imag.is_integer():
                raise ArithmeticError('exact doubled current left Gaussian integers')
            norm+=F(int(z.real)**2+int(z.imag)**2,4*len(blocks[r][1]))
        if norm!=F(9,2) or len(values)!=2 or any(blocks[k[0]][0]!=6 for k in values):
            raise ArithmeticError('source one-link current witness differs')
        rows.append({'axis':axis,'physical_one_link_coordinates':len(values),
                     'projected_commutator_norm_squared_at_a_one':norm})
    return {'actual_modes_in_two_cell_chart':120,'charge_six_native_index':q6,
        'local_flip_preserves_Gauss_and_native_occupations':True,
        'other_cells_are_unchanged_even_filled_beta_singlets':True,
        'axes':rows,'corner_projected_commutator_norm_squared_coefficient':sum(r['projected_commutator_norm_squared_at_a_one'] for r in rows),
        'spacing_power':-4,'projection_is_within_actual_link_cutoff':True,
        'witness_is_the_original_preparation':False,
        'full_commutator_replaced_by_projected_commutator':False}


def q6_two_cell_control(axis):
    """Independent finite Gauss sector for testing the projected current."""
    hop=np.array(transport.spin_blocks()[axis],complex)
    states=[(mask,2-(mask&15).bit_count()) for mask in range(256)
            if mask.bit_count()==4 and abs(2-(mask&15).bit_count())<=1]
    address={state:i for i,state in enumerate(states)};d=len(states)
    T=np.zeros((d,d),complex);onsite=np.zeros(d);local=np.zeros(d);electric=np.zeros(d)
    beta=(1,1,-1,-1)
    for col,(mask,flux) in enumerate(states):
        local[col]=3*sum(beta[k] for k in range(4) if mask>>k&1)
        onsite[col]=local[col]+3*sum(beta[k] for k in range(4) if mask>>(4+k)&1)
        electric[col]=F(63,4)*flux*flux
        for a,b in product(range(4),repeat=2):
            for target,origin,z,shift in ((a,4+b,hop[a,b],-1),(4+b,a,hop[a,b].conjugate(),1)):
                if not z:continue
                moved=move(mask,target,origin)
                if moved and (moved[0],flux+shift) in address:
                    T[address[moved[0],flux+shift],col]+=moved[1]*z
    H=np.diag(onsite+electric)+T;A=np.diag(local+electric/2)+T/2
    phi=((1<<2)|(1<<3)|(1<<6)|(1<<7)|(1<<0))^(1<<2)
    v=np.zeros(d);v[address[phi,0]]=1
    J=1j*(H@A-A@H);out=J@v
    one=np.array([flux!=0 for mask,flux in states]);norm=float(np.vdot(out[one],out[one]).real)
    return H,A,J,v,{'dimension':d,'axis':axis,'one_word_current_norm_squared':norm,
        'domain':'two-cell charge-six control; not the full selected source'}


def factored(prefactor,exponent):
    p=F(prefactor)
    if p<=0 or exponent<0:raise ValueError('positive factored bound required')
    return {'rational_multiplier':p,'base':2,'negative_exponent':int(exponent)}


def bind_member(member):
    contract.verify_member(member)
    energy=covariant.cell_energy(member);K=energy['norm_upper']
    source=member['preparation_plan']['source'];R=max(abs(F(source['spectral_lower'])),abs(F(source['spectral_upper'])))
    C=2*R*K;G=R+K
    floor=regional.common_contract_payload(member,1)['physical_density_floor_factored'];N=floor['negative_exponent']
    slope=modular.selected_modular_slope(member);beta=slope['beta'];a=F(member['source']['spacing'])
    witness=current_witness();norm=F(witness['corner_projected_commutator_norm_squared_coefficient'])/a**4
    qpref=slope['modular_slope_lower']*norm
    if qpref.numerator.bit_length()>N or min(R,K,C,G)<1:
        raise ArithmeticError('the rational pulse proof requires q_lower<1 and norm caps at least one')
    taupref=qpref/(16*G**2*C)
    contrast=qpref*taupref/(4*C)
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'same_source_state_reference_clock':True,'source_covariant_corner_energy':energy,
        'observable_A':'cube-covariant corner energy h_prime_0',
        'observable_B':'i[H,A], the actual net incoming energy current',
        'energy_norm_upper':K,'source_H_norm_upper':R,'current_norm_upper':C,
        'controlled_H_plus_A_norm_upper':G,'physical_projected_current_witness':witness,
        'projected_current_norm_squared_lower':norm,'selected_modular_slope':slope,
        'source_density_floor':floor,'original_double_commutator_lower':factored(qpref,N),
        'pulse_duration_source_clock':factored(taupref,N),
        'actual_pulse':'one Cayley step of H+A, then complete binary B readout',
        'Cayley_unitary':'(I+i tau(H+A)/2)^(-1)(I-i tau(H+A)/2)',
        'binary_effects':'(I +/- B/C)/2','numeric_record_labels':[-C,C],
        'no_pulse_original_plus_probability':F(1,2),
        'negative_original_current_mean_magnitude_lower':factored(qpref*taupref/2,2*N),
        'complete_binary_total_variation_lower':factored(contrast,2*N),
        'cyclic_source_work_lower':factored(2*contrast**2/beta,4*N),
        'relative_source_derivative':'d_lambda[W_lambda(T)^* O W_lambda(T)]_0 = i[integral g(t) alpha_t(A) dt,O]',
        'relative_source_unitary':'W_lambda(T)=exp(i H T) U_lambda(T)',
        'finite_Ward_response_commutes_with_spatial_coarsening':True,
        'controlled_energy_balance':'d< f_x h_x >/dt + sum_y f_x f_y <J_xy> = dot(f_x)<h_x>',
        'applied_control_work_term_retained':True,
        'pulse_control_is_a_proved_proper_metric':False,
        'complete_lapse_control_and_current_records_constructed':True,
        'large_source_histogram_evaluated':False,
        'uniform_physical_refinement_response_proved':False,
        'Lorentz_stress_or_Einstein_dynamics_proved':False,
        'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_energy_response.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('source energy response, pulse, probability or domain changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    value={'schema':'oph.primary_energy_response_certificate.v1',
           'members':[bind_member(contract.member(n)) for n in (0,1,2,8)],'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('energy response source changed during replay')
    return value
