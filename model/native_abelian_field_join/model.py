"""Exact reference, charge and interval controls for the native finite join."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import json
import sympy as s
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_source.model import PAULI
from model.coherent_defect_wave.chiral_joint.oriented_links.rotations import cubic_rotations
ROOT=Path(__file__).resolve().parents[2]


def local_reference_controls():
    currents=[s.kronecker_product(s.eye(2),g) for g in representation()[2]]
    spin=[s.kronecker_product(p/2,s.eye(15)) for p in PAULI]
    raw=s.diag(*(1 if j%15 in (10,11) else -1 if j%15==9 else 0 for j in range(30)))
    J=s.zeros(63);J[0,0]=1;J[61,62]=J[62,61]=1
    for j in range(30):J[31+j,1+j]=1;J[1+j,31+j]=-1
    parity=s.diag(1,*([-1]*60),1,1)
    assert J*J.conjugate()==parity and J.H*J==s.eye(63)
    energy=s.diag(0,*([4]*60),0,0)
    assert J*energy.conjugate()*J.H==energy
    def extend(matrix,charge=0):return s.diag(s.zeros(1),matrix,-matrix.conjugate(),s.Matrix([[charge]]),s.Matrix([[-charge]]))
    for j,g in enumerate(currents):
        G=extend(g,int(j==11))
        assert J*G.conjugate()*J.H==-G
        assert g.extract((9,24),(9,24))==int(j==11)*s.eye(2)
    for j,operator in enumerate(spin):
        S=extend(operator)
        assert J*S.conjugate()*J.H==-S
        assert operator.extract((9,24),(9,24))==PAULI[j]/2
    L=extend(raw);assert J*L.conjugate()*J.H==-L
    assert raw.extract((9,24),(9,24))==-s.eye(2)
    Q=s.diag(*([0]*61),-1,1);assert J*Q*J.H==-Q
    for rotation,U in cubic_rotations():
        old=s.kronecker_product(U,s.eye(15))
        whole=s.diag(s.eye(1),old,old.conjugate(),s.eye(2))
        assert (J*whole.conjugate()*J.H-whole).applyfunc(s.simplify)==s.zeros(63)
    rows=[]
    for n in range(1,31):
        # The current/number representation preserves every Fock sector.
        # Exchange phases are constant there and square to its parity.
        K=s.Matrix([[0,(-1)**n],[1,0]])
        assert K*K.conjugate()==(-1)**n*s.eye(2)
        assert K*(4*n*s.eye(2))*K.H==4*n*s.eye(2)
        rows.append({'N':n,'old_and_dual_dimension_each':comb(30,n),'energy':4*n,'parity':(-1)**n})
    dimension=3+2*sum(comb(30,n) for n in range(1,31))
    assert dimension==2**31+1
    return {'complete_reference_dimension':dimension,'old_full_Fock_dimension_retained':2**30,
      'vacuum_shared_and_unique':True,'port_charges':[1,-1],'port_energies':[0,0],
      'actual_ec_one_reference_gauge_factorization_checks':12,
      'actual_ec_one_reference_spin_factorization_checks':3,
      'complete_Fock_sector_exchange_controls':rows,
      'proper_cube_covariance_controls':24,
      'full_reference_energy_invariance':'number-preserving old/dual exchange; not particle-hole'}


def global_reference_controls():
    n=26;I=s.eye(n);Z=s.zeros(n)
    J=s.BlockMatrix([[Z,-I],[I,Z]]).as_explicit()
    Q=s.diag(I,-I)
    assert J*J.conjugate()==-s.eye(52)
    assert J*Q.conjugate()*J.H==-Q
    assert J*(4*s.eye(52))*J.H==4*s.eye(52)
    for p in PAULI:
        plus=s.kronecker_product(s.eye(13),p/2)
        spin=s.diag(plus,-plus.conjugate())
        assert J*spin.conjugate()*J.H==-spin
    return {'full_compensator_modes':52,'full_compensator_Fock_dimension':2**52,
      'opposite_charge_spinor_flavors_each':13,'Hamiltonian':'4N_ref',
      'charge_exchange_preserves_energy_on_all_number_sectors':True,
      'fixed_compensator_vertex':1}


def scalar_interval_controls():
    native=json.loads((ROOT/'model/coherent_defect_wave/chiral_joint/interacting_charged_response/results.json').read_text())
    row=native['same_H_interacting_channels']['ec']
    L,U=map(F,row['full_interacting_energy_readout']['source_excitation_energy_bin_after_subtracting_declared_reference_cost'])
    a=F(row['reference_sector_comparison']['zero_reference_cost_full_lattice_one_particle_sector_gap_lower'])
    beta=F(row['reference_sector_comparison']['positive_cost_full_Fock_occupied_reference_sector_gap_lower'])
    affine=json.loads((ROOT/'model/affine_gauss_spectrum/results.json').read_text())
    low,high=map(F,affine['spectral']['low_charged_levels'][0]['actual_H_energy_interval'])
    upper=U+2+high; triple=4+3*(L-F(1,2));inactive=6+a
    assert L>F(1,2) and a>0 and beta>0
    assert triple>upper and inactive>upper
    return {'native_selected_E_interval':list(map(str,(L,U))),
      'complete_port_source_gap_lower':str(a),'complete_old_and_dual_reference_gap_lower':str(beta),
      'unpenalized_actual_neutral_ground_active_charge_margin':str(L-F(1,2)),
      'exact_selected_one_charge_intertwiner_shift':'E_*+2',
      'selected_ground_total_energy_outer_interval':list(map(str,(L+2+low,upper))),
      'selected_source_energy_after_reference4_subtraction':list(map(str,(L-2+low,upper-4))),
      'unpenalized_full_charged_ground_identified':False,
      'alternative_penalty_explicitly_added':2,
      'alternative_three_active_charge_margin':str(triple-upper),
      'alternative_inactive_port_margin':str(inactive-upper),
      'continuous_physical_spin_or_Dirac_dispersion_proved':False,
      'full_native_selected_eigenvector_numerically_executed':False}


def whole_bin_bounds():
    native=json.loads((ROOT/'model/coherent_defect_wave/chiral_joint/interacting_charged_response/results.json').read_text())
    row=native['same_H_interacting_channels']['ec']['full_interacting_energy_readout']
    L,U=map(F,row['source_excitation_energy_bin_after_subtracting_declared_reference_cost'])
    wlo,whi=map(F,row['conditional_probabilities']['line'])
    assert 0<wlo<whi<F(1,2)
    epsilon=(U-L)/2
    prep=json.loads((ROOT/'model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json').read_text())['preparation']['finite_energy_agreement_filter']
    def down(value,digits):
        value=F(value);scale=10**digits
        return F((value*scale).numerator//(value*scale).denominator,scale)
    q=down(prep['success_probability_lower'],18)
    f=down(prep['conditional_ground_pair_fidelity_lower'],24)
    Q=q**13;p=f**13
    assert 0<Q<=1 and 0<p<=1
    mobile=json.loads((ROOT/'model/affine_gauss_mobile/results.json').read_text())
    universal=mobile['uniform_exact_bounds']['neighbor_arrival']
    time=F(universal['time']);amplitude=F(universal['actual_arrival_vector_norm_lower'])
    assert amplitude**2/4==F(universal['unconditional_probability_lower'])
    reduced=amplitude-epsilon*time;assert reduced>0
    value=Q*p*wlo*reduced**2
    assert value>F('0.0000098')
    records=[]
    for case in mobile['finite_Peierls_first_integral_controls']:
        time=F(case['time']);remainder=F(case['all_orders_remainder_interval'][1])
        amplitudes=[F(case[key][0])-remainder for key in ('first_arrival_norm_interval','first_arrival_nonvacuum_norm_interval')]
        lowered=[a-epsilon*time for a in amplitudes]
        assert all(a>0 for a in lowered)
        bounds=[Q*p*wlo*a*a for a in lowered]
        records.append({'source_vertex':case['source_vertex'],'arrival_vertex':case['arrival_vertex'],'time':str(time),
          'native_unconditional_arrival_lower':str(bounds[0]),
          'native_unconditional_arrival_and_nonvacuum_lower':str(bounds[1]),
          'native_unconditional_arrival_lower_decimal':format(float(bounds[0]),'.17g'),
          'native_unconditional_arrival_and_nonvacuum_lower_decimal':format(float(bounds[1]),'.17g')})
    return {'actual_native_bin_half_width':str(epsilon),'bin_weight_unnormalized_interval':list(map(str,(wlo,whi))),
      'literal_native_CAR_pump_success_probability':'1/2','conditional_native_bin_probability':'2w_J',
      'native_physical_preparation_filter_success_lower':str(Q),'additional_P0_success_lower':str(p),
      'P0_is_additional_finite_spectral_service_with_pointer_and_control_cost':True,
      'source_copies_retained':26,'P0_projects_all26_copies_including_partners':True,
      'native_spin_reference_modes_retained':26*16,
      'retained_filter_counter_dimension_each':prep['explicit_integer_M'],
      'complete_kinematic_storage_formula':'D^26 (2^31+1)^13 2^52 3^12 5456 2^416 M^13 2^17',
      'D':'2^240 5^8 119^24','retained_partner_Hamiltonian':'sum of13 unchanged native h; energy13E_g in successful baseline',
      'uniform_neighbor_arrival_time':'1/4','uniform_unconditional_native_arrival_lower':str(value),
      'uniform_unconditional_native_arrival_lower_decimal':format(float(value),'.17g'),
      'complete_six_record_law':['1-Q','Q(1-p)','Qp/2','Qp(1/2-w_J)','Qp w_J alpha','Qp w_J(1-alpha)'],
      'sharper_actual_field_amplitude_transfers':records,
      'extra_factor_one_half_or_one_quarter_inserted':False,
      'actual_native_H_replaced_by_center_energy':False,
      'comparison_error_scope':'one full bin charge, all bin states and all5456 photon states'}


def gauss_boundary_guards():
    from itertools import product
    # Three-site complete test carrier: g, four selected odd components,
    # and both parities of inactive opposite-charge states. Every ternary
    # longitudinal value is tested, including nonphysical inputs.
    charges=(0,-1,-1,1,1,-1,1,-1,1)
    parities=(0,1,1,1,1,1,1,0,0)
    def gauss(state,ell,reference):
        return (-sum(ell)-charges[state[0]]-reference,ell[0]-charges[state[1]],ell[1]-charges[state[2]])
    count=0;zeros=0;pumps=0;pumpzeros=0
    for state in product(range(9),repeat=3):
        for ell in product((-1,0,1),repeat=2):
            for reference in range(-3,4):
                before=gauss(state,ell,reference)
                for u,v in ((0,1),(1,0),(1,2),(2,1),(0,2),(2,0)):
                    if state[u] not in (1,2,3,4) or state[v]!=0:continue
                    q=charges[state[u]];target=list(state);target[v]=state[u];target[u]=0
                    shifted=list(ell)
                    if u:shifted[u-1]-=q
                    if v:shifted[v-1]+=q
                    if any(x not in (-1,0,1) for x in shifted):zeros+=1;continue
                    assert gauss(target,shifted,reference)==before
                    assert sum(parities[x] for x in target)%2==sum(parities[x] for x in state)%2
                    count+=1
                for y in range(3):
                    if state[y]!=0:continue
                    for target_label in (1,2,3,4,5,6):
                        q=charges[target_label];target=list(state);target[y]=target_label
                        shifted=list(ell)
                        if y:shifted[y-1]+=q
                        ref=reference-q
                        if any(x not in (-1,0,1) for x in shifted) or ref not in range(-3,4):pumpzeros+=1;continue
                        assert gauss(target,shifted,ref)==before
                        assert (sum(parities[x] for x in target)+1)%2==sum(parities[x] for x in state)%2
                        pumps+=1
    assert count and zeros and pumps and pumpzeros
    return {'full_kinematic_nonzero_hop_Gauss_checks':count,'retained_boundary_zero_hops':zeros,
      'full_kinematic_charge_pump_Gauss_checks':pumps,'retained_boundary_zero_pumps':pumpzeros,
      'fixed_reference_charge_unchanged_under_hopping':True,
      'source_odd_plus_reference_odd_pump_is_even':True,
      'physical_input_assumed_by_diagonal_charge_checks':False}


def calculate():
    return {'schema':'oph.native_abelian_join_structure.v1',
      'local_reference':local_reference_controls(),'global_reference':global_reference_controls(),
      'intervals':scalar_interval_controls(),'whole_bin_coupling':whole_bin_bounds(),
      'full_kinematic_charge_guards':gauss_boundary_guards(),
      'scope':'selected reducing charged sector of an explicit replicated-native plus finite-abelian-field extension'}
