"""Same compact-texture current/pairing join and exact branch obstructions."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import sympy as s

from model.compact_defect_pair import topology, spectral
from model.standard_model_source_map import koide_source
from model.common_source.chiral_candidate import action

ROOT=Path(spectral.__file__).resolve().parents[2]


def clean(a):
    return s.ImmutableMatrix(s.Matrix(a).applyfunc(lambda x:s.factor(s.radsimp(s.expand(x)))))


def rational(value):
    if type(value)not in (str,int,F):raise ValueError('exact rational required')
    return F(value.numerator,value.denominator)if type(value)is F else F(value)


def parity_schur(mass=64):
    m=rational(mass)
    if m<8:raise ValueError('the existing four-band theorem requires M>=8')
    inherited=topology.spectral_bounds(m,1,3)
    b2=F(inherited['core_to_exterior_norm_squared']);r=F(inherited['band_abs_energy_upper'])
    q=2/(m-r)
    if not 0<q<1:raise ArithmeticError('source Neumann expansion is not admitted')
    error=b2*r/(m*(m-r))+b2/(m-r)*q*q/(1-q*q)
    alpha=1/(4*m)
    if not error<alpha:raise ValueError('this mass does not certify sign-separated clusters')
    return {'mass':str(m),'old_abs_energy_upper':str(r),'coupling_norm_squared':str(b2),
            'leading_positive_energy':str(alpha),'parity_resolvent_error':str(error),
            'negative_cluster':[str(-alpha-error),str(-alpha+error)],
            'positive_cluster':[str(alpha-error),str(alpha+error)],
            'cluster_count_each_sign':2,'count_inherited_from_complete_four_band_theorem':True,
            'odd_hopping_terms_vanish':True,'source_hopping_or_texture_changed':False,
            'core_probability_per_eigenvector_lower':inherited['core_probability_lower']}


def phase_zero_branches():
    bound=parity_schur();error=F(bound['parity_resolvent_error']);alpha=F(bound['leading_positive_energy'])
    lower,upper=F(707106,500000),F(707107,500000)
    if not lower*lower<2<upper*upper or not error<F(1,25000):
        raise ArithmeticError('rational source margins failed')
    light=((15-10*upper)/216,(15-10*lower)/216)
    heavy=((15+10*lower)/108,(15+10*upper)/108)
    if alpha-error<=F(1,40000) or heavy[1]>=1:raise ArithmeticError('spectator/exterior gap failed')
    rows=[]
    for name,interval,families in (('repeated_lower',light,2),('upper',heavy,1)):
        for sign in (-1,1):
            lo=interval[0]+sign*alpha-error;hi=interval[1]+sign*alpha+error
            if lo<=F(1,40000):raise ArithmeticError('the composed static gap was not proved')
            rows.append({'family_coefficient_group':name,'spatial_energy_sign':sign,
                         'positive_frequency_interval':[str(lo),str(hi)],
                         'spatial_rows':2,'family_multiplicity':families,
                         'conjugate_internal_charges':['-1','1'],
                         'positive_Nambu_rows':4*families})
    if any(F(a['positive_frequency_interval'][1])>=F(b['positive_frequency_interval'][0])
           for a,b in zip(rows,rows[1:])):raise ArithmeticError('the four observable groups overlap')
    return {'program_phase':'0','Higgs_norm':'1','texture_mass':'64','half_separation':1,
            'rows':rows,'complete_positive_charged_Nambu_rows':24,
            'negative_charged_Nambu_rows':24,'static_gap_floor':'1/40000',
            'localized_positive_spectator_Nambu_rows':156,
            'full_static_Nambu_gap_floor':'1/40000','exterior_abs_energy_floor':'61',
            'distinct_positive_frequency_groups_at_least':4,
            'within_group_distinct_splittings_proved':False,
            'three_observed_particle_species_assigned':False}


def localization_readout():
    # Each of four orthonormal band columns has the inherited core leakage.
    # Hilbert--Schmidt domination and equal rank control every core probe.
    bound=parity_schur();per=F(bound['core_probability_per_eigenvector_lower'])
    outside=min(F(1),4*(1-per))
    return {'normalized_core_probe_original_band_weight_lower':str(1-outside),
            'all_four_band_rows_retained':True,
            'complete_ground_energy_window':['0','1'],
            'native_spectral_readout_executed':False,'other_modes_complement_retained':True,
            'probability_is_conditioned_on_band':False,
            'input_domain':'normalized core coefficient in the full retarded CAR probe; energy interpretation uses the chosen paired quadratic ground',
            'separate_addition_or_removal_weight_bound':False,
            'old_unpaired_filled_sea_asserted_to_be_new_ground':False}


def source_receipt():
    p=ROOT/'model/compact_defect_pair/results.json';value=json.loads(p.read_text())
    required=('PROOF.md','spectral.py','topology.py','certification.py')
    pins={}
    for name in required:
        key='model/compact_defect_pair/'+name;digest=hashlib.sha256((ROOT/key).read_bytes()).hexdigest()
        if value['source_bindings'].get(key)!=digest:raise ValueError('compact source proof/input changed')
        pins[key]=digest
    for key in action.SOURCE_FILES:
        digest=hashlib.sha256((ROOT/key).read_bytes()).hexdigest()
        if value['source_bindings'].get(key)!=digest:raise ValueError('the inherited action fixture changed')
        pins[key]=digest
    if value['complete_exact_certificate']['exact_core_band_dimension']!=4 or value['infinite_volume']['nonzero_isolated_band_dimension']!=4:
        raise ValueError('owning complete cluster count changed')
    return {'receipt_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'direct_source_pins':pins,
            'retained_complete_spectrum_dimension':value['complete_exact_certificate']['dimension'],
            'cold_spectrum_recomputed':False,'floating_splittings_used_as_proof':False}


def local_pairing_identities():
    sigma=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))
    j=s.Matrix([[0,1],[-1,0]])
    for p in sigma:
        hop=s.I*p/2
        if hop*j-j*hop.conjugate()!=s.zeros(2):raise ArithmeticError('source hopping conjugation failed')
        if p*j+j*p.conjugate()!=s.zeros(2):raise ArithmeticError('source texture conjugation failed')
    # An actual core-neighbor source block: core(0,0,1), neighbor(0,0,2).
    h=s.zeros(4);h[:2,2:]=s.I*sigma[2]/2;h[2:,:2]=-s.I*sigma[2]/2;h[2:,2:]=64*sigma[2]
    c=s.diag(-j,j);j0=s.diag(j,j)
    if h*c+c*h.conjugate()!=s.zeros(4):raise ArithmeticError('staggered source PH identity failed')
    defect=clean(h*j0-j0*h.conjugate())
    if defect==s.zeros(4):raise ArithmeticError('uniform pairing incorrectly restored Clifford square')
    k=s.Matrix([[s.Rational(1,5),s.I/100],[s.I/100,s.Rational(3,100)]])
    hk=s.kronecker_product(h,s.eye(2));delta=s.kronecker_product(c,k.conjugate())
    if hk*delta+delta*hk.T!=s.zeros(8):raise ArithmeticError('the static pairing does not preserve unpaired energy')
    bdg=hk.row_join(delta).col_join(delta.H.row_join(-hk.conjugate()))
    u=s.diag(s.eye(8),s.kronecker_product(c.H,s.eye(2)))
    target=hk.row_join(s.kronecker_product(s.eye(4),k.conjugate())).col_join(
        s.kronecker_product(s.eye(4),k).row_join(hk))
    if clean(u.H*bdg*u-target)!=s.zeros(16):raise ArithmeticError('additive spectrum dictionary failed')
    return {'exact_Pauli_checks':6,'actual_source_core_neighbor_dimension':4,
            'exact_complex_symmetric_Nambu_control_dimension':16,
            'staggered_pairing_additive_dictionary':True,
            'bare_one_body_and_pairing_Fock_generators_commute':True,
            'uniform_pairing_source_texture_Clifford_defect_nonzero':True,
            'uniform_pairing_defect_Frobenius_squared':str(s.trace(defect.H*defect)),
            'same_uniform_free_action_claimed':False}


@lru_cache(maxsize=1)
def _higgs_json():
    subsets,_,currents,_,channels=action._representation()
    tensor=channels['charged_lepton']['components'][1]
    if tensor!=tensor.T or clean(tensor**3-tensor)!=s.zeros(15) or s.trace(tensor.H*tensor)!=2:
        raise ArithmeticError('the inherited unit charged-Higgs tensor changed')
    y=koide_source.exact_source(1)['action_coefficient']
    k=s.kronecker_product(tensor,y.conjugate());q=action.electric_charge()
    if clean(q*k+k*q.T)!=s.zeros(45):raise ArithmeticError('source pairing lost electric charge')
    supports=[i for i in range(45) if any(k[i,j]!=0 for j in range(45))]
    charges={str(q[i,i]):sum(q[j,j]==q[i,i]for j in supports)for i in supports}
    if charges!={'-1':3,'1':3}:raise ArithmeticError('the actual charged source support changed')
    if len(currents)!=12:raise ArithmeticError('native current count changed')
    if y!=koide_source.exact_source(1)['action_coefficient']:raise ArithmeticError('coherent coefficient changed')
    a,b=s.sqrt(koide_source.WEIGHTS[0]),s.sqrt(koide_source.WEIGHTS[1])
    light=s.expand((a-b)**2/3);heavy=s.expand((a+2*b)**2/3)
    ph=s.ones(3)/3;pl=s.eye(3)-ph
    if clean(y-light*pl-heavy*ph)!=s.zeros(3) or ph*ph!=ph or s.trace(pl)!=2:
        raise ArithmeticError('phase0 source family projector identity changed')
    expected={light:2,heavy:1}
    return json.dumps({'operator_frame':'authenticated existing compact-action basis; native P15 dictionary inherited',
                       'native_current_architecture_cold_rebuilt':False,'internal_dimension':15,'declared_family_copies':3,'native_currents':12,
                       'charged_internal_family_supports':charges,'all_other_internal_family_rows':39,
                       'localized_CAR_modes':180,'localized_Nambu_rows':360,
                       'exact_source_family_spectrum':{str(x):n for x,n in expected.items()},
                       'electric_pairing_identity_exact':True,
                       'full_gauge_invariance_without_transforming_Higgs_claimed':False},sort_keys=True)


def certificate():
    return {'source':source_receipt(),'parity_refined_band':parity_schur(),
            'exact_source_local_algebra':local_pairing_identities(),'same_current_Higgs':json.loads(_higgs_json()),
            'phase0_complete_charged_branches':phase_zero_branches(),'complete_core_readout':localization_readout(),
            'full_Fock_reference_countercontrol':ground_countercontrol(),
            'fixed_compact_texture_is_declared':True,
            'texture_selects_internal_representation':False,'coherent_Y_is_selected_by_texture':False,
            'CAR_creation_changes_texture_degree':False,'autonomous_moving_defect_dispersion_proved':False,
            'physical_particle_assignments':[],'physical_mass_predictions':[]}


def ground_countercontrol():
    # Actual invariant pairing blocks have unpaired one-body energies e,-e.
    # A rational block checks the full four-state occupation algebra independently.
    e,m=s.Rational(1,256),s.Rational(1,100)
    n0=s.diag(0,1,0,1);n1=s.diag(0,0,1,1)
    pairing=s.zeros(4);pairing[3,0]=pairing[0,3]=m
    h=e*(n0-n1)+pairing;q=-n0+n1
    if h*q-q*h!=s.zeros(4):raise ArithmeticError('pairing changed total electric charge')
    if h.eigenvals()!={-m:1,-e:1,e:1,m:1}:raise ArithmeticError('full Fock source levels changed')
    sea=s.eye(4)[:,2];ground=(s.eye(4)[:,0]-s.eye(4)[:,3])/s.sqrt(2)
    if h*sea!=-e*sea or h*ground!=-m*ground or sea.H*ground!=s.zeros(1):
        raise ArithmeticError('the source-state obstruction changed')
    annihilator=s.zeros(4);annihilator[0,1]=annihilator[2,3]=1
    energy_window=ground*ground.H+sea*sea.H
    maps=tuple(clean(p*a)for p in (energy_window,s.eye(4)-energy_window)
               for a in (annihilator,annihilator.H))
    if sum((a.H*a for a in maps),s.zeros(4))!=s.eye(4):
        raise ArithmeticError('complete source energy/addition/removal Gram failed')
    weights=[str(s.simplify((ground.H*a.H*a*ground)[0]))for a in maps]
    return {'complete_energy_records':weights,'all_input_effect_sum_identity':True,
            'energy_window_relative_to_paired_ground':['0',str(m)],
            'Fock_dimension':4,'e':str(e),'pairing_strength':str(m),
            'complete_energies':list(map(str,(-m,-e,e,m))),
            'old_filled_sea_stationary':True,'old_filled_sea_is_paired_ground':False,
            'electric_charge_conserved':True,'this_control_is_the_full_compact_spectrum':False}
