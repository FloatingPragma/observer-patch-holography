"""Native pole residues and full-H energy-filtered arrival bounds."""
from functools import lru_cache
from pathlib import Path
import hashlib
import json

import sympy as s

from model.coherent_defect_wave.chiral_joint.interacting_charged_response.transfer import data as energy_data
from model.coherent_defect_wave.chiral_joint.interacting_charged_transport.algebra import routes
from model.coherent_defect_wave.chiral_joint.interacting_charged_transport.transport import certificate as transport_certificate
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.resolvent import field_data,field_value
from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import source
from model.coherent_defect_wave.chiral_joint.atomic_spin_geometry.cube import spectrum
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import outward,sqrt_interval,mul,add
from model.coherent_defect_wave.chiral_source.model import PAULI
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def rational(value):
    value=s.sympify(value)
    if value.is_Rational is not True:raise ValueError('exact rational input required')
    return value


def encoded_interval(values,bits=100):
    return outward(tuple(map(rational,values)),bits)


def upper(value,bits=120):
    return s.Rational(outward((0,rational(value)),bits)[1])


def lower(value,bits=120):
    return s.Rational(outward((rational(value),rational(value)),bits)[0])


@lru_cache(maxsize=1)
def native_pole_residue():
    """Prove the actual lowest e^c residue is a spin-triplet tensor.

    The two original CAR transfer maps act on the complete native local
    blocks. Their pole residues are compared algebraically, not fitted.
    """
    data=energy_data('ec');local=data['local']
    h=source(1)['H'];x=s.Symbol('x');z=s.Symbol('z')
    p=s.Poly(h.charpoly(x).as_expr(),x).sqf_part()
    q=local['q'];_,_,field=field_data()
    coefficients=[];powers=[s.eye(h.rows)]
    for i in range(p.degree()):
        poly=s.Poly(sum(p.nth(j)*x**(j-i-1) for j in range(i+1,p.degree()+1)),x)
        coefficients.append(field_value(poly)/field_value(p.diff()))
        if i+1<p.degree():powers.append(h*powers[-1])
    numerators=[]
    for block in local['blocks']:
        image=block['map'];moments=[]
        for _ in range(q.degree()):
            pullback=block['adjoint']*image
            moments.append(sum((coefficient*field.convert(s.trace(power*pullback))
                               for coefficient,power in zip(coefficients,powers)),field.zero))
            image=block['H']*image
        entries={}
        for i,moment in enumerate(moments):
            for j in range(i+1,q.degree()+1):
                entries[j-i-1]=entries.get(j-i-1,field.zero)+field.convert(q.nth(j))*moment
        numerators.append(s.Poly.from_dict({(k,):v for k,v in entries.items()},z,domain=field))
    lo,hi=local['rows'][0]['energy_interval']
    factors=[factor for factor,_ in s.factor_list(q)[1] if factor.count_roots(lo,hi)==1]
    assert len(factors)==1
    selected=factors[0]
    selected_field=s.Poly(selected.as_expr().subs(x,z),z,domain=field)
    difference=numerators[1]-2*numerators[0]
    assert difference.rem(selected_field).is_zero
    assert numerators[0]+numerators[1]==local['numerator']
    # Lowest visible root belongs to the actual spin-one highest-weight
    # block, and not to the raw-L=2 spin-zero highest-weight block.
    raw_l_two=[row for row in data['sector']['blocks'] if row['raw_L']==2]
    spin_one=[row for row in raw_l_two if row['twice_source_spin']==2]
    spin_zero=[row for row in raw_l_two if row['twice_source_spin']==0]
    assert len(spin_one)==len(spin_zero)==1
    one=s.Poly(s.sympify(spin_one[0]['complete_characteristic_polynomial']),x)
    zero=s.Poly(s.sympify(spin_zero[0]['complete_characteristic_polynomial']),x)
    assert one.rem(selected).is_zero and s.gcd(zero,selected).degree()==0
    return {'selected_atomic_pole_polynomial':str(selected.as_expr()),
            'native_removed_spin_target_dimensions':[b['H'].rows for b in local['blocks']],
            'actual_native_residue_identity':'residue_down=2 residue_up at the selected pole',
            'exact_residue_difference_divisible_by_selected_polynomial':True,
            'actual_selected_target_spin':1,
            'raw_source_lepton_number':2,
            'source_ground_spin':s.Rational(1,2).__str__(),
            'pole_tensor':'P_g c_s^dagger S_a c_t P_g = w delta_st I -(2w/3) sum_i sigma_i,ts S_i',
            'pump_pole_weight_interval':encoded_interval(local['rows'][0]['weight_interval'])}


def tagged_spin_effects():
    """Exact same-pole tag at both ends of the actual elementary hop."""
    w=s.symbols('w',positive=True)
    def tensor(i,j):
        return w*int(i==j)*s.eye(2)-2*w/s.Integer(3)*sum(
            (sigma[j,i]*sigma/2 for sigma in PAULI),s.zeros(2))
    effects=[]
    for axis,sigma in enumerate(PAULI):
        effect=s.zeros(4)
        for spin in range(2):
            for a,b in sigma.todok():
                for ap,bp in sigma.todok():
                    effect+=s.conjugate(sigma[ap,bp])*sigma[a,b]*s.kronecker_product(
                        tensor(ap,spin).H*tensor(a,spin),tensor(bp,b))/32
        exchange=sum(((1 if i==axis else -1)*s.kronecker_product(v/2,v/2)
                      for i,v in enumerate(PAULI)),s.zeros(4))
        expected=w**3*(s.eye(4)/12+s.Rational(2,27)*exchange)
        assert (effect-expected).applyfunc(s.expand)==s.zeros(4)
        effects.append(effect)
    return {'all_three_axes_exact':True,
            'same_pole_first_hop_effect':'w^3 [I/12+(2/27) S_origin dot R_axis S_arrival]',
            'same_cube_trial_moment':'w^3 (1/12+x_cube/648), x_cube=4 epsilon_cube',
            'normalization':'same actual |edge coefficient|^2=1/16, pump factor1/2, PW119 e^c blank-link Gram1'}


def projector_tail_bounds(*,coupling,tau,edge_norm,weight_intervals,distances,
                          tagged_energy_separation,pole_to_other_energy_max,
                          trial_error,moment_interval):
    """Full-H tagged projector tail without a global charged-band gap.

    F is an exact H0 eigenprojector at a+1; E is the actual full-H bin.
    Each atomic component S_j B P0 has squared norm w_j on a unit input.
    The leading F T S_0 B P0 norm is compared to the native cube trial.
    """
    lam,tau,k,gap,eta=map(rational,(coupling,tau,edge_norm,tagged_energy_separation,trial_error))
    ds=list(map(rational,distances));delta=list(map(rational,pole_to_other_energy_max))
    weights=[tuple(map(rational,row)) for row in weight_intervals]
    if not (lam>0 and tau>=k>0 and gap>0 and eta>=0 and len(weights)==len(ds)==len(delta)):
        raise ValueError('complete positive source controls required')
    if not all(0<lo<=hi<=1 and distance>0 for (lo,hi),distance in zip(weights,ds)):
        raise ValueError('positive complete atomic weights and bin distances required')
    roots=[sqrt_interval(row,110) for row in weights]
    v=lam*tau;beta=upper(v/gap)
    selected=lam*tau*(beta+v/ds[0])*roots[0][1]
    outside=beta*v*sum(roots[j][1]/ds[j] for j in range(1,len(weights)))
    remainder=upper(selected+outside)
    reference_error=upper(lam*k*eta)
    moment_low,moment_high=sqrt_interval(tuple(map(rational,moment_interval)),110)
    total_error=remainder+reference_error
    low=lower(lam*moment_low-total_error);high=upper(lam*moment_high+total_error)
    # Every P0-localized input satisfies ||B P0 z||<=1/sqrt(2).
    pump_norm=sqrt_interval((s.Rational(1,2),)*2,110)[1]
    energy_residual=upper(v*(pump_norm+sum(delta[j]*roots[j][1]/ds[j] for j in range(1,len(weights)))))
    return {'v':v,'tagged_projector_norm_upper':beta,
            'projection_remainder_upper':remainder,'ground_and_preparation_amplitude_error':reference_error,
            'leading_moment_root_interval':(moment_low,moment_high),
            'initial_tagged_amplitude_interval':(max(0,low),high),
            'strict_initial_amplitude_margin':low,
            'bin_state_energy_residual_upper':energy_residual,
            'tagged_phase_adjusted_amplitude_change_per_time_upper':upper(beta*energy_residual),
            'tagged_phase_adjusted_amplitude_change_absolute_upper':upper(2*beta*pump_norm)}


@lru_cache(maxsize=1)
def source_bounds():
    data=energy_data('ec');local=data['local'];neutral=data['neutral']
    prep=data['preparation']['finite_energy_agreement_filter']
    trace=s.Rational(prep['conditional_one_source_trace_distance_upper'])
    fidelity=s.Rational(prep['conditional_ground_pair_fidelity_lower'])
    assert trace**2>=1-fidelity
    eta=data['eta']+2*trace
    weights=[row['weight_interval'] for row in local['rows']]
    w=weights[0];x=tuple(4*s.Rational(value) for value in spectrum()['ground_energy_interval'])
    factor=add((s.Rational(1,12),)*2,mul((s.Rational(1,648),)*2,x))
    moment=mul(mul(mul(w,w),w),factor)
    assert moment[0]>0
    # The tag uses the same local pole after a native hop, plus one actual
    # link-energy unit. Ground and static-reference constants are retained.
    ground=local['e_interval'];atom=local['rows'][0]['energy_interval']
    a=(atom[0]+7*ground[0]+4,atom[1]+7*ground[1]+4)
    gap=a[0]+1-data['window'][1]
    delta=[max(abs(row['energy_interval'][0]-atom[1]),
               abs(row['energy_interval'][1]-atom[0])) for row in local['rows']]
    bounds=projector_tail_bounds(coupling=neutral['coupling'],tau=neutral['tau'],edge_norm=15,
        weight_intervals=weights,distances=data['distances'],tagged_energy_separation=gap,
        pole_to_other_energy_max=delta,trial_error=eta,moment_interval=moment)
    assert bounds['strict_initial_amplitude_margin']>0
    transport=transport_certificate('ec')
    line=tuple(map(s.Rational,transport['same_preparation_energy_setting']['conditional_line_interval']))
    return {'data':data,'moment':moment,'a':a,'bounds':bounds,'eta':eta,'gap':gap,
            'preparation_success_lower':lower(s.Rational(prep['success_probability_lower'])),
            'P0_success_lower':s.Rational(transport['localization_and_probe']['atomic_localization_success_lower']),
            'line_probability_interval':line}


def record_law():
    q,a,b=s.symbols('q a b')
    tagged=s.symbols('r0:6');untagged=s.symbols('s0:6')
    n=b-sum(tagged)-sum(untagged)
    probabilities=[1-q,q*(1-a),q*a/2,q*(a/2-b)]
    probabilities.extend(q*r for r in tagged)
    probabilities.extend(q*r for r in untagged)
    probabilities.append(q*n)
    assert s.expand(sum(probabilities))==1
    route_rows=routes()
    records=[{'record':'physical_preparation_failure','probability':'1-q'},
             {'record':'atomic_localization_failure','probability':'q(1-a)'},
             {'record':'native_pump_failure','probability':'q a/2'},
             {'record':'complete_full_H_energy_bin_complement','probability':'q(a/2-b)'}]
    records.extend({'record':'arrival_with_same_atomic_pole_tag','probability':f'q r{i}',**route}
                   for i,route in enumerate(route_rows))
    records.extend({'record':'arrival_with_other_atomic_energy','probability':f'q s{i}',**route}
                   for i,route in enumerate(route_rows))
    records.append({'record':'complete_nonarrival_after_energy_bin',
                    'probability':'q(b-sum_i r_i-sum_i s_i)'})
    assert len(records)==17
    return {'order':'physical preparation; P0 localization; B_ec pump; actual full-H energy bin; same-H evolution; complete tagged arrival instrument',
            'definitions':'q=preparation success; a=P0 success conditional on q; b=joint pump/bin weight conditional on q; r_i and s_i retain their original unnormalized weights',
            'records':records,'coarse_record_count':17,
            'all_preparation_failures_and_full_source_reference_environments_retained':True,
            'complete_probability_sum':'1',
            'energy_and_arrival_projectors_assumed_to_commute':False}


def timeline():
    values=source_bounds();bounds=values['bounds'];q=values['preparation_success_lower']
    low,high=bounds['initial_tagged_amplitude_interval']
    line_low,line_high=values['line_probability_interval']
    beta=bounds['tagged_projector_norm_upper']
    rows=[]
    for t in (s.Integer(0),s.Rational(1,2**24),s.Integer(1),s.Integer(100),s.Integer(1000),s.Integer(10**6)):
        change=min(bounds['tagged_phase_adjusted_amplitude_change_absolute_upper'],
                   abs(t)*bounds['tagged_phase_adjusted_amplitude_change_per_time_upper'])
        lo=max(0,low-change)**2
        hi=min((high+change)**2,beta**2*line_high)
        difference=min(beta**2*line_high,change*(2*high+change))
        assert 0<=lo<=hi<=line_high and 6*hi<line_low
        rows.append({'source_time':str(t),
            'per_route_joint_bin_and_tag_interval_conditional_only_on_preparation':encoded_interval((lo,hi),200),
            'per_route_original_input_probability_interval':encoded_interval((q*lo,hi),200),
            'per_route_probability_conditioned_on_energy_bin_interval':encoded_interval((lo/line_high,hi/line_low),200),
            'per_route_probability_change_from_t0_absolute_upper':str(difference),
            'decimal_original_input_interval':[str(s.N(q*lo,16)),str(s.N(hi,16))],
            'decimal_probability_change_from_t0_upper':str(s.N(difference,16)),
            'six_tag_union_original_input_interval':encoded_interval((6*q*lo,6*hi),200),
            'complete_line_tag_complement_original_input_probability_lower':str(q*(line_low-6*hi)),
            'positive_time_change_or_propagation_certified':False})
    assert s.Rational(rows[0]['per_route_original_input_probability_interval'][0])>0
    return rows


def source_pins():
    seeds=[str(path.relative_to(ROOT)) for path in HERE.iterdir() if path.suffix in ('.py','.md','.tex')]
    inherited=['model/coherent_defect_wave/chiral_joint/'+name+'/results.json' for name in
               ('interacting_charged_response','interacting_charged_transport','interacting_ground_preparation',
                'interacting_charge_localization','collective_spin_response')]
    paths=set(seeds+inherited)
    for name in inherited:
        receipt=json.loads((ROOT/name).read_text())
        for path,expected in receipt.get('input_sha256',{}).items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=expected:
                raise ValueError('changed inherited scientific source: '+path)
            paths.add(path)
    paths.update(import_closure(sorted(paths)))
    return {path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def certificate():
    pins=source_pins();residue=native_pole_residue();spin=tagged_spin_effects();values=source_bounds()
    bounds=values['bounds'];data=values['data']
    result={'schema':'oph.energy_resolved_native_charge_arrival.v1',
        'source_sha256':pins,
        'source':{'Hamiltonian':'unchanged N2 PW119 H0+2^-34 T+4N_ref',
                  'reference':'same full30-mode static CAR reference at cell0, one e^c particle after pump',
                  'ground_preparation':'same physical two-copy algebraic ground preparation with all failure outputs',
                  'full_H_energy_bin':list(map(str,data['window'])),
                  'atomic_selected_total_energy_interval':encoded_interval(values['a']),
                  'atomic_tagged_arrival_total_energy_interval':encoded_interval((values['a'][0]+1,values['a'][1]+1)),
                  'reference_energy_removed_from_H':False,
                  'joint_free_particle_mass_shell_assumed':False},
        'atomic_pole':residue,'tagged_first_hop':spin,
        'detector':{'P_e':'same six actual charged-arrival effects: selected link nonblank, other links blank, all nonarrival cells in true local ground',
                    'F_e':'P_e times the actual selected local charged-energy projector at arrival',
                    'full_detector':['six F_e','six P_e-F_e','I-sum_e P_e'],
                    'exact_atomic_relation':'H0_total F_e=(a+1)F_e',
                    'exact_compressed_full_H_relation':'F_e H_total F_e=(a+1) F_e',
                    'same_pole_tag_changes_or_resolves_the_measured_state':True,
                    'post_tag_state_remains_in_initial_full_H_bin':False,
                    'no_energy_free_position_measurement_claimed':True},
        'certified_bounds':{key:(list(map(str,value)) if isinstance(value,tuple) else str(value)) for key,value in bounds.items()},
        'same_cube_trial_moment_interval':encoded_interval(values['moment']),
        'complete_records':record_law(),'time_rows':timeline(),
        'controls':{'no_energy_filter_t0_tag_probability':'0, because B_ec P0 has all links blank',
                    'zero_hopping_t0_and_all_time_tag_probability':'0',
                    'actual_energy_filter_t0_tag_probability':'strictly positive by the retained projector-tail bound',
                    'global_preparation_probability_error_subtracted':False,
                    'global_isolated_charged_band_assumed':False,
                    'positive_arrival_at_t0_is_called_propagation':False,
                    'screened_and_other_link_complements_retained':True},
        'scope':'A complete energy-bin and native charge-arrival experiment. The positive same-pole arrival is present at t0 after the energy filter, hence is virtual dressing and measurement response, not established free transport. No physical mass shell, autonomous electron or continuum identification follows.',
        'verdict':'PASS_NATIVE_ENERGY_FILTERED_ARRIVAL_WITH_POSITIVE_T0_DRESSING'}
    if source_pins()!=pins:raise ValueError('source changed during energy-resolved replay')
    return result


def validate_pins(saved):
    for path,expected in saved['source_sha256'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=expected:
            raise ValueError('changed source binding: '+path)
