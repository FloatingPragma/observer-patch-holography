"""Reproducible non-Gaussian source field and full detector-law comparisons."""
from __future__ import annotations
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from .lattice import WeylLattice, TwoFermions, node_table
from . import source, bounds
from ..common_source.chiral_candidate import gaussian_source as gaussian
from ..local_repair_family import local_bias

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
INTERVAL = F(1, 16)
SPEED = F(1, 8)
INTERACTION = F(1)


def prepare(lattice):
    p, q = lattice.origin, lattice.neighbor
    return source.hop(lattice, lattice.basis_state(p, p+1), p+1, q+1, -1j, math.pi/4)


def protocol(lattice, *, steps=16, exact=False, feedback=True):
    state = prepare(lattice)
    evolve = (lambda x: lattice.exact_evolve(x, float(INTERVAL))) if exact else (
        lambda x: source.evolve(lattice, x, float(INTERVAL), steps))
    state = evolve(state)
    first = source.read_ramsey(lattice, state, lattice.origin+1, lattice.neighbor+1)
    outputs, probabilities = [], []
    for bit, branch in enumerate(first):
        if bit and feedback:
            branch = source.phase_mode(lattice, branch, lattice.neighbor+1, math.pi/2)
        branch = evolve(branch)
        for continued in source.read_ramsey(lattice, branch, lattice.origin+1, lattice.neighbor+1):
            outputs.append(continued)
            probabilities.append(float(np.vdot(continued, continued).real))
    return np.array(probabilities), outputs


def pure_cq_distance(first, second):
    """Exact rank-two formula per unnormalized cq pure output block."""
    total = 0.
    for a, b in zip(first, second, strict=True):
        na, nb = float(np.vdot(a, a).real), float(np.vdot(b, b).real)
        if na == 0 or nb == 0:
            total += (na+nb)/2
            continue
        # Stable Gram determinant. Subtracting nearly equal norm products
        # before taking a square root would fabricate an O(sqrt(eps)) error.
        perpendicular = b-(np.vdot(a,b)/na)*a
        determinant = na*float(np.vdot(perpendicular,perpendicular).real)
        total += math.sqrt((na-nb)**2+4*determinant)/2
    return total


def number_tomography(space, state):
    """Reconstruct the one-body density from complete source number records.

    Each setting prepares a fresh copy of the declared output state. It is
    not a claim that incompatible settings can be read from one undisturbed
    specimen. Every binary complement remains in the receipt.
    """
    settings = []
    diagonal = []
    for p in range(space.modes):
        branches = source.read_number(space, state, p)
        probs = [float(np.vdot(v,v).real) for v in branches]
        diagonal.append(probs[1])
        settings.append({'kind':'number','mode':p,'probabilities':probs})
    gamma = np.diag(diagonal).astype(complex)
    for p in range(space.modes):
        for q in range(p+1, space.modes):
            values=[]
            for name, coefficient in (('X',-1j),('Y',-1.)):
                reverse = gaussian.local_gate('hopping',coefficient,-math.pi/4)
                changed = space.apply_pair(state,p,q,reverse.fock_unitary)
                branches = source.read_number(space,changed,p)
                probs = [float(np.vdot(v,v).real) for v in branches]
                values.append(probs[1]-(diagonal[p]+diagonal[q])/2)
                settings.append({'kind':name,'modes':[p,q],'probabilities':probs})
            gamma[p,q] = values[0]-1j*values[1]
            gamma[q,p] = gamma[p,q].conjugate()
    return gamma, settings


def slater_witness():
    space = TwoFermions(4)
    state = space.basis_state(0, 1)
    for p, q in ((0, 2), (1, 3)):
        state = source.hop(space, state, p, q, -1j, math.pi/4)
    initial = state.copy()
    t = .25
    gate, info = source.contact(t)
    for pair in ((0, 1), (2, 3)):
        state = space.apply_pair(state, *pair, gate)
    gamma, settings = number_tomography(space,state)
    direct_gamma = space.one_body(state)
    tomography_residual = float(np.linalg.norm(gamma-direct_gamma,2))
    if tomography_residual > 2e-10:
        raise ArithmeticError('source number records failed independent one-body tomography')
    witness = float(2-np.trace(gamma@gamma).real)
    initial_gamma = space.one_body(initial)
    free = float(2-np.trace(initial_gamma@initial_gamma).real)
    # Pair occupancy is directly readable by complete number pointers. Wick
    # predicts gamma00 gamma11-|gamma01|^2; number-changing contractions vanish.
    actual = float(abs(state[space.index(0, 1)])**2)
    wick = float((gamma[0,0]*gamma[1,1]-abs(gamma[0,1])**2).real)
    # For this particular on-site pair both remain 1/4. A cross-site SAME-SPIN
    # pair is always empty, whereas its Wick prediction is sin^2(t)/4.
    pair_records=[]
    for first,branch in enumerate(source.read_number(space,state,0)):
        for second,out in enumerate(source.read_number(space,branch,2)):
            pair_records.append({'first':first,'second':second,'probability':float(np.vdot(out,out).real)})
    same_actual = pair_records[3]['probability']
    same_wick = float((gamma[0,0]*gamma[2,2]-abs(gamma[0,2])**2).real)
    lower = (F(1, 4)-F(1, 4)**3/6)**2
    copies=len(settings)+1
    counts={key:0 for key in source.COUNT_KEYS}
    swap=gaussian.local_gate('fswap')
    prep=gaussian.local_gate('hopping',-1j,math.pi/4)
    _,pointer=source.occupation()
    residual_sum=copies*(2*source.local_error(prep)+4*source.local_error(swap)+2*info['source_target_residual'])
    for key in counts:
        counts[key]=copies*(2*prep.resources[key]+4*swap.resources[key]+2*info['resources'][key])
        counts[key]+=18*pointer['resources'][key]
    residual_sum += 36*pointer['source_target_residual']
    for p in range(4):
        for q in range(p+1,4):
            for coefficient in (-1j,-1.):
                gate=gaussian.local_gate('hopping',coefficient,-math.pi/4)
                routes=2*(q-p-1)
                residual_sum+=source.local_error(gate)+routes*source.local_error(swap)
                for key in counts:counts[key]+=gate.resources[key]+routes*swap.resources[key]
    return {'model_time': '1/4', 'full_Fock_contact_columns_executed': 4,
            'detector_witness_resources':{**counts,'prepared_state_copies':copies,
                'fresh_pointer_blanks':18,'committed_detector_bits':18,
                'floating_local_residual_sum_with_routing':residual_sum,
                'floating_residual_estimate_outward_certified':False},
            'two_particle_dimension': 6, 'initial_slater_defect': free,
            'source_final_slater_defect': witness, 'ideal_slater_defect': math.sin(t)**2,
            'ideal_slater_defect_lower_exact': str(lower),
            'source_one_body_tomography_settings': settings,
            'source_pair_number_records': pair_records,
            'source_one_body_tomography_residual': tomography_residual,
            'independent_prepared_copies_for_settings': len(settings)+1,
            'same_spin_pair_probability': same_actual,
            'same_spin_Wick_prediction': same_wick,
            'connected_four_point': same_actual-same_wick,
            'ideal_connected_four_point': -math.sin(t)**2/4,
            'source_norm': float(np.vdot(state, state).real),
            'source_contact_target_residual': info['source_target_residual'],
            'Gaussian_covariance_evolution_used': False,
            'phase_intervention_and_geometry_supplied': True,
            'all_detector_outcomes_retained': True,
            'certificate_scope': 'ideal completed source words; floating replay is a separate diagnostic'}


def clock_obstruction():
    """Actual PR1151 reference has only three source modular energies.

    This spectral multiplicity is an exact finite obstruction to identifying
    its full H66 modular flow with arbitrary field flows. Tensoring an
    independent field factor is permitted, but needs its own joint reference.
    """
    p = local_bias.Parameters(port_weights=(25,25)+(8,)*10)
    reference = local_bias.reference(p)
    multiplicities = {str(x): reference.count(x) for x in sorted(set(reference))}
    retained=json.loads((ROOT/'model/local_repair_family/record_clock_results.json').read_text())
    inferred=retained['admission']['inference']['relative_weights']
    space=TwoFermions(4);initial=space.basis_state(0,1)
    for a,b in ((0,2),(1,3)):
        initial=source.hop(space,initial,a,b,-1j,math.pi/4)
    controls=[]
    for label,weights in (('stationary_target',[F(x) for x in (25,25,8,8)]),
                          ('actual_record_inferred',[F(x) for x in inferred[:4]])):
        vector=initial.copy()
        for mode,weight in enumerate(weights):
            vector=source.phase_mode(space,vector,mode,math.log(float(weight))/4)
        gamma=space.one_body(vector)
        controls.append({'reference':label, 'four_embedded_port_weights':list(map(str,weights)),
                         'source_executed_modular_Slater_defect':float(2-np.trace(gamma@gamma).real)})
    return {'source': 'model/local_repair_family/record_clock_campaign.json',
            'reference_probabilities_and_multiplicities': multiplicities,
            'stationary_target_modular_energy_levels': 3, 'dimension': 66,
            'record_inferred_product_reference_can_have_more_than_three_levels': True,
            'occupation_preserving_CAR_dictionary': {
                'generator_identity':'K = -sum_i log(theta_i) n_i + log(Z) I on N=2',
                'all_positive_port_weights_covered': True,
                'clock_is_quadratic_under_this_dictionary': True,
                'Slater_inputs_remain_Slater_under_clock_unitary': True,
                'native_repair_interleavings_covered': False,
                'controls':controls,
                'contact_interaction_generated_by_this_clock':False,
                'nonlinear_CAR_embeddings_excluded':False},
            'field_N2_dimensions_executed': [1431,31125,234955],
            'unrestricted_unitary_algebra_identification_possible': False,
            'reason': 'dimension and modular spectral multiplicities differ; no unitary algebra isomorphism',
            'embedding_or_compression_excluded': False,
            'field_action_time_is_recovered_native_modular_time': False,
            'remaining': 'joint state/algebra embedding, intervention covariance and physical normalization'}


def apparatus_resources(lattice):
    p, q = lattice.origin+1, lattice.neighbor+1
    swap = gaussian.local_gate('fswap')
    calls = [(gaussian.local_gate('hopping', -1j, math.pi/4), 1),
             (gaussian.local_gate('hopping', -1., math.pi/4), 2),
             (gaussian.local_gate('hopping', -1., -math.pi/4), 2)]
    counts = {key: 0 for key in source.COUNT_KEYS}
    residual_sum=0.
    for gate, multiplicity in calls:
        residual_sum += multiplicity*(source.local_error(gate)+2*(q-p-1)*source.local_error(swap))
        for key in counts:
            counts[key] += multiplicity*(gate.resources[key]+2*(q-p-1)*swap.resources[key])
    _, pointer = source.occupation()
    feedback = gaussian.local_gate('onsite', -math.pi/2, 1.)
    for key in counts:
        counts[key] += 2*pointer['resources'][key]+feedback.resources[key]
    residual_sum += 4*pointer['source_target_residual']+source.local_error(feedback)
    return {**counts, 'floating_local_residual_sum_with_routing':residual_sum,
            'floating_residual_telescoping_estimate':math.expm1(residual_sum),
            'floating_residual_estimate_outward_certified':False,
            'fresh_neutral_pointer_blanks': 2, 'committed_detector_bits': 2,
            'initial_filled_modes_supplied': 2, 'field_qubits': lattice.modes,
            'native_biased_repair_attempts_claimed': 0, 'record_reset_between_reads': False,
            'initialization_storage_readout_costs_derived': False,
            'resource_scope': 'one preparation, two full reads, worst-branch feedback'}


def level(radius, steps=16):
    lattice = WeylLattice(radius)
    probabilities, outputs = protocol(lattice, steps=steps)
    reference, reference_outputs = protocol(lattice, steps=steps, exact=True)
    nofeedback, _ = protocol(lattice, steps=steps, feedback=False)
    free = WeylLattice(radius, interaction=0.)
    free_probs, _ = protocol(free, steps=steps)
    witness_bound = bounds.interaction_witness(steps=steps)
    t = float(F(witness_bound['model_time']))
    def single_read(obj):
        evolved = source.evolve(obj, prepare(obj), t, steps)
        yes = source.read_ramsey(obj, evolved, obj.origin+1, obj.neighbor+1)[1]
        return float(np.vdot(yes, yes).real)
    contrast = single_read(lattice)-single_read(free)
    temporal = bounds.history(((INTERVAL,steps),(INTERVAL,steps)), SPEED, INTERACTION)
    spatial = bounds.boundary(radius, 2*INTERVAL, SPEED, INTERACTION)
    joint = min(F(1), F(temporal['complete_history_half_trace_upper'])+F(spatial['complete_history_half_trace_upper']))
    return {'radius': radius, 'sites': lattice.volume, 'CAR_modes': lattice.modes,
            'two_particle_dimension': lattice.dimension, 'split_steps_per_interval': steps,
            'complete_two_read_probabilities': probabilities.tolist(),
            'independent_sparse_probabilities': reference.tolist(),
            'free_complete_two_read_probabilities': free_probs.tolist(),
            'history_normalization_residual': abs(float(sum(probabilities))-1),
            'source_reference_history_TV': float(np.sum(abs(probabilities-reference))/2),
            'source_reference_complete_quantum_cq_distance': pure_cq_distance(outputs,reference_outputs),
            'interaction_deletion_history_TV': float(np.sum(abs(probabilities-free_probs))/2),
            'feedback_deletion_history_TV': float(np.sum(abs(probabilities-nofeedback))/2),
            'source_short_interaction_probability_contrast': contrast,
            'short_interaction_certificate': witness_bound,
            'temporal_certificate': temporal, 'spatial_certificate': spatial,
            'ideal_source_to_infinite_lattice_complete_history_bound': str(joint),
            'evolution_resources': source.resources(lattice,float(INTERVAL),steps),
            'apparatus_resources': apparatus_resources(lattice)}


def source_bindings():
    paths = sorted([p.relative_to(ROOT).as_posix() for p in HERE.glob('*.py')]+[
        'model/interacting_field/README.md', 'model/interacting_field/DERIVATION.md',
        'model/interacting_field/campaign.json',
        'model/common_source/chiral_candidate/gaussian_source.py',
        'model/common_source/chiral_candidate/kinetic.py',
        'model/observer/source_scalar/compiler.py',
        'model/observer/scalar_response/experiment.py',
        'model/local_repair_family/local_bias.py',
        'model/local_repair_family/record_clock_campaign.json',
        'model/local_repair_family/record_clock_results.json'])
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def run(radii=(1,2,3), steps=16):
    rows = [level(r,steps) for r in radii]
    slater = slater_witness()
    okay = abs(slater['source_final_slater_defect']-slater['ideal_slater_defect']) < 2e-10
    okay &= abs(slater['connected_four_point']-slater['ideal_connected_four_point']) < 2e-10
    for row in rows:
        okay &= row['history_normalization_residual'] < 2e-8
        okay &= row['source_reference_history_TV'] < 2e-7
        okay &= row['source_short_interaction_probability_contrast'] > float(F(
            row['short_interaction_certificate']['ideal_executed_interaction_minus_free_probability_lower']))
    return {'schema': 'oph.interacting_weyl_field.v1', 'all_pass': bool(okay),
            'scope': 'source-compiled interacting 3D lattice CAR field; finite N=2 histories; no chiral continuum acceptance',
            'base_advance': {'PR':1151, 'head':'67aa7b6f4a904a1a90ea6d77bf2f66c741f2ba0c',
                             'treated_as_available_for_development': True},
            'levels': rows, 'non_Gaussian_witness': slater,
            'naive_Weyl_nodes':node_table(), 'net_node_chirality':sum(x['jacobian_chirality'] for x in node_table()),
            'free_low_momentum_certificate': bounds.low_momentum(F(1,4),2*INTERVAL),
            'native_clock_compatibility': clock_obstruction(),
            'source_bindings': source_bindings()}
