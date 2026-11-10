"""Actual current-map reuse, finite scalar obstruction and copy symmetry."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from model.coherent_defect_wave.chiral_particle_map import model as fields
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
from model.native_higgs_spectral_join import model as join
from model.particle_common_contract.provenance import import_closure

JOIN = 'model/native_higgs_spectral_join/results.json'
JOIN_SHA = '5749af81e1612f0658733b51acd52c2d08a1a69d174a3707989094d66d149973'
MAP = 'model/coherent_defect_wave/chiral_particle_map/results.json'
RESEARCH = HERE/'research_reference.json'
RESEARCH_SHA = '63bc7d02fbd19066ed6afb09e79713d6ba330073838a6b874568f39172a1bb81'


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def scalar_controls(lower):
    phi, currents, _, _ = higgs_memory()
    for field in phi:
        assert field**3 == s.zeros(5) and field**2 != s.zeros(5)
        assert (field.H*field)**2 == field.H*field
        assert currents[-1]*field-field*currents[-1] == -field/2
    stack = s.Matrix.vstack(phi[0], phi[0].H, phi[1], phi[1].H)
    assert stack.rank() == 5
    # A nonzero expectation is possible despite the strict eigenbackground
    # obstruction. The collective average has a quantified large-copy control.
    eta = s.Matrix([s.sqrt(2), 1, 0, 1, 0])/2
    assert (eta.H*eta)[0] == 1
    means, variances = [], []
    for field in (*phi, *(x.H for x in phi)):
        mean = s.simplify((eta.H*field*eta)[0])
        residual = field*eta-mean*eta
        means.append(str(mean)); variances.append(str(s.simplify((residual.H*residual)[0])))
    assert means == ['sqrt(2)/2', '0', 'sqrt(2)/2', '0']
    assert variances == ['1/4', '1/2', '1/4', '1/2']
    low = F(lower); bound = low**3/(1+low+low*low)
    return {'actual_source_scalar_dimension': 5,
            'actual_phi_matrices': [[[str(value) for value in row] for row in field.tolist()] for field in phi],
            'both_phi_nilpotency_orders': [3, 3], 'both_full_operator_norms': ['1', '1'],
            'common_phi_and_adjoint_kernel_dimension': 5-stack.rank(),
            'nonzero_scalar_eigenbackground_exists': False,
            'arbitrary_entangled_reference_extension_has_same_obstruction': True,
            'normalized_vector_residual_lower': '|v|^3/(1+|v|+|v|^2)',
            'native_unit_Yukawa_shared_energy_as_background_residual_lower_exact': str(bound),
            'scalar_expectation_or_collective_limit_excluded': False,
            'explicit_expectation_control': {'normalized_scalar_state': ['sqrt(2)/2', '1/2', '0', '1/2', '0'],
                'means_phi0_phi1_and_adjoints': means, 'one_copy_variances': variances,
                'independent_n_copy_average_residual_squared': ['1/(4n)', '1/(2n)', '1/(4n)', '1/(2n)'],
                'native_interacting_ground_or_stationary_background_asserted': False}}


def family_controls():
    swap = s.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    cycle = s.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    basis = []
    for i in range(3):
        for j in range(3):
            E = s.zeros(3); E[i,j] = 1
            basis.append(s.Matrix.vstack((E*swap-swap*E).reshape(9,1), (E*cycle-cycle*E).reshape(9,1)))
    constraints = s.Matrix.hstack(*basis)
    assert constraints.rank() == 7
    d, b, x = s.symbols('d b x', real=True)
    matrix = (d-b)*s.eye(3)+b*s.ones(3)
    assert matrix*swap == swap*matrix and matrix*cycle == cycle*matrix
    assert s.expand((x*s.eye(3)-matrix).det()-(x-d-2*b)*(x-d+b)**2) == 0
    return {'three_copy_permutation_commutant_dimension': 2,
            'general_Hermitian_family_matrix': '(d-b) I3+b ones(3,3), d,b real',
            'one_particle_family_eigenvalues': ['d+2b', 'd-b', 'd-b'],
            'identical_independent_copy_case': 'd=m, b=0: three equal energies',
            'spatial_replication_is_three_observed_generations': False,
            'family_breaking_operator_selected_by_copy_count': False}


def dressed_weak_polynomial():
    """Original native CAR/scalar W, D, S, Ka and NL; no source compression."""
    car = join.poly; phi = higgs_memory()[0]
    K = [[phi[i]*phi[j].H for j in range(2)] for i in range(2)]
    scalar = K[0][0]+K[1][1]
    sigma = fields.pauli()
    Ka = [sum((sigma[a][j,i]*K[i][j] for i in range(2) for j in range(2)), s.zeros(5)) for a in range(3)]
    currents = []
    for a in range(3):
        one = s.zeros(30)
        for spin in range(2):
            for i in range(2):
                for j in range(2):one[15*spin+10+i,15*spin+10+j] = sigma[a][i,j]/2
        currents.append(car.one_body(one))
    nL = car.one_body(s.diag(*[int(j%15 in (10,11)) for j in range(30)]))
    D = car.add(*[car.multiply(car.adjoint(car.composite('charged_lepton',9,spin)),
                              car.composite('charged_lepton',9,spin)) for spin in range(2)])
    W = car.add(*[car.multiply(car.monomial(matrix=Ka[a]), currents[a]) for a in range(3)])
    assert D == car.add(car.scale(car.multiply(car.monomial(matrix=scalar),nL),s.Rational(1,2)),W)
    assert W == car.adjoint(W) and W
    for a in range(12):assert not car.commutator(car.gauge_current(a),W)
    for a in range(3):assert not car.commutator(car.spin_current(a),W)
    raw = car.one_body(s.diag(*[1 if j%15 in (10,11) else -1 if j%15==9 else 0 for j in range(30)]))
    assert not car.commutator(raw,W)
    assert all(len(word)%2 == 0 for word,_,_ in W)
    return W, D, scalar, Ka, nL


def dressed_weak_matrix():
    """Literal80-state scalar/four-L-mode matrix in scalar-major order."""
    _, _, _, Ka, _ = dressed_weak_polynomial()
    sigma = fields.pauli()
    # Full operator norm: the other26 native fermion modes are identity
    # spectators. Execute every one of16 four-L-mode occupations and5 scalars.
    aa = []
    for j in range(4):
        c=s.zeros(16)
        for mask in range(16):
            if mask>>j&1:c[mask^(1<<j),mask]=(-1)**((mask&((1<<j)-1)).bit_count())
        aa.append(c)
    jl = [sum((sigma[a][i,j]*aa[2*spin+i].H*aa[2*spin+j]/2
               for spin in range(2) for i in range(2) for j in range(2)),s.zeros(16)) for a in range(3)]
    matrix=sum((s.kronecker_product(Ka[a],jl[a]) for a in range(3)),s.zeros(80))
    return matrix


def dressed_weak_density():
    matrix = dressed_weak_matrix()
    eigenvalues=(s.Integer(-1),s.Rational(-1,2),s.Integer(0),s.Rational(3,2),s.Integer(2))
    annihilator=s.eye(80)
    for value in eigenvalues:annihilator=annihilator*(matrix-value*s.eye(80))
    assert annihilator==s.zeros(80) and matrix.H==matrix
    top=s.eye(80)
    for value in eigenvalues[:-1]:top=top*(matrix-value*s.eye(80))/(2-value)
    assert top*top==top and s.trace(top)==2
    return {'actual_full_domain_identity': 'sum_beta J_beta^dagger J_beta=(tr K/2) N_L+sum_a K_a j_a',
            'K_ij': 'phi_i phi_j^dagger', 'K_a': 'tr(sigma_a K)',
            'j_a': 'sum_spin c_L^dagger sigma_a c_L/2',
            'dressed_weak_current_operator': 'W=sum_a K_a j_a',
            'all_twelve_Gauss_current_commutators_zero': True,
            'spin_scalar_and_even_parity': True,
            'raw_lepton_number_change': '0',
            'executed_full_active_dimension': 80, 'other_native_fermion_modes_are_identity_spectators': 26,
            'complete_possible_eigenvalues': list(map(str,eigenvalues)),
            'actual_full_operator_norm': '2', 'top_eigenprojector_rank_on_active_carrier': 2,
            'actual_eight_cell_ground_expectation_or_form_factor_calculated': False}


def calculate():
    if digest(ROOT/JOIN) != JOIN_SHA: raise ValueError('the same-native-H Higgs source changed')
    parent = json.loads((ROOT/JOIN).read_text()); join.validate(parent['input_sha256'])
    if join.certificate() != parent: raise ValueError('fresh same-source Higgs result differs')
    saved_map = json.loads((ROOT/MAP).read_text())
    join.validate(saved_map['input_sha256'])
    current_map = fields.certificate()
    if current_map != {k:v for k,v in saved_map.items() if k != 'input_sha256'}:
        raise ValueError('the actual complete current-map producer differs')
    lower = parent['same_H_cross_residue_and_records']['positive_cross_selector']['source_only_energy_interval'][0]
    return {'exact_source_current_interface': {
                'actual_package': 'model/coherent_defect_wave/chiral_particle_map',
                'same_full_CAR_automorphism_independent_of_scalar_state': True,
                'all_thirty_native_modes_and_occupation_sectors_retained': True,
                'actual_current_map_controls': current_map['physical_currents'],
                'complete_four_mode_pair_CAR_controls': current_map['full_Fock_pair_control'],
                'separate_conditional_background_control': current_map['background_spectrum'],
                'state_or_Hamiltonian_replaced_by_field_variable_change': False,
                'physical_Lorentz_chirality_of_interacting_eight_cell_source_proved': False,
                'constant_background_BdG_intertwiner_applies_to_quantum_scalar_without_error': False},
            'fixed_scalar_carrier_obstruction': scalar_controls(lower),
            'exact_Higgs_dressed_weak_current_interface': dressed_weak_density(),
            'family_copy_symmetry': family_controls()}


def certificate():
    result = calculate(); parent = json.loads((ROOT/JOIN).read_text()); mapping = json.loads((ROOT/MAP).read_text())
    seeds = set(parent['input_sha256'])|set(mapping['input_sha256'])|{JOIN, MAP,
        'model/coherent_defect_wave/chiral_particle_map/PROOF.md',
        'model/coherent_defect_wave/chiral_source/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/higgs_composites/PROOF.md',
        'model/interacting_extension/higgs_join/PROOF.md',
        'model/particle_common_contract/provenance.py'}
    seeds.update(str(path.relative_to(ROOT)) for path in HERE.iterdir() if path.suffix in ('.py','.md'))
    seeds.add(str(RESEARCH.relative_to(ROOT)))
    result.update({'schema': 'oph.native_weak_current_interface.v1',
        'input_sha256': {name:digest(ROOT/name) for name in sorted(import_closure(seeds))},
        'local_sha256': {p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}})
    if digest(RESEARCH)!=RESEARCH_SHA:raise ValueError('the immutable weak-current source calculation changed')
    ignored={'schema','input_sha256','local_sha256'}
    science=lambda value:{key:item for key,item in value.items() if key not in ignored}
    if science(result)!=science(json.loads(RESEARCH.read_text())):
        raise ValueError('public weak-current science differs from the frozen calculation')
    result['independent_research_comparison']={
        'reference_sha256':RESEARCH_SHA,'all_scientific_fields_equal':True,
        'ignored_only_provenance_fields':sorted(ignored),'fresh_calculation_was_replayed':True}
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__); mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args(); path=HERE/'results.json'
    if args.check:
        saved=json.loads(path.read_text());join.validate(saved['input_sha256'])
        for name,value in saved['local_sha256'].items():
            if digest(HERE/name)!=value:raise ValueError('local source custody differs')
    actual=certificate()
    if args.write:path.write_text(json.dumps(actual,indent=2,sort_keys=True)+'\n')
    elif actual!=saved:raise ValueError('weak-current/scalar audit differs')
    print('PASS: exact same-CAR current map, scalar eigenbackground obstruction and family-copy control')


if __name__=='__main__':main()
