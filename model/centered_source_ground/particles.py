"""Exact charge-conjugate one-fermion response of the centered cube.

The exact rank-one, variance and spin certificate gives a unique N=8, S=0 ground g of
K=H-N/2.  Since K commutes with eta_+ and eta_-, uniqueness forces both
ladders to annihilate g; otherwise they give a ground at N=10 or N=6.
Thus g has (S,eta)=(0,0).  The checked CAR tensor identities imply
S^2 O g=eta^2 O g=3 O g/4 for every linear creator or annihilator O.
Commutation with K extends this to their entire K-cyclic subspaces.
This is an inclusion, not a claim that those subspaces fill every copy.

Let C be the actual signed bit-complement unitary in the hidden frame.
The full-source identities below give C K C*=K and C d_i C*=eta_i d_i*.
Uniqueness gives Cg=phase*g.  For every Borel energy set B, with
P_B=1_B(K-Eg), the two positive matrix measures

    A_ij(B)=<g|d_i P_B d_j*|g>, B_ij(B)=<g|d_i* P_B d_j|g>

therefore obey A=D B D, D_ii=eta_i.  In the physical frame c=R* d,
the relation is A=F B F*, F=R* D conjugate(R).  The usual removal
contribution to a retarded matrix uses the transpose of this B.  The
statement retains off-diagonal complex entries, all energies and dark
weights; it does not assume isolated poles or compute their residues.

The same C reverses each centered local number and oriented bond current.
More generally C K(A) C*=K(-A) for external Peierls links, since it
interchanges each directed hopping with its adjoint.  This fixes the
relative CAR charge/current sign; electromagnetic normalization and a
Lorentz/CPT identification are not supplied.

Every existing cube repair support column has N=6.  Consequently all
120 native raw maps vanish on the entire N=7 direct-sum N=9 domain and
all ten outside maps are identity there.  K preserves each number, so
every finite interleaving retains the unitary charged response and an
all-outside repair history, including arbitrary reference cross terms.
There is no active native defect record, local protecting invariant,
new preparation, quantified gap, physical mass or species theorem here.
"""
from functools import lru_cache
import math

import numpy as np
from scipy import sparse

from ..cube_eta_protection import cube as source, repair
from ..field_stress.energy import GaussianInteger as GI, commutator, sum_exact
from ..field_symmetry import symmetry
from ..interacting_vacuum import ground
from . import gap


DIMENSION = 1 << 16


def _zero(value, message):
    if not value.is_zero():
        raise ArithmeticError(message)


def _times(value, coefficient):
    return symmetry._scale(value, coefficient)


@lru_cache(maxsize=1)
def _particle_hole_bytes():
    """Cache only immutable payloads, never mutable ndarray metadata."""
    masks = np.arange(DIMENSION, dtype=np.int64)
    parity_mask = sum(1 << (2*x+s) for x, site in enumerate(source.GEOMETRY.sites)
                      for s in (0, 1) if (sum(site)+s) % 2)
    complement = (DIMENSION-1)^masks
    phase = np.array([(-1)**((int(m)&parity_mask).bit_count()) for m in masks], dtype=np.int64)
    if not np.array_equal(phase[complement]*phase, np.ones(DIMENSION, dtype=np.int64)):
        raise ArithmeticError('particle-hole address map is not an involution')
    return complement.tobytes(), phase.tobytes()


def particle_hole_data():
    """Fresh views of immutable bytes for the full hidden-frame involution."""
    return tuple(np.frombuffer(raw, dtype=np.int64) for raw in _particle_hole_bytes())


def conjugate(value):
    """C value C*, without complex conjugating a unitary symmetry's input."""
    if type(value) is not GI or value.real.shape != (DIMENSION, DIMENSION):
        raise ValueError('a full-cube Gaussian-integer operator is required')
    complement, phase = particle_hole_data()
    def part(a):
        out = a[complement][:, complement].multiply(phase[:, None]).multiply(phase[None, :]).tocsr()
        out.eliminate_zeros()
        return out
    return GI(part(value.real), part(value.imag))


@lru_cache(maxsize=1)
def _fields():
    return tuple(symmetry._annihilator(16, j) for j in range(16))


def _centered_action():
    ops = source.operators()
    if ops.scale != 16 or ops.hubbard_H.imag.nnz:
        raise ArithmeticError('the exact real numerator16H changed')
    return ops.hubbard_H-_times(ops.N, 8)


def selection_certificate():
    """Actual full-source Lie, tensor and action identities, before any state."""
    ops = source.operators(); fields = _fields(); k = _centered_action()
    identity = GI.from_exact(sparse.eye(DIMENSION, dtype=np.int64, format='csr'))
    spin_plus = sum_exact([fields[2*x].adjoint()@fields[2*x+1] for x in range(8)], DIMENSION)
    spin_z_twice = sum_exact([fields[2*x].adjoint()@fields[2*x]-fields[2*x+1].adjoint()@fields[2*x+1]
                             for x in range(8)], DIMENSION)
    eta_plus = ops.eta_plus
    eta_z_twice = ops.N-_times(identity, 8)
    _zero(commutator(spin_plus, spin_plus.adjoint())-spin_z_twice, 'spin algebra changed')
    _zero(commutator(eta_plus, eta_plus.adjoint())-eta_z_twice, 'eta algebra changed')
    _zero(commutator(spin_z_twice, spin_plus)-_times(spin_plus, 2), 'spin ladder changed')
    _zero(commutator(eta_z_twice, eta_plus)-_times(eta_plus, 2), 'eta ladder changed')
    for spin in (spin_plus, spin_z_twice):
        for eta in (eta_plus, eta_z_twice):
            _zero(commutator(spin, eta), 'spin and eta do not commute')
    for generator in (spin_plus, spin_z_twice, eta_plus, ops.N):
        _zero(commutator(k, generator), 'centered action lost a required symmetry')
    _zero(k-k.adjoint(), 'centered action is not Hermitian')
    checked = 0
    for x, site in enumerate(source.GEOMETRY.sites):
        up, down = fields[2*x:2*x+2]; sign = (-1)**sum(site)
        relations = (
            (spin_z_twice, up.adjoint(), up.adjoint()),
            (spin_z_twice, down.adjoint(), _times(down.adjoint(), -1)),
            (spin_plus, up.adjoint(), GI.zero(DIMENSION)),
            (spin_plus, down.adjoint(), up.adjoint()),
            (spin_plus.adjoint(), up.adjoint(), down.adjoint()),
            (spin_plus.adjoint(), down.adjoint(), GI.zero(DIMENSION)),
            (eta_plus, up, _times(down.adjoint(), -sign)),
            (eta_plus, down, _times(up.adjoint(), sign)),
            (eta_plus, up.adjoint(), GI.zero(DIMENSION)),
            (eta_plus, down.adjoint(), GI.zero(DIMENSION)),
            (ops.N, up.adjoint(), up.adjoint()),
            (ops.N, down.adjoint(), down.adjoint()),
        )
        for generator, field, expected in relations:
            _zero(commutator(generator, field)-expected, 'full-CAR spin/eta tensor identity changed')
            checked += 1
    # Bind the hidden spin action to the actual physical onsite frame.
    _zero(ops.frame@ops.frame.adjoint()-identity, 'physical Fock frame is not unitary')
    _zero(ops.frame@ops.K@ops.frame.adjoint()-k, 'centered physical action/frame mismatch')
    hidden_spins = (spin_plus+spin_plus.adjoint(),
                    (spin_plus.adjoint()-spin_plus).times_i(), spin_z_twice)
    for actual, hidden in zip(ops.spin_twice, hidden_spins, strict=True):
        _zero(ops.frame@actual@ops.frame.adjoint()-hidden, 'physical spin frame mismatch')
    reference = gap.certificate()
    if (reference['full_ground_multiplicity'] != 1 or reference['ground_hidden_spin'] != 0
            or reference['full_source_SU2_and_spin_frame_checked'] is not True
            or reference['imported_Lieb_theorem_used_for_gap_or_uniqueness'] is not False):
        raise ArithmeticError('the exact actual-ground certificate changed')
    rows = {(row['twice_spin'], row['twice_eta']): row['multiplicity']
            for row in symmetry.irrep_multiplicities(8)}
    return {'full_Fock_dimension': DIMENSION, 'exact_CAR_tensor_identities': checked,
            'ground_twice_spin_twice_eta': [0, 0], 'joint_singlet_isotypic_dimension': rows[0, 0],
            'charged_twice_spin_twice_eta': [1, 1],
            'charged_joint_irrep_multiplicity': rows[1, 1],
            'charged_joint_isotypic_dimension': 4*rows[1, 1],
            'charge_branch_isotypic_dimension': 2*rows[1, 1],
            'full_K_cyclic_response_inclusion': True,
            'every_isotypic_copy_reachable_claimed': False,
            'numerical_eigenvalue_or_residual_used_for_admission': False,
            'exact_ground_certificate_used': True,
            'ground_uniqueness_input': 'exact balanced rank-one, variance and spin certificate with the full-charge ground bound',
            'internal_spin_is_physical_Lorentz_spin': False}


def current_operators():
    """Oriented x<y bond currents, as exact numerators16J in both frames."""
    ops = source.operators(); fields = _fields(); a = ground.band_data()[0]
    index = {site: x for x, site in enumerate(source.GEOMETRY.sites)}
    bonds = []
    for x, site in enumerate(source.GEOMETRY.sites):
        for axis in range(3):
            neighbor = list(site); neighbor[axis] += 1
            y = index.get(tuple(neighbor))
            if y is not None: bonds.append((x, y))
    if len(bonds) != len(ops.bonds): raise ArithmeticError('physical bond catalogue changed')
    result = []
    for (x, y), actual_bond in zip(bonds, ops.bonds, strict=True):
        forward = _times(sum_exact([fields[2*x+s].adjoint()@fields[2*y+s]
                                   for s in (0, 1)], DIMENSION), int(a[x, y]))
        _zero(conjugate(forward)-forward.adjoint(), 'directed hopping does not reverse under charge conjugation')
        hidden = (forward-forward.adjoint()).times_i()
        physical = ops.frame.adjoint()@hidden@ops.frame
        number_x = sum_exact([fields[2*x+s].adjoint()@fields[2*x+s] for s in (0, 1)], DIMENSION)
        independent = commutator(number_x, actual_bond).times_i()
        _zero(physical-independent, 'current differs from the original physical bond orientation')
        result.append(((x, y), hidden, physical))
    return tuple(result)


def response_certificate():
    ops = source.operators(); k = _centered_action(); fields = _fields()
    identity = GI.from_exact(sparse.eye(DIMENSION, dtype=np.int64, format='csr'))
    _zero(conjugate(k)-k, 'charge conjugation does not preserve the centered clock')
    for j, field in enumerate(fields):
        _zero(conjugate(field)-_times(field.adjoint(), (-1)**sum(source.GEOMETRY.sites[j//2])),
              'actual particle-hole CAR map changed')
    currents = current_operators(); divergences = [GI.zero(DIMENSION) for _ in range(8)]
    for (x, y), hidden, physical in currents:
        _zero(conjugate(hidden)+hidden, 'charge-conjugate current sign changed')
        _zero(physical-physical.adjoint(), 'physical current is not Hermitian')
        divergences[x] = divergences[x]+hidden; divergences[y] = divergences[y]-hidden
    for x in range(8):
        number = sum_exact([fields[2*x+s].adjoint()@fields[2*x+s] for s in (0, 1)], DIMENSION)
        centered = number-identity
        _zero(conjugate(centered)+centered, 'local excess charge does not reverse')
        _zero(commutator(k, number).times_i()+divergences[x], 'local number-current continuity changed')
    physical_signs = []
    for site in source.GEOMETRY.sites:
        frame = symmetry.onsite_frame(site)
        local = (-1)**sum(site)*frame.conj().T@frame.conj()
        if not any(np.array_equal(local, sign*np.eye(2)) for sign in (-1, 1)):
            raise ArithmeticError('the actual frame does not preserve each local charge probe')
        physical_signs.append(int(local[0, 0].real))
    return {'full_Fock_matrix_entries_checked': True, 'CAR_particle_hole_maps_checked': 16,
            'oriented_bond_currents_checked': len(currents), 'local_charge_and_continuity_equations_checked': 8,
            'current_numerator_denominator': 16,
            'charge_conjugation': 'unitary C; C K C*=K; C d_i C*=eta_i d_i*',
            'physical_field_conjugation_matrix': 'F=R* D conjugate(R), c=R* d, D_ii=eta_i',
            'physical_C_operator': 'Fock_R* C Fock_R',
            'physical_site_conjugation_signs': physical_signs,
            'positive_matrix_measures': 'A_ij=<g|c_i P_B c_j*|g>; B_ij=<g|c_i* P_B c_j|g>',
            'complete_matrix_response_identity': 'A(B)=F B(B) F* for every Borel excitation-energy set',
            'retarded_removal_matrix_convention': 'transpose of B',
            'all_local_charge_measures_equal': True,
            'current_relation': 'C J_xy C*=-J_xy in the same physical frame',
            'external_link_relation': 'C K(A) C*=K(-A); external Peierls links are supplied',
            'physical_electromagnetic_normalization_or_CPT_claimed': False}


def repair_certificate():
    """Check all odd columns of the actual low-rank repair, without a dense Q."""
    a, b = repair.exact_embedding(); actual = repair.embedding()
    if a.shape != (DIMENSION, 24) or b.shape != a.shape or actual.shape != a.shape:
        raise ArithmeticError('complete repair embedding shape changed')
    for matrix in (a, b, actual):
        if any(int(row).bit_count() != 6 for row in matrix.nonzero()[0]):
            raise ArithmeticError('repair support escaped its exact N6 sector')
        for number in (7, 9):
            if matrix[source.weights(number)].nnz:
                raise ArithmeticError('an odd response column reaches native repair support')
    native = outside = 0
    for law in ('C', 'E'):
        for seam in repair.SEAMS:
            slots = repair.matrices(law, seam)
            labels = [('native', coin, raw) for coin in (0, 1) for raw in range(4 if law == 'C' else 8)]
            if [label for label, _ in slots] != labels:
                raise ArithmeticError('a complete native raw slot was omitted or relabelled')
            for _, matrix in slots:
                if (np.ma.isMaskedArray(matrix) or type(matrix) is not np.ndarray
                        or matrix.shape != (24, 24) or matrix.dtype != np.dtype(np.int64)
                        or min(int(matrix.min()), 0) < -1 or max(int(matrix.max()), 0) > 1):
                    raise ArithmeticError('unadmitted raw native partial injection')
            if not np.array_equal(sum(m.T@m for _, m in slots), 2*np.eye(24, dtype=np.int64)):
                raise ArithmeticError('complete native effects changed')
            native += len(slots); outside += 1
    return {'full_number_sectors': [7, 9], 'columns_per_number_sector': math.comb(16, 7),
            'joint_input_dimension': 2*math.comb(16, 7), 'native_raw_slots_retained': native,
            'outside_slots_retained': outside, 'repair_support_number': 6,
            'all_native_maps_zero_on_entire_domain': True, 'every_outside_map_identity_on_entire_domain': True,
            'arbitrary_reference_cross_terms_and_finite_histories_preserved': True,
            'nontrivial_native_defect_record_obtained': False,
            'new_native_preparation_or_detector_executed': False}


def mass_order_certificate():
    """Bind the four actual onsite mass channels to vanishing ground orders.

    The scalar channel is particle-hole odd.  The three other channels
    form a vector under the hidden global spin.  Invariance of the unique
    ground therefore kills their expectations, independently of any
    numerical ground vector.  This does not exclude an interacting gap,
    mass-order fluctuations, or symmetry breaking in another volume limit.
    """
    from ..native_mass_channels import algebra as mass
    ops = source.operators(); fields = _fields()
    hidden = [GI.zero(DIMENSION) for _ in range(4)]
    total_spin = [GI.zero(DIMENSION) for _ in range(3)]
    for x, site in enumerate(source.GEOMETRY.sites):
        up, down = fields[2*x:2*x+2]
        plus = up.adjoint()@down
        local = (up.adjoint()@up+down.adjoint()@down,
                 plus+plus.adjoint(), (plus.adjoint()-plus).times_i(),
                 up.adjoint()@up-down.adjoint()@down)
        for j in range(4): hidden[j] = hidden[j]+_times(local[j], (-1)**sum(site))
        for j in range(3): total_spin[j] = total_spin[j]+local[j+1]
    for j in range(4):
        blocks = [sparse.csr_matrix(np.asarray(mass.onsite_basis(site)[j]).astype(complex))
                  for site in source.GEOMETRY.sites]
        physical = symmetry._lift(sparse.block_diag(blocks, format='csr'), 16)
        _zero(ops.frame@physical@ops.frame.adjoint()-hidden[j], 'classified mass channel has the wrong physical frame')
    _zero(conjugate(hidden[0])+hidden[0], 'scalar staggered mass order is not particle-hole odd')
    for a, b, c in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        _zero(commutator(total_spin[a], hidden[b+1])-_times(hidden[c+1], 2j),
              'staggered mass spin vector does not transform under the actual spin')
    return {'classified_onsite_mass_channels_bound': 4,
            'physical_channels': ['(-1)^(x0+x1+x2) I', '(-1)^x0 sigma0', '(-1)^x1 sigma1', '(-1)^x2 sigma2'],
            'hidden_frame_channels': 'eta_x times (I,sigma0,sigma1,sigma2)',
            'centered_local_order_definition': 'c_x* B_x c_x-Tr(B_x)/2; the scalar is eta_x(N_x-1)',
            'uncentered_scalar_local_expectation': 'eta_x; its full balanced-cube sum is zero',
            'actual_unique_ground_order_expectations': ['0', '0', '0', '0'],
            'scalar_order_zero_reason': 'particle-hole invariance and the balanced sublattice sum',
            'spin_orders_zero_reason': 'the actual unique ground is a hidden-spin singlet',
            'finite_many_body_gap_is_evidence_for_onsite_mass_order': False,
            'mass_order_fluctuations_claimed_zero': False,
            'interacting_gap_or_thermodynamic_symmetry_breaking_excluded': False}


def certificate():
    return {'selection': selection_certificate(), 'response_and_current': response_certificate(),
            'complete_repair': repair_certificate(), 'onsite_mass_order': mass_order_certificate(),
            'new_gap_or_particle_mass_claimed': False, 'physical_species_identified': False}
