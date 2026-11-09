"""Joint source/reference translations and complete finite position records.

Crystal momentum belongs to the complete translated system. No individual
particle, physical momentum normalization or finite canonical commutator is
inferred. Small controls are explicitly separate from the full native field.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb

import sympy as s

from model.particle_local_gauge_state import state
from model.particle_local_wilson_density import density


SOURCE_INPUTS = (
    'model/particle_local_gauge_state/state.py',
    'model/particle_local_gauge_state/PROOF.md',
    'model/particle_local_history_family/PROOF.md',
    'model/particle_local_wilson_density/density.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/particle_observer_outputs/READOUT_CONTRACT.md',
    'model/periodic_momentum/README.md',
    'model/carrier_position_translation/PROOF.md',
    'Lean/QFT/PeriodicMomentum.lean',
    'Lean/EventAlgebra/Robertson.lean',
)


def clean(a):
    return s.ImmutableSparseMatrix(a).applyfunc(s.simplify)


def signed_fock_permutation(permutation):
    """Full exterior-power lift, including vacuum and filled occupations."""
    p = tuple(permutation)
    n = len(p)
    if sorted(p) != list(range(n)):
        raise ValueError('a permutation of every mode is required')
    entries = {}
    for mask in range(1 << n):
        image = [p[j] for j in range(n) if mask >> j & 1]
        inversions = sum(a > b for i, a in enumerate(image) for b in image[i+1:])
        target = sum(1 << j for j in image)
        entries[target, mask] = (-1)**inversions
    return s.ImmutableSparseMatrix(1 << n, 1 << n, entries)


def position_effects(sites, modes_per_site):
    """Complete effects n_x/N and the retained N=0 outcome for a small control."""
    v = state.integer(sites, 'sites', 1)
    d = state.integer(modes_per_site, 'modes per site', 1)
    if v*d > 12:
        raise ValueError('explicit control exceeds twelve modes; use native_law')
    diagonals = [[s.S.Zero]*(1 << (v*d)) for _ in range(v+1)]
    for mask in range(1 << (v*d)):
        N = mask.bit_count()
        if N == 0:
            diagonals[-1][mask] = 1
        else:
            for x in range(v):
                nx = sum(mask >> (d*x+j) & 1 for j in range(d))
                diagonals[x][mask] = s.Rational(nx, N)
    return tuple(s.ImmutableSparseMatrix.diag(*row) for row in diagonals)


def fourier_projectors(translation, period):
    """Complete exact cyclic Fourier menu for a declared finite unitary."""
    n=state.integer(period,'period',1)
    t=clean(translation)
    if t.rows != t.cols or clean(t.H*t)!=s.eye(t.rows) or clean(t**n)!=s.eye(t.rows):
        raise ValueError('a unitary with the declared finite period is required')
    omega=s.expand_complex(s.exp(2*s.pi*s.I/n))
    return tuple(clean(sum((omega**(-k*r)*t**r for r in range(n)),
                           s.zeros(t.rows))/n) for k in range(n))


def initial_overlap(sites, modes_per_site):
    """Exact binomial expression for E[sqrt(n_x n_y)/N], x!=y.

    This control evaluates small inputs only. The full native expression
    is exported symbolically, without enumerating the native Fock carrier.
    """
    v = state.integer(sites, 'sites', 2)
    d = state.integer(modes_per_site, 'modes per site', 1)
    if v*d > 18:
        raise ValueError('the exact sum is a small-control evaluator')
    outside = d*(v-2)
    total = s.S.Zero
    for i in range(1, d+1):
        for j in range(1, d+1):
            harmonic = sum((s.Rational(comb(outside, k), i+j+k)
                            for k in range(outside+1)), s.S.Zero)
            total += comb(d, i)*comb(d, j)*s.sqrt(i*j)*harmonic
    return s.simplify(total/2**(d*v))


def native_law(side=2):
    """All original native60-mode labels, represented by exact finite formulas."""
    length = state.integer(side, 'periodic side', 2)
    V = length**3
    return {
        'periodic_side': length, 'cells': V, 'modes_per_cell': 60,
        'source_state': 'site-gauge-twirled finite link product tensor I_Fock/2^(60 V)',
        'joint_state': 'canonical purification with conjugate gauge and translation reference',
        'source_H': 'P_K [H_E+dGamma(h_m)-kappa sum_p ReTr R(U_p)/15] P_K',
        'same_clock': 'exp(-it H_source) tensor I_reference; H_reference=0',
        'joint_translation': 'T_source(a) tensor conjugate(T_source(a))',
        'position_effects': ['F_x=n_x N_inverse_on_positive_N', 'F_vac=1_(N=0)'],
        'joint_momentum_projectors': 'Pi_k=V^-1 sum_a exp(-2pi i k.a/side) T_joint(a)',
        'vacuum_probability': {'base': 2, 'inverse_exponent': 60*V},
        'A': '(1-2^(-60V))/V',
        'C': '2^(-60V) sum_i=1^60 sum_j=1^60 binom(60,i) binom(60,j) sqrt(i*j) sum_l=0^(60(V-2)) binom(60(V-2),l)/(i+j+l)',
        'original_joint_position_momentum_law': {
            'vacuum_then_k0': '2^(-60V)', 'vacuum_then_nonzero_k': '0',
            'each_x_then_k0': '[A+(V-1)C]/V',
            'each_x_then_each_nonzero_k': '(A-C)/V'},
        'A_minus_C_positive_lower': {'numerator': 60, 'base': 2, 'inverse_exponent': 60*V},
        'position_labels': V+1, 'momentum_labels': V,
        'all_joint_labels': (V+1)*V,
        'initial_joint_momentum_is_exactly_zero': True,
        'full_law_uses_original_unconditioned_weights': True,
        'position_or_momentum_is_individual_particle_property': False,
        'joint_gauge_zero_domain_preserved': True,
        'full_ambient_projector_completeness': True,
        'full_ambient_gauge_commutation_claimed': False,
        'physical_Heisenberg_constant_or_GeV_assignment': False,
        'new_empirical_inputs': [],
    }


@lru_cache(maxsize=1)
def small_reference_control():
    """Complete two-site, one-mode Z2 gauge control of the measurement contract.

    It is a structural finite control, not the native SU3/SU2 field. Its
    full8-state field and64-state referenced carrier are retained.
    """
    f = signed_fock_permutation((1, 0))
    eye2 = s.eye(2)
    link_x = s.Matrix([[0, 1], [1, 0]])
    tf = clean(s.kronecker_product(f, eye2))
    g0 = clean(s.kronecker_product(s.diag(1, -1, 1, -1), link_x))
    g1 = clean(s.kronecker_product(s.diag(1, 1, -1, -1), link_x))
    tj = clean(s.kronecker_product(tf, tf.conjugate()))
    gj = tuple(clean(s.kronecker_product(g, g.conjugate())) for g in (g0, g1))
    omega = s.ImmutableSparseMatrix(64, 1, {(8*i+i, 0): 1/s.sqrt(8) for i in range(8)})
    effects = tuple(clean(s.kronecker_product(e, eye2)) for e in position_effects(2, 1))
    projectors = ((s.eye(64)+tj)/2, (s.eye(64)-tj)/2)
    rows = []
    for label, e in zip(('site0', 'site1', 'vacuum'), effects):
        k = s.diag(*[s.sqrt(e[i, i]) for i in range(8)])
        psi = s.kronecker_product(k, s.eye(8))*omega
        for momentum, projector in enumerate(projectors):
            probability = s.simplify((psi.H*projector*psi)[0])
            rows.append({'position': label, 'momentum_index': momentum,
                         'original_probability': probability})
    x = clean(s.kronecker_product(effects[1], s.eye(8)))
    phase = s.diag(*[s.exp(s.I*s.pi*effects[1][i, i]/2) for i in range(8)])
    psi = s.kronecker_product(phase, s.eye(8))*omega
    p = clean(s.pi*projectors[1])
    mean = lambda a: s.simplify((psi.H*a*psi)[0])
    varx = s.simplify(mean(x*x)-mean(x)**2)
    varp = s.simplify(mean(p*p)-mean(p)**2)
    comm = s.simplify(mean(x*p-p*x))
    detector_var = s.simplify(mean(x)-mean(x)**2)
    source_only = s.kronecker_product(tf, s.eye(8))*omega
    wrong_charge_residual = s.simplify(((gj[0]-s.eye(64))*source_only).norm()**2)
    if wrong_charge_residual == 0:
        raise ArithmeticError('source-only translation failed to expose reference charge')
    if any(g*omega != omega for g in gj) or tj*omega != omega:
        raise ArithmeticError('canonical referenced state lacks its advertised symmetries')
    return {'kind': 'distinct complete Z2 structural control',
            'field_dimension': 8, 'joint_dimension': 64,
            'source_translation': tf, 'joint_translation': tj,
            'joint_Gauss_generators': gj, 'canonical_state': omega,
            'joint_position_momentum_records': rows,
            'source_only_translation_Gauss_residual_squared': wrong_charge_residual,
            'X_variance': varx, 'record_coordinate_variance': detector_var,
            'P_variance': varp, 'commutator_expectation': comm,
            'Robertson_rhs': s.simplify(s.Abs(comm)**2/4),
            'Robertson_lhs': s.simplify(varx*varp),
            'unsharp_coordinate_noise': s.simplify(detector_var-varx)}


@lru_cache(maxsize=1)
def native_covariance_control():
    """All480 native one-body columns on a side-two cubic carrier, m=0.

    The comparison relabels every coordinate link together with matter.
    Electric/plaquette and arbitrary-order covariance follow by the proof.
    No exponential-size full Fock matrix is allocated.
    """
    side = 2
    cells = tuple(product(range(side), repeat=3)); index = {x:i for i,x in enumerate(cells)}
    links = {(x, mu):s.eye(15) for x in cells for mu in range(3)}
    links[cells[0], 0] = density.native_weak_link(0, s.Rational(3,5), s.Rational(4,5))
    links[cells[1], 1] = density.native_weak_link(1, s.Rational(5,13), s.Rational(12,13))
    gamma = density.gamma_data(); beta = gamma[3]
    def shift(x, a): return tuple((x[j]+a[j])%side for j in range(3))
    def kernel(config):
        h = s.kronecker_product(beta*s.Rational(25,17), s.eye(120))
        h = s.MutableSparseMatrix(h)
        for mu in range(3):
            ts = s.MutableSparseMatrix.zeros(120)
            a = tuple(int(j==mu) for j in range(3))
            for x in cells:
                y = shift(x,a)
                block=config[x,mu]
                for i,j in block.todok():
                    ts[15*index[y]+i,15*index[x]+j]=block[i,j]
            # h0=(beta+(4/17)beta(2I-sum W_mu)), with same source spin projectors.
            h -= s.Rational(4,17)*s.kronecker_product(beta*(s.eye(4)-gamma[mu])/2,ts)
            h -= s.Rational(4,17)*s.kronecker_product(beta*(s.eye(4)+gamma[mu])/2,ts.H)
        return clean(h)
    h = kernel(links)
    checks=[]
    for axis in range(3):
        a=tuple(int(j==axis) for j in range(3))
        perm=s.MutableSparseMatrix.zeros(120)
        for x in cells:
            for j in range(15):perm[15*index[shift(x,a)]+j,15*index[x]+j]=1
        t=s.kronecker_product(s.eye(4),perm)
        moved={(shift(x,a),mu):u for (x,mu),u in links.items()}
        residual=clean(t*h*t.H-kernel(moved))
        wrong=clean(t*h*t.H-h)
        if residual != s.zeros(480) or wrong == s.zeros(480):
            raise ArithmeticError('full link/matter covariance or negative control failed')
        checks.append({'axis':axis,'complete_columns':480,
                       'joint_covariance_residual_zero':True,
                       'matter_only_translation_with_fixed_links_fails':True})
    return {'native_internal_rows':15,'Dirac_rows':4,'cells':8,
            'one_body_dimension':480,'Legendre_degree':0,
            'all_axes':checks,'full_Fock_lift_is_algebraic_not_dense_execution':True}


def serial(x):
    if isinstance(x,(F,s.Basic)):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [serial(v) for v in x]
    return x


def certificate():
    control=small_reference_control()
    public={k:v for k,v in control.items() if k not in
            ('source_translation','joint_translation','joint_Gauss_generators','canonical_state')}
    return serial({'schema':'oph.joint_translation_readout.v1',
                   'source_inputs':SOURCE_INPUTS,
                   'native_complete_law':native_law(),
                   'native_covariance':native_covariance_control(),
                   'complete_reference_control':public,
                   'finite_CCR_i_identity_claimed':False,
                   'joint_crystal_momentum_is_physical_particle_momentum':False})
