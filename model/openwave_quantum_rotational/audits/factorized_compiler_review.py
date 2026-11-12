"""Independent complete-lattice / exterior-CAR review of the factorized compiler.

Run from the repository root with BLAS/OMP restricted to one thread.
Numerical controls supplement the symbolic reasoning recorded in this receipt;
they do not replace the producer's rational error theorem.
"""
from __future__ import annotations
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent


def maxabs(a):
    return float(np.max(np.abs(a)))


def complete_periodic_lattice():
    """Build in position space, then independently Fourier/spin transform."""
    side = 3
    xyz = np.array(list(np.ndindex(side, side, side)))
    sites = len(xyz)
    lookup = {tuple(x): i for i, x in enumerate(xyz)}
    sigma = [np.array([[0, 1], [1, 0]], complex),
             np.array([[0, -1j], [1j, 0]]), np.diag([1., -1.])]
    eps = np.array([[0, 1], [-1, 0]], complex)
    eye2 = np.eye(2)
    h = np.zeros((8*sites, 8*sites), complex)
    for axis in range(3):
        T = np.zeros((sites, sites), complex)
        for xidx, x in enumerate(xyz):
            y = x.copy(); y[axis] = (y[axis]+1) % side
            T[xidx, lookup[tuple(y)]] = 1
        h += np.kron(0.5j*(T-T.conj().T), np.kron(eye2, np.kron(sigma[axis], eye2)))
    effect = .31
    W = np.array([[effect, np.sqrt(effect*(1-effect))],
                  [np.sqrt(effect*(1-effect)), 1-effect]])
    pair = np.kron(np.eye(sites), np.kron(sigma[0], np.kron(eps, W)))
    assert maxabs(pair+pair.T) < 1e-15
    fourier = np.exp(2j*np.pi*(xyz@xyz.T)/side)/np.sqrt(sites)
    Ffull = np.kron(fourier, np.eye(8))
    spin = np.zeros_like(h)
    minus = [lookup[tuple((-p) % side)] for p in xyz]
    rotations = {}
    for p in range(sites):
        rep = min(p, minus[p])
        if rep not in rotations:
            vector = np.sin(2*np.pi*xyz[rep]/side)
            S = sum(v*s for v, s in zip(vector, sigma))
            if np.linalg.norm(vector) < 1e-14:
                U = np.eye(2, dtype=complex)
            else:
                _, U = np.linalg.eigh(S)
                U[:, 1] /= np.linalg.det(U)
            rotations[rep] = U
        U = rotations[rep]
        assert maxabs(U.conj().T@eps@U.conj()-eps) < 1e-14
        spin[8*p:8*p+8, 8*p:8*p+8] = np.kron(eye2, np.kron(U, eye2))
    change = Ffull@spin
    hp = change.conj().T@h@change
    dp = change.conj().T@pair@change.conj()
    def idx(p, leg, s, flag): return 8*p+4*leg+2*s+flag
    factors = [[idx(p, 0, s, f) for f in range(2)] +
               [idx(minus[p], 1, 1-s, f) for f in range(2)]
               for p in range(sites) for s in range(2)]
    assert sorted(i for g in factors for i in g) == list(range(8*sites))
    support_error = 0.; scalar_error = 0.; pair_error = 0.
    for p in range(sites):
        for s in range(2):
            g = factors[2*p+s]
            other = [i for i in range(8*sites) if i not in g]
            support_error = max(support_error, maxabs(hp[np.ix_(g, other)]), maxabs(dp[np.ix_(g, other)]))
            kappa = hp[g[0], g[0]].real
            scalar_error = max(scalar_error, maxabs(hp[np.ix_(g, g)]-kappa*np.eye(4)))
            expected = np.block([[np.zeros((2, 2)), (-1)**s*W],
                                 [-(-1)**s*W, np.zeros((2, 2))]])
            pair_error = max(pair_error, maxabs(dp[np.ix_(g, g)]-expected))
    assert max(support_error, scalar_error, pair_error) < 2e-14
    return {"side": side, "sites": sites, "complete_CAR_modes_one_character": 8*sites,
            "complete_Nambu_rows": 16*sites, "four_mode_factors_one_character": len(factors),
            "three_characters_factor_count": 3*len(factors),
            "zero_and_opposite_momenta_retained": True,
            "position_space_to_factor_support_max_error": support_error,
            "same_kappa_on_paired_modes_max_error": scalar_error,
            "epsilon_pairing_max_error": pair_error,
            "no_physical_Fourier_or_helicity_control_added": True}


def annihilator(n, mode):
    a = np.zeros((2**n, 2**n), complex)
    for mask in range(2**n):
        if mask & (1 << mode):
            a[mask ^ (1 << mode), mask] = (-1)**((mask & ((1 << mode)-1)).bit_count())
    return a


def own_block(n, modes, effect, kappa, phase):
    c = [annihilator(n, i) for i in modes]
    W = np.array([[effect, np.sqrt(effect*(1-effect))],
                  [np.sqrt(effect*(1-effect)), 1-effect]])
    H = kappa*sum(a.conj().T@a for a in c)
    V = sum(W[i, j]*c[i].conj().T@c[j+2].conj().T for i in range(2) for j in range(2))
    H = H+V+V.conj().T
    counts = np.array([((mask >> modes[1]) & 1)+((mask >> modes[3]) & 1) for mask in range(2**n)])
    r = (1j)**(phase*counts)
    return r.conj()[:, None]*H*r[None, :]


def exterior_factor_and_scalar_phase():
    params = [(.31, .71), (.67, -.9)]
    factor_errors = []; pulse_errors = []; centered_errors = []
    global_letters = []
    for phase in range(4):
        small = [own_block(4, range(4), e, k, phase) for e, k in params]
        whole = own_block(8, range(4), *params[0], phase)+own_block(8, range(4, 8), *params[1], phase)
        expected = np.kron(np.eye(16), small[0])+np.kron(small[1], np.eye(16))
        factor_errors.append(maxabs(whole-expected))
        dt = .013*(phase+1)
        pulse = expm(-1j*dt*whole)
        pulse_errors.append(maxabs(pulse-np.kron(expm(-1j*dt*small[1]), expm(-1j*dt*small[0]))))
        scalar = 2*sum(k for _, k in params)
        centered_errors.append(maxabs(pulse-np.exp(-1j*dt*scalar)*expm(-1j*dt*(whole-scalar*np.eye(256)))))
        global_letters.append(whole)
    assert max(factor_errors+pulse_errors+centered_errors) < 1e-13
    # A physical, even hopping between distinct factors invalidates the factorization.
    c = annihilator(8, 3); d = annihilator(8, 4)
    hopping = .2*(c.conj().T@d+d.conj().T@c)
    mutation = maxabs(expm(-.04j*(global_letters[0]+hopping))-expm(-.04j*global_letters[0]))
    assert mutation > .001
    # Each symmetric word has signed weights summing to one; composition 4p+q=1.
    rational_second_sum = sum([F(1,8), F(1,8), F(1,8), F(1,4), F(1,8), F(1,8), F(1,8)])
    assert rational_second_sum == 1
    return {"full_CAR_modes": 8, "full_Fock_rows": 256,
            "all_four_phase_source_factor_errors": factor_errors,
            "all_four_exact_pulse_factor_errors": pulse_errors,
            "all_four_scalar_phase_errors": centered_errors,
            "same_phase_reason": "second-order signed time is 1 and 4p+(1-4p)=1; scalar 2*kappa is subtracted only in both error operands",
            "even_interfactor_hopping_mutation_defect": mutation,
            "fermionic_exterior_signs_explicitly_constructed": True}


def energy_and_resources():
    rng = np.random.default_rng(217902)
    worst_ratio = 0.
    for _ in range(32):
        d, reference = 7, 3
        X = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
        Q, _ = np.linalg.qr(X)
        levels = rng.random(d)*5
        A = Q@np.diag(levels)@Q.conj().T
        Aref = np.kron(A, np.eye(reference))
        v = rng.normal(size=d*reference)+1j*rng.normal(size=d*reference); v /= np.linalg.norm(v)
        Y = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d)); Y = (Y+Y.conj().T)/2
        actual = np.kron(expm(-.015j*Y), np.eye(reference))@v
        w = actual-v
        delta = np.linalg.norm(w)
        mu = np.vdot(v, Aref@v).real
        limit = 2*np.sqrt(max(levels)*mu)*delta+max(levels)*delta**2
        error = abs(np.vdot(actual, Aref@actual).real-mu)
        assert error <= limit+1e-14
        worst_ratio = max(worst_ratio, error/limit)
    receipt = json.loads((ROOT/'model/openwave_factorized_compiler/results.json').read_text())
    r = F(35, 3)
    def eps(n, time, sites): return min(F(2), F(6*sites, 10)*r**5*time**5/n**4)
    resource = receipt['resources']
    controls = []
    for row in resource['finite_examples']:
        n = row['cycles']; s = row['sites']; target = F(row['target_trace_error'])
        assert eps(n, F(1), s) == F(row['complete_trace_bound']) <= target
        assert n >= r and eps(n-1, F(1), s) > target
        controls.append({"sites": s, "minimal_admissible_cycles": n})
    packet = resource['same_localized_packet']
    time = F(packet['source_time']); sites = packet['spatial_sites']; energy = F(282*sites)
    for name in ('same_localized_packet', 'state_sensitive_same_packet'):
        row = resource[name]; n = row['cycles']; actualeps = eps(n, time, sites)
        assert actualeps == F(row['complete_trace_bound']) <= F(1, 10**12)
        assert F(row['physical_energy_mean_error_upper']) <= F(row['physical_energy_error_budget'])
        if name == 'same_localized_packet':
            assert energy*actualeps == F(row['physical_energy_mean_error_upper'])
        else:
            # Verify energy inequality independently, without importing producer sqrt_bounds.
            mu = F(row['target_mean_energy_upper']); delta = actualeps/2
            upper = F(row['physical_energy_mean_error_upper'])
            residual = upper-energy*delta**2
            assert residual > 0 and residual**2 >= 4*energy*mu*delta**2
            previous = eps(n-1, time, sites)/2
            budget = F(row['physical_energy_error_budget'])
            residbudget = budget-energy*previous**2
            assert 2*previous > F(1,10**12) or residbudget < 0 or residbudget**2 < 4*energy*mu*previous**2
    return {"arbitrary_reference_positive_energy_random_controls": 32,
            "maximum_observed_fraction_of_state_sensitive_bound": worst_ratio,
            "finite_resource_binary_search_minimality": controls,
            "localized_packet_trace_and_energy_inequalities_exactly_rechecked": True,
            "state_sensitive_integer_optimality_rechecked_without_producer_square_root": True,
            "state_domain": "target physical excitation mean <= m_upper+d; complete ancilla/reference outputs retained; no ancillary energy guarantee"}


def main():
    data = {"schema": "oph.independent_factorized_compiler_review.v1",
            "reviewer_package": "openwave_quantum_rotational",
            "target": "openwave_factorized_compiler",
            "independent_position_space_Fourier_control": complete_periodic_lattice(),
            "independent_exterior_CAR_and_global_phase_control": exterior_factor_and_scalar_phase(),
            "independent_state_energy_and_exact_resources": energy_and_resources(),
            "review_findings": [
                "The homogeneous nearest-neighbor source has the stated 6*sites even four-CAR factors after a mathematical Fourier/SU(2) change of coordinates, including self-inverse momentum.",
                "The centered norm estimate preserves the original scalar phase because actual and ideal products have identical signed duration.",
                "The state-sensitive energy bound follows by expanding a purified vector and requires positivity, a target mean bound, and the unsaturated unitary error.",
                "The disjoint-factor hypothesis excludes extra inter-factor interactions; no large-carrier execution, local hardware efficiency, or source clock is inferred."
            ],
            "numerical_diagnostics_are_not_interval_proofs": True,
            "verdict": "PASS_INDEPENDENT_FACTOR_AND_STATE_ENERGY_REVIEW"}
    paths = [Path(__file__).resolve(), ROOT/'model/openwave_factorized_compiler/compiler.py',
             ROOT/'model/openwave_factorized_compiler/PROOF.md', ROOT/'model/openwave_factorized_compiler/results.json',
             ROOT/'model/openwave_kinetic_phase_cycle/compiler.py', ROOT/'model/openwave_kinetic_phase_cycle/PROOF.md']
    data['source_sha256'] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE/'factorized_compiler_review.json').write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
    print(data['verdict'])


if __name__ == '__main__':
    main()
