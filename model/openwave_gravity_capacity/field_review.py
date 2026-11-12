"""Independent small-matrix audit of the finite field-observable construction."""
from fractions import Fraction as F
from itertools import product
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import sympy as s

from model.openwave_field_observables import lattice, orbitals, mass_shell

ROOT = Path(__file__).resolve().parents[2]


def independent_periodic_operators(n):
    shift = np.roll(np.eye(n), -1, axis=1)-np.eye(n)
    identity = np.eye(n)
    derivatives = [np.kron(np.kron(shift if j == 0 else identity,
                                  shift if j == 1 else identity),
                           shift if j == 2 else identity) for j in range(3)]
    # Pick the forward convention by its action on one nonperiodic test word.
    test = np.arange(n**3).reshape(n, n, n)
    for j, derivative in enumerate(derivatives):
        if not np.array_equal((derivative@test.ravel()).reshape(test.shape),
                              np.roll(test, -1, axis=j)-test):
            derivatives[j] = derivative.T
    d = np.stack(derivatives, axis=1).reshape(3*n**3, n**3)
    zero = np.zeros((n**3, n**3))
    cblocks = [[-derivatives[1], derivatives[0], zero],
               [-derivatives[2], zero, derivatives[0]],
               [zero, -derivatives[2], derivatives[1]]]
    c = np.block(cblocks)
    order = np.arange(3*n**3).reshape(3, n**3).T.ravel()
    return d, c[np.ix_(order, order)]


def quantum_cutoff():
    states = [x for x in product(range(3), repeat=3) if sum(x) <= 2]
    lookup = {x:i for i, x in enumerate(states)}
    frequencies = (s.Integer(1), s.Integer(2), s.Integer(5))
    energies = [sum(a*b for a,b in zip(state, frequencies)) for state in states]
    h = s.diag(*energies)
    quadratic = s.zeros(len(states))
    for j, w in enumerate(frequencies):
        lowering = s.zeros(len(states))
        for i, x in enumerate(states):
            if x[j]:
                y = list(x); y[j] -= 1
                lowering[lookup[tuple(y)], i] = s.sqrt(x[j])
        q = (lowering+lowering.T)/s.sqrt(2*w)
        p = s.I*s.sqrt(w/2)*(lowering.T-lowering)
        assert (s.I*(h*q-q*h)-p).applyfunc(s.simplify) == s.zeros(len(states))
        assert (s.I*(h*p-p*h)+w*w*q).applyfunc(s.simplify) == s.zeros(len(states))
        quadratic += (p*p+w*w*q*q)/2
    top = s.diag(*(int(sum(x) == 2) for x in states))
    correction = top*(h+sum(frequencies)*s.eye(len(states)))/2
    expected = h+sum(frequencies)*s.eye(len(states))/2
    assert (quadratic+correction-expected).applyfunc(s.simplify) == s.zeros(len(states))
    assert quadratic != expected
    return {'independent_modes': 3, 'cutoff': 2, 'dimension': len(states),
            'complete_Heisenberg_equations_exact': True,
            'positive_top_boundary_correction_exact': True,
            'omitted_boundary_term_detected': True}


def run():
    files = [p for p in (ROOT/'model/openwave_field_observables').glob('*')
             if p.suffix in ('.py', '.md')]
    files.append(Path(__file__))
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    rng = np.random.default_rng(24119)
    n = 4
    d, c = independent_periodic_operators(n)
    assert np.linalg.norm(c@d) == 0
    x = rng.normal(size=n**3)
    e = rng.normal(size=3*n**3)
    b = rng.normal(size=3*n**3)
    errors = {
        'gradient': float(np.max(np.abs(d@x-lattice.gradient(x.reshape(n,n,n)).ravel()))),
        'curl': float(np.max(np.abs(c@e-lattice.curl(e.reshape(n,n,n,3)).ravel()))),
        'adjoint': float(np.max(np.abs(c.T@b-lattice.curl_adjoint(b.reshape(n,n,n,3)).ravel()))),
    }
    assert max(errors.values()) < 2e-15
    lap = d.T@d
    pinverse = np.linalg.pinv(lap, hermitian=True)
    green_error = float(np.max(np.abs(pinverse[:,0]-lattice.green(n).ravel())))
    assert green_error < 2e-15
    rho = np.zeros(n**3); rho[0]=1; rho[21]=-1
    phi = pinverse@rho
    assert np.max(np.abs(lap@phi-rho)) < 2e-14
    energy = float(np.linalg.norm(d@phi)**2/2)
    assert abs(energy-float(rho@phi/2)) < 1e-14
    assert np.max(np.abs(lap@(.5*phi)-rho)) > .49

    n=7
    lap1, potential, spacing = orbitals.dirichlet_source(n, 12.)
    eye=sparse.eye(n)
    lap3=(sparse.kron(sparse.kron(lap1,eye),eye)
          +sparse.kron(sparse.kron(eye,lap1),eye)
          +sparse.kron(sparse.kron(eye,eye),lap1)).tocsr()
    delta=np.zeros(n**3); delta[(n**3)//2]=1
    direct_green=spsolve(lap3,delta)
    dirichlet_error=float(np.max(np.abs(-4*np.pi*direct_green/spacing-potential.ravel())))
    assert dirichlet_error < 5e-15
    full=(lap3/(2*spacing**2)+sparse.diags(potential.ravel())).toarray()
    all_eigenvalues=[]
    max_intertwining=0.
    for parity in product((False,True), repeat=3):
        frames=[]
        for odd in parity:
            vectors=[]
            for j in range(n//2+int(odd),n):
                v=np.zeros(n)
                if j==n//2:
                    v[j]=1
                else:
                    v[j]=1/np.sqrt(2);v[n-1-j]=(-1 if odd else 1)/np.sqrt(2)
                vectors.append(v)
            frames.append(np.column_stack(vectors))
        iso=np.kron(np.kron(frames[0],frames[1]),frames[2])
        reduced,_=orbitals.reduced_hamiltonian(n,12.,parity)
        error=float(np.max(np.abs(full@iso-iso@reduced.toarray())))
        max_intertwining=max(max_intertwining,error)
        all_eigenvalues.extend(np.linalg.eigvalsh(reduced.toarray()))
    spectrum_error=float(np.max(np.abs(np.linalg.eigvalsh(full)-np.sort(all_eigenvalues))))
    assert max_intertwining < 2e-15 and spectrum_error < 2e-14
    opposite=lap3/(2*spacing**2)-sparse.diags(potential.ravel())
    assert np.linalg.eigvalsh(opposite.toarray())[0]>0

    worst_shell=0.; worst_bound_ratio=0.
    masses=(.8,.3,.4,0.)
    for index in range(24):
        spacing=2.**(-1-index%4)
        p=rng.uniform(-1,1,3)
        h=mass_shell.hamiltonian(p,spacing,masses)
        kinetic=np.sum((np.sin(spacing*p)/spacing)**2)
        expected=sorted([sign*np.sqrt(kinetic+mass*mass)
                         for mass in (.3,1.3) for sign in (-1,1) for _ in range(4)])
        err=float(np.max(np.abs(np.linalg.eigvalsh(h)-expected)))
        worst_shell=max(worst_shell,err)
        delta=np.dot(p,p)-kinetic
        bound=spacing**2*np.sum(p**4)/3
        assert -1e-14 <= delta <= bound+1e-14
        worst_bound_ratio=max(worst_bound_ratio,float(delta/bound))
    assert worst_shell<3e-14
    h=mass_shell.hamiltonian((.3,.2,-.1),.25,masses)
    assert abs(np.trace(h+.025*np.eye(16)))>.39
    if pins != {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}:
        raise ValueError('field source changed during independent audit')
    return {'schema':'oph.independent_field_observables_audit.v1','verdict':'PASS_SCOPED_CONSEQUENCES',
            'source_sha256':dict(sorted(pins.items())),
            'quantum_cutoff':quantum_cutoff(),'incidence_errors':errors,
            'periodic_Green_dense_inverse_error':green_error,
            'Dirichlet_Green_sparse_inverse_error':dirichlet_error,
            'all_eight_parity_sectors_compared':True,'complete_small_grid_dimension':n**3,
            'maximum_parity_intertwining_error':max_intertwining,
            'full_spectrum_vs_eight_reduced_spectra_error':spectrum_error,
            'alternative_mass_profile':list(masses),'independent_mass_momenta':24,
            'maximum_signed_spectrum_error':worst_shell,
            'maximum_remainder_to_analytic_bound_ratio':worst_bound_ratio,
            'mutations':{'Gauss_amplitude_change_detected':True,'repulsive_atomic_sign_detected':True,
                         'cutoff_completion_deletion_detected':True,'unreported_energy_origin_detected':True},
            'scope':{'free_field_and_declared_bound_state_results_confirmed':True,
                     'source_selected_physical_alpha_or_atomic_mass_derived':False,
                     'complete_interacting_particle_source_join_established':False,
                     'finite_lattice_exact_Lorentz_invariance':False,
                     'finest_95_grid_atomic_solve_repeated_by_this_auditor':False},
            'review_notes':['Time reversal and energy conservation fix the reciprocal Ampere block once Faraday curl is supplied.',
                            'The positive top-cutoff energy completion is essential; canonical commutation at the finite boundary is not assumed.',
                            'Periodic field and Dirichlet atom use the same interior stencil but distinct boundary conditions; atom also supplies a kinetic mass and g squared equal to 4pi.',
                            'Finite cubic split levels and the observed 0.42040 versus 0.375 spacing remain a numerical refinement diagnostic, not a measured hydrogen match.']}


if __name__=='__main__':
    result=run()
    Path(__file__).with_name('field_review.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['verdict'])
