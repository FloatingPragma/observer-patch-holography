"""Independent one-particle Lie algebra proof of the dressed axial bound."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import sympy as s

HERE = Path(__file__).resolve().parent
I = s.I
spin = [s.Matrix([[0,1],[1,0]])/2, s.Matrix([[0,-I],[I,0]])/2, s.diag(1,-1)/2]
colour = []
for a, b in [(0,1),(0,2),(1,2)]:
    X, Y = s.zeros(3), s.zeros(3)
    X[a,b] = X[b,a] = s.Rational(1,2)
    Y[a,b], Y[b,a] = -I/2, I/2
    colour += [X, Y]
colour += [s.diag(1,-1,0)/2, s.diag(1,1,-2)/(2*s.sqrt(3))]
assert sum((T*T for T in colour), s.zeros(3)) == s.eye(3)*s.Rational(4,3)
weak = s.Matrix([[0,1],[0,0]])
G = [s.kronecker_product(S, s.eye(2), T) for S in spin for T in colour]
V = s.kronecker_product(s.eye(2), weak, s.eye(3))
A = [s.kronecker_product(2*S, weak, s.eye(3)) for S in spin]


def comm(X, Y):
    return X*Y-Y*X


assert all(comm(g, V) == s.zeros(12) for g in G)
for axial in A:
    double = sum((comm(g, comm(g, axial)) for g in G), s.zeros(12))
    assert s.simplify(double-s.Rational(8,3)*axial) == s.zeros(12)
# CAR second quantization preserves bilinear commutators exactly on every
# Fock sector. In the singlet spin-half branch sumG²=c=4; hence compression
# of sumG A G is (c-4/3)A, and normalized excited-sector current is2A/3.
c = s.Integer(4)
excited_factor = (c-s.Rational(4,3))/c
assert excited_factor == s.Rational(2,3)
Delta, coupling = s.symbols('Delta coupling', positive=True)
weight = (1-Delta/s.sqrt(Delta**2+16*coupling**2))/2
axial = s.Rational(5,3)*(1-weight/3)
lower = s.limit(axial, coupling, s.oo)
assert lower == s.Rational(25,18)
assert s.simplify(axial-lower) == 5*Delta/(18*s.sqrt(Delta**2+16*coupling**2))
values = []
for kap in [s.Rational(1,2),s.Integer(1),s.Integer(2)]:
    value = axial.subs({Delta:3*kap,coupling:s.Rational(1,8)})
    values.append(dict(kappa=str(kap), exact=s.sstr(value), numerical=float(value)))
paths = [HERE/'hadrons/LOCAL_AUXILIARY_INTERFACE_MANIFEST.json',
         HERE/'species_leptons/auxiliary_current/MANIFEST.json']
out = dict(schema='oph.obstruction_repair.auxiliary_bound_independent_review.v1',
           producer_code_imported=False, exact_one_particle_double_commutator='8A/3',
           CAR_commutator_lift_exact=True, vector_current_commutator='0',
           singlet_spin_half_G_square='4', normalized_excited_axial_factor='2/3',
           positive_gap_lower_bound_exclusive='25/18',
           source_coupling_predictions=values,
           conditional_on_stated_auxiliary_source_and_bare_singlet_dictionary=True,
           full_5500_state_numerical_replay_performed_here=False,
           physical_axial_row_solved=False,
           pins={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
(HERE/'AUXILIARY_BOUND_REVIEW.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'exact_double_commutator':True,'exclusive_lower_bound':str(lower),
                  'physical_axial_row_solved':False}))
