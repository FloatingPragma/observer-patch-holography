# Complete quadratic conservation for the Weyl/contact action

Let a finite open rectangular box have V sites and two fermionic modes at
each site. Every present nearest-neighbor hopping coefficient is nonzero.
Use the unchanged source action

```
H_mu = T + U D - mu N,
T = sum_<x,y=x+ej> c_x^* (i c sigma_j/2) c_y + adjoint,
D = sum_x n_x,up n_x,down,       c>0, U>0, mu real.
```

The original energy convention is H=H_0; the centered convention is
K=H_(U/2). The chemical term changes relative energies between number
sectors. It is not an additive scalar on the full Fock carrier, and the
two conventions are not identified as one clock. There is no restriction
to an eigenstate, a Gibbs reference or one particle-number sector.

An even self-adjoint polynomial of degree at most two has the unique form

```
Q = lambda I + dGamma(q) + A^* + A,
q=q^*,        A^*=sum_(p<q) Delta_pq c_p^* c_q^*.
```

The repeated q in the ordered-mode sum is an index, independent of the
one-particle matrix q. There are `(4V)(4V-1)/2` traceless real coefficients,
equivalently the coefficients of `i gamma_p gamma_q`, p<q, in the complete
Majorana basis. The scalar identity is a separate dimension.

**Theorem.** For every such connected box, the real vector space of
self-adjoint conserved even quadratics, modulo scalar identity, is

```
span_R {N-V, S_1, S_2, S_3}                         if mu != U/2,
span_R {N-V, S_1, S_2, S_3, eta_x, eta_y}          if mu = U/2,

S_a = (1/2) sum_x (-1)^(sum_(j!=a) x_j) c_x^* sigma_a c_x,
eta_+ = sum_x (-1)^(x_1+x_2+x_3) c_x,up^* c_x,down^*,
eta_x=(eta_+ + eta_-)/2,   eta_y=(eta_+ - eta_-)/(2i),
eta_-=(eta_+)^*,          eta_z=(N-V)/2.
```

Thus the traceless conserved Lie algebra is hidden `su(2)+u(1)` for H
and hidden `su(2)+su(2)` for K. Dimensions are four and six, or five and
seven after including the scalar identity. These are internal algebras;
no Lorentz-spin or measured-species identification is part of the theorem.

## Number-preserving part

Since H_mu commutes with N, its commutator preserves the number change of
each term of Q. Terms with number change zero, plus two and minus two
cannot cancel. First consider `dGamma(q)`.

The one-particle restriction of H_mu is h-mu I because D vanishes there.
Consequently conservation implies `[h,q]=0`. The CAR identity
`[dGamma(h),dGamma(q)]=dGamma([h,q])` gives kinetic commutation on the
entire Fock algebra, hence `[D,dGamma(q)]=0` separately.

Fix distinct sites x,y and mode indices p=(x,s), q=(y,t). Take the input
with particles at `(x,opposite(s))` and q, and the output doublon at x.
Only `c_p^* c_q` connects this input to this output. Its matrix element
is the corresponding one-body coefficient times a nonzero fermionic sign.
The D eigenvalues differ by one. Therefore the coefficient must vanish.
All intersite coefficients vanish, without assuming a particular spin form.

Every onsite bilinear commutes with its local doublon projection. Kinetic
commutation on a nonzero axis-j edge is precisely

```
q_x sigma_j = sigma_j q_y,       q_y = sigma_j q_x sigma_j.
```

The three conjugation maps commute and are involutions. An arbitrary
Hermitian root matrix `q_0=a I+sum_a b_a sigma_a` therefore extends
uniquely over the connected box, giving the four stated coefficients.
The transported matrices have equal norm. No nonzero member can vanish
on even one site while remaining supported elsewhere.

## Every intersite anomalous pair is excluded

The kinetic commutator of a quadratic hopping with a pair creator is
another pair creator, by direct CAR expansion. Write
`[T,A^*]=B^*=sum b_pq c_p^* c_q^*`.

For p and q on distinct sites, compare two entries of
`[H_mu,A^*]=0`. The vacuum-to-pair entry, divided by its CAR sign, is

```
b_pq - 2 mu Delta_pq = 0.
```

For the second entry, place one spectator fermion in the opposite spin
mode at the site of p. Create the same p,q pair. The input has no doublon
and the output has exactly one. Only this particular pair creator reaches
that output. A quadratic pair operator has the same coefficient in the
two entries after dividing by their individual CAR signs, so the second
equation is

```
b_pq + (U-2 mu) Delta_pq = 0.
```

Subtracting gives `U Delta_pq=0`. Since U>0, every intersite anomalous
coefficient vanishes. This argument retains the one-to-three transition;
a vacuum or one-number-sector calculation alone would miss its constraint.
`spectator_witnesses()` checks the signs and uniqueness for all 112 complex
intersite pair coefficients on the eight-site cube.

Only `A^*=sum_x f_x c_x,up^* c_x,down^*` can survive. Its interaction
commutator is exactly `U A^*`, and its number commutator is `2 A^*`.
The kinetic commutator has only intersite pair entries. Thus the onsite
entries of conservation give `(U-2 mu)f_x=0` at each site. When
mu differs from U/2 this eliminates every anomalous term.

At mu=U/2 only `[T,A^*]=0` remains. Write
`epsilon=[[0,1],[-1,0]]`. The elementary Pauli identity
`epsilon sigma_j^T=-sigma_j epsilon` gives the x,y pair block

```
(i c/2)(f_x+f_y) sigma_j epsilon.
```

It vanishes precisely when `f_y=-f_x`. Connectedness and bipartiteness of
the open box give the unique complex alternating profile
`f_x=f_0 (-1)^(x_1+x_2+x_3)`. Its real and imaginary parts supply eta_x
and eta_y. The complex conjugate equation supplies no additional real
dimension. On an odd-cycle real-hopping graph the alternating recursion
instead forces f=0. That control is a separate graph, not a periodic
identification of the stated open source.

The usual two-site CAR identities give the two commuting su(2) algebras;
equivalently onsite spin acts on the singly occupied doublet and onsite
eta on the empty/doubly occupied doublet. The signs in the displayed
generators are the exact source transport. This proves both necessity and
sufficiency at every positive c,U, including values with accidental energy
degeneracies. Degeneracies can enlarge the full Hamiltonian commutant; they
do not enlarge this degree-two subspace.

## Localization and the finite stress premise

Every nonzero conserved quadratic modulo scalar has a nonzero onsite
coefficient at every site. The one-body and anomalous coefficients cannot
cancel each other because they have different number change. Hence no
nonconstant conserved quadratic is supported in a proper site region of
this connected source. This excludes a strictly localized quadratic
protecting charge under this Hamiltonian, without excluding localized
states, metastability, higher-order invariants or a changing source law.

In particular, all conserved quadratics contain no intersite coefficient.
Allowing anomalous pair terms therefore supplies no omitted quadratic
charge with a kinetic momentum component. Such a charge cannot provide an
exact closed-boundary conserved momentum column for the finite interacting
source. The exact forced momentum equations in
[`lattice_momentum_balance`](../lattice_momentum_balance/README.md) retain
their contact forces. This is a statement about total quadratic charges,
not a proof that all possible local stress tensors fail. Higher-degree
dressings, nonquadratic translations, continuum limits and different
physical state/clock attachments require separate analysis.

## Independent exact finite certificates

The analytic argument covers arbitrary finite connected rectangular boxes.
The executable certificate has a declared resource domain of at most eight
sites and dyadic source coefficients with denominators at most 2^16. The
retained examples have c=1/8 and U=1. It reconstructs every coordinate link
from the Pauli specification and compares the full Clifford Hamiltonian
with the actual source's emitted groups. A changed hopping strength or
axis therefore fails even if the conserved algebra happens to be unchanged.

For each complete traceless Hermitian Majorana-quadratic basis, the exact
commutator is expanded in the linearly independent full Clifford algebra.
Real and imaginary rational coefficients give a real linear system. Exact
elimination modulo the trial-division-verified prime 65521 gives a lower
bound on its rational rank, because all admitted denominators invert. The
explicit four or six commuting generators are independent over that same
field, hence independent over the rationals and reals. They give the
matching upper rank bound. Agreement proves the exact rank and exhausts
all real coefficients; a rank modulo a prime by itself would not do so.

| Source box | Quadratic coefficients | Rank for H | Rank for K |
| --- | ---: | ---: | ---: |
| Two sites on each of the three axes | 28 | 24 | 22 |
| Three-site chain | 66 | 62 | 60 |
| Four-site plaquette | 120 | 116 | 114 |
| Eight-site cube | 496 | 492 | 490 |

An independent Jordan-Wigner occupation construction forms the complete
256-by-256 plaquette source matrices and all 120 quadratic generators.
Its matrix-entry commutator ranks agree at 116 and 114. All 65,536
input/output entries and all number sectors enter these checks. Integer
sparse products use the predecessor's checked overflow bounds. The cube
certificate is symbolic coverage of the full 65,536-state algebra, not
execution of every state amplitude or a native readout experiment.

Independent rational elimination for small controls gives traceless
kernel dimensions 16 for the free dimer, 12 for two uncoupled centered
sites, four for a chemically detuned centered dimer, and four for the
separate centered real-hopping triangle. A pair supported on just one
site fails conservation on the coupled dimer. These controls expose the
positive-contact, connectedness, centering and bipartite hypotheses.

The complete deterministic receipt is reconstructed in a fresh interpreter.
It checks exact values and every scope flag, distinguishes booleans from
integers, rejects duplicate JSON keys and pins all transitive repository
Python sources. This is finite Hamiltonian algebra; it selects neither the
source specification nor a physical continuum theory from observer closure.
