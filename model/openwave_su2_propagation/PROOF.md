# Propagation on the full-Haar Yang-Mills lattice

Use the same pure SU(2) lattice Hamiltonian
H=alpha sum_e C_e^2+beta sum_p(1-W_p), with
alpha=g^2/(2a), beta=2/(g^2 a), and all Haar link representations.
The additive magnetic scalar does not affect Heisenberg dynamics. Regard
links as sites: the Casimir is an unbounded on-site term and each remaining
interaction is -beta W_p, supported on four links with norm beta.

Each link belongs to four plaquettes. A plaquette can overlap at most sixteen
plaquettes, including itself in this deliberately loose bound. A chain of n
interactions starting at a finite link support X therefore has at most
4|X|16^(n-1) choices. Give the link graph an edge when two links belong to one
plaquette. A chain connecting disjoint supports at distance d needs at least
d interactions. The iterated commutator inequality in the proof of
[Nachtergaele and Sims, Theorem 3.1, equations 71-73](https://arxiv.org/abs/1410.8174)
then gives, for unit-norm observables,

```
||[tau_t(A),B]|| <= min(2, 2|X| sum_(n>=d) (32 beta |t|)^n/n!).
```

The unbounded electric evolution changes neither support nor the norm of a
bounded plaquette in its interaction picture. The cited argument applies to
separable infinite local spaces with strong on-site evolution. No finite
Peter-Weyl truncation or bounded electric-operator norm is used. The source
commutes with local gauge transformations, so the bound also holds on physical
observables and in their gauge-invariant representations.

For a local observable within a fixed central unit box, enclose it by a cube
with R interaction layers of buffer. Crossing plaquettes can be bounded by
12(2R+3)^3; a volume overcount suffices. Duhamel comparison of a larger box
with the decoupled interior/exterior source, followed by the preceding bound,
gives

```
||tau_t^large(A)-tau_t^interior(A)||
 <= min(2,24 |X| beta |t| (2R+3)^3
               sum_(n>=R) (32 beta |t|)^n/n!).
```

The factorial tail dominates the polynomial prefactor as R tends to infinity,
uniformly on bounded time intervals at fixed a and g. This proves norm
convergence of each local observable's finite-volume dynamics. It does not
claim norm continuity in time on the entire algebra of bounded local
operators with unbounded on-site generators. The resulting limit retains
local gauge invariance. Integer incidence and overlap counts on complete
periodic cubes independently check the constants; the proof of the general
counts follows from the four coordinate plaquettes incident on an oriented
link. Rational Taylor sums with geometric remainder enclose every printed
propagation and boundary error.

The strong-coupling ground-state theorem supplies its own existential
uniform gap and nonzero Wilson variance. In that domain the physical spectral
gap gives connected Euclidean-time correlations bounded by the corresponding
variances times exp(-Delta t). The propagation theorem supplies an explicit
volume error at every positive coupling even where that gap theorem does not
numerically admit the coupling. Those conclusions are distinct.

At fixed g the entire Hamiltonian scales exactly as 1/a. Any nonzero fixed
dimensionless excitation energy therefore diverges as a tends to zero with
the declared clock. The positive bound from the uniform strong-coupling
theorem has the same scaling. A finite physical mass requires a changing
coupling or clock/energy-unit trajectory for which the dimensionless gap scales
as a, together with the appropriate Euclidean or Lorentz continuum
properties. A volume limit and a time continuum do not establish that spatial
limit. A scalar reference-energy subtraction cannot change a gap; a different
vacuum sector needs a separate argument. No explicit strong-coupling threshold,
GeV glueball mass or Clay result is claimed.

## Compact gauge groups

The propagation estimate uses only the plaquette support and the bound
|W|<=1. It therefore holds for every compact connected Lie gauge group with
a bi-invariant Casimir and a bounded real normalized Wilson class function.
A finite-dimensional unitary representation supplies such a function by
W=Re Tr(U_loop)/dim. The full Haar space is separable, and the compact
Laplacian is self-adjoint. Its constant zero mode is unique for a connected
group, with a positive first Casimir eigenvalue c_G.

The existing strong-coupling product-state theorem also transfers directly.
Set h_x=(1/c_G) sum_(y in x+star) sum_i C_(y,i)^2. Its local gap is one,
and each link is counted four times. Dividing by alpha c_G/4 leaves local
perturbation norm at most12 beta/(alpha c_G)=48/(c_G (g^2)^2). These are the
same infinite-local-space stability premises. The constant Haar product
vacuum is the unique classical ground. Compact semisimple groups have no
nontrivial continuous one-dimensional characters, so the unique perturbed
ground is gauge invariant, as required for the physical restriction. For a
general connected compact group this last character argument requires its
own invariance check; the propagation theorem itself has no such restriction.

For SU(N) with the conventional trace normalization, the positive entries of
the inverse Cartan matrix make C2(lambda)=(lambda,lambda+2rho)/2 increasing
in every nonnegative Dynkin coordinate. A nonzero integral weight dominates
a fundamental weight. The fundamental values are
k(N-k)(N+1)/(2N), minimized at k=1 or N-1, hence
c_G=(N^2-1)/(2N). The exact program checks the complete finite list for
N=2,...,8; the quadratic k(N-k) proves the all-N minimum. Fundamental Haar
character orthogonality gives Var(W)=1/4 for SU(2) and1/(2N^2) for N>=3,
where the fundamental representation is not self-conjugate. Continuity of
local expectations in the same stability theorem preserves a positive
variance in an existential small-perturbation interval. This supplies
nonzero physical Wilson excitations for each SU(N), at fixed spacing.

For other compact connected semisimple groups a nonconstant Wilson class
function has strictly positive Haar variance; its numerical value and c_G
can remain group-specific inputs to the same proof. No numerical coupling
threshold or spatial-continuum gap follows from this group extension.
