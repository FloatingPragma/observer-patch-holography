# Complete Fourier records and the relative-observation momentum

## 1. One observation state and a regular physical seed

Use the exact open Wilson/Gauss source, polynomial density, idle canonical
reference and Hamiltonian clock in the
[primary contract](../particle_common_contract/PROOF.md). The
[retained-origin construction](../particle_primary_origin_momentum/PROOF.md)
adds a classical origin X uniform on its bulk cube F_n and keeps that same
origin throughout a complete history. Its relative regular algebra has a
common subsequence state eta and commuting continuous actions alpha_t and
beta_y. In its canonical representation

```
U(t)=exp(i t H), V(y)=exp(i y.P), H>=0, H Omega=P Omega=0.
```

The finite box is not translation invariant. These generators belong to
the relative-observation domain. The quantum reference and X are retained
at each finite level; no limit of their growing full joint density is used.

The seed is the actual neutral dimension-six collective contraction
`C_x=sum_(alpha=1,2;i=1..6) c†_(+,alpha i,x)c_(-,alpha i,x)/sqrt(42)`.
Its full-Fock norm is at most one, and it is even and Gauss invariant.
Its native block and norm are authenticated by
[the refinement response](../particle_primary_refinement_response/PROOF.md).
Let phi_R be the product of three normalized triangular densities with
support [-R,R], and tau_T the normalized time triangle on [-T,T].
With exact spatial cell-integral weights define

```
A_n,x = integral tau_T(u) alpha_u^box(sum_z w_n,z(phi_R) C_(x+z)) du.
```

All source evolutions have the unchanged source energy origin. Every
source/reference gate is physical, with the reference idle. At levels
whose common origin collar cannot contain this experiment, set the probe
to zero and retain its failure. This affects finitely many levels for fixed
R,T and aperture. The receipt marks comparison with the fixed regular
target as inapplicable at these levels and supplies no error bound for
that comparison. The limiting regular contraction A obeys

```
||beta_y A-beta_z A|| <= (2/R) sum_i |y_i-z_i|,
||alpha_t A-alpha_s A|| <= 2|t-s|/T.                         (1)
```

These are the spatial L1 and triangular time-filter bounds in the owning
algebra proofs. No norm continuity of the unsmeared point seed is assumed.

## 2. An exact complete coherent-offset experiment

Take spatial half-aperture Q and temporal half-aperture S, positive
rationals with Q/a an integer. Offsets are the Cartesian grid
`y_j=-Q+a j`, `0<=j_i<N=2Q/a`. Outer times are
`t_l=-S+h l`, `0<=l<Nt`, `h=2S/Nt`. Define frequencies
`q_i=pi k_i/Q`, `nu=pi m/S` for fixed integers k_i,m. The phases at these
grid points are exact roots of unity. Negative frequencies have their
signed integers; a finite Fourier row represents that integer modulo its
grid size and retains its known phase convention.

The inner spatial triangle has positive rational cell weights summing to
one. Approximate the inner time triangle by its positive exact cell
integrals at M even midpoint nodes u_r. Their weights also sum to one.
Each combined selector address j comprises an outer offset/time and an
inner spatial/time address. Its positive weight p_j is the product of
these weights and `1/(N^3 Nt)`. At the retained origin x, its physical
success operator is

```
D_(0,j,x)=exp(-i(q.y_j+nu t_j)) alpha_(u_j+t_j)^box(C_(x+y_j/a+z_j)).
```

Complete it with the conjugated positive root of `I-C†C`, with the same
scalar phase, as D_(1,j,x). Thus each address has
`sum_b D_(b,j,x)†D_(b,j,x)=I`. A coherent selector isometry sends
`psi -> sum_j sqrt(p_j)|j> D_(b,j,x)psi`. After a unitary analyzer W
whose distinguished row is `(sqrt(p_j))_j`, its Kraus operators are

```
K_(b,r,x)=sum_j W_(r,j) sqrt(p_j) D_(b,j,x),
sum_(b,r) K_(b,r,x)†K_(b,r,x)=I.                            (2)
```

The distinguished success is exactly the weighted coherent sum Btilde_n,x.
An exact algebraic Householder analyzer supplies this row. A separate
choice uses the inner weight analyzer and outer discrete Fourier analyzers; their
complete rows retain all spatial/time Fourier labels and inner failures.
Multiplying a row by its known scalar phase gives the stated signed
frequency convention. No source translation unitary is applied.
These analyzer choices share the distinguished success operator; their
other outcome maps need not agree. The source-bound plan uses the inner
weight analyzer and the Cartesian outer Fourier analyzer. The bin is
`(k_1 mod N,k_2 mod N,k_3 mod N,m mod Nt)`; its known row phase is
`(-1)^(k_1+k_2+k_3+m)` because all outer nodes start at the left aperture
endpoint. Small controls execute both analyzer constructions. All Fourier
rows survive, including every other inner analyzer outcome and local
failure. Neither construction is an exact finite momentum measurement.

The finite output retains X, every selector/analyzer label, both local
probe branches, the source and the full reference. The original success
probability is

```
p_n = |F_n|^-1 sum_x Tr(rho_n Btilde_n,x† Btilde_n,x).       (3)
```

It is not divided by success, selector or aperture probabilities. A
classical mixture of offsets erases the cross terms and gives a different
experiment. The positive-root and controlled matrix services are declared
architecture choices. The large selector, source matrix and complete
histogram are specified, not executed. No efficient compiler cost is proved.

## 3. Actual finite errors and one common limit

Let D_n=U_n-L_n be the certified spectral width and K_n=max(|L_n|,|U_n|)
the source norm cap. The positive inner triangle quadrature error is
`D_n T/M`; this consumes the growing source bound and is not replaced by
a uniform stress or Hamiltonian norm. Approximate every source evolution
of duration at most T+S by the inherited rational Cayley service with r
subdivisions. Its operator error is at most

```
e_U = K_n^3 (T+S)^3/(12 r^2).                              (4)
```

Each address's complete column changes by at most 2e_U, because both
source evolutions in its conjugation are approximated. The selector and
analyzer are unitary, so the complete all-record output half-trace error
is at most 2e_U for every quantum reference, with no address-count or
origin-count multiplier. Any extra preparation/compiler error adds
separately. This compares the finite node instruments; it does not
identify their failure law with a continuous filter's failure law.

The continuous target in the relative algebra is

```
B_(q,nu)=1/((2Q)^3 2S) integral_[-Q,Q]^3 integral_[-S,S]
          exp(-i(q.y+nu t)) alpha_t beta_y(A) dt dy.         (5)
```

Its contraction norm is at most one. The finite outer grid approximation
has the explicit regularity bound

```
e_grid <= a(6/R + pi sum_i |k_i|/Q)
                 + h(2/T + pi |m|/S).                    (6)
```

At each regulator, nonlattice y means the shifted spatial test function
phi_R(. - y), not a fictitious finite spatial automorphism. Its exact
cell-integral weights have the L1 bound (1). Therefore (6) holds at finite
n, and its limiting operator is (5). It needs no finite translation
symmetry. The rational receipt uses pi<22/7.

Compare infinite and finite source times using the actual retained-origin
bulk bound with support radius R+Q and time extent T+S. Its operator error
e_bulk tends to zero along the declared source family. The collar is the
same for every address and every origin. Positive weights have L1 total
one, so the bulk comparison introduces no selector factor.

Thus the distinguished success operator comparison is bounded by

```
e_amp <= D_n T/M + 2e_U + e_grid + e_bulk.                  (7)
```

Original success probabilities differ by at most 2e_amp. Choose M,r,Nt
so their terms vanish; the code supplies finite integer schedules using
the actual D_n,K_n. One common subsequence is chosen for all rational
R,T,Q,S, integer frequencies, translates and finite words. The resulting
probabilities converge to eta(B†B) on the same state and source clock.
The state subsequence has no supplied effective convergence modulus:
(7) controls apparatus and operator comparisons, not the distance from
a finite state to its limiting state.

For a fixed complete-history interface, a separately declared binary
instrument `Btilde, sqrt(I-Btilde†Btilde)` can be used. Its failure channel
differs from the full selector instrument (2). If two success contractions
differ by delta, the binary columns differ by at most
`sqrt(delta^2+2delta)`. This follows from
`||B†B-C†C||<=2delta` and square-root Holder continuity. It bounds complete
outputs including the idle reference, and telescopes over fixed histories.
Continuous functional calculus puts the limiting root in the same regular
algebra. No limit of the growing fine selector register is claimed.

## 4. The joint spectral dictionary

Fix the Fourier convention
`fhat(p)=integral f(y)exp(+i y.p)dy` and the same positive sign for time.
Since both canonical unitary groups fix Omega, Bochner integration gives

```
pi(integral f(y) beta_y(A)dy) Omega = fhat(P) pi(A)Omega,
pi(B_(q,nu))Omega = sinc(S(H-nu))
                    product_i sinc(Q(P_i-q_i)) pi(A)Omega,             (8)
```

where sinc(z)=sin(z)/z and sinc(0)=1. The strong commutation of H and P
supplies their joint spectral measure. If dmu_A is that measure on the
vector pi(A)Omega, the original limiting success law is exactly

```
eta(B†B)=integral sinc^2(S(E-nu)) product_i sinc^2(Q(p_i-q_i)) dmu_A.
                                                                    (9)
```

The measure has mass at most one and support E>=0. Equations (2)-(9)
connect actual complete finite coherent-offset records to canonical
momentum kernels on the declared relative-observation domain. A finite
Fourier outcome is an aperture filter, not a projection onto an exact
finite-box momentum. No fixed-origin state invariance is asserted.

For positive radii dE,dp_i let G be the closed joint rectangle centered
at (nu,q). Outside G, at least one factor in (9) is at most

```
leak = min(1, max(1/(S dE)^2, max_i 1/(Q dp_i)^2)).          (10)
```

Hence `mu_A(G)>=max(0,eta(B†B)-leak)`. This keeps original mass, including
zero-response cases. A positive uniform lower bound on measured success
exceeding the leakage and vanishing finite errors would therefore give
nonzero joint spectral weight. Such a primary lower bound is not supplied
here. The kernels can also be treated as scaled approximate identities;
that data-analysis scaling is not success conditioning or a probability
readout. No sharp spectral line, dispersion, Lorentz law, particle species
or GeV mass follows solely from this dictionary.

## 5. Executed controls

Small exact controls check a non-Hermitian contraction's complete weighted
selector instrument, an entangled idle reference, and the difference from
an incoherent offset mixture. A commuting finite H/P control with an
invariant vacuum checks both signs in (8), the joint kernel and finite
quadrature convergence. A noninvariant input countercontrol fails the
vacuum-fixed identity. The source-bound receipts regenerate the actual
member and geometry, authenticate the native seed, compute finite
selector counts and all error schedules, and execute no large primary
matrix. The common limiting state may give A Omega=0 or trivial momentum.
