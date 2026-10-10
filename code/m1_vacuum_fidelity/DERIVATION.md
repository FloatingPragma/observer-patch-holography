# Global vacuum fidelity of local read programs

## 1. The reference state cannot be changed during comparison

Let B be a real symmetric strictly positive N by N matrix, with frequencies
omega_k=sqrt(lambda_k). On L2(R^N), set hbar=1 unless displayed and

    H = (P^T P + Q^T B Q)/2,       [Q_i,P_j]=i hbar delta_ij.

The reference is its normalized ground state Omega. Its evolution under H
is stationary up to phase. All excitation numbers, energies and overlaps below
are relative to this same state and Hamiltonian. We compare actual products
of kinetic and potential unitaries. Quadratic identities hold on Schwartz
space and their metaplectic unitaries and bounded Weyl effects are defined
on the full finite-population Hilbert space.

A potential kick of signed duration s sends P to P-s B Q, leaving Q fixed;
a kinetic drift sends Q to Q+s P, leaving P fixed. The symmetric second-order
step is kick(tau/2), drift(tau), kick(tau/2). In the original vacuum's
dimensionless mode coordinates X=sqrt(omega) Q, Y=P/sqrt(omega), its matrix is

    S2(z) = [[1-z^2/2, z], [-z+z^3/4, 1-z^2/2]],   z=tau omega.       (1)

Put rho=cuberoot(2), w=1/(2-rho), v=-rho/(2-rho). The fourth-order step is

    S4(z)=S2(w z) S2(v z) S2(w z).                                  (2)

This is the standard symmetric composition of Yoshida, not a new integrator.
The signed middle step is retained. Inverse Hamiltonian gates are declared
operations, not negative elapsed transport time. Both matrices have the form
S=[[a,b],[c,a]] and determinant one. Direct multiplication in Q(rho) gives

    a4 = 1-z^2/2+z^4/24+A z^6,
    b4 = z-z^3/6-B0 z^5,
    c4 = -z+z^3/6+C z^5-D z^7,
    A = 1/48+5rho/288+rho^2/72,
    B0 = 1/36+rho/36+rho^2/48,
    C = (rho+rho^2)/144,
    D = 25/1728+5rho/432+rho^2/108.                                    (3)

In particular b4+c4=d4 z^5-D z^7, where
d4=-1/36-rho/48-rho^2/72 !=0. For (1), b2+c2=z^3/4 and d2=1/4.

For 0<z<=1/10, both steps satisfy b>0, c<0, 0<a<1 and
0.9 z^2 <= -bc <= 1.1 z^2. For S4,
z^5/16 <= |b+c| <= z^5/8; for S2 the defect is exactly z^3/4.
These bounds follow by substituting 1259/1000<rho<126/100 into (3),
bounding the remaining powers by z^2<=1/100. The exact certificate performs
these rational coefficient bounds. This domain is conservative, not a fitted
stability limit. Write cos(theta)=a, sin(theta)=sqrt(-bc)>0 and
chi=sqrt(-c/b). Then for every integer j>=0,

    S^j = [[cos(j theta), sin(j theta)/chi],
           [-chi sin(j theta), cos(j theta)]].                    (4)

## 2. Exact original-vacuum error

The initial covariance of (X,Y) is hbar I/2. After j steps it is
hbar S^j(S^j)^T/2. Hence the mean original particle number in that mode is

    n_j(z) = (chi-chi^-1)^2 sin^2(j theta)/4
           = (b+c)^2 sin^2(j theta)/(-4bc).                         (5)

At one step n_1=(b+c)^2/4. Thus S2 produces exactly z^6/64 particles in
one mode after one step, although its matrix is stable and symplectic.
Stability does not imply stationarity of the original vacuum. The entire
population has

    N_j = sum_k n_j(tau omega_k),
    Delta E_j = hbar sum_k omega_k n_j(tau omega_k),
    F_j = |<Omega,U^j Omega>|^2 = product_k (1+n_j)^(-1/2).          (6)

The overlap formula follows by integrating the normalized one-mode squeezed
Gaussian; diagonalizing the real matrix B gives N independent real modes.
Fourier degeneracies count with their full real multiplicities, not half
of the complex Fourier catalog. No ultraviolet zero-point subtraction is
being adjusted: Delta E is the nonnegative energy above the original vacuum.

### Evaluating the original-vacuum signal without cancellation

Cayley--Hamilton gives `S^2-2aS+I=0`, hence
`S^j=U_(j-1)(a)S-U_(j-2)(a)I`, where `U_-1=0`, `U_0=1` and
`U_(k+1)=2a U_k-U_(k-1)`. The diagonal entries of every power agree,
so its original-vacuum particle number is exactly

    n_j = (b+c)^2 U_(j-1)(a)^2/4.

With `s=lambda*tau^2`, equations (1)--(3) therefore become

    order 2: n_j = s^3 U_(j-1)(1-s/2)^2/64,
    order 4: n_j = s^5 (D0+D*s)^2 U_(j-1)(a4)^2/4,
    D0 = 1/36+rho/48+rho^2/72 = -d4.

These are identities, not asymptotic replacements. Forming the positive
defect from these coefficients avoids subtracting the nearly opposite
off-diagonal entries. It also exposes the separate conditioning problem at
a zero of `U_(j-1)`: a factored defect alone does not resolve a nearby revival.
The numerical evaluator encloses the recurrence and refines precision until
the returned observables have unambiguous 23-digit formatting; it refuses
unresolved results rather than replacing them with zero. The checker
reconstructs the original gate evolution independently. Both paths retain
the original exact `s` for the domain decision, all multiplicities and all
32 times. See the [numerical contract](README.md#numerical-contract).

For p=2 or 4, uniformly on the stated domain,

    c_p z^(2p) <= (b+c)^2/(-4bc) <= C_p z^(2p).                    (7)

The positive lower coefficient matters. A formal order bound alone would
not establish a necessary precision threshold. As z->0, that amplitude is
d_p^2 z^(2p)/4 times (1+O(z^2)).

Replacing Omega by the invariant squeezed vacuum does not fix this comparison.
That vacuum already has original-mode excess energy
hbar omega (chi+chi^-1-2)/4. Its change of reference must be declared if used.

## 3. A full-spectrum threshold, not a selected-mode estimate

Consider a sequence of these matrices. Let sigma^2=N^-1 Tr B, suppose
sigma->infinity and lambda_max(B)<=2 sigma^2. The concrete graph families
below satisfy this bound exactly. At least N/3 eigenvalues exceed sigma^2/2:
otherwise their mean would be strictly less than sigma^2 even if all the
remaining eigenvalues attained 2 sigma^2. Consequently, for each fixed s>0,

    (N/3)(sigma/sqrt(2))^s <= sum omega_k^s
                                      <= N(sqrt(2) sigma)^s.       (8)

Assume tau sqrt(lambda_max)<=1/10 and tau->0. Fix T>0 and
J=floor(T/tau). The average is over every j=1,...,J. The geometric sum gives

    | J^-1 sum_j sin^2(j theta)-1/2 | <= 1/(2J sin(theta)).          (9)

For the bulk modes in (8), sin(theta)>=sqrt(0.9) tau sigma/sqrt(2),
so the error in (9) is O(1/(T sigma)). It tends to zero. Applying (7)
to all modes for the upper bound and this bulk for the lower bound proves

    mean_j N_j           = Theta(N tau^(2p) sigma^(2p)),
    mean_j Delta E_j     = Theta(hbar N tau^(2p) sigma^(2p+1)),
    mean_j [-log F_j]    = Theta(N tau^(2p) sigma^(2p)).             (10)

The matching upper bounds hold for every j, even without restricting time.
For the last estimate use (6) and x/(1+x)<=log(1+x)<=x, with the uniform
bound on n_j provided by (7). All constants are independent of population.

It follows that, on every fixed nonzero time interval,

    sup_j (1-F_j)->0  iff  N tau^(2p) sigma^(2p)->0,
    sup_j Delta E_j->0 iff N tau^(2p) sigma^(2p+1)->0.              (11)

For the first necessity, sup_j(1-F_j)->0 implies sup_j(-log F_j)->0.
If either scale diverges, its corresponding error is large on a positive
fraction of the sampled times: its mean has a lower bound and its maximum
has the same-scale upper bound. In the fidelity case that fraction has
F_j bounded above by exp(-c N tau^(2p) sigma^(2p)). This is not a claim
about every predetermined isolated time; revivals are allowed.

When tau sigma->0, the weighted form of (9) additionally yields

    mean Delta E = (hbar d_p^2/8) tau^(2p)
                     sum omega_k^(2p+1) (1+o(1)),                 (12)

and the analogous formula for mean particle number. The phase-average error
is bounded by C tau^(2p) sum omega^(2p)/(T), which is smaller than the energy
main term by O(1/(T sigma)) using (8). The Taylor remainder is O((tau sigma)^2).

These thresholds distinguish convergence in state norm from convergence of
an unbounded energy observable. For example, tau chosen so that the energy
scale tends to a positive constant makes the number scale tend to zero:
the states approach the vacuum in norm while retaining excess energy in
rarer and increasingly energetic ultraviolet excitations.

## 4. A global bounded-effect comparison

For the vacuum, the trace distance of the two pure states is sqrt(1-F_j).
Thus (11) bounds every bounded effect 0<=E<=I, including collective effects
that act on all N modes. No restriction to smooth detectors is used here.

Let a coherent input have original excess energy at most E0, supported on
frequencies at most a fixed Omega0, with omega_k>=mu>0. In the normalized
(X,Y)/sqrt(hbar) coordinates its displacement d obeys
||d||^2<=2E0/(hbar mu).
For tau Omega0<=1/10, equations (1)--(4) give on that band and 0<=j tau<=T,

    ||S_p^j-R(j tau omega)|| <= tau^p Omega0^p(1+T Omega0).          (13)

where R is the exact oscillator rotation and the norm is the Euclidean
operator norm. Here are explicit constants. The sign and coefficient bounds
in section 1 imply |chi-1|<=z^p and |chi^-1-1|<=z^p. For S2,
|a-cos(z)|<=z^4/24; for S4 the alternating cosine remainder and A<1/10
give |a-cos(z)|<0.102 z^6. Both sin(theta) and sin(z) exceed 0.94 z,
and both angles lie in (0,pi/2). The mean value theorem therefore gives
|theta-z|<=z^(p+1). The difference of (4) from R(j theta) is off-diagonal
with operator norm at most z^p, while two rotations differ by at most their
angle difference. This proves (13) with coefficient one.
Weyl covariance and the triangle inequality for trace distance
now give a bound, uniform over every bounded effect and these coherent inputs,

    probability error <= sqrt(1-F_j)
          + tau^p Omega0^p(1+T Omega0) sqrt(E0/(hbar mu)).           (14)

The squared coherent overlap exp(-||delta d||^2/2) bounds its trace distance
by ||delta d||/sqrt(2), supplying exactly the second term.
This statement
does not assert an unconstrained diamond-norm bound for arbitrary input
states or identify finite matrix registers with oscillator Hilbert spaces.

## 5. The actual ultraviolet scale of the repaired sparse reads

Use a unit periodic cube, N=q^3 sites, c=mu=hbar=1, and a symmetric nonzero
integer stencil V without wrap-identifications. Set

    d=|V|, M2=sum_(v in V)|v|^2, alpha=6q^2/M2,
    A f(x)=alpha sum_v [f(x)-f(x+v)], B=I+A,
    W=alpha d, sigma^2=1+W.                                      (15)

The positive graph Laplacian has trace N W and norm at most 2W. Thus
lambda_max(B)<=1+2W<=2 sigma^2, supplying (8) for these actual matrices.
It is the total outgoing action weight W, not the longest edge alone,
that controls the bulk ultraviolet spectrum.

For integer t>=1 define

    q=2^(16t), m=2^(8t), K=2^(6t), r=2^(7t), R=mK,
    U=(m Z^3 intersect closed B_R) union (Z^3 intersect closed B_r)
             union {+/-2^j e_i: 0<=j<8t}, with zero removed.       (16)

These sets define every read, independently of any test outcome. The two
balls meet only at zero since r<m. The dyadic points not already in the
short ball number O(t). Ordinary cube Riemann sums on a ball give

    d ~ (4pi/3)r^3 = Theta(q^(21/16)),
    M2 ~ (4pi/5)m^2 K^5 = Theta(q^(23/8)),
    W ~ 10 q^(7/16), sigma = Theta(q^(7/32)), a=R/q=q^(-1/8).      (17)

The short ball dominates the degree; the spaced ball dominates the second
moment. In particular sigma is not Theta(1/a). Losing this distinction
would undercount both stability and vacuum-fidelity work.

For this family, (11) requires

                 vacuum fidelity              energy fidelity
    S2           tau=o(q^(-31/32))             tau=o(q^(-131/128))
    S4           tau=o(q^(-19/32))             tau=o(q^(-159/256)). (18)

Here is an explicit repair with no success filter:

    use S4 and tau=(1/64) q^(-5/8).                               (19)

Then sup_j Delta E_j=O(q^(-1/32)) and sup_j(1-F_j)=O(q^(-1/4)).
Together with (14), this gives the full bounded-effect comparison for the
declared preparations over each fixed time interval.

For completeness, (19) lies in the proved stability domain at every t>=1.
The integer ball has at most 27s^3 points for s>=1. Retain the box
K/4<=v_1<=K/2 and |v_2|,|v_3|<=K/4 inside B_K, for K divisible by four.
It has at least K^3/16 points with |v|^2>=K^2/16, proving the conservative
bound M2>=m^2 K^5/256. Also d<=60r^3 for t>=1.
Thus W<=92160 q^(7/16), and tau sqrt(1+2W)
<=sqrt(184321)/64 * q^(-13/32)<1/10 at q>=2^16.

The same tick on S2 has mean energy Theta(q^(51/32)), so a stable second-order
implementation fails decisively. The fourth-order cancellation is substantive.

## 6. Count every operation, and keep the clock honest

For J repeated S2 steps, merging adjacent potential kicks leaves J+1 kicks
and J drifts. For S4, the explicit seven-stage word has coefficients

    K(w/2), D(w), K((w+v)/2), D(v),
                K((v+w)/2), D(w), K(w/2).

Repeating it and merging only adjacent kicks leaves 3J+1 kicks and 3J drifts.
Each reference kick reads the actual current Q at its d spatial neighbors
and locally; each reference drift reads the current local P. The ledger
retains writer identifiers and classical reference values. It does not copy
unknown quantum fields. The quantum potential unitary factors into Nd/2
commuting two-mode phases exp[-i s alpha (Q_i-Q_j)^2/(2 hbar)] and N onsite
mass phases; the drift is N onsite kinetic unitaries. These are explicit
unitary implementations of the same canonical map. Thus the reference
nonlocal ordered incidence count is exactly
(J+1)Nd or (3J+1)Nd for these specified programs. All negative coefficients,
intermediate writes and inverse replays are retained in the finite evidence.
This is an ideal arithmetic/read count, not a bit, energy or gate-synthesis
complexity claim.

At fixed T and vanishing energy tolerance epsilon, (10)--(11) give matching
orders for the allowed step and resulting incidence count:

    tau = Theta((epsilon/(N sigma^(2p+1)))^(1/(2p))),
    reads = Theta(T d N^(1+1/(2p)) sigma^(1+1/(2p)) epsilon^(-1/(2p))). (20)

The necessary side concerns uniform error on the full interval, using the
time average; it is not an isolated-time lower bound or a lower bound for
every possible algorithm. For the dense radius q^(-1/2) stencil,
d=Theta(q^(3/2)), sigma=Theta(q^(1/2)). At the same q and energy tolerance:

    program       sparse U exponent       dense exponent      epsilon factor
    S2            683/128                  47/8                epsilon^(-1/4)
    S4            1263/256                 87/16               epsilon^(-1/8).

The explicit schedule (19) uses Theta(q^(79/16)) incidences and achieves its
stated energy decay. This is below the dense unintegrated q^5 incidence
benchmark, but that benchmark is not an equal-fidelity simulation. Formula
(20) is the matched comparison. Sparse and dense spatial errors also differ;
no equal-total-continuum-accuracy claim follows from equal q.

Gate locality is exact: each kick expands canonical support by at most one
stencil edge and each drift is onsite. After J macros the support bound is
(J+1)a for S2 or (3J+1)a for S4. The signed coefficients do not make time
negative. These program durations are simulation parameters: a/ tau diverges
in (19). Native transport at a fixed physical speed has not been supplied.
Renaming the extra work as faster physical reads would violate the clock
accounting. This result repairs global quantum approximation and its read
cost; it does not repair that distinct native timing identification.

### The support bound is attained, so the clock mismatch is an obstruction

For the stated U and the dense ball, v=(R,0,0) is the unique displacement
with first coordinate R. Its coefficient in B is -alpha, alpha>0. Let j>=1
and exclude periodic wrapping by requiring (3j+1)R<q/2; the S2 version only
needs (j+1)R<q/2. The coefficient from initial Q to final P at displacement
(j+1)R e1 for S2, or (3j+1)R e1 for S4, is exactly

    S2:  (-1)^(j-1) tau^(2j+1) (-alpha)^(j+1) / 4,
    S4:  -D (2A)^(j-1) tau^(6j+1) (-alpha)^(3j+1).                 (22)

Both are nonzero; A and D are the positive coefficients in (3).
Cayley--Hamilton gives the off-diagonal entry of S^j as U_(j-1)(a) times
that of S, where the Chebyshev polynomial U_(j-1) has leading coefficient
2^(j-1). In unscaled (Q,P) coordinates, the top lower-left term of S2 is
tau^3 B^2/4 and the top diagonal term is -tau^2 B/2. For S4 those terms
are -D tau^7 B^4 and A tau^6 B^3. Multiplication gives (22). In B^d, the
extreme displacement d R e1 can only be reached by taking R e1 at every
factor; all lower powers fail to reach it. Its coefficient is (-alpha)^d.
The non-wrapping condition prevents an alternative image of the same site
from contributing. The position support radii ja and 3ja are attained by
the leading-diagonal argument. Thus the canonical support radii (j+1)a and
(3j+1)a are exact, not only upper estimates.

These are operational signals. A local Q displacement generated by P_y
changes the final P_x mean by the nonzero coefficient times its displacement.
Against a finite Gaussian baseline, the bounded effect (I+sin(eta P_x))/2
has a nonzero response for a suitable finite displacement and eta. The
corresponding canonical commutator is nonzero. This uses the declared full
oscillator operation family; the classical reference program also
distinguishes two initial scalar values at that writer.

Any faithful realization subject to a complete influence-speed bound c must
charge elapsed time at least (3j+1)a/c for S4 or (j+1)a/c for S2 in this
non-wrapping regime. It cannot assign elapsed time j tau when that is smaller.
The exact native cone bound proved from LC1--LC3 in
`code/rg_principle/DERIVATION.md` applies to every admitted influence, so
this particular simulation-clock identification is excluded.

On the repaired schedule take j=floor(1/(8a)). At sufficiently large levels
it is positive and (3j+1)a<1/2. The attained distance tends to 3/8, while
j tau=O(q^(-1/2))->0. No fixed finite c can identify these simulation
durations with a complete native influence clock. Small amplitude does not
make an exact coefficient zero. Conversely, the global fidelity theorem
concerns specified preparations and bounded errors; it does not promote
every ultraviolet local intervention to a fixed-band input. The results are
compatible. The repair closes global quantum approximation, and (22) proves
the finite timing attachment it cannot supply. Other operation grammars are
outside this obstruction.

Six independent 32-site ray executions check the entire support and the
nonzero extremal coefficient for both methods at j=1,2,3. They use B=I plus
the nearest-neighbor line Laplacian, alpha=1 and tau=1/64. These are ray
controls for the algebraic argument, not executions of a large 3D source.

## 7. Why a fixed local scalar split cannot preserve the vacuum exactly

Consider a scalar update with real polynomial entries a(B),b(B),c(B),d(B),
of bounded degree, acting on (Q,P). Suppose it is symplectic and preserves
the original positive energy, equivalently its pure original Gaussian vacuum.
For a scalar eigenvalue lambda>0, these two identities imply d=a, c=-lambda b
and

    a(lambda)^2+lambda b(lambda)^2=1.                              (21)

This follows by equating the symplectic inverse [[d,-b],[-c,a]] to the
energy inverse diag(lambda,1)^(-1) S^T diag(lambda,1). If the number of
distinct eigenvalues exceeds the degree of the polynomial in (21), it is
a polynomial identity. But the degrees of a^2 and lambda b^2 have opposite
parity, so their top nonzero coefficients cannot cancel. Hence b=0 and
a=+1 or -1. A near-identity consistent step is therefore the identity.

A fixed number of kicks and drifts has bounded polynomial degree and cannot
give nontrivial exact preservation on all refinements. The families (15)--(16)
have an unbounded number of distinct eigenvalues: along Fourier indices
(k,0,0), 0<=k<=floor(q/(4R)), every nonzero x-displacement contributes a
strictly increasing 1-cos(2pi k v_1/q), since its angle lies in [0,pi/2].
At least one such displacement occurs. Thus there are at least
floor(q/(4R))+1 distinct values, and q/R->infinity. The dense case has the
same property. This proves the cardinality premise for the stated families.

This obstruction is for scalar polynomial action updates, including fixed
kick/drift compositions. It does not prohibit nonpolynomial exact flows,
population-dependent large-degree interpolation, spatial translations,
extra internal fields, or different operation grammars. The fourth-order
repair supplies asymptotically exact fidelity without claiming exact finite
vacuum preservation or evading the obstruction.

## 8. Relation to the existing corpus and imported mathematics

The existing finite golden scalar execution and its original-vacuum bounded
probability theorem are in `paper/tex_fragments/SOURCE_SCALAR_QUANTUM.tex`.
They correctly distinguish the squeezed split baseline from the stationary
reference. The present result studies a growing entire population, derives
the necessary global thresholds, and constructs a sparse repair. It does
not retract those finite probabilities or substitute a different vacuum.

The higher-order composition is due to H. Yoshida, *Construction of higher
order symplectic integrators*, Physics Letters A 150 (1990), 262--268,
https://doi.org/10.1016/0375-9601(90)90092-3;
an accessible original paper is
https://tlakoba.w3.uvm.edu/math337/for_final_topics/SSM_1990_Yoshida.pdf.
Gaussian oscillator algebra, finite Fourier diagonalization, polynomial root
counting and ball Riemann sums are standard inputs. The specific contribution
is (10)--(20) for these read programs and populations, with independent
all-mode and actual-gate evidence. Universal asymptotic claims are analytic;
the finite evidence and Lean algebra reductions have their narrower scopes.
