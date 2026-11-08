# Aligned holes, source dispersion and position

The signed Hubbard action has an exact all-volume free fermion sector.
Its holes have a complete native position detector, a fixed translation
algebra and a controlled nonrelativistic continuum. The sector, centered
reference and scaling are specified inputs. They are not consequences of
the unread repair programme or an identification with measured particles.

## 1. The action and the repair use one chart

Use the same declared composition as `DEFECT_BUDGET.md`: the two-flavour
full-domain repair and the signed Hubbard Hamiltonian act on the same
site-outer occupation basis. Write its fermions as d_(x,s), and set

    H = sum_(x,j,s) (c/2) eta_j(x)
          [d_(x,s)* d_(x+ej,s) + d_(x+ej,s)* d_(x,s)]
        + U sum_x n_(x,up)n_(x,down),
    K = H - (U/2)N,
    eta_1(x)=(-1)^(x2+x3), eta_2(x)=(-1)^x3, eta_3(x)=1.  (1)

These are the actual `field_symmetry.symmetry.link_sign` values, not
chosen hopping signs. Every elementary plaquette has product -1. The
readable `full_fock_field.FieldGeometry` defaults are c=1/8 and U=1.
The action, coefficient bundle and coherent controls are declared source
data; this calculation does not derive their selection from repair.

The original source Weyl links are i c sigma_j/2. Its exact onsite frame
R_x=(i sigma_1)^x1(i sigma_2)^x2(i sigma_3)^x3 gives d_x=R_x c_x and
transports that full action to (1). A transport of this entire experiment
also conjugates every repair Kraus operator, preparation and detector.
The bare original-flavour repair need not preserve hidden alignment.
For the two-site axis-1 source, d-up at site 1 is c-down. Its original
mask 8 is moved by E, coin 1, raw 6 to original mask 2, which is d-down
at site 0. Its aligned output weight is zero. `frame_audit()` verifies
both the full Hamiltonian identity and this precise repair countercontrol.

## 2. Exact hole Fock space on the full source

Fix all d-down modes empty, allowing every d-up occupation. This whole
2^V-dimensional subspace is invariant under H, K and the same-chart
repair. Contact vanishes identically. Let Omega_F be the filled d-up
vector and define

    f_x* = d_(x,up),    f_x = d_(x,up)*.

These satisfy the canonical anticommutation relations exactly, with
f_x Omega_F=0. Their q-hole states are the qth exterior power of the
site Hilbert space, including their CAR signs. If h is the real symmetric
one-body hopping matrix of (1), tr h=0 and

    (K-K(Omega_F))|aligned = dGamma(U/2 - h^T),
    K(Omega_F)=-UV/2.                                     (2)

Indeed d_x* d_y=delta_xy-f_y* f_x. This fixes the transpose and the minus
sign; guessing a particle hopping sign from hole motion gives the wrong
answer. The producer compares every column of the full fixed-number
source operator with (2), including q=0 and q=V. It checks the contact
and all off-sector entries, rather than projecting before applying H.

For a finite one-hole propagator u(t)=exp[-it(U/2-h^T)], the q-hole
propagator is wedge^q u(t). Its amplitude between ordered configurations
X and Y is det[u(t)_(Y_i,X_j)]. Exchange of two labels changes the sign.
Number and CAR statistics are attached to this actual sector; no
Lorentzian spin or topological particle type is inferred. The interaction
U has no scattering effect inside this sector. Its many-hole scattering
relative to its own free dynamics is the free comparison.

## 3. Open boxes and a reference gap

Let B_L be half the adjacency matrix of an open chain of length L and
Z_L=diag((-1)^x). Then Z_L B_L=-B_L Z_L and, on the actual rectangular box,

    h/c = B_L1 tensor Z_L2 tensor Z_L3
           + I tensor B_L2 tensor Z_L3 + I tensor I tensor B_L3,
    h^2/c^2 = B_L1^2 tensor I tensor I
              + I tensor B_L2^2 tensor I + I tensor I tensor B_L3^2. (3)

The three summands of h anticommute. The chain sine vectors have
eigenvalues cos(pi r/(L+1)). The largest eigenvalue of (3), together
with the bipartite anticommutation of h, therefore gives the exact norm

    ||h|| = c sqrt(sum_j cos^2(pi/(L_j+1))).

The lowest one-hole energy relative to Omega_F on the centered clock is

    Delta_L = U/2 - c sqrt(sum_j cos^2(pi/(L_j+1))),
    Delta_infinity = U/2 - sqrt(3)c.                        (4)

If Delta_infinity>0, Omega_F is the ground within this aligned Fock
space; every q-hole energy is at least q Delta_infinity. This is not a
claim about the ground of the full interacting spinful source. The exact
finite-size control is

    0 <= Delta_L-Delta_infinity
       <= (c/sqrt(3)) sum_j pi^2/(L_j+1)^2.                 (5)

It follows by writing the difference of the two square roots with its
conjugate denominator. The executable enclosure uses rational outward
bounds for sqrt(3) and pi<22/7. At the readable defaults the bulk gap is
(4-sqrt(3))/8, approximately 0.28349365 in the declared source time units.
This is a reference-dependent sector gap. The physical mass field in the
certificate is null.

## 4. Operational position and source translations

The hole position is the vacancy projector, not a repair current or an
assigned continuum label. On q=1, X_j=sum_x x_j |x><x|. For general q
the position record is the full unordered vacancy configuration.

A complete E scan along a rooted spanning tree reads this configuration
in every fixed-q sector. If n_a,n_b are the binary particle occupations,
raw 0 says n_a-n_b=0. An odd E record determines both endpoint bits:
with current j and recorded coin b, n_a=j+b and n_b=1-j-b.
Subtract the recorded cumulative currents to express every such
difference in the initial occupations. Along a tree they determine all
initial occupations up to one common constant. The declared total V-q
fixes that constant. This covers the all-empty and all-filled sectors.

Each fine record branch is consequently

    M_w = 2^(-(V-1)/2) s_w |terminal_w><initial_w|,

or zero, with the exact signed CAR transport s_w. Every initial
configuration has all 2^(V-1) coin histories. The configuration effects
sum to its exact position projection-valued measure. Arbitrary reference
columns retain their correct original probabilities and conditional
back-action. There is no success normalization. The allocated alphabet
is 32^(V-1); absent sparse dictionary entries denote zero maps. The
detector writes the full archive. It is a destructive position
measurement, and it does not assert that its committed classical record
retains coherent superpositions of different positions.

On the infinite signed lattice define the actual magnetic translations

    (U_j psi)(x)=zeta_j(x) psi(x-e_j),
    zeta_j(x)=(-1)^(sum_(i<j) x_i).

Direct substitution in (1) gives

    [h,U_j]=0,       U_j* X_i U_j=X_i+delta_ij,
    U_i U_j=-U_j U_i (i!=j),
    U_j^2 psi(x)=psi(x-2e_j).                              (6)

The three two-step translations commute and provide a joint
quasimomentum. Ordinary one-step translations cannot all replace U_j:
they change the actual pi-flux coefficients. On open boxes shifts have
boundary losses and are partial isometries. The executable test checks
(6) only on genuinely padded source columns and retains that boundary.
No finite-dimensional canonical commutator is asserted.

Use an eight-parity cell x=2n+s and the physical Bloch convention

    (F psi)_s(k)=sum_n exp[-ik.(2n+s)] psi(2n+s),
    k in [-pi/2,pi/2)^3,        measure dk/pi^3.

Then X=i partial_k on its quasiperiodic domain and U_j^2 multiplies by
exp(-2ik_j). The source matrix and hole energies are

    h(k)=c sum_j cos(k_j) Gamma_j,
    Gamma_j^2=I, {Gamma_i,Gamma_j}=0 (i!=j),
    epsilon_+/- (k)=U/2 +/- c sqrt(sum_j cos^2 k_j).        (7)

The producer constructs Gamma directly from eta and the parity flips.
All eight cell components are retained. Each sign has rank four away
from a crossing. Cell folding and this multiplicity are not a physical
species count. The conical band crossings occur at energy U/2 relative
to Omega_F, not at its gapped hole vacuum threshold.

Finite coherent momentum or band analyzers require the declared
coherent source services. They are not supplied by the unread E scan.
Any finite computable one-body unitary has a Givens decomposition with
at most V(V-1)/2 two-mode rotations and V phases. The existing
`gaussian_source.local_word` provides the complex hopping rotations and
onsite phases, with every nonadjacent mode routed and restored by FSWAP.
Its full-Fock lift acts on all q, so the analyzer followed by the exact
configuration scan realizes the corresponding complete observable.
If the complete local dilations have norm errors delta_l, their sum
bounds the ideal-versus-native full-dilation error, including arbitrary
references and all export complements. A unitary or isometric target
then has trace error at most 2 sum_l delta_l; if numerical outputs are
not presumed exactly isometric, use 2eta+eta^2 instead.

`native_dimer_analyzer()` executes 78 actual source pulses on all sixteen
input columns, retains 288 first-export failure slots and all 32 E
records, and distinguishes the two signed dimer hole energy modes. Its
floating residual is a numerical admission diagnostic. The analytic
compiler error implication is separate from a numerical interval or
hardware certificate. The finite Fourier analyzers are accessible
supplied controls, not autonomous selections by the repair schedule.

## 5. Complete repair histories for every hole number

For the aligned q-hole space let C_e be the projection that exactly one
endpoint of seam e is a hole. Let V implement one uniform global repair
attempt, with E-law probability theta, including seam, law, coin and all
raw records. Let V0 write the same seam/law/coin with raw zero and act
as identity on the source. The literal full-domain columns give

    (V-V0)* (V-V0) = (1+theta)/|E| sum_e C_e
                   <= (1+theta)d_max q/|E| I.              (8)

For zero or two holes at a seam the source records raw zero and holds.
For one hole the one-hole calculation applies, with its actual CAR sign.
Different raw records are orthogonal; a moving branch never cancels
an identity branch. The producer verifies the entire Gram matrix,
including off-diagonal entries, for every q of its bounded witnesses.

Both the exact source Hubbard propagation and these repairs preserve the
aligned q-hole space. Interleave m repairs with the same exact propagations
and the same record-adapted controls preserving this sector. A dilation
telescope and trace contractivity give

    ||rho_actual-rho_identity||_1
      <= min(2,2m sqrt((1+theta)d_max q/|E|)).               (9)

This includes the entire source/archive and any reference. It compares
identity repair with the same propagation, rather than omitting the
propagation. On open cubes of side n^2+1 and m<=n^2, (9) is at most
min(2,4sqrt(q)/n). Thus fixed-q free hole histories survive this
specified global-attempt scaling. The full position detector can follow
both histories without increasing their difference; its measurement
disturbance is common to both. A fixed rate per seam is a different
clock. No relation between a repair attempt and model Hamiltonian time,
or any capacity trajectory, is selected by (9).

## 6. A controlled continuum of these observables

Assume c>0 and Delta_infinity>0, and keep these source coefficients fixed.
Prepare wavepackets in the lower hole band, with k=a xi and continuum
Fourier support |xi|<=R. This is an explicit state/scale prescription.
Let n(k)=(cos k_j)_j/|(cos k_j)_j|, n0=(1,1,1)/sqrt(3), and define

    P(k)=[I+n(k).Gamma]/2,
    W(k)=[I+(n(k).Gamma)(n0.Gamma)]/sqrt[2(1+n(k).n0)].

For |k|<=1/2, W is a smooth unitary satisfying W P(0) W*=P(k).
It is fixed by the source matrices and W(0)=I. For a continuum function
f(xi) valued in Ran P(0), the isometry E_a is specified in Bloch space by

    (F E_a f)(k)=pi^(3/2) a^(-3/2) W(k) f(k/a).            (10)

It retains the whole rank-four low-band fibre. No momentum or band
outcome is conditioned on; (10) specifies the prepared state.

Put epsilon(k)=U/2-c sqrt(sum cos^2 k_j), Delta=Delta_infinity.
Writing S=sum sin^2 k_j and using
|S-|k|^2|<=|k|^4/3 gives, for |k|<=1/2,

    |epsilon(k)-Delta-c|k|^2/(2sqrt(3))| <= c|k|^4/3.       (11)

Indeed epsilon-Delta=cS/(sqrt(3)+sqrt(3-S)); the difference from
cS/(2sqrt(3)) is cS^2/[2sqrt(3)(sqrt(3)+sqrt(3-S))^2].
On the model duration T/a^2, remove the fixed-q scalar phase qDelta.
For q=1,

    || exp[-iT(K_rel-Delta)/a^2] E_a f
          - E_a exp[-iT c|xi|^2/(2sqrt(3))] f ||
      <= |T| c a^2 R^4 ||f||/3.                           (12)

Exterior powers multiply this bound by q for q normalized wavepacket
inputs in the specified one-particle band. Density extends the strong
limit to the free q-hole continuum. The kinetic generator is
c P^2/(2sqrt(3)). Its curvature parameter sqrt(3)/c and the offset
Delta are distinct; no relativistic relation between them has been used.
The limit has the free CAR exchange structure, not interacting scattering.

The physical coordinate is the same measured vacancy position, aX.
Its same-source quasimomentum is k/a. On smooth compactly supported f,

    (k/a) E_a f=E_a(xi f),
    ||a X_j E_a f-E_a(i partial_(xi_j)f)||
         <= 2a^2 R ||f||.                                (13)

To prove the second estimate, differentiate (10). The difference is
i a (partial_(k_j)W) f. For |k|<=r<=1/2,
|partial_j n|<=sin(r)/(sqrt(3)cos(r)) and
||partial_j W||<=2|partial_j n|<=2r. The physical-phase Bloch convention
includes the original cell offsets s, so (13) has no omitted internal
position. Momentum and position moments converge on this common core.
Their bounded Weyl operators converge strongly there and then by density;
for position this follows directly by translating the compact Fourier
argument in (10) and using W(a xi)->I. Consequently the complete original
position laws converge weakly to the Fourier position law, with the usual
continuity-set qualification for position bins. All postmeasurement
source components are retained; no low-band projection after readout is
inserted.

The continuum readout therefore has X_j=i partial_(xi_j), P_j=xi_j on
their usual common core and

    [X_i,P_j]=i delta_ij,      Delta X_j Delta P_j >= 1/2.  (14)

This proves the uncertainty relation for the same operational position
and translation momentum in this specified nonrelativistic limit. It
does not obtain an exact finite-dimensional commutator or identify the
spin-response field's position with these hole positions. On q holes the
configuration is antisymmetric; X_center=(sum X)/q and P_total=sum P
have this canonical relation for q>0.

Infinite-lattice formulas in (10)-(14) are limits of the actual open boxes.
For an initially supported one-hole vector, let d be the number of
nearest-neighbour hops needed to reach an omitted boundary. The power
series agree through degree d-1, and ||h||<=3c gives the explicit
propagation difference

    2 sum_(l>=d) (3c|t|)^l/l! <= 2^(1-d) exp(6c|t|).      (15)

For q holes use 3cq in the same estimate. A wavepacket truncation adds
twice its initial discarded norm after normalizing the prescribed finite
preparation; this is a state-approximation error, not an outcome success
division. Compact Fourier smooth functions have computable position-tail
bounds from their derivatives. Thus finite boxes, propagation times,
complete analyzer errors and (9) can be chosen in order to realize each
fixed continuum experiment with any specified tolerance. Model-time
T/a^2, source-box growth, wavepacket preparation and uniform repair
attempts are declared compatible choices. The theorem does not claim
that the native capacity law selects them.

The exact source sector closes an observable position/translation/energy
map and its free continuum. Full-ground relativistic matter, a common
interacting gauge field, laboratory calibration and measured species
assignments do not follow from this aligned construction.
