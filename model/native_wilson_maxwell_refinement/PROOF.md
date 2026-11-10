# Total-photon cutoff on the fixed native Wilson cube

The source is the finite construction in
`model/native_wilson_maxwell`, with its actual native constituents, selected
positive-cross energy `m`, complete four-mode occupation code at each of eight
vertices, references, charge density, Gauss embedding, preparation and clock.
The only refinement parameter in this result is the total transverse photon
cap `K`. The finite member with `K=3`, `g=1/64`, `mT=1/64` is the parent
source. This proof compares its complete records with an unbounded-oscillator
mathematical approximation on the same finite cube. It supplies no spatial
continuum, interacting relativistic limit or particle-identity transfer.

## Source, domain and literal finite link matrices

Let `F` be the five-mode oscillator Fock space, `N` total number, and
`P_K=1_{N<=K}`. The dimension is `binom(K+5,5)`. Frequencies in units of `m`
are `nu=(2,2,2,sqrt(6),sqrt(6))`. With the parent orthonormal transverse
frame, write

`Q_e=a(c_e)+a†(c_e)`,
`||c_e||²=v=1/16+1/(12sqrt(6))`.

The mathematical comparison has
`H_ph=m[dGamma(nu)+Tr(nu)/2]`. Its finite compression is exactly the parent
completed field Hamiltonian at cap `K`: adding
`P_{N=K}[dGamma(nu)+Tr(nu)]m/2` to the finite quadratic expression gives
`P_K H_ph P_K`. Every boundary sector is retained. In particular the free
photon evolution commutes with every `P_K` and intertwines its finite
compression without error.

The actual finite link is `W_{e,K}=exp(i g P_K Q_e P_K)` on `P_K F`.
The unbounded comparison link is `W_e=exp(i g Q_e)`. These are different
operations: `W_{e,K}` is not identified with `P_K W_e P_K`. The finite
commutator and noncommuting-link boundary controls in the parent apply.
This proof bounds the literal exponentials themselves.

On the full logical native/reference carrier the Hamiltonian can be written

`H_K=H_M + H_ph,K + sum_{j=1}^{24} A_j tensor W_{j,K}`.

Each oriented term is accompanied by its adjoint, `||A_j||<=2m`, and its
phase is either sign of a parent `Q_e`. `H_M` includes the unchanged native
sum, the free code correction with its zero-phase hopping removed, Coulomb,
the complete `4|R|` reference cost and every retained native complement.
It is finite and self-adjoint, independent of photon number and of `K`.
The oriented hopping terms vanish on the native code complement.
Consequently `H_infinity` with `W_j` is self-adjoint on `D(H_ph)` by the
bounded perturbation theorem; the perturbation norm is at most `48m`.
Coulomb is bounded on this fixed finite source and is included in `H_M`.
There is no infinite physical register or unbounded control primitive in
any finite member.

The same outer integer-Gauss isometry `Y` acts trivially on the photon
factor. Each physical Hamiltonian is `Y H_K Y†` plus zero on the outer
Gauss-orthogonal kinematic complement, as in the parent. The comparison
uses the identical extension. Inner native Gauss-violating states belong
to the retained native complement, with zero probability under the stated
preparation. No charge sector or reference cost is removed.

## A displacement-word estimate

In the interaction picture of `H_M+H_ph`, each photon factor is a displacement
`D(z_j)` with `||z_j||=r=|g|sqrt(v)`. Free rotation changes its direction,
not its length. A product of `n` such factors applied to vacuum is, up to
a phase, a normalized coherent vector of radius at most `R=nr`. The
finite product uses precisely the compressed, freely rotated generators.

For a coherent vector of radius at most `R`, its total number is Poisson
with mean at most `R²`. The falling-factorial moment gives

`||(1-P_K)|z>|| <= R^(K+1)/sqrt((K+1)!)`.

The exact layer bound is
`||1_{N=K+1}|z>|| <= R^(K+1)/sqrt((K+1)!)`.
There are no independent-mode or product-cutoff assumptions.
The only projected generator residual is

`P_K Q_e (1-P_K)=P_K a(c_e) 1_{N=K+1}`,

whose norm is `sqrt(v(K+1))`. Projected Duhamel along all `n` displacement
segments, of total coefficient length `nr`, therefore bounds the
projected difference by
`sqrt(K+1) R^(K+2)/sqrt((K+1)!)`.
Adding the final orthogonal tail gives the full embedded-word bound

`b_K(R)=min(2, R^(K+1)/sqrt((K+1)!) [1+sqrt(K+1)R])`.       (1)

The initial vector is vacuum in both programs. All matter coefficient
products factor from the photon word; their norms multiply this estimate.
No estimate on the boundary operator in the full Hilbert-space norm is
claimed to vanish with `K`.

## Complete autonomous evolution

Put `tau=mT`, `x=48tau`, and choose any integer `J` with `x<J+2`.
Both interaction-picture Dyson series converge in operator norm, since
the interaction is bounded. Their terms through order `J`, evaluated on
an arbitrary native/reference vector tensored with photon vacuum, obey
(1). The two omitted absolute series have sum at most

`R_J(x)=2 x^(J+1)/(J+1)! / [1-x/(J+2)]`.

Thus, including every virtual fermion occupation and every native sector,

`||(U_K(T)-U_infinity(T)) psi|| <= delta_K`,
`delta_K=min(2, sum_{n=0}^J x^n/n! b_K(nr)+R_J(x))`.       (2)

Embedding of `U_K` and the finite initial vector is understood in (2).
The equality of free evolutions supplies the same bound in the source
picture. It holds with arbitrary retained initial quantum references.

For each fixed `n`, `b_K(nr)` tends to zero. Since `b_K<=2`, dominated
convergence of the factorial series proves cutoff convergence at fixed
bounded `g,T`. A constructive diagonal choice is `J=floor(sqrt(K+1)/8)`
when `r<=1`, or division by `8 max(1,r)` in general. For sufficiently
large `K` its maximal radius is at most `sqrt(K+1)/8`, while the omitted
Dyson tail tends to zero. The factorial bound
`(K+1)! >= [(K+1)/e]^(K+1)` proves exponential decay of the retained
word bound. This argument does not use an unsplit high-number-moment
series that could cease to be a useful bound at large `K`.

The unbounded comparison is also a well-defined vector evolution on the
oscillator domain. A coherent word satisfies
`||N^q |z>||²=Touchard_{2q}(||z||²)`, a fixed polynomial in `||z||²`.
Multiplying its square root by `x^n/n!` is summable for every fixed `q`.
The vacuum Dyson series therefore converges in each `D(N^q)` graph norm.
This establishes the domain for the electric/magnetic expectation and
commutator statements below; no numerical occupation tail is substituted
for a domain argument.

## Actual preparation, all records and the cap-three result

The deterministic parent trial retains all original spin references and
all added `R_loc` vacua. Its exact support identity is

`P_code Q_all_Rloc_vac=|g^32,refvac><g^32,refvac|`.

It is unaffected by changing `K`, by `Y`, or by free photon embedding.
The actual complete count law at every cap, and in the comparison, is
`p_code` times the code-vacuum law, with the same source probability
`p_code>=p0^32`, plus `1-p_code` on native-complement/photon-zero.
The complement never interacts with the photons on this input. No
postselection, initial exact-ground projection, discarded trial, or
normalization by a rare response probability is used.

For a common finite outcome alphabet, count `N_a,N_b=0,...,16`, native
code/complement, and `min(N_ph,K+1)`, with `K+1` explicitly labelled
overflow. Retain the distinct outer-Gauss failure effect. There are
`290(K+2)+1` comparison labels. The finite `K` law embeds the parent's
`290(K+1)+1` labels by zero probability on the `290` added overflow
labels. This is an embedding of record laws, not a claim that the
parent executed a larger pointer. Binary hit/complement is a common
coarse graining without any additional label.

The same complete count Stinespring isometry acts on the embedded finite
and infinite vectors. Trace distance, complete classical total variation,
and any fixed subsequent common detector processing are bounded by
`min(1,delta_K)`. The unbounded overflow probability is at most `delta_K²`
because the embedded finite vector has no overflow component. This
includes terminal quantum branches immediately after this common count
isometry. It does not identify terminal source states after two different
Hamiltonians have acted for an additional uncounted duration.

The parent readout pulse is `exp(-i H_K) U_count`: its final common source
evolution leaves every recorded pointer probability unchanged. Hence (2)
at the parent's autonomous wait `T` bounds its actual complete pointer
law after the full elapsed source time `2+T`. If that terminal source is
used in another experiment, its additional duration must be included in
the history theorem below. Every finite preparation/readout control is
the parent's declared Hermitian-log finite-matrix service with its actual
`H_K` cancellation and elapsed duration. No uniformly bounded infinite
controller or efficient passive native compilation is asserted.

At `K=3`, `g=tau=1/64`, `J=24`, rational outward calculation gives

`delta_3 < 1.863e-9`.

Let `a_*` be the parent's positive ideal target amplitude lower bound,
and `b=tau*nu`, where `nu=4667/16384`. The exact count target is
`(N_a,N_b,N_ph)=(1,1,1)`. Its unbounded-comparison actual probability is
at least `p0^32(a_*-delta_3)²`; each ideal excitation marginal is at most
`(b+delta_3)²`. The actual binary connected response is at least

`p0^32[(a_*-delta_3)²-(b+delta_3)^4] > 0`.                  (3)

These are amplitude estimates before squaring. The unconditioned
comparison uses the same imperfect native trial and exact support weight.
Its positive response is a finite-time dressing of the prepared free
vacuum, not radiation from the coupled ground or an asymptotic scattering
amplitude. `results.json` gives rational bounds and decimal diagnostics.

For an arbitrary finite cap `L`, its comparison to the cap-three record
law is bounded by `delta_3+delta_L` through the same mathematical source.
No empirical datum or parameter is introduced by this comparison.

## Bounded electric/magnetic histories

For a fixed complete history, allow arbitrary number-preserving
contractions/isometries, including matter/reference operations, coherent
pointer controls, retained number measurements, feedback controlled by
stored orthogonal records, and free photon rotation. Also allow finitely
many literal finite displacement generators with total radius bound `s`.
The same generators define the infinite comparison. Their branch controls
can be coherent; no successful-record conditioning is used.

For a state or retained branch vector write `f_l=||1_{N=l}psi||`. A linear
photon generator of coefficient norm `r` is dominated in these sector
norms by `r(b+b†)` on the one-mode number sequence. Number-preserving
contractions do not increase any `f_l`. The positive majorant after
total radius `R` is

`exp[R(b+b†)]|0> = exp(R²/2) sum_l R^l/sqrt(l!) |l>`.

This follows directly by expanding each exponential in an absolutely
convergent number-sector power series on finite analytic vectors and
then taking its norm limit. It also applies to a bounded-control multiple
of the generator, because its block norm has the same ladder bound.
Consequently the tail is at most
`exp(R²) R^(K+1)/sqrt((K+1)!)`, and the layer is at most
`exp(R²/2) R^(K+1)/sqrt((K+1)!)`.
Projected Duhamel therefore gives the history-word estimate

`b^hist_K(R)=min(2, R^(K+1)/sqrt((K+1)!)`
`                         [exp(R²)+sqrt(K+1) R exp(R²/2)])`. (4)

Replace `b_K(nr)` in (2) by `b^hist_K(nr+s)`. Here `T` is the sum of
all autonomous source waits, including free source factors after
intermediate readout pulses. The multilinear Dyson simplex volumes give
the same total factor `T^n/n!` when waits are separated by retained
instruments. This proves complete joint-history and terminal-branch
trace-distance control, for every fixed finite such program.

The explicit diagonal schedule still proves convergence: for sufficiently
large `K`, `Jr+s<=sqrt(K+1)/4`. In (4) the exponential factor contributes
at most `exp[(K+1)/16]`, whereas the factorial tail contributes
`(sqrt(e)/4)^(K+1)` times a polynomial. Their product tends to zero.
The Dyson tail also vanishes. For an explicit graph-norm control, let
`A_j psi` be the tensor of all ordered `j` annihilators applied to `psi`.
Its squared norm is the falling-factorial moment `(N)_j`. The norm
`G_q(psi)=sum_{j=0}^{2q} ||A_j psi||` is equivalent to the `D(N^q)` graph
norm. A displacement of radius `r` obeys
`G_q(D(z)psi)<=(2q+1)(1+r)^(2q) G_q(psi)` by the binomial expansion of
`a+z`. Number-preserving contractions decrease every summand. Thus the
history Dyson graph-norm bound is a factorial series with only an
additional fixed exponential factor per word. It converges. Arbitrary instruments
that can inject an uncontrolled photon tail are outside this theorem.

## Local gauge-invariant comparison algebra

On the exact global occupation code let `Psi_x` be the parent's local
canonical Wilson spinor, `rho_x=Psi_x†Psi_x-2`, and
`rho_ref=rho+R delta_r`. The actual Foldy-Wouthuysen frame convention,
full-code charge intertwiner and coherent Gauss label map are unchanged.
In the unbounded comparison define

`E=P_T + g D G rho_ref`, `B=C Q`,
`F_e=(Psi_head† B_e† Psi_tail) exp(i g Q_e)`.

The incidence convention is `D_e,tail=-1`, `D_e,head=+1`, with every
edge directed along its positive coordinate axis. A face indexed by
`(x,i,j)`, `i<j`, has boundary
`+(x,i)+(x+e_i,j)-(x+e_j,i)-(x,j)`. The electric edge and magnetic face
signs in the following identities use exactly these parent conventions.

Here `P_T` denotes the transverse electric operator vector in this formula;
the transverse projector is written `Pi_T` below. The exact cochain
relations give `Pi_T+D G D†=I`, `C Pi_T=C`, `D†Pi_T=0`.
The two contributions to the electric-hop commutator are

`[P_T,e,F_f]=g Pi_T,ef F_f`,
`[g(D G rho_ref)_e,F_f]=g(D G D†)_ef F_f`.

Thus `[E_e,F_f]=g delta_ef F_f`. Also
`[B_f,E_e]=i C_fe`, `[rho_x,F_e]=D_ex F_e`,
and `[E_e,rho_x]=0`. Electric operators commute mutually, magnetic
operators commute mutually, and disjoint even fermion hops commute.
These are local incidence relations in the gauge-invariant algebra.
The source's actual transported current is
`I_e=-i(F_e-F_e†)=-J_phase,e`, so `dot(rho)=D† I`.
The full field energy is exactly
`m/2 sum_e E_e² + m/2 sum_f B_f²`; its longitudinal term is the parent
Coulomb energy. On this code the matter energy is the local Wilson
onsite and dressed-edge sum, minus its fixed filled-free-band constant.
This is an equality on the admitted code, not a replacement of the
globally projected native Hamiltonian on its full complement.

The same commutators imply, on the common finite-particle core,
`dot(B)=m C E` and `dot(E)=-m C†B+g I`.
These are finite-graph Maxwell/current identities. They do not require
commuting finite-cap link coordinates: the literal cap operators have
boundary defects. Their bounded observational comparison is supplied
by (4), rather than by asserting these identities exactly at the cap.

An observer can test bounded `exp(i u B_f)` or `exp(i v E_e)` by a
controlled displacement and a complete two-outcome Hadamard readout.
The electric longitudinal factor `exp(i v g(D G rho_ref)_e)` acts only
on the retained finite matter/reference space and preserves photon
number. The transverse factors are literal finite compressed linear
generators. A magnetic test of strength `u` adds radius
`|u| sqrt(<B_f²>_vac)`; an electric test of strength `v` adds
`|v| sqrt(<P_T,e²>_vac)`. All control outcomes are retained, and their
finite sequences are covered by (4), including state disturbance.
The exact covariances and incidence cancellation are machine checked.

There is a quantitative observer-order consequence on this same source.
Take electric edge `0=(000,100)` and magnetic face `3`, the `xy` face at
`z=1`, in the parent orientation. Their incidence is zero. The infinite
electric and magnetic Weyl tests commute, including their separate
coherent Hadamard controls. Perform the two test factors in either order
with all control records retained, and with no autonomous wait inserted
between the factors. Each finite history is within (4) of the same
comparison history. Their complete joint record laws are therefore
within twice the history bound. At strengths `1/64`, following the parent
wait, the resulting bound is below `1.856e-6`.

The actual finite linear generators do not commute in operator norm on
this nonincident pair. Let `B=q(s)`, `E_trans=p(t)` with real coefficient
vectors, `s dot t=0`, `||s||²=var_B`, `||t||²=var_E`. Their truncated
commutator is `-i P_{N=K} dGamma(s t^T+t s^T)`. The one-particle matrix
has eigenvalues `+sqrt(var_B var_E)`, its negative, and zeros, so its
norm on the top sector is exactly `K sqrt(var_B var_E)`. The longitudinal
electric part commutes with `B`. This explicit boundary countercontrol
excludes operator-norm commutator convergence on the entire finite
carrier, while the complete prepared-record order bound tends to zero.
The order comparison is an equal-time algebra statement; no
separated-time propagation speed is inferred from it.

This gives a shared algebra of bounded observer tests and a controlled
cutoff interpretation on one finite cube. It does not prove a uniform
spatial Lieb-Robinson cone or Lorentz cone. The full native program uses
global spectral code compression and a spatially extended frame; its
total photon cap is also global. The unbounded constrained comparison
has the local gauge-invariant relations proved here, while a growing
volume, shrinking spacing, ultraviolet control and uniform propagation
bound require their own estimates. The fixed-cube photon refinement
does not supply those estimates.

## Capacity and controller boundary

The finite source capacity replaces the parent's `ln(56)` by
`ln binom(K+5,5)` and its executed count pointer `ln(1161)` by
`ln[290(K+1)+1]`. The native/reference capacity is unchanged. The
mathematical overflow comparison is an isometric padding of the law;
if implemented as an extra finite pointer it instead costs
`ln[290(K+2)+1]`. Additional bounded-field tests carry their actual
retained control/result registers. `K` is a declared register capacity,
not a scheduler prediction or a physical-universe entropy estimate.

The free photon norm grows linearly with `K`; a safe bound is
`m[3+sqrt(6)+sqrt(6)K]`. The oriented Peierls sum is uniformly bounded
by `48m`, while the photon-free native/reference/Coulomb part has the
parent's finite full-carrier bound. Finite Hermitian-log control costs
therefore have explicit cap dependence. The comparison requires exact
finite services as in the parent. For a separately established complete
implemented-pulse isometry error `eta`, add `eta` to the state-vector
comparison before converting to probabilities; rare-event lower bounds
use amplitudes before squaring. No unbounded service is hidden in a
finite member's capacity or clock.
