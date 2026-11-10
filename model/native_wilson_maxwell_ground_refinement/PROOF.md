# Interacting ground records as photon capacity grows

This theorem fixes the eight-vertex open cube, all 32 native constituents,
their complete old/dual/port and spin references, the outer integer-Gauss
isometry, the Wilson current, the source-selected mass `m`, spacing `a=1/m`,
and coupling `g=1/64`. The only autonomous-source refinement parameter is
the total occupation cap `K>=3` of the same five transverse oscillators.
The Hamiltonian is the positive macrocell completion of
`model/native_wilson_maxwell_ground/PROOF.md`:

`H_K=H_original,K+m c(g) sum_x(I-P_x)`.

The rational `c(g)`, all code projectors and all physical references are
identical at every cap. This is the declared completion, not a consequence
that the primitive repair law selects its coefficient. The global code
guards on the kinetic/current terms remain part of the source. Every
finite member has finite photon dimension `binom(K+5,5)`. Five unbounded
oscillators occur only in the mathematical comparison `H_infinity`.
No volume, spacing, ultraviolet or relativistic limit is taken.

The raw baseline is retained:

`E_base=32 E_g+m(3+sqrt(6))`.

The positive field boundary completion makes the finite photon Hamiltonian
exactly `m[3+sqrt(6)+dGamma(2,2,2,sqrt(6),sqrt(6))]` on every cap. It is
not the uncompleted quadratic expression at the top layer. The link is
the literal `W_e,K=exp(i g P_K Q_e P_K)`. No equality with
`P_K exp(i g Q_e) P_K` is used.

All ground assertions concern the full physical Gauss carrier. The active
source ground excludes degeneracy from retained zero-energy preparation
spectators. Native Gauss-violating sectors and the outer constraint-failure
sector remain in the complete apparatus. The outer failure Hamiltonian
has its declared zero extension; it is not assigned a physical ground.

## A common lower form bound and isolated ground

On the complete code write `h_K=(H_K-E_base)/m=H_0+V_hop+V_C`.
The completion vanishes there. The free operator satisfies `H_0>=2N`,
has vacuum `|0>`, and on `Q=I-|0><0|` has lower bound 2. This includes
all charge sectors: a nonzero total charge carries the unchanged `4|R|`
reference cost. Coulomb is nonnegative on the complete carrier.

Set

`v=1/16+1/(12sqrt(6))`, `A=96 g sqrt(v)`.

The 24 oriented matter hopping coefficients have norm at most 2 in
units of `m`. The coordinate estimate
`||Q_e psi||<=2sqrt(v)||sqrt(N+1)psi||` holds on the infinite oscillator
form domain and on every literal truncation. Spectral calculus gives
`||(exp(i g Q_e)-I)psi||<=g||Q_e psi||`. Thus

`||V_hop psi||<=A||sqrt(N+1)psi||`.                    (1)

For a normalized vector, put `u=<H_0>`. Cauchy-Schwarz and
`N<=H_0/2` give `|<V_hop>|<=A sqrt(u/2+1)`. Young's inequality yields,
for any positive `epsilon`,

`V_hop >= -epsilon H_0-b I`,
`b=A^2/(8epsilon)+2epsilon`                           (2)

as a quadratic-form bound. This argument uses no reversed Jensen
inequality. With `epsilon=1/10`,

`b<0.471464444`, `D=0.9 H_0-b`,
`D|_Q>=d_0=1.8-b>1.328535556`.

The vacuum moments of orders 2, 4 and 6 are `v,3v^2,15v^3` at every
`K>=3` and at infinity. Each closed length-six ladder path from vacuum
reaches occupation at most three. Consequently the parent's exact
vacuum expectation bounds hold uniformly:

`mu_L<=<0,V_hop+V_C 0><=mu_U`,
`0.0001002309438139758<mu_L<=mu_U<0.0001002309438154073`.

The exact Coulomb contraction is

`C_0=4055/4368-sqrt(21)/14-3sqrt(39)/208-27sqrt(91)/728`.

The exact opposite hopping expectation is
`h=-sqrt(3)+6sqrt(7)/7+17sqrt(13)/13`. Explicitly,

`mu_L=g^2 C_0+h(g^2 v/2-g^4 v^2/8)`,
`mu_U=mu_L+h g^6 v^3/48`.

The actual first current vector has `||V_1|0>||=sigma`, with the parent's
independently executed exact current coefficient. It lies in the
one-particle, one-antiparticle, one-photon sector and has free energy
at least `2sqrt(3)+1`. The complete vacuum off-diagonal vector obeys

`Q(V_hop+V_C)|0>=g V_1|0>+r`,
`||r||<=g^2 R_v`, `R_v=24sqrt(3)v+sqrt(8C_0)`.

The Coulomb norm 8 used here is restricted to the neutral code, where
the vacuum and its image lie. It is not a bound on every charged sector.
Applying (2) on `Q` and the monotonicity of positive inverse forms,
the Schur complement at `c_*=10^-6` is positive if

`mu_L-c_* > [g sigma/sqrt(0.9(2sqrt(3)+1)-b-c_*)`
`             +g^2 R_v/sqrt(d_0-c_*)]^2`.

Outward rational evaluation gives a strict margin greater than
`7.774765e-6`. Therefore `h_K>=10^-6` on the entire code, uniformly
in `K>=3` and in the comparison. Min-max on the codimension-one space
`Q` gives `lambda_2(h_K)>=d_0`, while the vacuum trial gives
`lambda_1(h_K)<=mu_U`. There is a unique active code ground `G_K`, with

`gamma_code>=d_0-mu_U>1.3284353259` in units of `m`.    (3)

Writing `G_K=a|0>+q` with `a>0`, its projected eigenvalue equation gives

`||q||<=kappa a`,
`kappa=(g sigma+g^2 R_v)/(d_0-mu_U)<0.012958608`.

Thus `a^2>=1/(1+kappa^2)`. Charge conservation and nonzero overlap
with the neutral vacuum make this unique ground neutral.

The native complementary sectors are unchanged by photon refinement.
Their minimum is `E_base+Delta_native+m c(g)`, attained by one of the
96 active neutral triplet states with photon vacuum. Arbitrary retained
spectator states multiply this active multiplicity. The parent's exact
native complement proof includes all charged ports and old/dual reference
sectors. Since the same coefficient satisfies

`c(g)-mu_U>3.64456284e-10`,

the completed full physical ground is `G_K`, and its full physical gap is

`Delta>=min[m_L(d_0-mu_U),Delta_native,L+m_L(c(g)-mu_U)]`
`     >2.5207860907995e-10`.                           (4)

The code gap (3) and the full gap (4) are distinct. The uncompleted
Hamiltonian has its active physical ground in the rank-96 native
complement at every cap, because its code lies at least `10^-6 m`
above the common baseline. Photon refinement does not repair that source.

## Stationary occupation and actual finite preparation

Let `T=T_111` be the complete-code count of one positive quasiparticle,
one antiparticle and one photon. Its free energies are below 11.
Photon parity and the exact sixth moment give the uniform insertion bound

`||T V|0>-g V_1|0>||<=4sqrt(15) v^(3/2) g^3`.

On the range of `T`, (1) gives
`||Q V T||<=A sqrt(2)+8g^2`; by adjointness this also bounds `T V Q`.
The ground equation therefore proves

`||T G_K|| >= a_G`
`=[(g sigma_L-4sqrt(15)v^(3/2)g^3)/sqrt(1+kappa^2)`
`  -(A sqrt(2)+8g^2)kappa]/(11+mu_U)`
`>0.0006720456715`.

The probability is greater than `4.51645384e-7`. This is stationary
virtual occupation in the declared free Wilson/photon frame. It is not
radiation from the ground or an asymptotic scattering observable.

The original deterministic native trial has exact support
`P rho_trial P=p_code |0><0|`, with `p_code>=p0^32`. Its source ground
probability is at least `p=p0^32/(1+kappa^2)>.9998311854`. No exact
ground projection or additional conditional preparation is assumed.

Choose the fixed integer `B=19397`. This bounds the absolute raw ground
energy at every cap, since

`|E_G|<=32(480+360*2^-34)+m_U(3+sqrt(6)+mu_U)<B`.

It does not bound all excited energies as `K` grows. The Cayley unitary
`C_K=(I-i H_K/B)(I+i H_K/B)^-1` is nevertheless well-defined and unitary
for every finite source and for the self-adjoint comparison. If `E` is
an excited physical energy, `d=E-E_G>=Delta`, and
`A_0=sqrt(B^2+E_G^2)<=sqrt(2)B`, then

`|C_K(E)-C_K(E_G)|=2Bd/[A_0 sqrt(B^2+E^2)]`
` >=2Bd/[A_0(A_0+d)]>=Delta/(2B)`.

The final inequality uses `Delta<=B`. It holds for arbitrarily high
excited energies and avoids phase wrap at an uncontrolled energy cutoff.
Take

`M=ceil[4B/(epsilon_target Delta_L)]`
` =315179904752653412`, `epsilon_target=1/1024`.

The parent's complete energy-agreement instrument uses the fixed
M-state pointer and controlled powers of `C_K tensor C_K^dagger`.
Every pointer basis outcome is retained. Its zero outcome has mixed
ground/excited amplitude at most `epsilon=4B/(M Delta_L)`.
This is a separate apparatus setting from the parent's smaller K=3
pointer; the autonomous source family and the frozen parent receipt
are unchanged.

For ground odds `r`, one successful two-copy filter has odds at most
`r^2+2epsilon^2 r`. The depth-two tree uses four independent physical
trials and three filters. With `r_0=(1-p)/p`,

`r_{j+1}=r_j^2+2epsilon^2 r_j`,
`q_tree>=(1+r_0)^-4(1+r_1)^-2>.99932485508`.

The retained main source has conditional fidelity
`F>=1/(1+r_2)>1-5.583e-14`, and purification vector error
`e_p=sqrt(2-2sqrt(F))<2.362628e-7`. All main/partner sources and their
purifying references, all failures, every pointer outcome and all
conditional no-op records are retained. Contractivity transfers this
conditional state estimate to arbitrary fixed subsequent complete
instruments. The original-attempt stationary count lower bound is

`q_tree (a_G-e_p)^2>4.51023170e-7`.                    (5)

This uses the actual full physical ground of the completed source at
the same cap, rather than replacing a prepared state by an exact ground.

## Positive heat baseline and a quantitative ground comparison

Heat evolution in this section is a mathematical proof device, not the
physical clock or an apparatus primitive. Remove all zero-phase hopping
from the free matter Hamiltonian. The resulting matter operator is

`H_M=Psi^dagger (5 beta/2) Psi+S_0+V_C+4|R|/m`,
`S_0=4sqrt(3)+8sqrt(7)+4sqrt(13)`.

The fermion term has lower bound -40, and `S_0-40>2.5164188`.
Therefore `H_M>=0`. This requires no commutativity between the onsite
matrix and Coulomb: each summand after the displayed constant is
nonnegative. The positive photon number operator `H_N` commutes with
matter. The remaining hopping has 24 oriented terms `A_j W_j` with
`sum||A_j||<=48`. In this split there are no subtraction words. The
Schrodinger heat-Dyson expansion alternates number-preserving
contractions and literal displacements. It never conjugates bounded
terms by an inverse heat semigroup.

The parent's general-history cutoff bound applies to such contractions.
For total displacement radius `R`, cap `K`, it gives the sufficient bound

`b_hist(K,R)<=R^(K+1)/sqrt((K+1)!)`
`             *exp(R^2)*(1+sqrt(K+1)R)`.

For heat parameter `s`, subsequent real waiting time `tau` and fixed
bounded field tests of total displacement radius `r_app`, the absolute
Dyson weight is `x=48(s+tau)`. Truncating at order `J`, use
`R=J r+r_app`, `r=g sqrt(v)`. If `3R^2<=(K+1)/4`, the elementary
factorial bound `n!>=(n/e)^n`, `e<3`, proves

`delta_heat <= exp(x+R^2)(1+sqrt(K+1)R)2^(-K-1)`
`  +2x^(J+1)/(J+1)!/[1-x/(J+2)]`.                    (6)

The normalized heat vector differs from the ground by at most
`kappa exp(-gamma_code s)`. Its unnormalized norm is at least
`exp(-mu_U s)/sqrt(1+kappa^2)`. Normalizing the two heat vectors and
then applying their complete histories gives the vector bound

`delta_ground_history <=2kappa exp(-gamma_code s)`
` +2sqrt(1+kappa^2)exp(mu_U s)delta_heat`.              (7)

Every complete record law has total-variation distance at most (7).
For `s=12`, `tau=1/64`, two electric/magnetic tests of strengths `1/64`,
`r_app=(sqrt(1/2)+sqrt(2))/64`, `K=1024`, `J=1700`, (7) is less than
`3.093163e-9`. The photon dimension is `binom(1029,5)`; no matrix of
that dimension is numerically enumerated. This is the ideal interacting
ground comparison, distinct from the actual finite preparation below.

Each field test is a complete controlled unitary instrument using the
same total electric field, magnetic face flux and Gauss/reference
coordinates as `model/native_wilson_maxwell_refinement`. The longitudinal
electric part is finite matter/reference data. The transverse part is
a literal finite Weyl exponential. All control and pointer outcomes
are included, so the general-history estimate applies before readout.
For the displayed two-test settings, their controlled factors and both
pointer records form one duration-one final readout service, with no
autonomous wait between the factors. Its final common source evolution
does not affect either pointer law. They are an alternative to the final
count setting, rather than tests after an uncounted disturbed count state.

## Norm resolvents and a fixed complete preparation alphabet

Embed each cap into the common mathematical Fock space. Extend its
coordinate by zero outside the cap and its exponential by identity.
Use the number-preserving zero-phase-removed baseline `a=H_M+H_N`
on the code, and the unchanged native/Coulomb/completion operator on
its complement. The added photon sectors then have the free matter
operator after the identity link is reinserted. This is a mathematical
comparison extension; the finite program contains no omitted sectors.
On the code the extended operator is `a+b_K`, with `a` fixed,
`b_K=sum_j A_j W_j,K` bounded by 48 and
`b_K -> b_infinity` strongly. Strong convergence follows first on
finite-occupation vectors by the displacement estimate, then by density.
The resolvent of `a` is compact because matter is finite and there are
five positive oscillator frequencies. A uniformly bounded strongly
convergent family times a compact operator converges in norm. The
resolvent identity consequently proves norm-resolvent convergence.
Outside the code the difference of hopping operators is zero. The
identical outer zero extension also contributes no difference, even
though that extension itself has no compact resolvent. The isolated
rank-one active ground projections therefore converge in norm.

For a quantitative estimate, let `L=floor(K/2)`, `d=K+1-L`, and
`r_Q` be an individual displacement radius. The ladder estimate

`||Q^j P_{<=L}||<=(2sqrt(v))^j sqrt((L+j)!/L!)`

and `j!>=(j/3)^j` bound the first escaping Taylor term by
`rho^d`, where `rho=6 r_Q sqrt(K+1)/d`. Subsequent terms have ratio at
most `rho`. The full-Fock overflow is at most `rho^d/(1-rho)`.
For the retained part, boundary Duhamel has a layer-K term starting at
Taylor degree `d-1`. Integration contributes `1/(j+1)` and the boundary
ladder coefficient contributes `r_Q sqrt(K+1)`. Its first term is at
most one half the escaping degree-d term. Thus

`||(W_K-W_infinity)P_{<=L}||<=2rho^d/(1-rho)`.         (8)

For `rho<=1/16` and `K>=7`, (8) is at most `(K+1)^-2`.
This follows directly at `K=7,8` and by increasing K by two: d increases
by one, the exponential bound decreases by 16, and the squared
polynomial factor increases by less than 4. The link radius is below
`1/200`; the two field tests individually have radius below `1/32`.
All explicit budgets below satisfy the stronger guards in (8).

For the raw Cayley transform set `B=19397` as above and assume
`K>=ceil(4B/m_L)`. On low photon occupation the fixed baseline resolvent
has norm at most `m_U/B`; on the higher sector it has norm at most
`2/(K+1)`. The latter uses the retained raw baseline bound and the
positive number energy. The hopping difference has global norm at
most 96 and low-sector norm at most `48/(K+1)^2`. The resolvent
identity, with Neumann factor `1/(1-48m_U/B)`, gives

`||C_K-C_infinity||`
` <=2/(1-48m_U/B)[48m_U/(B(K+1)^2)+192/(K+1)]`
` <400/(K+1)=d_C`.                                   (9)

The identical complement blocks make (9) a full extended-carrier bound.
For the two-source Cayley comparison W, the error is at most `2d_C`.
Controlled powers through `M-1` therefore have norm error at most
`2(M-1)d_C`; three complete filters contribute at most
`6(M-1)d_C`. This is a comparison of the complete isometries with the
same finite pointer alphabet, before conditioning on any outcome.

The duration-one controller for each filter includes the common source
free factor. Every controlled Cayley gate and failure/no-op control
commutes with each individual source Hamiltonian. Hence the three
common free factors may be moved before the filter tree. They act on
the original four trials with photon vacuum, with dimensionless time
`tau=3m`. The infinite-source commutator `[a_j,W_e]` is a bounded
multiple of W, and the vector annihilator estimate gives

`||a_vector exp(-i h t)psi_vac||<=48 r t`.

This can be proved on finite-particle cores by Duhamel and then by
closedness of the annihilator vector. The photon-number-preserving
baseline rotates that vector unitarily; no exponential growth factor
is introduced. Splitting the evolved vector at photon L, using (8),
Markov's inequality and ordinary unitary Duhamel yields

`d_U(K,tau)<=48tau/(K+1)^2`
`              +2304 r tau^2 sqrt(2/(K+1))`.          (10)

Four copies and the complete filter tree therefore obey

`delta_preparation_and_count<=4d_U(K,3m_U)`
`                              +6(M-1)400/(K+1)`.    (11)

The count projectors agree under the common embedding, with all photon
overflow and all native/outer complement outcomes retained. For finite
K the count alphabet is `290(K+1)+1`; the mathematical comparison has
the corresponding countable photon alphabet. Finite zero-probability
overflow records are padded for comparison. The final common source
factor after readout acts after the record is stored and cannot change
its law. At `K=2^96`, (11) is less than `9.548479e-9`.

For additional bounded field tests, the complete filter tree commutes
with each source energy, so the unconditioned energy expectation is
unchanged, including all failed branches. The exact trial support and
its reference vacuum give

`<h><=mu_U+(1-p_code)[64B_native/m_L+8c(g)]`,
`B_native=480+360*2^-34`.

On the trial complement the port references and all charges are vacuum;
the bound therefore needs the native norm, not an occupied reference
cost. Formula (2), and the positive complement energies, imply the
full-branch mean bound

`<N><=(<h>+b)/1.8<0.284617466`.

On the infinite comparison branch, applying a displacement increases
`||sqrt(N)psi||` by at most its radius. No exact finite-cap Weyl
commutator is asserted. Splitting that comparison branch at L bounds
each test's finite/infinite discrepancy
by `(K+1)^-2+2sqrt(<N>/(L+1))`. With the two-test accumulated radius
included before both tests, (11) plus both discrepancies is less than
`9.548491e-9` at the same `K=2^96`. These are unconditioned complete
record comparisons. Conditional-law bounds require the explicit
success denominator, such as `q_tree` above; no failure is discarded.

## Physical clock, controls, finite capacity and scope

Every finite member runs its actual raw `H_K`, with its native and
photon zero points and every reference/completion term. Trial preparation
occupies one source-time unit, the three energy-agreement slots occupy
three more, and complete count readout occupies one: elapsed time is
five on every success and failure branch. Each matrix service is
realized by the same specified Hermitian-logarithm driven pulse as in
the ground parent, including its final common free factor where stated.
No exact spectral projection or unknown ground vector is a coefficient.
The fixed Cayley transform is a declared exact finite matrix service;
this is not an efficient autonomous local compilation.

A sufficient full raw norm bound for every finite cap is

`B_raw(K)=32(600+360*2^-34)`
` +m_U[185+(K+1)sqrt(6)+224g^2+8c(g)]+128`.

It includes the full native references, all charged Coulomb sectors,
the unchanged zero point, positive completion and all retained kinematic
extensions. The four-copy norm is at most `4B_raw(K)`; each duration-one
control has norm at most `4B_raw(K)+pi`. Net source work is at most
`8B_raw(K)`; gross switching work of the four-pulse preparation tree is
at most `8[4B_raw(K)+pi]`. A subsequent readout has its own control and
work budget. These resources grow with K even though the pointer M is
fixed. The infinite comparison has no finite control-norm claim.

The full preparation log dimension in nats is

`4[32(256ln2+8ln5+24ln119+ln(2^31+1))`
`  +7ln9+ln65+ln binom(K+5,5)]+3ln M`.

The complete count pointer adds `ln[290(K+1)+1]`. Every retained source,
partner and spin reference is counted. Coefficient-storage workspace
and exact-matrix service implementation resources are separate. At
the displayed huge cap the receipt evaluates this capacity outward;
the source matrix is not enumerated and feasibility is not asserted.
Separately certified physical service errors must be added at the
complete-isometry level before any success conditioning.

The result connects actual finite preparations, the interacting ground
and complete gauge-invariant field/count records to a fixed-cube
oscillator comparison. The global code guards, finite volume, boundary
Gauss reference, Foldy-Wouthuysen coordinates and apparatus choices are
preserved. A spatial propagation cone, local relativistic field limit,
interacting ultraviolet limit and observed particle identification do
not follow from photon-capacity convergence.
