# Actual ground of the finite native Wilson–Maxwell source

The parent is `model/native_wilson_maxwell`, with `K=3`, `g=1/64`,
`a=1/m`, the same 32 native constituents, full old/dual/port references,
outer integer-Gauss isometry `Y`, positive field boundary completion and
complete complement. Its finite nonstationary response theorem is
unchanged. This proof determines the active full physical ground of that
Hamiltonian and of an explicitly declared positive macrocell completion. Original
spin preparation references are retained zero-energy spectators; ground
multiplicities below exclude their arbitrary states.

## Constants and exact free-vacuum contractions

Write the common raw energy baseline as

`E_base=32 E_g+m(3+sqrt(6))`.

This retains both the native ground constants and the photon zero point.
Let `P` be the complete global occupation-code projector, and `|0>` its
free filled-Wilson-band vacuum with photon vacuum. On the code,

`(H_original-E_base)/m = H_0+V_hop+V_C`,

where `H_0|0>=0`, `H_0` contains all positive Wilson quasiparticle energies,
completed photon number energies and `4|R|/m`, and `V_C>=0` is the
dimensionless Coulomb operator. The free gap is at least 2. Neutral
fermion excitations require a pair with energy at least
`2(sqrt(3)-1/2)>2`. The least photon energy is 2. A nonzero total charge
has the reference cost at least `4/m` in addition to a positive
quasiparticle energy; the certified `m<1` makes its energy greater than 2.

The source's exact 32-dimensional Wilson matrix has
`sign(H_W/m)=sum_{r=3,7,13} M_r/sqrt(r)`, with rational Hermitian `M_r`,
`M_r M_s=delta_rs r E_r` and `sum E_r=I`. These are obtained by exact
CRT from the executed degree-twelve annihilator, with ranks `8,16,8`.
No numerical eigenvector supplies these projectors. The filled projector
is `P_-=(I-sign(H_W/m))/2`; the trace of each onsite sign block is zero,
so all eight free-vacuum density means vanish exactly.

For onsite blocks `P_x`, the full CAR Wick contraction gives

`<rho_x rho_y>=Tr(P_- P_x P_+ P_y)`.

With the same cube Green matrix `G`, exact contraction yields

`C_0=<rho^T G rho>/2`
`   =4055/4368-sqrt(21)/14-3sqrt(39)/208-27sqrt(91)/728`
`   =0.15714743657949333...`.

The sum of directed hopping expectations has positive opposite sign

`h=-<H_hop(0)>/m`
` =-sqrt(3)+6sqrt(7)/7+17sqrt(13)/13`
` =5.250687698400934...`.

All complete magnetic/spin components and all site blocks are included.
The exact rational Coulomb CRT contraction matrix and hopping coefficients
are retained in the receipt. The total charge of `|0>` is zero, hence
`R=0` and `rho_ref=rho` in these expectations. There is no substitution
of a classical Higgs expectation or target mass.

Every edge coordinate has variance
`v=1/16+1/(12sqrt(6))`. Total-cap rotational invariance makes
`w=<0_ph,exp(i g Q_e)0_ph>` independent of the edge. Photon parity makes
it real. At cap three the exact fourth and sixth moments are `3v²` and
`15v³`. Taylor inequalities for the scalar cosine, applied by spectral
calculus, therefore give

`g²v/2-g^4v²/8 <= 1-w <= g²v/2-g^4v²/8+g^6v³/48`.

The exact free-vacuum energy shift `mu=<0,V_hop+V_C 0>` satisfies

`mu_L=g² C_0+h(g²v/2-g^4v²/8) <= mu`
`mu <= mu_U=g² C_0+h(g²v/2-g^4v²/8+g^6v³/48)`.

All algebraic quantities in the following inequalities are evaluated
with outward rational intervals. In particular
`0.0001002309438139758 < mu < 0.0001002309438154073`.

## A lower bound on the entire interacting code

The parent full-Fock Peierls bound on every code charge sector is

`||V_hop||<=v_h=18g+27g²/4=0.28289794921875`.

Put `Q=I-|0><0|` within the code. Positivity of Coulomb gives
`Q(H_0+V_hop+V_C)Q >= Q H_0 Q-v_h Q >=(2-v_h)Q`.
This bound includes all code fermion/photon occupations and all charge
references. No neutral-sector Coulomb norm is applied to charged states.

The vacuum off-diagonal vector is `w_0=Q(V_hop+V_C)|0>`.
Its first-order term is `g V_1|0>`. The actual spinor current annihilates
the filled-vacuum expectation, so this vector lies entirely in the
one-particle, one-antiparticle, one-photon sector. Its squared norm is
the parent's exact positive coefficient `sigma²`, and its free energy
is at least `2sqrt(3)+1`.

For the remainder, the literal exponential integral gives
`||(exp(i g Q_e)-I-i g Q_e)|0_ph>||<=g² sqrt(3)v/2`.
The 24 oriented matter coefficients have norm at most 2, hence the
vacuum Peierls remainder is at most `24sqrt(3)v g²`.
On the neutral code `0<=V_C<=8g²`; positivity implies
`||V_C|0>||²<=8g²<0,V_C0>=8 C_0 g^4`.
Thus

`w_0=g V_1|0>+r`, `||r||<=g² R_v`,
`R_v=24sqrt(3)v+sqrt(8 C_0)`.

Choose `c_*=1/50000`. The weighted Schur complement at energy `c_*`
is bounded using the free spectral support of the first term:

`||[QH_0Q-v_h-c_*]^-1/2 w_0||`
` <= g sigma/sqrt(2sqrt(3)+1-v_h-c_*)`
`    +g² R_v/sqrt(2-v_h-c_*) = S`.

The full interacting complementary resolvent is no larger, in the
positive-operator order. Exact rational evaluation proves

`mu_L-c_*-S² > 3.4530e-6 >0`.

The block Schur criterion therefore proves the operator inequality

`P(H_original-E_base)P >= m c_* P`.                       (1)

The second code eigenvalue is at least `m(2-v_h)` by the min-max
principle applied to the codimension-one subspace `Q`. The first is at
most `m mu_U` by the free-vacuum trial. Hence the code itself has a unique
lowest state `G`, with internal gap at least

`m(2-v_h-mu_U)>1.7170018 m`.                              (2)

Its original energy relative to the common baseline is `m lambda_0`,
with `c_*<=lambda_0<=mu_U`. This code gap does not determine the full
physical ground.

## The unchanged full source ground is in its native complement

The actual native constituent has a unique neutral ground and an exact
first triplet at its three native X momenta. The atomic-spin geometry
theorem applies to the unchanged `lambda_native=2^-34`. Its full
complementary resolvent bound gives

`Delta_native in [2 lambda_native² J_L delta_L-2r_4,`
`                 2 lambda_native² J_U delta_U+2r_4]`,
`r_4=25(360 lambda_native)^4/g_atomic³`.

This exact interval is approximately
`[8.207584263816053e-22,8.210347190546957e-22]`.
The complete charged-port and occupied old/dual reference lower gaps
are checked to exceed its upper endpoint. Consequently the first full
physical constituent excitations are exactly this neutral triplet,
with all added `R_loc` references in their vacuum.

Place one such excitation in any of the 32 constituents and put every
other constituent in `g`. These 96 orthogonal active source states
are outside the complete occupation code: a native neutral triplet
has no code excitation port. On this complement `U_FW=I` and
`Q_native=0`, so all outer density labels, the outer compensator and
Coulomb energy vanish. Kinetic/current couplers are zero by the
declared global code guard. With photon vacuum these are exact
eigenstates of energy `E_base+Delta_native`.

Conversely, every active physical source vector orthogonal to the native
product ground costs at least `Delta_native` in the complete uncoupled
constituent sum, and the code complement adds only nonnegative field,
reference and Coulomb energies. The native product ground belongs to
the code. Thus the complement ground is precisely this rank-96
eigenspace. Equation (1), the source lower bound on `m`, and the gap
upper endpoint prove it lies strictly below the whole interacting code.
The margin is greater than `1.38331e-5` in source energy units.

This is a statement on the complete physical Gauss carrier. Native
Gauss-violating kinematic states are retained for complete instruments
but are not assigned a physical ground interpretation. The outer
Gauss-orthogonal extension and every original zero-energy preparation
reference retain their parent definitions. Spectator reference
degeneracies multiply the quoted 96, which counts active source states.

## An explicitly changed, source-only macrocell Hamiltonian

Define the exact rational coefficient

`c_2=2^-48 ceil[2^48 (C_0+h v/2)]`,
`c(g)=g² c_2`,
`H_repaired=H_original+m c(g) sum_x(I-P_x)` on the logical source.

Here `P_x` is the complete sixteen-occupation code projector at macrovertex
`x` and `P=product_x P_x`. Each macrovertex contains four complete native
constituents. The term is onsite in this declared macrocell decomposition.
The preexisting global guards on kinetic/current terms are retained.

Each `P_x` commutes with the unchanged native/reference constituent sum,
because its code is an invariant spectral subspace. The full port charge
preserves that subspace and its orthogonal complement. This is checked
on the actual selected 81-state four-constituent port carrier and follows
on every additional native state from spectral invariance and fixed
port charge. Since `U_FW` acts only within `P` and is identity outside it,
it commutes with every `P_x`. Hence `rho=U_FW Q_native U_FW†`, its joint
spectral projectors, Coulomb and `Y` commute or intertwine every `P_x`.
An exact coherent pair-rotation and integer-Gauss lift independently
checks this identity with retained charge-carrying complements. All
globally guarded kinetic/current terms are supported where every
`P_x=I`, so they also commute. Thus every occupation-code pattern is a
reducing sector of the full source and of the complete parent counts.

The physical operator is transported by the same `Y`, with the same zero
extension on the outer constraint complement. Its coefficient uses only
exact source contractions and a declared dyadic precision. It vanishes
at `g=0`. At the stated member, `c_2=0.410547438680851...` and
`c(g)=0.00010023130827169213...`. This is a different Hamiltonian.
The unchanged source and receipt keep their original ground and energies.
No implicit vacuum subtraction is made in an old energy measurement.
This choice of source-derived completion is a declared architecture;
the primitive repair law is not proved to select it.

The strict quartic cosine correction gives

`d_sel=c(g)-mu_U>3.64456284e-10`.

The code Hamiltonian and its absolute energies are unchanged. A
complement with `q>=1` defective macrovertices is raised by `q m c(g)`.
Its minimum is exactly `E_base+Delta_native+m c(g)`, attained by a
single native triplet in one macrovertex. The code trial has energy at
most `E_base+m mu_U`, strictly below that complement minimum.
Therefore the full active physical ground of `H_repaired` is the unique
code state `G`. Its gap obeys

`Delta_repaired >= min[m(2-v_h-mu_U),`
`                      Delta_native,L+m d_sel]`
`                 > 2.52078609e-10`.                     (3)

Equation (2) and (3) are different gaps. Every spectator reference is
retained; uniqueness concerns the active source factor. No common raw
energy constant is subtracted. The actual ground energy is
`E_base+m lambda_0`, with `c_*<=lambda_0<=mu_U`.

All code eigenvectors, absolute code energies and internal differences
are unchanged. A q-defect sector acquires phase `exp(-i q m c(g)t)`.
Absolute complement energies and transitions between different patterns
are changed. For any fixed complete history whose Kraus maps commute
with every `P_x`, the extra phase commutes through every intervention
and cancels in each recorded probability. Thus the entire parent
code/count experiment, with its deterministic imperfect preparation,
all complements and all 1161 labels, has exactly the same probabilities
under the two autonomous Hamiltonians. Intended preparation/readout
unitaries are realized with their actual changed Hamiltonian in the
declared controller; raw work and phase accounting use `H_repaired`.
Arbitrary coherent interference between different code patterns is
outside this count-law invariance claim.

## Ground overlap and finite energy-agreement preparation

Choose the phase of `G` so `<0|G>>0`, write `G=a|0>+q`, and use its
eigenvalue equation on `Q`. The vacuum remainder bound gives

`||q|| <= kappa |a|`,
`kappa=[g sigma+g²R_v]/[2-v_h-mu_U] <0.01002601`.

Thus `|<0|G>|²>=1/(1+kappa²)`. The ground is charge zero, since it is
unique, the Hamiltonian conserves total charge, and its overlap with
the charge-zero free vacuum is nonzero.

The parent's deterministic physical trial, with every spin reference
and `R_loc` vacuum retained, has its exact support decomposition
`P rho_trial P=p_code |0><0|`, `p_code>=p0^32` on the active code.
Consequently its actual repaired-ground probability is at least

`p_init=p0^32/(1+kappa²)>0.9998985719755`.

This preparation does not use `G`, its unknown exact energy, or an exact
ground projector. For the full raw Hamiltonian choose the integer bound
`B=19397`, bounding the parent raw norm plus the complete `8m c(g)`
completion norm, and the algebraic Cayley unitary

`C=(I-i H_repaired/B)(I+i H_repaired/B)^-1`.

It is unitary on the full finite kinematic carrier. On the physical
ground/excited pairs the phase separation is at least
`Delta_repaired/B`, since both raw energies are in `[-B,B]` and

`|C(a)-C(b)|=2|a-b|/[B sqrt((1+a²/B²)(1+b²/B²))]`.

Set `epsilon_target=1/1024` and
`M=ceil[2B/(epsilon_target Delta_repaired,L)]`.
The exact integer in the receipt is `157589952376326706`.
Use an M-state uniform pointer, controlled powers of
`W=C tensor C†`, inverse uniform preparation, and a measurement of
every pointer basis state. Pointer zero realizes
`F_M(W)=M^-1 sum_{j=0}^{M-1} W^j`.
All other outcomes are retained. Mixed ground/excited amplitudes are
bounded by the exact rational `epsilon=2B/(M Delta_repaired,L)`;
equal ground energies pass with amplitude one.

For two independent source marginals with ground probability at least
`p=1/(1+r)`, the success probability is at least `p²`. The conditional
ground-pair fidelity is at least
`1/[1+r²+2epsilon² r]`. The same lower bound applies to either retained
main source. This follows from the orthogonal ground/ground, mixed and
excited/excited energy sectors and holds for mixed source states with
arbitrary retained purifying references.

The complete two-level tree uses four original source trials and three
filters. With `r_0=(1-p_init)/p_init`,
`r_{j+1}=r_j²+2epsilon²r_j`, its all-pass probability is at least
`(1+r_0)^-4(1+r_1)^-2>0.9995943286660`.
The final main source has ground fidelity at least `1/(1+r_2)`, and
trace distance below `1.417923e-7`. All four source copies, their
references, every pointer outcome and every failed branch are retained.
If a preceding filter fails, its descendant slot executes a declared
no-op instrument and retains the earlier failure label; its elapsed
source time is identical. No fresh trial is silently substituted.

For any fixed subsequent complete instrument on the main source, the
conditional output law is within this trace-distance bound of the
actual stationary ground law. Its free evolution under the actual
`H_repaired` preserves the bound. With all purifying partners retained,
the corresponding vector error is bounded by
`sqrt(2-2sqrt(F))`, and the comparison partners undergo their own actual
source evolution. They need not be asserted to be exactly in their
ground states.

## A complete stationary dressed-ground record

Let `T_111` project onto the declared free-basis count `(1,1,1)` inside
the exact code. Its free energies lie below `11m`, and the current's
first-order vacuum vector lies entirely in it. Photon parity removes
the even Peierls terms and Coulomb from its vacuum insertion. The
parent odd remainder bound gives

`||T_111(V_hop+V_C)|0>-g V_1|0>||<=27g³/16`.

On the entire neutral code
`||V_hop+V_C||<=nu=18g+59g²/4`. Projecting the exact code-ground
equation and using `||q||<=kappa`,

`||T_111 G|| >=`
` [(g sigma_L-27g³/16)/sqrt(1+kappa²)-nu kappa]`
` /(11+mu_U) = a_G >0`.                                 (4)

The receipt gives its exact lower probability and the actual physical
preparation/readout joint lower bound
`q_tree,L (a_G-epsilon_purification)²`.
The final detector is the parent's complete 1161-count instrument, and
all filter failures and partner states remain in the joint record.
Equation (4) is stationary virtual occupation in the specified free
quasiparticle/photon frame. No radiation from a ground state, asymptotic
scattering, observed photon production or SM mass is inferred.

## Source clock, finite controls and capacity

Every trial copy evolves under the same repaired raw source Hamiltonian,
including zero points, all reference costs and positive macrocell terms. The
four-copy Hamiltonian has norm at most `4B`. The deterministic preparation
uses one duration-one pulse `i Log(S)-H_copies`. A filter is a declared
finite matrix service `U_filter`; a duration-one pulse
`i Log(exp(-i H_copies) U_filter)-H_copies` implements it with a final
common source evolution. All Cayley powers and uniform pointer gates
commute with `H_copies`, including their conditional failure controls.
The final common evolution cannot change their energy probabilities or
the fidelity estimates. All branches have elapsed source time four
after the three filter slots. An additional duration-one complete count
readout includes its own final common evolution and all labels.

These are explicit programmable finite-matrix controls, with a specified
Hermitian logarithm branch in `(-pi,pi]`; they are not an efficient
autonomous local implementation of the passive repair law. The inverse,
controlled powers and M-state Householder preparation are finite
algebraic services. Pointer dimension, coefficient description/workspace
and Hamiltonian-control strength are separate resources. No unknown
ground vector is an apparatus coefficient.

Each control norm is at most `4B+pi`; total net bare-source work is at
most `8B` and gross four-pulse switching work at most `8(4B+pi)`.
The complete preparation capacity is

`4[32(256ln2+8ln5+24ln119+ln(2^31+1))+7ln9+ln65+ln56]`
` +3ln M` nats.

This counts all four native sources and their original preparation
references, longitudinal labels, static compensators, photon caps and
the three M-state filter registers. A subsequent count readout adds
`ln1161`. Controller coefficient/workspace costs are excluded and
explicitly declared. A separately certified error of the whole physical
preparation/readout isometry must be subtracted from the unconditional
target amplitude before squaring; it cannot be placed inside a
conditional probability multiplied by an unchanged success weight.

## The guarded alternative and locality boundary

For comparison, the distinct Hamiltonian
`H_guarded=H_original-m c(g)P` also selects the same code ground vector.
It leaves complements unchanged and shifts all absolute code energies.
It is related to the primary positive completion by
`H_repaired-H_guarded=m c(g)[sum_x(I-P_x)+P]`.
This difference is the constant `m c(g)` on code and single-defect
sectors and is `q m c(g)` on a q-defect sector. The two operators are
not identified by one global energy subtraction on the entire source.
The primary receipt, raw norm, filter clock and ground energies all use
the positive macrocell completion without any scalar subtraction.

The parent kinetic/current Hamiltonian continues to use its global
code guard and extended Foldy-Wouthuysen/Gauss coordinates. A local
completion of the energy penalty alone does not prove a full-source
local propagation cone or an interacting continuum limit. All ground
estimates here concern the stated cap-three finite member. The cutoff
theorem for the parent's nonstationary record laws does not supply
uniform interacting-ground control as cap, volume or spacing changes.
