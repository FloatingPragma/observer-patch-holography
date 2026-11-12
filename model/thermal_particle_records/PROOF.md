# Complete finite thermal particle-count records

The construction consumes the calibrated quadratic source of
`particle_joint_mass_identity` and the finite Gibbs preparations of
`particle_thermal_history`. Both receipts must bind the same
`conditional_particle_spectrum/results.json` bytes. The separately
compiled lepton source and the continuum ideal gas are not used here.
All empirical Yukawa, gauge, Higgs and optional Weinberg inputs retain
their ancestry. Source temperatures are supplied parameters.

## Group number projectors on every occupation

The parent joint observables partition the 180 canonical fermionic
quasiparticle modes into 84 groups, and the 56 bosonic modes into 28
groups. Each group is a simultaneous one-particle eigenspace of energy,
electric charge, colour Casimir, helicity and momentum. Its comparison
species label is attached after this spectral grouping. The original
canonical quasiparticle fields, rather than an unrelated number basis,
define

`N_g=sum_(j in g) beta_j† beta_j`.

For fermions of group rank `r`, the spectrum is `0,...,r`, with
projector ranks `binom(r,n)`. For bosons the global occupation cap two
gives spectrum `0,1,2`. Define the literal spectral projections by

`P_(g,n)=product_(ell!=n) (N_g-ell I)/(n-ell)`.

Their orthogonality and sum-to-identity follow by polynomial evaluation
on that complete finite spectrum. Different group number operators
commute by the canonical relations and the orthogonal one-particle
partition. Products of their projectors give a complete joint projective
instrument, with output `P_n rho P_n`. Tensoring each projector with the
identity preserves the whole original 474-coordinate reference and any
additional quantum reference. This is a declared number-measurement
service; a narrower primitive pulse implementation is not claimed.

There are 48 rank-three fermion groups and 36 rank-one groups. Their
joint count alphabet has `4^48*2^36=2^132` outcomes. It is stored exactly
as 84 finite factors, not sampled or truncated. The rank sum is
`product_g sum_n binom(r_g,n)=2^180`. The bosonic count alphabet has
435 nonzero-rank outcomes: zero, one, or two quanta distributed among
28 groups. A count vector has rank
`product_g binom(r_g+n_g-1,n_g)`. Their sum is 1,653. All other vectors
in `{0,1,2}^28` have the zero operator because their total count exceeds
two; the receipt retains their exact symbolic zero-effect rule and count.
The total rank is the complete source dimension `2^180*1653`.

## Gauge invariance and the actual polarization frame

An independent partition of the original mode inventory gives 48 full
quark colour triples and four full gluon colour octets. Every number
group contains or excludes each whole colour fiber. Thus its one-particle
projector commutes with every unitary on those fibers, including the
actual global colour representation. Second quantization gives
`[N_g,dGamma(T_a)]=0` and then `[P_(g,n),dGamma(T_a)]=0`.
The source energy, excitation charge and momentum also commute with
the number projectors. Local colour Gauss constraints, confinement and
the total many-body colour Casimir are not resolved by these counts.
In particular one-particle Casimir labels are not added as if they were
the total many-body Casimir.

The boson groups use circular-polarization columns. The linear-polarization
occupation words of the Gibbs producer cannot be assigned those helicity
labels directly. The executable constructs the actual block rotation U
and its full capped symmetric-Fock lift on all 1,653 states. For input
two-boson words `i<=j` and output `a<=b`, its matrix entry is

`(U_ai U_bj+U_aj U_bi)/sqrt((1+delta_ij)(1+delta_ab))`.

Vacuum and one-boson blocks complete this isometry. Every column is
checked for unitarity and energy/charge intertwining. The rotation mixes
only exactly equal-frequency and equal-charge modes and preserves the
global cap. Consequently `Gamma(U) rho_Gibbs Gamma(U)†=rho_Gibbs` exactly;
all entries of this identity are checked for each finite Gibbs control.
Only after this attachment are the original word probabilities read as
occupation probabilities in the joint polarization basis.

A group-count projection leaves all coherences within a given count
block. An entangled one-quantum source/reference state in a colour fiber
is unchanged by its count-one outcome. Measuring the fine colour instead
changes that conditional density by unhalved trace norm `2(1-1/r)`.
Independent full density controls at `r=3,8` give `4/3` and `7/4`.
This distinguishes a gauge-invariant group count from a classical colour
label measurement.

## Complete finite Gibbs output laws

For the supplied source temperatures `T=0.1,1,10,100 GeV`, let
`q_g=exp(-E_g/T)` and `f_g=q_g/(1+q_g)`. The complete fermion group law is

`Pr(N_g=n)=binom(r_g,n) f_g^n (1-f_g)^(r_g-n)`.

The full joint law is the product of all 84 factors. Every count from
zero through the rank is retained, with finite log probabilities even
when binary probabilities underflow. The bosonic producer's 1,653
occupation words are assigned to their 435 group-count vectors, with
every original word index retained. Log-sum-exp sums each original-weight
fiber. Its combinatorial rank agrees with the independent symmetric-Fock
formula, and all group count probabilities sum to one.

The resulting group and species means are measurable finite counts.
Reported fractions are `E[N_g]/E[N_total]`, obtained by pooling complete
trials including empty trials. They are not the expectation of a
per-trial ratio conditioned on nonempty output. Independent checks
reconstruct total number, energy, excitation charge and momentum from
the group means and compare them with the original Gibbs calculation.
Binary quadrature or probability evaluations are not certified hardware
accuracy or rigorous arithmetic intervals. Exact spectral identities
and numerical diagnostic precision are kept distinct.

## The vacuum/one-excitation detector's initial complement

The parent's vacuum/one-excitation projectors are complete only with their
many-excitation complement. On the thermal preparation, before applying
any creation, removal or no-operation probe, write

`Z_F=product_(i=1)^180(1+exp(-E_i/T))`,
`Z_B=1+sum_j q_j+sum_(i<=j)q_i q_j`.

Then

`Pr(N_total<=1)=(1+sum_F exp(-E_i/T)+sum_B exp(-E_j/T))/(Z_F Z_B)`.

The complementary probability is retained, together with the log of
the small accepted-sector probability. This is an initial-state
diagnostic. It is not the parent's 570-record law after its additional
particle probes. The full number instrument resolves every occupation
while preserving colour coherence, so it supplies a finite thermal
measurement without requiring that thermal states lie in a one-particle
sector. No native temperature, interacting particle identity or
early-Universe abundance is derived.
