# Complete massive neutral current construction

## Declared source

The complete native contact of `model/openwave_weak_decay` has fifteen
internal coordinates, three families, two spins and three momenta. Add one
gauge-neutral Weyl coordinate per family, giving 288 fermionic modes. The
31464 native quartic words are relabelled into these modes unchanged. All
quark, charged-lepton and other spectator words occur; the added singlets
have zero native weak-current entries. The interaction coefficient sum is
exactly 41952, hence its full Fock-space operator norm is at most 41952.

Charged masses are `(m,1,M)`, the native analyzer ratios in middle-charged
units. Its rational eigenvalue enclosures give `0<m<1` and `M<17`.
Let `R=F† V_(2,3,10) F`, with the complete native family word and the ordered
charged-mass frame. Declare neutral masses `mu_i=epsilon(m,1,M)_i` and
tensor `R diag(mu)`. The primary member takes `epsilon=2^-32`; a control
takes `epsilon=1/4`. The fixed Higgs is `(0,1)` and the outgoing box has

`q=sqrt((delta²-(mu_0+mu_1)²)(delta²-(mu_0-mu_1)²))/(2delta)`,

where `delta=1-m`. This box is selected from the source masses and is a
declared architecture input. It places the `(0,1)` neutral pair on the
parent/recoil free threshold; it does not impose a lifetime or remove other
channels. The control retains channels whose free thresholds are closed.

The native Higgs tensor `L` satisfies `t_a^T L+L S_a=0` for all three actual
weak generators, with `S_a=-(L†L)^-1 L†t_a^T L`. The actual hypercharge
matrix gives `Y^T L+L/2=0`. Its fixed neutral component has electric charge
zero. The added Weyl coordinate is a colour, weak and hypercharge singlet,
so it contributes zero to every gauge and mixed gauge anomaly polynomial.
Assign Dirac lepton number +1 to the active lepton doublet and -1 to its
charged and neutral conjugate singlets. Every Yukawa pair and native current
conserves it. A Majorana term is absent by this declared symmetry. Neither
the added module nor its coefficient is claimed unique or forced by OPH.

## Complete free diagonalization

For spin sign `sigma=1-2s`, momentum `k=pq` and positive mass `mu`, put
`E=sqrt(k²+mu²)`, `z=-sigma k`,

`u=sqrt((E+z)/(2E))`, `v=sigma sqrt((E-z)/(2E))`.

Algebraically equivalent formulas in the code avoid cancellation when one
helicity weight is very small. The charged transformation is
`c_L=u a_L-v a_R†`, `c_R=u a_R+v a_L†`. For neutrals rotate the active
row by `R` and keep the sterile mass frame. These maps preserve every
canonical anticommutator and diagonalize the complete doubled operator,
whose pair block is `Delta_(nu_j,N_i)=sigma R_ji mu_i` and whose transpose
is antisymmetric. All massless spectator coordinates are retained with
their appropriate filled negative bands and vacant zero bands. Massive
neutral zero-momentum reference modes are filled in the negative Dirac
band. This reference differs from a massless zero-vacancy preparation.

The free spectrum is positive after its explicit reference constant is
removed. Each mass is a Dirac particle/antiparticle pair with two helicities;
an active and a sterile chirality do not constitute two mass branches.
The complete 576-dimensional doubled-operator calculation checks all signs
and family rotations against this declared source.

## Actual native vertex and exact isometry

Prepare a middle charged particle at rest with spin `p`. Retain every output
with a light charged rest recoil of spin `s_e`, neutral mass label `i` at
momentum `o q`, and antineutral mass label `j` at momentum `-o q`, for
`o=+1,-1`, every mass label and every spin. There are 144 output states.
Their matrix is obtained from the complete native quartic current, using
the actual Bogoliubov transformation and canonical fermionic signs.

The direct exact spin contraction is zero unless
`s_antinu=p` and `s_e!=s_nu`. On that support it gives

`W_(oij,s_e,s_nu,s_antinu;p)`
`= -o (-1)^(p+s_e) u_i(o q,s_nu) v_j(o q,1-p) conjugate(R_1i) R_0j`.

The producer checks the exact rational native CAR coefficient for every
spin and orientation and compares this formula against every direct
matrix element in both massive members. The charged rest coefficients
contribute the exact product one half; the native spin contraction supplies
the corresponding coefficient. No effective transition matrix is installed.

For each fixed momentum, the two neutral-spin `u²` weights sum to one.
For each fixed antineutral spin, its two opposite-momentum `v²` weights sum
to one. Both relevant rows of the unitary family matrix have squared norm
one. Thus each parent column has norm one. Different parent spins have
disjoint antineutral-spin supports, so their overlap vanishes. Therefore
`K=P_all144 V P_parent` satisfies `K†K=P_parent` exactly. The Kraus pair
`K, I-P_parent`, with separate flags, is complete on the entire source and
preserves arbitrary external quantum references by the finite instrument
theorem. This is a specified current apparatus, not autonomous decay.
All output energies and the apparatus mean work relative to the same free
Hamiltonian are retained, including the off-shell outputs of the control.

## Unfiltered short-time source bound

The autonomous source is the complete `H0+GV`, with no vertex projection
in its Hamiltonian. Its interaction-picture Duhamel remainder after first
order is at most `G² ||V||² t²/2`, by unitarity. In the primary member all
144 target free energies differ from the parent energy by less than one:
`m+2q<=1`, `mu_max<17/2^32` and `E_i<=q+mu_i` give an upper energy below
two, while every energy is positive. The first-order phase correction
therefore has norm at most `G t²/2` on the exact unit-norm target vertex.
Consequently, with `r=(G+G² 41952²)t²/2`,

`(Gt-r)² <= P_all144(t) <= (Gt+r)²`.

At `G=1/64`, `t=10^-8`, one has `0<r<Gt`. The rational receipt stores both
probability bounds and retains the full complementary Fock outcome.
This short-time theorem is independent of the Fock dimension and applies
with idle external references. It supplies no long-time decay law.

## Flight instruments and work

A complete controlled experiment selects the spin-up light recoil and one
momentum/spin orientation from the current apparatus while retaining the
coherent three-by-three neutral/antineutral mass matrix. Every other
production outcome belongs to its retained failure branch. The contact is
switched off during propagation; this is an explicitly controlled protocol
on the same free 288-mode source. Relative phases are calculated with the
stable identity `(mu_i²-mu_0²)/(E_i+E_0)`.

For each neutral, the three active-port rows are its actual native
annihilation amplitudes in the free positive band. The remaining diagonal
Bogoliubov weight is a quantum sterile-failure Kraus map. Their four
Kraus squares sum to the identity. Applying both complete instruments gives
sixteen active/failure outcome pairs, plus the original production failure.
Each probability is an original-input probability; no production-success
renormalization occurs.

An active-port success leaves the vacuum in that particle factor and
transfers its conserved Dirac lepton number to a retained reference. The
reference has declared zero energy; an absorbed neutrino contributes +1
and an absorbed antineutrino -1. A sterile failure retains its mass-state
quantum factor. Summing the actual output free energies, including the
charged recoil and every failure, gives the recorded mean detector work.
These are declared native-coordinate ports. No physical charged-current
detector implementation or measured apparatus energy is inferred.

## Energy calibration and scope

The existing electron anchor is the sole empirical energy scale. The
neutral coefficient, family word, fixed Higgs, reference, box, current
apparatus and detection protocol are declared inputs. The inherited
charged-lepton phase and retrospective comparison ancestry are retained.
The conditional neutral masses are about `0.000118976`, `0.0246007` and
`0.413737` eV. There is no neutrino-data fit, physical neutrino identity,
measured mixing matrix, sterile-particle discovery or physical weak lifetime
claim. Source time orders the specified experiment; no universal clock is
identified. Exact identities and rational short-time bounds are separated
from floating source/record diagnostics and their recorded tolerances.
