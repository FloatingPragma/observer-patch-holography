# Conserved field spin and a native spin-current continuity equation

## Source and carrier

Use the existing full-domain CAR extension with two modes at each site. The
ordered basis is site outer, flavor inner. The actual endpoint arithmetic is
the binary restriction of `model/observer/native.py`: an occupied endpoint
moves into an empty endpoint or stays, according to the fair coin. The retained
`full_domain.py` law flags a hold when either endpoint has occupation two.
This is an explicit extension of the binary source, not a claim that the scalar
source already specified collisions or quantum flavor. Its complete coin/raw
instrument and source compiler already exist. We retain that flagged policy;
the separate coherent-collision policy is not substituted during this test.

For each branch h=(coin,raw), write K_h=V_h/sqrt(2). V_h is a signed partial
permutation, with the CAR sign determined by occupied modes crossed by a hop.
Every basis state appears exactly once for each coin. The inverse of a moving
branch is fixed by its recorded signed current. Consequently
`sum_h V_h^dagger V_h=2I`. The reverse count gives
`sum_h V_h V_h^dagger=2I`. Empty raw slots and collision slots are retained.

## Field Noether charge

Let `S_x^a=(1/2) sum_(alpha,beta) c†_(x,alpha) (sigma_a)_(alpha,beta)
c_(x,beta)` and `S^a=sum_x S_x^a`. These are operators of the CAR field, not
an externally attached scalar record counter. On an empty or doubly occupied
site all three vanish; on single occupation they act as sigma/2. They generate
the common flavor rotation `Gamma(I_sites tensor exp(-i theta.sigma/2))`.

The occupation projectors commute with every S^a. A moving V_h has the form
`sum_alpha c†_(b,alpha)c_(a,alpha)` restricted to its one-to-zero domain, or
the reverse hop. It contracts the flavor index with the identity and therefore
commutes with each S^a. The hold branches are occupation-domain projectors and
also commute. This proves, including CAR signs and the collision complement,

`[V_h,S^a]=0` and `[V_h,N]=0` for every complete branch.

Products of such branch maps obey the same equations. The claim therefore
holds on arbitrary untouched references and under adaptive seam choices based
on past records. It does not extend automatically to spin-changing feedback,
nor to a final spin measurement that destroys incompatible spin coherence.

The field generators satisfy `[S^a,S^b]=i epsilon_(abc) S^c`. On the full
one-particle sector, of dimension twice the number of sites,
`sum_a (S^a)^2=(3/4)I` and `spec(S^z)={-1/2,+1/2}`. A full 2pi flavor
rotation is exactly `(-1)^N` by second quantization. Thus field occupation,
fermionic exchange and the Spin2 central sign live on the same chosen CAR
object. Selection of this fundamental flavor representation is an architecture
input; its conservation is a new consequence of the executed repair law.

## Local transport is carried by the same recorded hop

Let e=(a,b), let j_h be the signed current in raw h, and let
`E_h=V_h^dagger V_h`. Define the doubled spin `T_x=2S_x`. Then for each
component,

`V_h^dagger T_a V_h - T_a E_h = -j_h E_h(T_a+T_b)`,

`V_h^dagger T_b V_h - T_b E_h = +j_h E_h(T_a+T_b)`.

All spectator differences vanish. On a moving branch precisely one endpoint
is occupied, and its spin is carried unchanged to the other endpoint; this
proves the formula. On any hold j_h=0 and the relevant spin commutes with the
domain projector. Multiplying by one half converts V effects to normalized
Kraus effects; multiplying by another half converts T to S. This identifies
the same branch's local spin current as `j_h E_h(S_a+S_b)/2` in the normalized
instrument. Spin current is operator-valued: a signed charge-current bit alone
does not reveal a noncommuting spin component. The complete pointer experiment
implements a specified component read and retains its measurement disturbance.

## Stationary field state

On N=1, `rho=(I_sites/sites) tensor |up><up|` is stationary under each fixed
seam channel, since the spatial channel is unital and the flavor factor is
unchanged. It has sharp S_z=1/2, Casimir3/4 and zero mean occupation current
through every seam. Its preparation is specified as a field input marginal;
this packet does not claim a newly executed uniform spatial preparation.
Arbitrary record-adaptive feedback need not preserve this stationary marginal.
The result is a stationary polarized field, not yet a stationary magnetic
dipole. In particular zero occupation circulation is reported rather than
renamed as a nonzero magnetic response.

## Geometric covariance and the physical boundary

The actual finite port-spin bridge is stronger than an arbitrary two-state
label. `code/a5_closure/sl2f5_port_spin_bridge_certificate.py` reconstructs
the source SL(2,F5), its committed12-port quotient with kernel +/-I, and its
faithful complex Spin2 representation. It checks all14,400 products and the
commuting cover square over Q(sqrt(5),i). The port-spin lift is built from the
declared sigma-dot-vector current representation; it does not derive that
representation from scalar probabilities. The source-executed independent
Spin action is also implemented in `model/observer_particles/spin_matter_join.py`.

Attach that existing lift to this field by
`c†_(x,alpha) -> sum_beta (u_g)_(beta,alpha)c†_(g x,beta)`.
Second quantization gives the correct reordering signs automatically. Flavor
contraction in a hop and its endpoint occupation domains are invariant, so a
geometric relabel sends its instrument to the relabeled seam instrument.
For reversed endpoint order the coin is complemented, current is negated and
the signed raw slots are relabeled. This is an exact finite covariance proof
for the declared attachment; no new all-group numerical run is claimed here.

The established conserved quantity is the field's global Spin2 Noether charge.
An identification with physical spatial angular momentum must still relate
this action to the observer's spatial rotations and orbital contribution.
No Maxwell Poynting angular momentum, physical electron, continuum rotation
group or hbar unit calibration follows from the finite conservation proof.
This limitation does not negate the finite port-spin theorem or assert a
scalar-memory obstruction for this quantum field architecture.

## Magnetic interface still to construct

`Lean/EventAlgebra/PauliMinimalCoupling.lean` proves the exact Pauli square.
To consume it here, one needs a kinetic triple X on this same field carrier,
whose commutator equals the curvature measured by an independently specified
gauge/field detector. Charge must enter that kinetic connection through the
same transported number current above. One then needs a stationary field
response and the orbital rate in the same operational clock before reporting
a gyromagnetic ratio or Larmor law. Merely defining B from [X_i,X_j], adding
an independently programmed Pauli rotation, or importing the installed
infinite-photon Hamiltonian would not supply that join. This packet makes no
new g claim and does not mark rows5 or6 physically complete.

This is a same-field composition gap, not an absence of native gauge or loop
constructions. `native_electric_gauss` and `native_observer_program` already
execute a declared integer/Gauss extension on tagged two-load states;
`observer/holonomy_reconstruction.py` already supplies the complete conditional
face-phase readout. `native_circulation` has nonzero cycle-current fluctuations
and an exact axial covariance. NEXT_INTERFACE.md identifies these providers,
their actual assumptions, and the field/carrier map still needed here.

## Executable checks

The preregistered integer audit checks all32 branches on the full two-site16
and separated-endpoint four-site256 Fock spaces, and all66 unordered pairs
on the twelve-port24-dimensional one-particle sector for both C/E laws.
The Y generator is represented as i times an integer matrix, so the algebra
and continuity checks use exact integer arithmetic. The full source replay
uses the existing six-level pulse compiler and complete512-address dilation,
then copies occupation bits to four fresh native pointer bits. It retains
all8192 joint addresses for all16 input columns in compressed arrays.

The changed-law operation is a source-compiled global fermionic flavor swap.
It leaves every native branch effect, current and occupancy unchanged, but
reverses S_z. It therefore tests sensitivity to the transported spin law
while keeping the charge-record object fixed. Complete pointer laws distinguish
planted up/down inputs and reverse under that control. Binary64 source errors
are recorded as diagnostics, separately from exact algebra.
