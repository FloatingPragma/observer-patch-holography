# A source-specific colour/weak current connection

This additive package reuses the finite SU(3) chain native source e222d339 and its
independent full review d301c89b. It does not replace that source's validation,
identify it with the fixed-Q193-mode Dirac source, or award a physical row.
The source, preparation, clock and terminal instrument are specified below.
All coefficients remain architectural source parameters; no empirical mass,
Fermi constant, lifetime or axial target has been inserted.

## The existing baryon band has an exact weak-current intertwiner

On the full colour-epsilon three-quark domain Sym³(C^(4L)), let
J_a=sum_(x,s) q†_(x,s) tau_a/2 q_(x,s). The chain hopping, colour Casimirs and
local spin Casimirs are weak-label invariant, so [H,J_a]=0 exactly. They also
preserve total spin. This reuses the existing60680-CAR/2840-singlet and tree-Gauss
isometries; those expensive enumerations are not repeated.

The co-located symmetric spin/weak triple decomposes as
(S,I)=(3/2,3/2) plus(1/2,1/2). In the fixed Sz=1/2 middle weak sectors, the old
isolated rank-two band lies in[55/16,57/16]. A co-located spin3/2 vector has floor
69/16; the old P/Q square completion subtracts at most1/8, leaving67/16 above
the old61/16 exclusion threshold. An I3/2,S1/2 vector has no co-located component
and is above the separated floor281/48. Thus the entire low band has S=I=1/2.
This identifies the band representation; a numerical eigensolver does not
supply the representation assumption.

Let F flip every quark weak label. On three fundamental doublets,
exp(-i*pi*Jx)=iF. On I=1/2 the same rotation is-i(J_++J_-). Therefore the neutral
to charged low-band map is J_+=-F, and J_-J_+=I there. The p/n spectra agree
exactly. Use the phase convention |p,s>=-F|n,s>; the encoded ladder is then a
unit tag raising operation. The explicit L2..5 source checks preserve this sign.
Outside the band, the I3/2 middle ladder has amplitude2, not1. The local current
at the endpoint also has off-band matrix elements. Neither fact is erased.

## A native measurable current, with an explicit charge reference

A bare n/p tag rotation changes charge. The reader instead adjoins a one-lepton
spin-up reference with labels nu,e and uses

    C = J_+ tensor e†nu + J_- tensor nu†e
      = (X_b X_l - Y_b Y_l)/2

on the two low bands. Here b=0,1 is n,p and l=0,1 is nu,e; total charge is b-l.
The reader exchanges n nu and p e; the two unmatched charge sectors are fixed.
The four-bit old data code has14 physical labels and2 padding labels. Multiply
C by the exact valid-data projector, making all padding identity. A Ramsey probe
controls exp(-i*pi*C/2). All Pauli strings commute; there is no splitting error
in this reader. Its emitted rational native controls still have their explicit
rounding/export error, so exact ideal encoded charge conservation is not a
claim that a finite-error native logical matrix is algebraically exact.

The new raw executor exponentiates the literal six-level response generators
for every pulse and accumulates every first-export effect. It checks all128
input columns, including false controls, unmatched charges and invalid data.
The charge-zero native Bell preparation has C eigenvalue+1. Its controlled
current gives X/Y probe contrasts0,-1. The same clock exp(-i*tau*Z/16) gives
clock contrasts cos(tau/8),sin(tau/8), with tau=1/8 in the fresh controls.
The floating matrix replays are diagnostics; their roundoff is not an interval
proof. The exact rational control intervals and reused six-level response theorem
supply the rigorous native bounds independently.
All four baryon/lepton/probe/clock pointers are actual native eight-pulse
instruments. Every16-choice C/E combination uses the joint four-pointer density,
not a product of marginals. Fine coin/current or coin/energy records, pointer
exports, source exports and16 finite blank-reset timeouts are retained.

These fresh zero-code executions test the apparatus. The full band experiment
is the separately executable lazy program in band_protocol.py: run the old
native Bell purification and31-round QPE filter, keep every rejection and phase
record, and **do not measure the baryon data** before the new current reader.
The full paired n/p code has equal H blocks because of the exact F symmetry.
Evolve that common source and apply its clock once. The complete bound adds
the old exact preparation/evolution certificate and the new reader errors;
BAND_PROTOCOL.json states the success lower bound and conditional comparison
explicitly while retaining unconditional failure outcomes. The old huge filter
and complete new joint word are unexecuted. No new ground-state or full-word
execution is inferred from component tests.

The global band reader preserves total electric charge but need not preserve
each cell density when its compensating lepton reference is fixed at endpoint0.
A future Maxwell-coupled version of this global apparatus needs the corresponding
Gauss/Weyl transport dressing or a co-moving local reference. The separate local
V_x below already preserves each cell density. Those are distinct statements.

The reader is a source observable and controlled apparatus service. It is not
an assertion that the complete physical weak Hamiltonian equals a two-level
ladder, nor does equal-energy p/n conversion give beta decay.

## The full local interaction is a genuine source amendment

The L2 source contains the complete120-dimensional three-quark colour-singlet
orbital domain and the complete16-dimensional four-mode lepton CAR Fock space.
The quark ordering is(spin,weak), and the lepton ordering is(weak,spin). Define
Aq=bar_sigma tensor(tau/2), Al=(tau/2)tensor bar_sigma, with
bar_sigma=(I,-sigma), eta=(1,-1,-1,-1). At endpoint0 retain the complete
normal-ordered current

    V=sum_(a,mu) eta_mu [Jq²-dGamma_q(Aq²)
                         +Jl²-dGamma_l(Al²)+2Jq Jl].

All self, cross and neutral terms remain. Nl=0..4 are present. The free lepton
register generator is exactly0; it is not a massive Dirac or outgoing lepton
model. The new source is H_joint=H_chain tensor I+(1/64)V. Its weak-off limit
recovers H_chain tensor I exactly; [V,H_chain tensor I] is nonzero, so the new
interaction is not a decoupled spectator.

The exact representation uses raw symmetric monomials with metric
D_nn=prod occupation factorials. H_chain=H48/48 and V=W4/4 in that representation,
where both numerators are integer matrices. Thus

    H_joint_raw=(16(H48 tensor I)+3W4)/768,
    H_joint_normalized=D^(1/2)H_joint_raw D^(-1/2).

Weighted Hermiticity, local/total charge, lepton number and weak-off recovery
are checked by exact integer sparse arithmetic. Pad120 quark codes to128,
retaining128 extra joint codes with energy16 and zero weak interaction. The
source has2048 total code states. The reused pair-current norm is6, so the complete selected Fock domain has
||V||<=6 choose(7,2)=126. It is a bounded perturbation; adding it to a future
self-adjoint positive-field source preserves that source domain. The sharper
vacuum P/Q floors below do not require the occupied-sector bound to be negative.
Radical coefficients are determined solely
by integer entries and factorial weights. No numerical eigenvalue defines it.
The lazy native compiler uses exactly bounded dyadic radical intervals,
complete diagonal and pair rotations, a Lie-product bound, and explicit
control/export/reset errors. Its huge finite pulse budget is not an execution.
Its known-code preparation and complete occupation/clock terminal instrument
include every padded code, all pointer records and all failures.

Full2048-state source controls evolve an incoming neutrino with the neutral
quark source occupation, retain every output probability, and compare G=0.
The contact-only small block is selected by exact connected components, not by
choosing a receiver. It retains every complementary state in that invariant
contact block. It is not invariant under the colour hopping; the larger source
and its native budget remain necessary for that dynamics. These are finite
scattering/current-oscillation controls, not a lifetime or asymptotic decay law.

## Endpoint binding and field compatibility

The independently isolated ENDPOINT_MANIFEST.json is an immutable nested annex.
Its proof establishes on the invariant empty-lepton sector

    V_x=3S_x²-(3/4)N_x(N_x+2),
    P V_x P >= -9P,  Q V_x Q >= -6Q, [V_x,P]=0.

P includes both spin1/2 and spin3/2: the first bound is a LOWER bound, not a
negative uniform upper bound. The neutral endpoint spin-half trial has exact
expectation-9. Consequently adding G V improves the old P/Q energy gap by3G
without consuming the positive trial allowance. Occupied-lepton floors are
not obtained from this vacuum statement. The contact is not positive.

For rho_y=(2/3)n_u,y-(1/3)n_d,y-n_e,y, exact [V_x,rho_y]=0 at every cell.
The contact also commutes with colour Gauss generators and each quark cell
count. A future field dictionary based on this **total** density therefore
requires no extra Weyl kick for the contact; the particular Gauss field,
Coulomb density functional and free photon field commute with it. Dressed
hopping does not generally commute with it. This is compatibility, not a
claim that the external new colour/Maxwell producer has already implemented
or audited the joined operator. The frozen source for that next transfer must
be bound separately. The endpoint lemma can be consumed in its lepton-vacuum
sector without inventing a dynamical lepton field.

## Constructive remaining transfers

ROW16 still needs the existing six-quark/two-baryon sector, binding source and
quadrupole readout joined to this current. ROW18 needs the separately owned
strange/antiquark source and its open energy channels. ROW24 needs the existing
outgoing relaxation source and current spectral measure transferred to this
specific Hamiltonian; sharing G=1/64 does not make them equal. ROW25 needs both
an energy-releasing channel and the existing physical-position/velocity parity
statistic, with the Wilson/doubler control preserved. Fable's calibrated
first-emission source is distinct from this degenerate p/n chain and from a
full-Fock lepton ring. Fable-B's Euclidean/native bridge also needs its own
real-time source intertwiner. These mechanisms exist elsewhere; their exact
source/preparation/time/units transfers remain the obligations.
