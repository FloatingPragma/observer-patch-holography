# Joint selection remains a separate physical input

The responding quartic source fixes an entropy-area obstruction. It does not
lock an independent matter excitation to a geometry response. The same-graph
Gaussian family currently admits a particle/hole mixture at fixed coframe;
its entropy changes by the positive modular coefficient measured in the
curved-ball packet while its area does not change. Therefore this family does
not yet satisfy generalized-entropy stationarity.

## A necessary restriction on admissible variations

For a faithful finite joint state rho, suppose a fixed-volume state tangent X
can vary the matter state within one unchanged geometry block. Then

    delta[S(rho)+eta Tr(rho A)] = Tr[X(-log rho)]

because `Tr X=0`, `delta A=0` and `delta V=0`. If every traceless Hermitian
tangent in that block is admissible, stationarity forces `-log rho` to be
scalar on the block, hence the conditional matter state maximally mixed.
That cannot also have a nontrivial geometrically normalized modular energy.

This is a statement about unrestricted tangent spaces, not a no-go for
constrained gravity. A successful joint source must specify why the independent
matter-only tangent is not an allowed physical variation, or how an actual
source process necessarily supplies the correlated geometry change. A
stationary density under time evolution does not establish this variational
restriction. Nor does the Gaussian first law.

## Existing joint constructions checked

1. `coherent_defect_wave` supplies an autonomous reciprocal link/wave unitary,
   complete readouts and a faithful stationary reference. Its links encode
   protected finite fluxes. Its own proof states that it has no geometric
   field action, spatial refinement or physical energy map. Therefore it is
   reusable as interaction architecture, but it does not attach the needed
   geometric area, bulk stress, or admissibility constraint.
2. `reciprocal_area_dynamics` supplies genuine thermal matter/geometry feedback,
   positive area coefficient, memory and complete bath records. It is the
   strongest existing joint-selection example found here. `replay.py`
   independently reconstructs all sixteen energy levels and the transition
   from its action matrices, with 512 complete ideal record labels.
3. `particle_geometric_response` supplies gauge-covariant coframe strain and
   work operators. The coframe is declared; no source law selects it.
4. `higgs_curvature_response` supplies classical identities for an explicitly
   declared curved action, including `F(X)R`. It does not provide a microscopic
   entropy state. Using that action to assert the missing Einstein response
   would insert the requested gravity dynamics.

## What the actual thermal source does

The existing source has two geometry blocks at r=.9,1.1, both proper volume
one, with area `20/9,20/11`. In dimensionless thermal energies,

    M(r)=[r/4+1/(7r)+3r^2/28] D -(XXI+IXX)/(7r),
    D=diag(2 popcount(mask)+3).

Its Gibbs weights and complete cross-geometry transition reproduce
`eta in (.578,.579)`, memory coefficient `k0=.02170170003...`, and initial
area response `-.0006270566...` under the published field randomization.

The trace-preserving tangent that transfers weight from the lowest to highest
energy in geometry zero has exactly zero area and volume variation. Its
entropy derivative equals the nonzero thermal energy difference. Hence
`S+eta A` is not stationary even on this interacting source's admitted state
space. The required correction is exactly `-beta E_ex`, where E_ex is matter
energy above the conditional thermal energy. The independently executed
first five updates retain that term and reproduce the complete heat/entropy
balance. Dropping it to claim the Einstein premise changes the result.

The source's restricted conditional-Gibbs family does have an entropy-area
maximum. Its actual dynamics does not preserve that family: k0 is nonzero.
Consequently there is no justified transfer of that restricted maximum to
arbitrary local excitations of the quartic Dirac vacuum. A new local joint
preparation/constraint/dynamics construction remains necessary.

## Strength and one common clock

The quartic candidate has approximately `eta=c/Rphys^2` with c near .073,
so the conditional Einstein relation would be `G=Rphys^2/(4c)` in natural
units. A correct sign and a nonzero area response do not fix Rphys.
Its length must follow from the same physical clock, energy anchor and limiting
speed as the matter source, independently of measured G.

`UNIT_DIAGNOSTIC.json` conditionally identifies the coordinate clock with the
separate 432-register electron anchor, U=8.4741660083 GeV. That identification
has not been proved: the register coefficient is a compiled mass slot, not
automatically this regulator's physical energy or elapsed time. Without the
transfer, no physical G comparison can be assigned to the actual candidate.

If the common clock and source speed one are assumed, length is fixed as
`hbar c/U`; it is not a second adjustable parameter. For a source area
coefficient eta in coordinate units the inferred value is
`G_model=1/(4 eta U^2)`. The script computes its discrepancy using the
repository's external G comparison input; it does not choose R to fit G.
It also reports the enormous required ultraviolet/matter hierarchy as a
counterfactual requirement, not as a solution. The old source's exact
length/clock elimination theorem is reused: a length-only rescaling cannot
match gravity while keeping both the energy anchor and limiting speed fixed.

## Next constructive task

The current viable object is the responding quartic field, not a completed
joint source. The next architecture must provide a physical restricted
matter/geometry state family or a local reciprocal dynamics whose actual
prepared/evolved records enforce the common-ball stationarity. It must retain
the source's regulator entropy, bulk energy, reference and clock, and account
for all feedback/preparation work. A rule defined by solving Poisson with a
chosen G, or by resetting each ball independently, does not derive this link.

No such local source with these properties was found in the audited existing
packages. The positive thermal response and its exact failure are preserved
for reuse. This is the precise outstanding construction, not a proof that it
cannot exist.
