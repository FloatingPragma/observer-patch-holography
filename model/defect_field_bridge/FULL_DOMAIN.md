# A complete repair law after particle creation

The original colored native instrument acts on two occupied ports with at
most one particle per port. Particle creation and coherent propagation can
leave that domain. `full_domain.py` supplies an explicit extension on the
entire finite spatial CAR space, including all particle numbers and
same-port collisions. The occupation basis is the static quasiparticle
basis at each site. Bare-field operations must first be transformed to
this basis, or the repair maps transformed back; the two are not identified.

For a selected seam, count its two endpoint occupations. If both are zero
or one, apply the original C/E floor/ceiling rule with its fair recorded
coin. A particle moves with its internal state, and the full CAR ordering
fixes the sign across every occupied spectator. If either endpoint has
more than one particle, hold the state and write the added collision record
`raw=8`. There are sixteen padded raw slots per coin. All original zero
and abort slots, the added collision slot and every unused slot remain.

This is a chosen repair extension, not a conclusion about how the original
binary rule treated collisions. It uses the same action and amplitudes on
the old domain. Exhaustive tests compare every original colored C/E column
on all thirty native seams, for internal dimensions two and three.

Within a recorded branch, the resulting map is a signed partial injection.
A hold is its own inverse; a hop is inverted using its recorded direction,
seam and unchanged internal state. Each input occurs exactly once for
each fair coin. Therefore, on the complete Fock space,

\[
 K_h^\dagger K_h=\tfrac12 P_{D_h},\qquad
 \sum_h K_h^\dagger K_h=I.
\]

This proves complete positivity and trace preservation with any untouched
reference. Composition with any unitary propagation or particle-changing
operation on that same Fock space is consequently well-defined. Total
particle number and each total diagonal internal charge are conserved by
the repair itself; the supplied unitary can change particle number. Common
internal unitaries also commute with this repair rule, since endpoint
particle counts and identical internal transport are basis independent.
Site-dependent gauge covariance would additionally require link transport.

The fair-coin repair channel is also unital: every output has exactly two
preimages across the recorded branches. Thus the normalized identity is
a faithful stationary state for each fixed-seam repair channel. Its
modular flow is trivial. Arbitrary record-conditioned field operations or
feedback schedules need a separate stationarity argument.

The extra record has a physical consequence within this model. A coherent
superposition of vacuum and an opposite-charge pair at one site becomes
an incoherent mixture if the collision record is forgotten. The equal
superposition has output purity one half. Keeping every record gives the
corresponding classical-quantum state, not an undisturbed pure field state.
An inverse pair-drive echo can detect this disturbance. Completeness of
the law does not establish compatibility with measured particle dynamics.

The source witnesses execute all input columns for a two-site/two-internal-
state collision case and a three-site/one-internal-state crossing case,
for both C and E. A source-prepared fair coin, a full signed Gray
permutation and all four record bits give a dilation on the padded
512-address register. Executed source matrices retain numerical residuals;
the code admits a Frobenius error below `2e-7`, rather than reporting a
floating calculation as an interval proof. Proper-code export complements
have zero probability on the exact invariant code and are explicitly
counted. The construction on 1,080 modes is an exact algorithmic extension;
its complete source dilation is not claimed to have been executed.

```sh
python3 -m model.defect_field_bridge.full_domain --check
python3 -m pytest -q -W error model/defect_field_bridge/test_full_domain.py
```
