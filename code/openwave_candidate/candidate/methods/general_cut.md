# General-cut charge and the apparatus contract

## Exact result on the canonical ring

`oph_openwave.general_cut` uses the actual portable `flat_ring(1,1)` geometry: 24 ports,
63 seams, and 87 port/transit vertices. Its oriented integer half-incidence
matrix A has shape 87 × 126. Each column is +1 at the port and -1 at the transit
vertex. Let s select any region, n be the occupation vector, and G=n+AF. Then

    Q_S = s n = s G - sum_h (s_source(h)-s_target(h)) F_h.

Symbolic coefficients prove that every internal half-edge cancels, for an
arbitrary region and arbitrary records. This is an operator identity for the
diagonal record/occupation observables, so it extends to any density matrix
supported in the declared Hilbert space. No classical path is inferred. With
only modular registers specified, the corresponding identity holds modulo N.
The complete-boundary cut need not be a unique bridge.

The result is source-bound to the actual ring geometry and construction files.
There is no quantum evolution or enumeration of the exponentially large modular
ring sector. Three rejecting tests cover orientation, altered current,
arbitrary signed records, and modular aliasing. The tiny algebraic controls
consume no reserved acceptance configuration.

## What integer information can be retained

Two distinct ambiguities must be separated.

An integer circulation can wind N times around a triangle inside one carrier. Its
half-record vector is N times a nonzero cycle z with Az=0. It has the same
modular records, same position and same integer Gauss charge as zero winding.
No apparatus reading only the modular state can distinguish those integer
seam windings. Appending a classical history of native commit outcomes does
not recover coherent windings: coherent hops are not measured trajectories.
One must add an actual integer winding register, an intervention measuring the
relevant current, or a proved quantum-state extension. Such additions can
change interference and require their own dynamics.

However, the net current through a COMPLETE boundary is not the circulation
winding. For a one-load state known to arise from zero integer records in a
fixed integer Gauss sector, it is exactly the initial regional occupation minus
the final regional occupation, in {-1,0,1}. Integer cycle ambiguities cancel
on the complete boundary. This is an endpoint charge identity, not a reconstructed
history. A modular source can therefore support this restricted integer net-cut
readout if its fixed integer-lift ancestry is independently justified.

An unrestricted modular sector also aliases different integer Gauss sectors:
adding N to one half-record changes integer G at its endpoints while leaving
all residues unchanged. The exact control constructs this case. Pure modular
Gauss constancy does not justify selecting one integer sector or unbounded
integer charge on arbitrary states. Modulus, occupancy restrictions and the
integer-lift ancestry belong in every apparatus claim.

## Minimal row-1 physical apparatus contract

The following is a proposed contract; the current one-load package does not
implement its signed-particle and pair sectors.

| Part | Required definition and observable |
|---|---|
| Charge | Specify a calibrated elementary charge unit e and the operational regional measurement Q_S/e, with its relation to signed occupation and complete-boundary record flux. Do not identify fractional expectations with fractional individual outcomes. |
| Signed species | Admit positive and negative unit excitations with q=n_+−n_−. State the statistics and occupancy restrictions relevant to the test without assuming spin or rest mass results. |
| Charge preparation | Name the admissible vacuum, one-charge and multi-charge sectors, record arithmetic, any integer-lift ancestry and the apparatus's reference/background charge. |
| Local propagation | For every allowed coherent and committed move prove the signed Gauss/cut identity and retain complete outcome probabilities. The coherent generator must operate on the same charge source as the measurement. |
| Pair creation | Define a local operation producing + and − with zero net charge; if the members separate, include the compensating record/current update. A missing-record mutation must be rejected. Charge conservation alone is not a pair-creation mechanism. |
| Annihilation | Define the inverse charge-removal channel and all retained outcomes. Radiation identification is a separate row; row 1 needs its charge ledger to close. |
| Interventions | Specify region cuts, attempted separation/splitting, external apparatus operations, allowed local record changes and preparations. Test the whole admitted charge domain, rather than only a one-particle basis that excludes all other outcomes. |
| Literal exclusion | Show that no admissible intervention produces two separated excitations with individual charge ±e/2 while preserving the stated topological/Gauss sector. Identify the relevant topological or quantized observable on the actual construction. |

A declared integer charge lattice may provide a conditional exact theorem. A
claim that the carrier derives that lattice or reproduces electric charge needs
the source-to-apparatus identification. The contract does not demand electron
mass, magnetic moment or spin, which belong to other criteria.

