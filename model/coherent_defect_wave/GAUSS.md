# Interior gauge constraint for the joint defect and wave program

The [joint program](PROOF.md) fixes the boundary and exterior links of the
five-face region around carrier vertex 5. Its ten link words contain one
defect in that region and its compensating exterior defect. The permitted
interior gauge transformation flips the five links incident on vertex 5.
This transformation preserves every face holonomy and the fixed boundary.
Gauge transformations that change the boundary belong to a different domain.

Let `s` be that five-link bit mask and let `J|g>=|g xor s>`. The ten-word
class is closed under J; J has five two-element orbits and no fixed word.
On the wave register define `D` as minus one on the doublet at vertex 5
and plus one on every other basis vector. Thus `D` is the central code
frame at the interior vertex, including its matter charge. The joint gauge
operator is

`Gamma = J tensor D`, `Gamma*=Gamma`, `Gamma^2=I`.

This is an interior Gauss constraint on the actual link and matter registers.
The gauge-invariant subspace is `H_phys=ker(Gamma-I)`. The complementary
sector is not discarded in the preparation instrument.

## Preservation by the same joint evolution

For `B=direct_sum_g W_g`, the executed central-link law obeys

`W_(g xor s)=D W_g D`.

Indeed, changing the frame at vertex 5 reverses the central sign of exactly
its incident links; the address reflection commutes with D. Consequently
`B Gamma=Gamma B`.

Each local defect permutation `T_e` depends only on the two adjacent face
fluxes. Those fluxes are unchanged by J, and bit flips commute, giving
`T_e J=J T_e`. The wave occupation projector `Pi_e` commutes with D. Hence
every controlled repair

`C_e=I+(T_e-I) tensor Pi_e`

commutes with Gamma. Their fixed ordered product, including B, satisfies

`U Gamma=Gamma U`.

The same evolution therefore preserves the interior Gauss constraint, odd
regional flux and the declared internal charge. No averaging over physical
histories or modification of the wave clock is needed for this statement.

## Complete preparation and readout

The two Kraus operators `P_+=(I+Gamma)/2` and `P_-=(I-Gamma)/2` form a
complete projective instrument. Both outcomes are retained. For any normalized
wave state psi and definite link word g, each outcome has probability one
half; the normalized positive branch is

`(|g> tensor psi + |g xor s> tensor D psi)/sqrt(2)`.

Gamma has zero trace, so `dim(H_phys)=72000/2=36000`. The success branch is
stationary under U only when its input preparation makes it so; the projection
alone does not turn an arbitrary input into a stationary state.

The [same-U heralded reference](PREPARATION.md) commutes with Gamma because
it is a polynomial in U and U*. Restricting it with P_+ and normalizing
therefore gives a faithful stationary state on H_phys. For this ordered
repair cycle, no conditional link permutation sends g to `g xor s`. Thus
`Tr(Gamma U)=0`, `Tr(U|H_phys)=-23925`, and its reference is

`rho_phys=[(5/4)I_phys+(U_phys+U_phys*)/2]/21075`.

Its eigenvalues are at least `1/84300`. From the maximally mixed physical
input the complete heralded filter succeeds with probability `281/1080`.
These values use the constrained dimension and trace explicitly. A complete
way to prepare that input is to apply the Gauss instrument to the maximally
mixed ambient state, retaining both outcomes, then filter the positive branch.
The positive Gauss outcome has probability one half, so the joint success
probability per such fresh ambient trial is `281/2160`.

Face-syndrome readouts and wave register projectors commute with Gamma.
For a definite-g input, its Gamma transform has a distinct link word. One U
maps these two words into disjoint gauge partners for each output wave
coordinate. The interference term vanishes for either of these diagonal
readouts. Gauge covariance makes the two diagonal terms equal. Consequently
the normalized positive preparation has exactly the same one-tick
face-syndrome and wave-register probabilities as the corresponding
unprojected input. This particular equality is not asserted for an arbitrary
long protocol: later dynamics can make the gauge-paired amplitudes overlap.

This transfers the joint program's reciprocal readout tests to the constrained
sector. Entanglement across the chosen link/wave tensor split is an additional
register diagnostic; the gauge-invariant marginal response tests establish
interaction without relying on that split as a physical tensor factorization.

## Gauge-invariant recombination readout

The complete apparatus `covariant_preparation_return` in [gauss.py](gauss.py)
tests coherence without leaving the positive constraint sector. For each
two-element orbit choose the smaller link word c. For a normalized owned
wave vector `psi=A(|5,3>+|6,3>)/sqrt(2)`, set

`P_c=|psi><psi|`, `P_(c xor s)=D P_c D`,
`Pi=sum_g |g><g| tensor P_g`.

The implemented wrapper applies D on noncanonical link words, performs the
owned wave-return instrument, then restores D. This is coherent conditional
control; it does not measure the link word. Its two Kraus operators are
`Pi` and `I-Pi`. They are orthogonal projections summing to identity, and
`P_(g xor s)=D P_g D` proves `[Pi,Gamma]=0`. Thus both actual outcomes
preserve either Gauss sector. The chosen frame and return apparatus are
declared services, not consequences of idle dynamics alone.

Prepare the normalized positive projection of
`|1> tensor W_1* psi`, evolve for one U, and undo just the conditional
transport product `C=C4 C3 C2 C1 C0`. Without an intervening measurement,
the source is the positive projection of `|1> tensor psi`, in the range
of Pi. The complete return probabilities are therefore `(1,0)`.

Alternatively, insert the complete five-outcome regional-face instrument
between U and `C*`. It separates the two orthogonal face branches associated
with wave ports 5 and 6. Each has probability one half, and each has squared
overlap one half with its appropriate conditional-frame return vector.
The ten face/return histories therefore contain four probabilities `1/4`
and six zeros; summing the complete records gives return/complement
probabilities `(1/2,1/2)`. Every branch remains in the positive Gauss sector.
The [executed receipt](results.json) records these values under
`interaction.constrained_coherence` and checks the operator commutation on
an input that is not already gauge invariant.

This is an invariant coherence and recording-disturbance witness for the
specified finite interaction. The return service does not change U, and
the argument asserts no arbitrary effect compiler or physical particle.

## Scope

This proves one interior central gauge constraint with fixed external boundary
data. It does not impose all independent gauge constraints of a growing
federation, derive continuous Standard Model gauge dynamics, supply additional
matter occupations, or identify relativistic particles and rest masses.
