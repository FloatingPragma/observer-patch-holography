# Complete planar wave observations on a native lapse profile

## Source and domain

The exact dyadic q profile is produced by `openwave_energy_lapse.prepare`
on the full31³ or47³ Dirichlet lattice with one central native level1
source. The complete native8-level energy matrix and its ground subtraction
are retained in that producer. This experiment declares a planar scalar
extension using q(x,y,z=center); it does not identify this plane with an
invariant subspace of a three-dimensional field. The q source is static
during the wave observation and distinct from the two-cell autonomous
quantum-geometry source. No transfer between those two states is assumed.

For the planar Dirichlet incidence matrix D, with outside q=0, let q_e
be the average at an edge's endpoints. For γ in {−1,0,1}, put

    S_x=exp((1+3γ)q_x/2), A_e=exp((1−γ)q_e),
    K=S D* A D S, H=0 direct-sum sqrt(K).

This is the canonical scalar wave equation for the declared discrete
static metric weights. There is no transverse-polarization or Maxwell
transfer claim. The quantum field carrier consists of a vacuum and all
n² one-quantum modes. K is positive definite, H is Hermitian, and U(t)
is exactly unitary with conserved source energy. No mode or detector
selection enters the dynamics. Adding a spectator reference gives
U(t) tensor I without changing these statements.

## Complete apparatus and work

Prepare the equal mixture of vacuum and the normalized full-grid vector

    ψ(x,y) ∝ exp(−((x−(c−10))²+(y−(c+6))²)/18 + 0.9 i x), c=floor(n/2).

The curved and flat experiments use this same vector and source time.
Every site and the vacuum are retained in an orthogonal position instrument.
Probabilities are 1/2 and |Uψ|²/2. Record copying is an isometry on these
orthogonal pointer values and does not clone the incident quantum state.

For Y=diag(y), the velocity operator is V=i[H,Y] on the excited sector.
Since ||H||≤3 and ||Y−(n−1)I/2||=(n−1)/2,
||V||≤B=3(n−1). The complete velocity POVM consists of the vacuum projector
and (I±V/B)/2 on the excited sector. Their positive square roots define a
Lüders instrument, including arbitrary references. On the prepared mixture
the probabilities are 1/2 and (1±〈V〉/B)/4. The velocity formula
〈V〉=2 Im〈ψ,YHψ〉 follows directly by differentiating 〈Y〉. These velocity
and position instruments describe separate experiments. Both are
trace-preserving, and energy before or after any such measurement lies
in [0,3], so the absolute apparatus work is at most3 source-energy units.
Conservation during source evolution does not imply zero measurement work.

## Finite computation and error scope

For n∈{31,47}, |q|<0.001 is checked. The Dirichlet square-lattice Laplacian
has λ_min≥8/(n+1)² and λ_max≤8. The S,A weights give
spec(K)⊂[7/(n+1)²,9] for all three γ choices. The executable applies
degree4096 Chebyshev functions of the entire sparse K, using16384 midpoint
DCT samples; it retains every mode.

The Bernstein ellipse parameter ρ=1+1/(2(n+1)) lies in the analytic
positive-square-root domain. Writing s=n+1, its image under the spectral
affine map obeys Re(z)≥103/(16s²)>0 and |Im(z)|≤9/(4s). Thus
|Im sqrt(z)|≤9/(2√103)<9/20. The analytic function bounds are
M=exp(9|t|/20) for exp(−it sqrt(z)) and M=4 for sqrt(z). The discarded
series and midpoint-DCT aliasing obey the conservative operator bound

    2 M ρ^(-4096)/(ρ−1) + 8 M ρ^(-2·16384+4096)/(ρ−1).

This is a mathematical truncation/aliasing bound, not a roundoff enclosure.
The printed decimal is a numerical rendering of that exact formula.
Floating DCT, sparse products and exponential arithmetic have separate
roundoff. The experiment independently evolves the complete second-order
local equation with matrix [[0,I],[−K,0]], initial data (ψ,−i sqrt(K)ψ),
and sparse exponential action. The agreement is a numerical diagnostic,
not a rigorous bound on all arithmetic error. Exact instrument positivity,
normalization and reference extension follow from the source theorem;
record numbers are floating evaluations with explicit residuals.

## Deflection, controls and the ray relation

At n=31, γ=1 and t=16, the conditional transverse displacement relative
to flat propagation is approximately −1.9417333·10⁻⁵ sites and the
transverse velocity change is −3.6701444·10⁻⁶ sites per source time.
Both the full record distribution and the velocity POVM record this
response. γ=0 produces approximately half the response. The side47
control preserves its sign and comparable magnitude. Reflecting the
incident packet across the source reflects the deflection. At q=0 the
three source matrices are exactly identical.

For locally constant q, the dispersion is

    ω(x,p)=exp((1+γ)q(x)) sqrt(4sin²(px/2)+4sin²(py/2)).

Its continuum principal speed is exp((1+γ)q), giving ray index
exp(−(1+γ)q). At the packet's p_x=0.9,p_y=0, the lattice group speed is
cos(0.45). To first order in q, along the straight central trajectory,

    δv_y(t)=−(1+γ) integral_0^t ∂_y q(x0+v_x s,y0) ds,
    δy(t)=−(1+γ) integral_0^t (t−s)∂_y q(x0+v_x s,y0) ds.

The producer compares both v_x=1 continuum and v_x=cos(0.45) lattice
rays with the full finite-wave records, using an explicitly declared
bicubic interpolation of the same q slice. These rays discard packet
width, diffraction, boundaries and variable-coefficient lower-order
terms. No small wave/ray difference theorem is supplied. In particular
γ=−1 gives exactly zero principal-ray deflection while the full finite
wave produces a small nonzero residual. This mismatch is retained.

The construction supplies finite planar wave bending on a source-derived
static lapse profile. γ, geometry stiffness, scalar extension and clock
units are declared. It does not derive physical G, a solar lensing angle,
the γ=1 selection, an observed equivalence principle, a continuum Einstein
closure, three-dimensional photon dynamics or a quantum-geometry join.
