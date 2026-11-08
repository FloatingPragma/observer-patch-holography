# Uniform source response on the joint clock

The uniform staggered response has a source-derived concentration estimate.
Its inverse-energy-weighted centre and energy-square variance both approach
zero on the prescribed joint clock. This concerns a volume-wide uniform
mode; it does not identify a nonzero-momentum dispersion branch, a particle,
or a limiting local field pole. The unweighted original preparation measure
and its complete complementary records remain distinct from the weighted
probabilities used in the proof.

## 1. Exact source moments

On a finite bipartite cubic box with N sites and edge set E, use precisely

    H_B=J sum_(xy in E)(S_x dot S_y-1/4)-B sum_x epsilon_x S_x^1,
    K=H_B-E_B, L=K/(Jp), p=n^-2, B=J n^-5.

The positive-field ground is unique and has zero transverse means. Write

    A=-sum_x S_x^2, M=sum_x epsilon_x S_x^3,
    O=sum_x epsilon_x S_x^1, m=<O>/N, q=M Omega.

The exact Ward identity from `WARD_DYNAMICS.md` gives

    K A Omega=iB q,   chi_K=<q,K^-1q>=Nm/(2B).          (1)

The same source double commutator supplies the first moment, without a
spin-wave or Gaussian replacement:

    h_K=<q,Kq>
       =-2J sum_(xy in E)<S_x^1S_y^1+S_x^2S_y^2>
         +(B/2)<O>
       <=J |E|+(B/2)Nm.                               (2)

Each transverse two-spin sum has norm 1/2. Since |E|<=3N, equations
(1)-(2) imply, for m>0,

    R_L=h_L/chi_L <= 6B/(J p^2 m)+B^2/(J^2 p^2)
                         =6/(mn)+n^-6.                (3)

Here h_L=h_K/(Jp) and chi_L=Jp chi_K. The source field momentum A is
not the generator of spatial translations. No momentum label is assigned
to the following scalar energy measures.

## 2. Two quantitative concentration statements

Let nu be the original spectral measure of q for L. On every finite
positive-field box its support is in [g_L,infinity), where
`g_L=B/(8Jp)=1/(8n^3)`. Define auxiliary probabilities

    dpi_1=E^-1 dnu/chi_L,
    dpi_3=E^-3 dnu/M_-3,  M_j=int E^j dnu.

They are mathematical reweightings of the same source measure, not a
claim that an experimental success probability can be divided out.
Equations (1)-(3) and Markov's inequality give

    pi_1(E>=omega) <= R_L/omega^2.                     (4)

Cauchy--Schwarz gives `M_-1^2<=M_1 M_-3`. Thus, with
`w^2=M_-1/M_-3`,

    0<w^2<=R_L,
    Var_pi3(E^2)=w^2(R_L-w^2)<=R_L^2/4.                (5)

Equation (5) is the actual source weighted-variance estimate left by the
uniform Ward identity. A uniform order lower bound m>=m0>0 makes its
centre O(1/n) in energy squared and variance O(1/n^2) on the joint clock.
Both concentrate at zero energy. Neither a finite positive rest energy nor
a nonzero-momentum shell follows. Constant-in-space modes over growing
boxes are not fixed compactly supported field test functions.

## 3. A finite order floor from the existing source construction

The ordered-state theorem proves maximal order greater than 1/10.
The infinite-volume positive-field order satisfies m_B>=1/10: the ground
energy density is concave in B, hence its negative right derivative is
nondecreasing in B; its B-down-to-zero value is that maximal selected
order. At fixed positive B the source gap and locality construction give
its unique local ground marginal and an effective finite-box approximation.
No convergence rate of m_B to its zero-field value is needed here.

Choose a source cube of radius r=384d, and a bulk whose sites are at least
d bonds from its boundary. Use the existing one-site marginal bound with
crossing-edge upper count `6(2r+1)^2` until its half-trace error is at most
1/40. The bulk signed order is then at least 3/40. The excluded fraction
is at most `6d/(2r+1)<=1/64`; an excluded site's signed order is at least
-1/2. Therefore the exact finite-source ground satisfies

    m >= (63/64)(3/40)-(1/64)(1/2)
      =169/2560 >1/20.                                (6)

The count-only doubling search in `uniform_mode.py` supplies such a box
for every n. The bound is large, but uses the actual fixed-B locality
modulus and source graph. It neither diagonalizes that large box nor
asserts that the eight-site control has this order. Equations
(3)-(6) yield the explicit uniform bound `R_L<=120/n+n^-6`.

## 4. Original preparation weight and retained outcomes

The contraction `M/(N/2)` gives a complete two-outcome preparation, with
its complementary Kraus operator retained. Follow it with the same-source
energy PVM, retaining every energy outcome. For any omega with
`omega^2>R_L`, (4) and the finite source gap give

    Pr(preparation and g_L<=E<=omega)
       >= m/(4N) (1-R_L/omega^2).                      (7)

Indeed nu on this interval is at least `g_L chi_L pi_1(interval)` and
`g_L chi_L=Nm/16`. At `omega^2=2(120/n+n^-6)` the finite-box floor is
`1/(160N)`. This is an original probability, including preparation failure;
no conditioning or normalized spectral count is used. Its N-dependence
means the estimate does not supply a nonvanishing original continuum
residue. It also does not prove that the actual residue vanishes.
Finite preparation, energy readout and clock errors must be added to (7)
under the existing declared coherent-service interface. No large pulse
tape is executed by this argument.

## 5. Finite spectral controls and the order of limits

The complete eight-site control uses the actual source at each specified
B=J n^-5 and p=n^-2, its own ground and the same scaled clock. Its local
charge preparation has both original branches and every even-cube
reflection outcome. These reflection characters are finite source
symmetries, not continuum translation momenta. Full-spectrum rational
certificates bound the energy intervals and original readout weights;
no solver output is accepted without those bounds.

The finite zero-field cube has a separately certified positive gap g0.
Since `||H_B-H_0||<=4B`, its positive-field gap is at least `g0-8B`.
On the joint clock every non-ground energy in this fixed box is therefore
at least

    n^2 g0-8/n^3.                                     (8)

This diverges as n increases. It is an exact finite-box obstruction to
using that box's B-down-to-zero table as the infinite-volume joint field.
It does not obstruct the growing-volume construction in (6) and does not
exclude a different original-weight spectral subbranch of that field.

## 6. Exact acceptance of the finite spectral table

The integer proposals in `uniform_mode_witness.npz` have 48 fractional
bits and retain all 256 vectors and eigenvalues at each stated n. They
are untrusted. The verifier reconstructs the actual integer source H_B.
For the proposed matrix V, diagonal D, and norm bound S, it bounds

    delta=||V*V-I||,
    eps=||H_B V-V D|| <= its exact Frobenius bound,
    radius=eps+2S delta.

The polar unitary W differs from V by at most delta. Therefore
`||H_B-W D W*||<=radius`, certifying all sorted eigenvalue intervals.
The separate ground residual and isolated first gap give the displayed
normalized ground-vector error e_g. All arithmetic after proposal is
integer or rational, with outward rational square-root bounds.

For an energy cutoff c disjoint from every certified excitation interval,
let d be its minimum distance from the proposed excitation levels, minus
their common error eta. Each W-column's leakage into the wrong exact
spectral half-line is at most eta/d. There are 256 columns; their
Frobenius bound and equality of the two projection ranks bound the
half-line projection error by `32 eta/d`. A finite energy bin uses both
cutoff errors. The approximate V weights add `2 delta+delta^2`, and the
prepared ground adds `2 e_g`. These are original probability errors,
not conditional errors after selecting a charge or reflection outcome.
Empty energy bins and ground-only bins have exactly zero created-charge
weight and remain allocated records.

For inverse moments, add the exact ground projector to K before inversion.
This changes no moment of the transverse created state. The corresponding
proposed ground eigenvalue is set to one. The positive operators differ
by at most `eta+2e_g+4delta`. If a is a common positive spectral floor,
the inverse and third-inverse operator errors are at most this difference
times `a^-2` and `3a^-4`, respectively. The same V/ground errors are then
added with their score norms. This certifies the weighted variance from
all original moments without omitting the spectral complement.

The uniform preparation's opposite charge measures are exactly equal:
the cube automorphism reversing the internal charge sends M to -M and
preserves the unique ground and H_B. They have the same positive weighted
variance in the finite table. Thus charge resolution alone does not
collapse this actual cube response to the dimer's one-energy example.
Every even-reflection label of the local two-branch preparation is also
retained. No continuum translation generator, finite-to-continuum rate,
or physical species is extracted from those cube labels.
