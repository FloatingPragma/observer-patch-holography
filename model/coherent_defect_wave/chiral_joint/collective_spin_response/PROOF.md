# An interacting spin line with complete finite records

## Source, reference and intervention

The source is the full side-two Peter–Weyl carrier with thirty native
fermion modes and five scalar states at each of eight cells, and twenty-four
119-state links. Its dimension is `D=(5*2^30)^8*119^24`. The physical
subspace satisfies all twelve Gauss constraints per cell. The Hamiltonian
is `H_lambda=H0+lambda T`, with all onsite Yukawa, scalar and link
coefficients one, `lambda=2^-34`, and source clock `exp(-it H_lambda)`.
The [ground theorem](../atomic_spin_geometry/PROOF.md) proves a unique
physical ground `g_lambda`, zero ground momentum, and a first isolated
rank-three excitation at the three X momenta. This construction uses that
ground as its stationary reference. It does not transfer to unit hopping.

For axis a define the actual full-carrier current

`S_x,a = c_x^dagger (sigma_a/2 tensor I_15) c_x`,
`O_a = sum_x (-1)^x_a S_x,a`.

The one-particle matrix of O_a has 120 eigenvalues +1/2 and 120 eigenvalues
-1/2. Its full fermionic second quantization has norm 60, including every
quark and lepton occupation. The operators commute with the local gauge
generators and preserve fermion parity. They commute with the onsite
Hamiltonian, preserve its ground doublet, and compress to the actual
spin-half current of that doublet. These facts are verified with the
native CAR polynomials and the frozen local ground representation.

The apparatus retains a four-dimensional reference: a scalar vacuum and
three vector states |a>. The reference clock is identity, and all its
states have even fermion parity. Its translations act on |a> with the
inverse X_a character. Its rotations are the actual proper cubic vector
representation. The successful map out of the reference vacuum is

`B=sum_a O_a tensor |a>/(60 sqrt(3))`.

Consequently `B^dagger B=sum_a O_a^2/10800<=I`. The complete finite
apparatus is the Halmos unitary

`[[sqrt(I-B^dagger B), -B^dagger], [B, sqrt(I-B B^dagger)]]`.

It acts on the direct sum of the scalar and vector reference sectors.
Both square roots are exact finite algebraic services. The failure map
is `sqrt(I-B^dagger B)`. Gauge symmetry, parity, joint translations and
proper cubic rotations are preserved. Reference compensation makes joint
momentum zero; each successful source branch itself has momentum X_a.
The code checks all 24 cubic rotations, all three Fourier components,
and actual graded antiperiodic translations on occupied full-source
columns with nonblank scalar and link memory. The general translation
identity follows by relabeling the complete local sum; the seam phase
cancels within each even fermion bilinear.

## Exact residue in the atomic comparison

The canonical graph theorem compares the source to the full 256-state
atomic ground space. Its second-order generator is
`-24 A_edge I+2 J_edge U H_cube U^dagger`; A_edge and J_edge are the
complete charged pair-resolvent coefficients, with the scalar term
retained. The coordinate change U obeys

`U^dagger O_a U = N_a = sum_x (-1)^(x1+x2+x3) S_x,a`.

The code verifies this identity as a complete 256-dimensional matrix,
along with actual source translations. In cube coordinates let
`x=4 epsilon0` be the least real root of
`f(x)=x^3+20x^2-16x-576`. The magnetic-zero block of `4 H_cube` has
dimension 70. The kernel of f on this block has dimension three, and the
first-triplet magnetic-zero eigenspace at -16 has dimension one.

Writing A for the restriction of `4 H_cube` to the cubic kernel, its
ground projector in the algebraic field is

`P(x)=[A^2+(20+x)A+(x^2+20x-16)I]/f'(x)`.

This identity follows from f(A)=0 and the three simple roots of f. Exact
Gram matrices, the complete triplet vector and N_z give the three moments
of the line-residue operator as `(23,-436,8368)/3`, and the three moments
of the total-residue operator as `(32,-448,8704)/3`. Thus the unnormalized
one-component line residue and total squared norm are respectively

`L(x)=(23x^2+24x-720)/(3(3x^2+40x-16))`,
`S(x)=32(x^2+6x-24)/(3(3x^2+40x-16))`.

The cube ground is a singlet. Its three component weights are equal by
the cube spin symmetry used only in this exact comparison. The covariant
pump therefore has total line probability `L/3600`, total success
probability `S/3600`, and outside probability `(S-L)/3600`. Rational
root isolation and interval arithmetic prove each is strictly positive.
Approximately L=7.48730716 and S=7.54675104. No floating eigensolver is
needed to calculate or bound these residues.

## Full-H transfer of the observable

Use the frozen full-source atomic gap g, hopping norm tau=360, exact
exchange interval and even-remainder theorem. Put `v=lambda tau`,
`r=2v/(g-v)` and `h=r+r^2/2`. If I_0 denotes the atomic embedding and
V_lambda the canonical graph isometry, its normalization gives

`||V_lambda-I_0||<=h`.

Indeed the graph contribution has norm at most r and
`||(1+X^dagger X)^(-1/2)-1||<=r^2/2`. Every complementary source state
is included in X. The effective generator obeys

`||K_lambda-E_atomic I-lambda^2 K2||<=R4`,
`R4=25(lambda tau)^4/g^3`.

Set `eps=R4/lambda^2`, let d_g be a lower bound on the ground gap of K2,
and let d_star be the smaller separation of its first triplet from both
neighboring spectral sets. The code uses `d_g=2 J_lower delta_lower`
and `d_star=2 J_lower min(delta_lower, next_lower)`. Here
`eps<min(d_g,d_star)/4`.

For the perturbed unit ground vector q, the residual equation and the
min-max eigenvalue bound give
`||(I-P_g)q||<=eps/(d_g-eps)`. Since both ground projectors have rank
one, their operator norm difference is at most this value. A phase can
be chosen so the two ground vectors differ by at most twice this value.

For the perturbed rank-three first-cluster projector Q, all three
eigenvalues lie within eps of the unperturbed triplet energy E. Hence
`||(K2-E)Q||<=2eps`. On the complement of the unperturbed triplet P_t,
`|K2-E|>=d_star`, so
`||(I-P_t)Q||<=2eps/d_star`. Equal finite ranks imply
`||P_t-Q||=||(I-P_t)Q||`. This residual proof requires no continuous
spin symmetry of the full source.

Embedding both bounds into the full carrier gives

`eta_g=h+2eps/(d_g-eps)`,
`eta_t=2h+2eps/d_star`.

Because B is a contraction, line-record and outside-record amplitudes
differ from the atomic comparison by at most `eta_g+eta_t`. Total-success
and failure amplitudes differ by at most eta_g. For an atomic probability
in [a,b] and the corresponding amplitude error eta, the transferred
probability lies in

`[max(0,sqrt(a)-eta)^2, min(1,sqrt(b)+eta)^2]`.

All square-root endpoints are enclosed with outward rational arithmetic.
The reported source coupling satisfies the stronger detector condition
`(eta_g+eta_t)^2<=min(a_line,a_outside)/16`. It yields, conditioned on
successful ground preparation, approximately

| Record | Lower probability | Upper probability |
| --- | ---: | ---: |
| First excitation line | 0.00204894419 | 0.00211090162 |
| Entire outside complement | 0.00001386725 | 0.00001938783 |
| Pump success | 0.00208084364 | 0.00211185317 |
| Pump failure | 0.99788814683 | 0.99791915636 |

The excitation energy lies between `8.20758426381e-22` and
`8.21034719056e-22` source units. Its exact degeneracy and the equal
one-third axis weights follow from the actual source cubic symmetry,
unique ground, and covariant pump. The receipt stores outward rational
bounds; printed decimal endpoints in this proof are explanatory.

## Complete preparation and readout law

Prepare `I/D` on the complete kinematic carrier and the scalar reference
vacuum. Apply the commuting joint physical-Gauss and unique-ground-energy
projector. It succeeds with probability 1/D and fails with probability
1-1/D. Failure retains the source and reference and appends a terminal
failure record; subsequent intervention and readout are controlled off.
This initial filter and the energy readout are declared finite spectral
services. Their probability laws do not assert an executed or efficient
native gate schedule.

On ground success, apply the complete pump unitary and measure its
orthogonal scalar/vector sectors. On pump success, evolve under the same
H_lambda for any real time, and measure the physical spectral projector
of the first rank-three cluster versus its full complement, together
with actual source and reference momenta. These energy and momentum
projectors commute with H_lambda, so their probabilities hold at every
source time. No outside state is discarded. Pump failure has source
momentum zero because its effect commutes with source translations.

There are eight positive coarse records: initial preparation failure,
pump failure after preparation, and line/outside for each of three
reference axes. Each successful axis has probabilities `p_line/(3D)`
and `p_outside/(3D)`. Pump failure has probability
`(1-p_line-p_outside)/D`. All other flag/phase/reference outcomes are
retained with probability zero. Their exact sum is

`(1-1/D)+(p_line+p_outside+1-p_line-p_outside)/D=1`.

The interval endpoints are individual marginal enclosures, not free
independent parameters; the exact probabilities obey this identity.
The result is an actual finite interacting collective-spin response.
It supplies no zero-momentum rest mass, relativistic dispersion,
thermodynamic limit, empirical generation assignment or GeV conversion.
