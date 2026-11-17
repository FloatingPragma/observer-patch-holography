# What the fixed-lattice test excludes, and a constructive alternative

This packet independently audits the F087 claim. It does not establish an
Einstein realization, Newton's constant, or any full OpenWave gravity row.
The consumed producer files are copied by content hash into `inputs/`.

## 1. Exact spectral response of the existing source

Write the producer's Hermitian generator as

    H(k) = i[U(k)-U(k)†]/2,
    U(k) = exp(-ik_z alpha_z) exp(-ik_y alpha_y)
           exp(-ik_x alpha_x) exp(-i mu beta).

Its negative-energy projector is the reference. A lapse Fourier component has
vertex `V(k+q,k)=[H(k+q)+H(k)]/2`. For an empty state b and a filled state a,

    <b,k+q|V|a,k>
      = [E_b(k+q)+E_a(k)] <b,k+q|a,k>/2.

At a momentum with a filled/empty gap, differentiating the eigenvalue equation
gives the first-order vertex coefficient

    M_ba = [(E_b+E_a)/(2(E_b-E_a))] <b|partial_k H|a>.

This equation is basis-independent when summed over degenerate eigenspaces.
For a fixed gapped quadrature the static susceptibility has

    chi(q) = c2 q^2 + O(q^3),
    c2 = (2/N) sum_(k,b empty,a filled) |M_ba|^2/(E_b-E_a) >= 0.

Symmetric momentum integration removes odd terms. It is not legitimate to
replace `c2` by zero merely because a quartic term dominates the momenta sampled.
Opposite partners with `E_b=-E_a` contribute zero; other interband transitions
need not. `audit.py` uses the exact derivative of the ordered matrix product,
independent of the producer's response routine.

At `mu=.05` it gives `c2=2.7111070084e-5, 2.4884518180e-5,
2.4283655695e-5` on 12, 20, 28 cubed half-twisted quadratures. The corresponding
massless values are below `1.7e-32`. At fixed 20 cubed quadrature,
`chi(q)/q^2=2.5012104380e-5` for `q=.002`, tending toward the analytic
`2.4884518180e-5`. These non-grid q values are continuous Bloch quadrature,
not allowed Fourier modes of that finite torus. No certified thermodynamic
quadrature error is asserted.

The allowed L=32 Fourier mode independently reproduces the producer's
`chi=3.6121233730276285e-5` with absolute difference `2.1e-20`.
Thus the producer's data are reproduced; the correction concerns the exact-zero
interpretation. This energy susceptibility is not the required area coefficient.
It gives no sign, normalization, or Einstein identification for an entropy
equilibrium equation.

### Exact positivity witness in the same source family

`exact_witness.py` uses `sin(mu)=1/20`, `cos(mu)=sqrt(399)/20` and
`k=(pi/4,pi/4,pi/4)`. This mass is `asin(1/20)`, not exactly `.05`.
Put `D=partial_(k_x)H`, `K=HD+DH`, and `S=sign(H)`. Exact arithmetic gives

    H^4 - (1199/800) H^2 + (319201/640000) I = 0,
    tr(K^2) = 1/400,
    tr(H K H K) = 599/64000000,
    tr(H^3 K H K) = 1037401/102400000000,
    tr(H^3 K H^3 K) = 446561/40960000000.

The positive and negative eigenvalues are separated from zero at this momentum.
Writing `d=sqrt(319201)/800` and `A=1199/800+d`,

    S = (A H-H^3)/(d sqrt(1199/800+2d)).

Consequently

    I_inter = tr(K^2)-tr(S K S K)
            = 80199/32240000
              + 79401 sqrt(319201)/10291040240000 > 0.

This is four times the squared filled-to-empty block norm of K. Since
`||H||<=1`, every transition gap is at most 2, and the local susceptibility
coefficient is at least `I_inter/64 > 0`. Continuity gives a positive-measure
momentum neighborhood with positive coefficient. An exact absence of a q²
term throughout this finite-mass source family is therefore false. This
witness neither proves a finite integrated upper bound nor matches Newton G.

## 2. The omitted fixed-proper-volume variation

The producer's `stage2_SOLUTION.md` section 8 explicitly lists "Balls resized
to fixed proper volume" among the uncomputed items. The matched Einstein
theorem requires them. Its premise is not the entropy of a fixed set of sites
under a Hamiltonian rescaling.

For the explicit, fixed-topology spatial metric

    g_ij(x)=exp(2 a |x|^2) delta_ij,

the curvature at the center is `R3=-24a`. The coordinate radius that encloses
proper volume `4 pi ell^3/3` is

    r = ell - (3/5) a ell^3 + O(a^2 ell^5).

Its proper boundary area is

    A = 4 pi ell^2 + (16 pi/5) a ell^4 + O(a^2 ell^6)
      = 4 pi ell^2 - (2 pi/15) R3 ell^4 + O(a^2 ell^6).

Holding coordinate radius instead gives `(8 pi) a ell^4`, a factor 5/2
different. `audit.py` independently integrates proper volume and solves for
the radius for both signs of curvature. At `a=.04, ell=.15`, the correct
area ratio is `.9990498`; the fixed-coordinate control is `2.5022514`.

This solves the geometric bookkeeping, not the microscopic entropy-to-area
attachment. Merely multiplying the old entropy by a chosen area factor would
not repair the source. The actual state and its ultraviolet regulator must
respond consistently.

## 3. A fixed-topology source with a responding proper regulator

The following is a declared local candidate, not derived from the existing
repair law. Let `v=exp(-sigma)` and use the same four Clifford matrices:

    H_sigma(k) = v sum_j alpha_j sin(k_j)
               + beta R^3 v^4 [2 sum_j(1-cos(k_j))]^2.

The squared lattice Laplacian has finite range two. No topology change,
empirical gravitational value, or selected Newton coefficient is used.
With physical momentum `p=v k/h`, the low-momentum symbol, after the common
unit conversion, is

    alpha.p + beta R_phys^3 |p|^4,
    R_phys=R h.

Thus the proper ultraviolet scale is unchanged under the coframe variation.
The infrared Dirac principal symbol survives. The regulator violates Lorentz
symmetry in the ultraviolet; its removal from low-energy predictions requires
an actual error bound. The candidate is not the original ordered walk.

For a constant coframe and `|p|<=P`, the elementary sine/cosine remainders
give the explicit symbol estimate

    ||H_(h,sigma)(p)-[alpha.p+beta R_phys^3 |p|^4]||
       <= h^2 [P^3+R_phys^3 P^6]/(6 v^2).

The additional difference from the massless Dirac symbol is at most
`R_phys^3 P^4`. Thus the relative low-energy correction at `P~1/ell` is
of order `(R_phys/ell)^3`, plus the displayed lattice error. This operator
bound does not by itself control a reduced state's logarithm; it cannot be
substituted for the missing three-dimensional modular estimate.

An explicit finite-range nonuniform extension, for later tests, is

    H_sigma = sum_j alpha_j {v,D_j}/2 + beta R_phys^3 (v L_h v)^2,

where `D_j` is the centered Hermitian difference, `L_h` the positive lattice
Laplacian, and v multiplication by `exp(-sigma(x))`. A lapse can be attached
as `{N,H_sigma}/2`. The operators are Hermitian for real positive N,v and
reduce to the tested symbol when uniform. Their principal metric is
`ds^2=-N^2 dt^2+v^-2 dx^2`. This specifies an executable route for a moving
proper-volume cut; its state, local entropy response and stationarity are
not asserted to have been tested by the uniform-coframe slab calculation.

In the continuum regulated family, the uniform coframe is a spatial dilation.
Unitary coordinate transport implies `S_sigma(A)=S_0(A/v)`. If its large smooth
region entropy has leading area coefficient eta, then

    S_sigma(A) = eta e^(2 sigma) Area_coordinate(A) + subleading terms.

This is a response of the actual Gaussian state, rather than a weighted
relabeling of its entropy. For a slowly varying coframe, extending this local
area law with controlled derivative errors is a remaining obligation.

`regulator_candidate.py 4` computes the actual four-component source's slab
entropy by exact Clifford reduction and transverse quadrature. The area
exponents are `1.8706, 1.9441, 2.0206` for `R/h=1,2,4`, versus geometric
exponent 2. Physical area coefficients `R^2 S/A_coordinate` are approximately
`.06780, .07204, .07327` in units `R_phys^-2`. `candidate_checks.py` refines the
last transverse quadrature and verifies the Clifford reduction against a full
four-component calculation.

The ordinary quadratic Wilson regulator was tried first and retained:
its exponents are `.9812, 1.1818, 1.3325`. It is not a successful area
dictionary on those grids. The quartic variant was introduced to suppress the
large-momentum projector tail, not selected using a measured gravity value.

### What survives of the modular premise

On the same quartic candidate, a transverse-zero sector is a one-dimensional
massless Dirac limit. For a half-circle interval of L cells, the continuum
geometric modular prediction for the lowest particle/hole energy variation is
`delta S/delta E=L/(pi v)`. At `R=4`, numerical ratios to that value are
`.98330` for `L=128` and `.99559` for `L=256`. A coframe-shifted `sigma=.1`
case and modular eigenvalue-floor controls are retained in
`CANDIDATE_CHECKS.json`.

This is an actual compatibility diagnostic: a nonzero area response and a
converging geometric bulk response occur on one declared source. It does not
prove the three-dimensional ball modular estimate, all observer directions,
or generalized-entropy stationarity. Calling premise 2 fully solved here
would be incorrect.

## 4. Remaining realization obligations

1. Derive the regulator/coframe law and its reference preparation from an
   admissible OPH source. The finite local formula above is a new declaration.
2. Construct a slowly varying physical metric and genuine fixed-proper-volume
   region instruments on that same state; prove the local curvature area
   error. A frozen-site ball or area-weighted number is insufficient.
3. Establish the bulk modular/stress coefficient on three-dimensional balls
   for all ten observers, with one clock and controlled errors.
4. Derive generalized-entropy stationarity from the source's dynamics or
   preparation. Imposing it to solve for the metric is a hypothesis, not
   its microscopic derivation.
5. Fix the positive physical eta and the physical ultraviolet scale
   independently of G; retain stress/Ward conservation and the vacuum
   reference needed for Lambda.
6. Specify the limiting regime. Refining the coordinate mesh at fixed proper
   regulator R is compatible with this candidate. Taking physical ball radii
   to zero at fixed R is not its low-energy regime. Conversely, if
   `R_n/ell_n -> 0` and `ell_n -> 0`, the unrenormalized coefficient
   `eta_n=c/R_n^2` diverges instead of tending to the matched theorem's finite
   eta. A renormalized area law or an explicitly finite-error mesoscopic
   Einstein statement is needed; silently changing physical units is not.

The existing matched theorem already gives a finite residual bound and can
be used in a declared mesoscopic window `R << ell << curvature radius` once
the three actual errors have been supplied. This packet does not supply them.
