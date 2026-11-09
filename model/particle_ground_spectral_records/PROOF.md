# Complete charged records and a positive-energy bound on one source

## The source, preparation and clock

The source is the native quantum-flux control in
`particle_observer_outputs.field_route.exact_rotor(6)`. Its Hilbert space is
the complete four-mode CAR space tensored with flux labels `-6,...,6`, of
dimension208. With the native singlet charge `q6Y=6`,

```
H = I tensor E²/2 + dΓ(beta) tensor I
    + dΓ(hop) tensor S + dΓ(hop)† tensor S†,
hop = (-beta+beta gamma)/2,
G = 6(N-2),     U(t)=exp(-itH).
```

The source matrices satisfy `beta²=gamma²=I` and
`beta gamma+gamma beta=0`. The shift raises flux by6 and is compressed
at the declared cutoff. No occupation, flux, reference, or recorded
outcome is removed. The82 reducing blocks, of dimension at most12,
partition all208 basis states. They are a way to compute the unchanged
operator, rather than a replacement source.

`particle_spatial_ground_state.ground.prepared_control` supplies

```
Q = [sum_(j=0)^m (beta_T (22I-H)/(2s))^j/j!]^s,
rho = Q²/Tr(Q²),     |Psi_Q> = vec(Q)/sqrt(Tr(Q²)).
```

Here `beta_T` is inverse temperature, distinct from the Dirac matrix.
The source density commutes with H. The retained purification need not
be stationary under a source-only clock. Its canonical half-trace error
relative to the Gibbs purification is bounded by that package's exact
Taylor estimate. The spectral upper bound22 is a preparation parameter;
changing it changes the finite Taylor polynomial. Its exponential target
differs from `exp(-beta_T H/2)` only by a scalar that cancels at
normalization. The preparation changes neither H nor its clock.

## Exact ground and charge thresholds

Every reducing matrix is real rational in the inherited diagonal phase
chart. For rational x, the LDL test checks `H_block-xI>0` exactly: each
successive pivot must be positive. A nonpositive pivot proves failure of
positive definiteness because all preceding pivots were positive.
Therefore every block is positive definite exactly when x lies below
the minimum of their union. Rational bisection starting at `[-4,22]`
encloses each occupation-sector minimum. Its returned width is at most
`2^(-32)` in the default calculation.

The N=2 minimum lies strictly below every other occupation sector.
This identifies the charge sector of the global minimum without preparing
a state conditioned on that sector. The exact interval is

```
-74228867013/34359738368 <= E0 <= -148457734013/68719476736.
```

The minima in N=1 and N=3 have equal enclosures. Subtracting the full
ground enclosure gives, for both adjacent charges,

```
19030367577/17179869184 <= E_(N=1,3)-E0
                      <= 38060735167/34359738368,
```

or approximately `1.107713183` source energy units. This is a sector
threshold. An operator's spectral residue at that threshold is a separate
quantity; the displayed residue calculations are floating diagnostics.
No unit conversion, dispersion relation, Lorentz representation, or
observed particle mass follows from the threshold alone.

## Exact charge conjugation of the finite source

Let `B=product_(j=0)^3(c_j+c_j†)`, let `Gamma(gamma)` be the exterior
lift of the native gamma matrix, let F reverse the flux basis, and let
K mean complex conjugation in the native occupation/flux basis. Define

```
C = [Gamma(gamma) B tensor F] K.
```

The executable checks the full208-dimensional identities

```
C²=I,     CHC^(-1)=H,     CGC^(-1)=-G,     CEC^(-1)=-E.
```

Consequently C maps N to `4-N` without changing source energy. Their
complete spectra coincide. Every real polynomial Q(H), and therefore
its source density, is C-invariant. Addition and removal responses,
summed over the complete four-mode frame, have identical energy-transfer
measures. These statements describe a symmetry of this finite regulator.
Identification with physical particle/antiparticle multiplets requires
their spatial fields and transformation laws on the same construction.

## A charge- and parity-compensated complete probe

The eight Kraus operators are

```
K_(+,i)=c_i†/2,     K_(-,i)=c_i/2,     i=0,1,2,3.
```

Each acts trivially on the flux register. Their effects obey

```
sum_i K_(+,i)†K_(+,i)=I-N/4,
sum_i K_(-,i)†K_(-,i)=N/4,
sum_(sign,i) K_(sign,i)†K_(sign,i)=I.
```

Attach an orthogonal reference state to each outcome. The addition
reference carries charge-6 and odd fermion parity; the removal reference
carries charge+6 and odd parity. The total Stinespring map is therefore
charge conserving and even. A pointer can measure the sign and declared
Dirac-component label while the compensating quantum reference is
retained. Coherent component readout is also possible by retaining that
label quantum mechanically. No physical spin measurement is inferred
from the component label.

The joint preparation's dual reference carries the opposite source
charge and matching parity. Since Q commutes with G and parity,
`vec(Q)` has total charge zero and even total parity. All probability
and distance statements apply with any additional untouched reference.

There is no rejected branch in the complete eight-outcome probe. If an
experiment calls addition its success event, the removal outcomes are
retained failures with their original probabilities. Pauli-blocked
outcomes remain as zero-weight records. On the C-invariant prepared
source, the total addition and removal weights are each1/2. Neither is
renormalized to one.

## Lehmann weights and the finite theorem premises

Diagonalize commuting H and rho. Write their eigenvalues as E_n and p_n,
and define original weights

```
w_(n,m,a)=p_n |<m|K_a|n>|²,
mu(B)=sum_(n,m,a) w_(n,m,a) 1_(E_m-E_n in B).
```

Every weight is nonnegative, and instrument completeness gives
`sum_(m,a) w_(n,m,a)=p_n`. Thus `mu(R)=1`, including all initial
occupations, all final states and every probe outcome. For any eta>0,
`E_m-E_n<=-eta` implies `E_n-E0>=eta`. It follows that

```
eta mu((-infinity,-eta]) <= sum_n p_n(E_n-E0) = epsilon.
```

This is exactly the input domain of
`OPH.QFT.GroundRecordTail.downwardWeight_le_excess_div`, with
`w(n,m)=sum_a w_(n,m,a)`. The code supplies the same source energy floor,
the stationary normalized density and the complete contraction
instrument. The formal theorem supplies the finite weighted inequality.
Its zero-excess corollary excludes every strictly negative transfer.
For a single bounded A, the row sum is at most `||A||² p_n`, giving
`mu_A((-infinity,-eta])<=||A||² epsilon/eta` without normalizing mu_A.

The proof needs no spectral gap. In particular, a full-rank cooled state
has negative-frequency weight in general. Deleting that weight would
change its observable law. The preparation's exact mean energy and the
certified lower bound on E0 give a rational upper bound on epsilon.

## An implemented finite clock apparatus

The apparatus reads energy transfer using the same source clock on the
source and its conjugate replay on the retained preparation reference.
After probe outcome a, its unnormalized vector is `vec(K_a Q)/sqrt(Z)`.
For `W(t)=U(t) tensor conjugate(U(t))`, its energy-difference eigenphase
is `exp(-it(E_m-E_n))`. The reference replay is a declared controlled
apparatus operation on the prepared reference. It is not an inference
that the reference's idle Hamiltonian is negative, and it does not
change the source Hamiltonian or rescale time.

Prepare a uniform M-state clock pointer, apply controlled powers
`W(t)^j`, then a finite Fourier transform and record all M pointer
labels. The effects on energy transfer have kernels

```
k_b(omega)=|M^(-1) sum_(j=0)^(M-1)
                  exp(-ij(t omega-2pi b/M))|²,
sum_b k_b(omega)=1.
```

For the executed M=2 apparatus, replace U by the inherited exact Cayley
unitary with time1/64 and four substeps. In the real computational chart,
for every probe block from n to m, put

```
X=K_(m,n) Q_n,
Y=U_m X U_n†,
p_(a,0)=sum_blocks ||X+Y||_F²/(4Z),
p_(a,1)=sum_blocks ||X-Y||_F²/(4Z).
```

All16 probabilities are calculated with exact rational arithmetic.
The raw sum is exactly one. Output intervals are outward dyadic
enclosures of these exact values. Commuting Q with the source unitary
allows evaluation from `Q_n²` and small clock matrices, avoiding repeated
multiplication of the large preparation denominators. This is an algebraic
trace identity, with commutation checked before use.

If the one-copy Cayley operator error is delta, the replay operator
error is at most2delta. Contractivity bounds the complete output
half-trace error by `preparation_error+2delta`. All mode pointers,
compensation registers and untouched references are included. No
successful output normalization enters this error estimate.

Two coarse phase bins do not resolve the source lines or classify a
negative-frequency event. The separate finite-resolution construction
in `model/spectral_readout/kernels.py` supplies positive kernels,
alias bounds and interval-readout bounds. Here the entire source energy
transfer has `|omega|<=26`. If `0<Delta<eta`, its reused estimates give

```
P(reported energy<=-eta)
 <= min(1, epsilon/(eta-Delta)+filter_leakage+clock_error).
```

The code supplies an exact resolution plan for a requested mean excess
of1/100, with the number of bins and maximum observation time explicit.
The gap-free preparation family can meet that tolerance; the executed
beta_T=4 control has the larger excess about0.1564. The finer plan's full
probability vector is not claimed to have been executed. The two-bin
probability vector is executed. An additional plan applies to that same
executed preparation: with resolution1/16 and threshold1/4, all filter
and clock terms give an original negative-record probability below0.89.

## Positive energy in a compatible limit

For a bounded observable A the same-clock characteristic function is

```
F_A(t)=Tr rho A† exp(itH) A exp(-itH)
      = integral exp(it omega) dmu_A(omega).
```

Its positive measure is the Gram spectral measure of
`vec(A sqrt(rho))` for the energy-difference generator. At finite
temperature that generator can have negative spectrum. If a compatible
family has `epsilon_L->0`, convergent local states and strongly
continuous convergent local dynamics, every compactly supported
negative-energy test vanishes in the limit by the bound above. Density
of local GNS vectors then gives a nonnegative implementing generator.

This reuses the Fourier-test argument in
`model/thermodynamic_ground/README.md` and the complete finite-resolution
readout argument in `model/spectral_readout/README.md`. Their separately
constructed contact-interaction ground is not substituted for this
gauge-coordinate source. The algebraic GNS modules
`Lean/QFT/ColimitSelectedStateGNS.lean` and
`Lean/QFT/ColimitGNSRepresentation.lean` consume compatible positive
states; they do not supply local dynamics or the energy inequality.
The fixed-spacing, fixed-link-cutoff spatial ground construction is
specified in `particle_spatial_ground_state/PROOF.md`.

The native208 calculation is a finite quantum-flux control. It supplies
neither a three-dimensional momentum variable nor a relativistic mass
shell. A cutoff-refined interacting chiral field limit, observed species,
GeV masses and Einstein dynamics require independent source attachments.
