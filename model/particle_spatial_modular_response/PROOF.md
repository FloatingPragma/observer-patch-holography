# Spatial restriction and same-state modular response

## 1. Declared source and regional injection

The source is the finite Wilson coordinate field of
[the spatial preparation construction](../particle_spatial_ground_state/PROOF.md):
all sixty CAR modes per cell, three outgoing native link-word factors
per cell, product per-link cutoff, one Hamiltonian `H` and clock
`tau_t(A)=exp(itH)Aexp(-itH)`. The parameters `a,e0,kappa0`, the state
preparation and any source-Gauss projection are part of its specification.
This consumer does not change them. A scalar shift of H leaves this
clock unchanged.

`common_contract_bound` authenticates the selected member of
[the common construction](../particle_common_contract/PROOF.md), consuming
its exact source identifiers, physical projection, beta, Taylor order,
dyadic scale and spectral interval. Its certificate binds the resulting
support floor to that member. It does not insert the separate warm
208-state control or replace the physical support with the ambient space.

For a set R of cells, put its sixty-mode groups before the exterior
groups in the fermionic ordering. A basis vector acquires sign
`(-1)^I`, where I counts pairs of occupied modes whose order changes.
Outgoing link factors belong to their tail cell. This gives an explicit
unitary permutation J and the factorization

```
H_ambient = H_R tensor H_out,
R(rho)=Tr_out(J rho J*),
j(A)=J*(A tensor I)J.                               (1)
```

The code computes the exact labels and CAR sign on the full sixty-mode
blocks without allocating their exponential matrix. The injection is
a unital star homomorphism; restriction is its trace-dual completely
positive map. Tensoring with an untouched reference preserves these
identities. No preparation outcome is conditioned away. The tracial
expectation `j(Tr_out(X)/d_out)` realizes the trace-based restriction
and retraction fields of `Lean/QFT/CPRestrictionNet.lean`.
It need not preserve the interacting prepared state.

This is a spatial regional algebra for the declared lattice ownership.
It is not an operational proper ball, a Lorentz diamond, or a statement
that the global evolution preserves the regional algebra.

## 2. Faithfulness from the actual source bounds

Let `[L,U]` contain the spectrum of H and `W=U-L`. The source preparation
gives these bounds from its electric, magnetic and kinetic terms. For
the exact finite Gibbs state, with `D=d_R d_out`,

```
rho_beta=exp(-beta H)/Tr exp(-beta H)
       >= exp(-beta W) I/D,
sigma_R=R(rho_beta) >= exp(-beta W) I_R/d_R.           (2)
```

The same inequality holds directly for the positive polynomial
preparation `rho_Q=Q(H)^2/Tr Q(H)^2` used in the source construction.
Its positive Taylor amplitude satisfies
`1<=Q(H)<=exp(beta W/2)I`. Thus its full and regional lower bounds
follow without subtracting an approximation error from a tiny weight.
In particular the finite prepared density itself is faithful.

For `r` regional cells and per-link word degree ell,
`d_R<=2^(60r)451^(3r ell)<=2^((60+27ell)r)`. Since
`exp(x)<=2^(2x)` for nonnegative x, a fully rational factored floor is

```
sigma_R >= 2^[-ceil(2 beta W)-(60+27ell)r] I_R.       (3)
```

The certificate keeps the exponent factored. It does not allocate a
huge denominator. This finite-volume floor generally deteriorates with
volume or cooling; no faithful limiting vacuum follows from (3).

## 3. Physical Gauss support and the boundary algebra

If the contract requires a physical source projection P, use
`[P,H]=0` and the directly prepared state

```
rho_phys=P Q(H)^2 P/Tr(P Q(H)^2),
rho_phys >= exp(-beta W) P/rank(P).                  (4)
```

Equation(4) is faithful on the physical carrier. It is not faithful
on the full ambient space. Partial trace gives a lower bound
proportional to `Tr_out P`; replacing that operator by an identity
without proving its support would be an error.
Using the source dimension bound `rank(P)<=2^b` gives the bound
`rho_phys>=2^[-ceil(2 beta W)-b] P` and its actual regional image
with P replaced by `Tr_out P`. This is the selected common-contract
consumer. Its support is source determined. Full spatial boundary
multiplicities and the large reduced matrix are not allocated by
this certificate.

After imposing the interior Gauss constraints on each side of a cut,
decompose the actual boundary gauge action into irreducible sectors:

```
H_R = direct_sum_lambda (M_R,lambda tensor V_lambda),
H_out = direct_sum_lambda (M_out,lambda tensor V_lambda*),
H_phys = direct_sum_lambda (M_R,lambda tensor M_out,lambda
                            tensor span(Omega_lambda)).           (5)
```

Only matching conjugate representations contribute to H_phys.
Here `Omega_lambda` is the normalized invariant pair, `d_lambda=dim V_lambda`,
and the multiplicities are `m_lambda,n_lambda`. This decomposition must
come from the specified finite boundary representation; dimension alone
does not determine it. General boundary gauge groups include products of
the cut-vertex groups. Nonmatching sectors with `n_lambda=0` remain
zero-probability ambient outcomes.
Exterior irreducible sectors absent from H_R contribute no physical
support and can be omitted from this regional-sector indexing.

The physical projector has rank `rP=sum_lambda m_lambda n_lambda`, and

```
Tr_out P = direct_sum_lambda (n_lambda/d_lambda) I_(m_lambda d_lambda).
```

On every admitted sector the extended regional density has floor
`exp(-beta W)n_lambda/(rP d_lambda)`. Gauge-invariant regional
observables form `A_R=direct_sum_lambda B(M_R,lambda) tensor I_Vlambda`.
Writing their state as the unnormalized blocks sigma_lambda gives

```
sigma_extended = direct_sum_lambda sigma_lambda tensor I_dlambda/d_lambda,
sigma_lambda >= exp(-beta W)n_lambda I_mlambda/rP.   (6)
```

All center weights `p_lambda=Tr sigma_lambda` are retained. With the
unweighted block trace on the observable algebra,

```
S_algebra = -sum_lambda Tr sigma_lambda log sigma_lambda,
S_extended = S_algebra + sum_lambda p_lambda log d_lambda.          (7)
```

The modular Hamiltonians differ by the central term `log d_lambda`.
It affects entropy variations of the center weights, while commuting
with every observable in that algebra. The edge multiplicity is a
specified algebraic entropy convention; (7) supplies no geometric area
or inverse Newton constant.

The code tests the multiplicity factors using the actual native
M15/conjugate-M15 representation. Its five imported irreducible charge
blocks have dimensions3,2,6,3,1. The invariant-pair projector has rank5
on225 dimensions and regional eigenvalues `1/d_lambda` before state
normalization. This is a boundary-representation control, not a
substitution of its tracial state for the spatial Wilson preparation.

## 4. Modular operator and finite entropy remainder

Apply the existing matrix-entropy proof in
[prepared_regional_entropy](../prepared_regional_entropy/PROOF.md)
and the source adapter
[`first_law_remainder`](../native_gravity_identification/entropy.py)
to the actual faithful sigma_R, or separately to the supported blocks
in (6). For any admitted trace-zero Hermitian X,

```
K_R=-log sigma_R,
S(sigma_R+tX)-S(sigma_R)
 = t Tr(X K_R)-D(sigma_R+tX || sigma_R),
0<=D<=t² ||X||_HS²/[2(m-|t| ||X||)],                 (8)
```

provided the stated path floor is positive. These are matrix logarithms;
no commutation of X and sigma_R is required. The bound follows from the
resolvent derivative of log and is used here with the source-derived
floor. The center contributes its ordinary probability entropy, so
sectorwise application retains both classical and internal variations.

This instantiates the finite entropy/modular pairing premise. The
Shannon restriction agrees with `CapFirstLaw.lean` and
`EinsteinPremiseLink.lean`. The latter's additional split
`K=2pi B+Z` is not derived merely by defining `B=K/(2pi)`.
Identifying B with a physical diamond-energy kernel and Z with a
proper-area response requires separate source data.

For the executed rational matrices the code certifies
`q=||I-d sigma_R||_infinity<1`. Set Y=`I-d sigma_R` and

```
P_n=sum_(j=1)^n Y^j/j,
||K_R-log(d)I-P_n|| <= q^(n+1)/[(n+1)(1-q)].         (9)
```

All coefficients, matrix products and bounds in (9) are rational.
The scalar `log d` is retained in K and cancels in its pairing with
trace-zero tangents. This comparison is made to the logarithm of the
actual prepared marginal, not directly to the ideal Gibbs state.

## 5. A complete source-admitted intervention

Choose a regional unitary U commuting with the admitted gauge action
and parity. The two-outcome instrument

```
M0=sqrt(1-t) I,   M1=sqrt(t) U,   0<=t<=1           (10)
```

has complete effect sum I. It acts on the original state and leaves
the canonical reference untouched. Its original pointer branches
are retained. On the regional marginal, the aggregate response is
`sigma_t=(1-t)sigma+t U sigma U*`, so `X=U sigma U*-sigma`.
Both endpoints have floor m; the whole convex path has the same
floor. The remainder in (8) therefore improves to
`t²||X||_HS²/(2m)`. The general bound in (8) is also evaluated and
checked for the actual control.

To read the modular pairing with a finite apparatus, take
`B_n=sum_(j=1)^n q^j/j`, which bounds `||P_n||`, and use effects
`E_±=(I±P_n/B_n)/2`. They are positive and sum to I. The four
original preparation-coin/modular-pointer probabilities are

```
Pr(b,±) = w_b [1±Tr(U_b sigma U_b* P_n)/B_n]/2,
(w0,w1)=(1-t,t), (U0,U1)=(I,U).                     (11)
```

They are rational in the control and sum exactly to one. Their
positive square-root Kraus matrices define the complete finite
apparatus; no preferred outcome is selected. The effects are in the
same gauge-invariant algebra. Unsupported ambient outcomes may be
completed by `I/2` on the orthogonal support and have zero probability
for these preparations. Equation(11) preserves all source/reference
outputs under the declared finite instrument service.

The entropy in (8) is the regional marginal entropy after ignoring
the pointer for that observable. The complete program retains the
pointer and reference. Their joint entropy is not replaced by (8),
and the two conditional unitary branches individually preserve the
regional entropy. A deterministic source clock is not inferred from
an irreversible coarse observation.

## 6. Executed spinor and flux controls

The source is exactly the inherited208-state quantum-flux Hamiltonian:
all sixteen four-mode CAR occupations times thirteen fluxes. The
control state is the inherited positive polynomial with beta=`1/64`
and Taylor order12. It is a separate declared warm preparation,
not the canonical spatial low-energy family's state.
The source-phase chart and full rational density are imported intact.

The two factors are the spinor register and flux register, with
dimensions16 and13. They are internal factors of the directional
one-cell control. They are not two spatial regions. The full source
has `[-4,22]` spectral enclosure, and its positive amplitude is at
least I. If its raw squared-amplitude normalizer is Z_Q, exact
partial trace gives floors `13/Z_Q` and `16/Z_Q`.

A separate directly prepared source-Gauss state uses N=2. It has
dimension78; its spinor marginal has faithful6-dimensional support
and ten zero ambient spinor outcomes. Its flux marginal remains
13-dimensional. Its normalizer supplies floors `13/Z_phys` and
`6/Z_phys`. This state is constructed canonically on the physical
carrier, not obtained by keeping the successful outcome of a
measurement on the warm208-state preparation.

For the spinor intervention U is the second-quantized swap of modes0
and2, including every CAR sign. For the flux intervention U swaps
flux0 and flux6. Both preserve the actual source Gauss generator
`6(N-2)` and parity; they are admitted apparatus operations, not
asserted to be free source evolution. Take `t=1/8` in (10).

All four control cases have certified nonzero modular slopes and
complete probabilities (11). A trace-zero matrix-element witness,
including the error in (9), excludes
`K_R=beta Tr_out(H)/d_out+constant` in each case. Interacting regional
reduction must therefore precede any proposed geometric local-energy
identification. Exact tangents of the source energy and its declared
isotropic strain are also reported from this same global varied state.
The resulting mechanical pressure derivative uses the declared
volume convention; it does not identify the modular energy with
proper stress.

## 7. Clock and geometric boundaries

For a Gibbs state, `-log rho=beta H+log Z` and the convention
`sigma_s(A)=rho^(is) A rho^(-is)` equals `tau_(-beta s)(A)` globally.
The polynomial source state instead has
`-log rho_Q=log Z_Q-2log Q(H)`. It is stationary, but its modular
generator is generally a nonlinear function of H. Neither identity
makes a regional modular flow the restriction of global H evolution.

The reference has idle Hamiltonian zero. The source marginal is
stationary; its canonical joint purification need not be stationary.
An apparatus counter-evolution is an explicit intervention, not a
negative Hamiltonian assigned to that idle reference. Charged probes
outside the source-only physical carrier require their full dressing
or compensating-reference representation and total Gauss generator.

The inherited strain derivatives vary the declared coframe at fixed
bare regulator spacing and fixed prepared density. Re-preparing a
state for every strained action defines a different experiment.
This construction provides an actual spatial restriction, supported
modular operator and source-admitted entropy tangent. It supplies no
proper-ball localization, geometric diamond kernel, fixed-proper-volume
area variation, generalized-entropy stationarity or matched shrinking-
ball estimates. Those are the additional premises required by the
same-state Einstein consumer. The dimension decoder and tensor
reconstruction theorem are not repeated here.
