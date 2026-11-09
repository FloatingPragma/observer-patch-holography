# Positive-energy support from the selected cooling family

This is a typed instance of
[the existing compatible-limit theorem](../particle_ground_spectral_records/PROOF.md#positive-energy-in-a-compatible-limit).
Its finite inequality is
[`GroundRecordTail.downwardWeight_le_excess_div`](../../Lean/QFT/GroundRecordTail.lean),
implemented by `particle_ground_spectral_records.records.negative_transfer_bound`.
The finite cooling proof and theorem are reused without a separate gap
premise or identification of finite and limiting spectra.

## Exact source correlations and their finite tail

For primary member n, the stationary source density rho_n lies in the
physical Gauss carrier and has mean excess above its restricted ground
at most `epsilon_n=2^(-n-4)`. Let A_n preserve this physical carrier and
obey `||A_n||²<=C_n²`. Source time is unchanged:

```
alpha_n,t(A)=exp(itH_n) A exp(-itH_n),
F_n(t)=Tr rho_n A_n† alpha_n,t(A_n)
      = integral exp(it omega) dmu_n(omega),
mu_n=sum_(i,j) p_i |(A_n)_(ji)|² delta_(E_j-E_i).        (1)
```

The measure is positive, with total mass `F_n(0)<=C_n²`. Its weights are
original Gram spectral weights; they are not conditioned on a successful
probe. An ordered correlation is not a sequential pointer probability.
An experiment reconstructing (1) needs its actual coherent apparatus and
complete-output errors.

The finite Lean theorem applies with `p'_i=C_n² p_i` and
`w_ij=p_i |(A_n)_(ji)|²`. Its row inequality follows from
`A_n† A_n<=C_n² I`, including C_n=0. It gives, for every delta>0,

```
mu_n((−infinity,−delta])
 <= min(C_n², C_n² epsilon_n/delta).                    (2)
```

`selected_tail_bound` regenerates the actual primary member and inserts
its prescribed cooling tolerance into this bound. For contractions and
delta=1/4, levels0,2,6 give upper bounds1/4,1/16,1/256.
The exact three-level matrix control rotates H, rho and A into a
non-energy coordinate basis before evaluating all nine transition traces.
Its measure has weights `1-epsilon/2, epsilon/4, epsilon/4` at frequencies
`3,-1,-2`, and mean ground excess epsilon. This is an independent finite
identity check, with no substitution for the spatial source.

The theorem uses the stationary source algebra. The canonical purification
with idle-reference Hamiltonian zero generally is not stationary under
source-only evolution; its arbitrary joint correlations cannot replace
(1). A bare charged operator need not preserve the source Gauss carrier.
An enlarged compensated charge experiment must use its actual accessible
ground floor and excess, which can differ from epsilon_n.

## A continuous correlation limit has no negative spectral weight

Assume that F_n(t) converges for every real t to F(t), including t=0,
and F is continuous at0. Assume also `C_n² epsilon_n->0` and bounded
total masses `F_n(0)`. For uniformly bounded probes the selected cooling
schedule ensures the first condition. Arbitrarily rescaled fields need
their own norm/excess estimate: an `a_n^(-3/2)` amplitude scaling does
not follow from this schedule alone.

Every F_n is positive definite. Its pointwise limit is positive definite,
and continuity at0 makes it the Fourier transform of a finite positive
measure mu. The finite-measure continuity theorem gives weak convergence
mu_n to mu. To reduce this to the usual normalized statement, if
`F(0)>0`, divide by the eventually positive masses F_n(0), apply the
characteristic-function continuity theorem, then restore those masses.
If F(0)=0, all masses tend to zero and weak convergence to zero is immediate.

For a nonnegative continuous compactly supported test h contained in
`(−infinity,−delta]`, (2) gives

```
integral h dmu_n <= ||h||_infinity C_n² epsilon_n/delta ->0.
```

Weak convergence therefore gives `integral h dmu=0`. Such tests exhaust
the negative half-line, so `supp mu subset [0,infinity)`. No convergence
of commutators or finite spectral eigenvalues is required for this step.

For a nonnegative implementing Hamiltonian on an entire representation,
more than one chosen correlation is needed. A common stationary limiting
algebra and its canonical strongly continuous GNS implementing group
`U(t)=exp(itL)` must exist, with `U(t)Omega=Omega` and
`F_A(t)=<pi(A)Omega,U(t)pi(A)Omega>`. An arbitrary implementing group can
carry a scalar phase and does not fix the same vacuum energy origin.
The covered cyclic vectors must have dense linear span there.
Their negative spectral projections vanish by the preceding argument;
density then gives `1_(-infinity,0)(L)=0`. The finite tail theorem supplies
this spectral-sign implication. It supplies no common algebra, state
limit or strongly continuous dynamics by itself.

Continuity at0 is substantive. A separate finite ground-probe control
has equal original weight `1/n²` at energies `k/n`, k=1,...,n². Its
correlation equals1 at0 and tends to0 at every fixed nonzero time:
the geometric sum is bounded by `2/(n²|1-exp(it/n)|)`. The discontinuous
limit loses all finite-energy mass. Zero limiting correlations also
satisfy positive spectral support and contain no observed excitation.

## A fixed temporal filter supplies a continuous response subsequence

A bounded-response corollary avoids assuming continuity for a single
filtered family. Fix T>0 and a uniform bound `||A_n||<=C`. On each finite
source form the norm-convergent operator integral

```
f_T(s)=(1-|s|/T)_+/T,
A_n,f=integral f_T(s) alpha_n,s(A_n) ds.                (3)
```

This preserves the physical carrier, and Gauss invariance when the original
probe is gauge invariant; it also gives `||A_n,f||<=C`. Its measure is
`dnu_n=|fhat_T(omega)|² dmu_n`, with multiplier
`sinc(omega T/2)^4`, where `sinc(x)=sin(x)/x` and `sinc(0)=1`.
Since `||f_T||_1=1` and
`||f'_T||_1=2/T`, integration by parts bounds
`|fhat_T(omega)|<=min(1,2/(T|omega|))`. Consequently

```
nu_n(R)<=C²,
integral |omega| dnu_n <= 2 C²/T,
|F_n,f(t)-F_n,f(u)| <= (2 C²/T)|t-u|.                  (4)
```

For the middle bound, use
`x min(1,(2/(Tx))²)<=2/T` for x>=0. Thus the filtered correlations are
uniformly bounded and equicontinuous on every compact time interval.
Arzela-Ascoli and a diagonal subsequence give compact-time convergence
to a continuous positive-definite function. Their negative tail is
bounded by (2), since the multiplier is at most one. The limiting filtered
measure therefore has support in `[0,infinity)`. The same diagonal choice
works for a countable list of specified probes and fixed filters, provided
each list member has its stated uniform norm bound.

This is a response-subsequence statement without an effective choice or
convergence modulus. Filtering changes the temporal observation domain;
the operator may be spatially spread by source dynamics. No local detector
implementation of the exact integral is asserted. The response can vanish
or concentrate at zero energy, and the collection is not thereby a common
field algebra or a faithful physical vacuum. This corollary does not
establish joint refinement, nontrivial interactions, Lorentz covariance,
a positive gap, species identification or particle masses.
