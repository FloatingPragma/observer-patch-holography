# A finite native charge-energy centroid from the actual preparation

## Source, state and retained charge reference

The source consists of the two complete native eight-cell copies used by
the physical preparation. Each has the unchanged Hamiltonian
h=H0+2^-34 T. The active copy is read; the partner, all spin references,
filter counters and failures remain. Add the complete thirty-mode static
CAR reference with Hamiltonian4N_ref, initially in its vacuum. The total
autonomous Hamiltonian is

```
Htotal=h_active+h_partner+4N_ref,
```

with the preparation's zero-energy references and pointers included.
The outer affine field and the added matrix-unit hopping of the replicated
native join are not part of this readout source.

Run the original physical two-copy preparation once. It succeeds with
probability q>=q_min. Its actual conditional active-source density rho
satisfies D(rho,|g><g|)<=delta, where D is half the trace norm and the
frozen receipt gives delta<2.355923e-7. The partner and reference state
may remain correlated with the active source. Every readout operation
acts only on the active source and the added charge reference, so the
reduced-state bound controls the relevant statistics without discarding
those other outputs. Subsequent common free evolution preserves delta
because g is an actual eigenstate of the same h.

There is **no additional P0 ground projection** in this experiment. The
actual preparation state is not asserted to be exact or stationary.
Its error enters the final centroid bound below. The compensated native
e^c CAR contraction B is the original charged-response pump, including
its full reference-vacuum guard and Halmos completion. Source and
reference change charge and parity together. All Gauss constraints and
every complementary kinematic state remain present. Reference energy4
on successful pump branches is retained explicitly.

## A narrow actual-H admission and an exact-ground comparison

Use the radius r=128v window I from `model/native_charge_resolution`,
where v=360*2^-34. Write A=h_active+4N_ref and c for the rational midpoint
of the exact-root isolator for a0=7e+epsilon0+4. Then
I=[c-r,c+r], and P_I is the complete actual A spectral projector on the
one-reference physical sector. This bin test remains a declared exact
spectral service; the final finite clock does not replace it.

For comparison only, let chi=P_I B g, w=||chi||^2. The native resolution
certificate proves w>=w_lo>0.4197064484, and gives all nine atomic support
roots and weights. The comparison atomic pump splits into orthogonal
vectors psi_j, with norms squared w_j and energies a_j. Its selected
component is psi0. All these vectors have every PW119 link blank.
The code independently replays all59 actual link coordinates and their
zero blank-to-blank entries. Consequently

```
P_allblank T P_allblank=0,
<psi0,(A-a0)psi0>=0,
||(A-a0)psi0||<=v sqrt(w0).
```

This first-moment identity holds for every matter/scalar vector on the
blank-link carrier. It does not require psi0 to be an actual interacting
eigenvector.

Let d0 be the certified distance from a0 to the complement of I, and
R=r+(a0_hi-a0_lo)/2. Put a=P_I psi0. Spectral orthogonality and the zero
full first moment give

```
<a,(A-a0)a>=-<(I-P_I)psi0,(A-a0)(I-P_I)psi0>.
```

Outside I, |E-a0|>=d0. Hence the absolute selected first moment inside I
is at most v^2 w0/d0, by its complete second spectral moment. This bound
includes every high-energy source state without multiplying by an
unrestricted Hamiltonian norm.

The full actual-ground vector and all other atomic components give
chi=a+e_err with

```
||e_err||<=rho_err=eta+ell_other.
```

Here eta is the frozen actual-ground/atomic-trial vector error, and
ell_other is the complete weighted leakage of the other eight visible
atoms into I. Both are reused from the refined-window proof. Since both
a and e_err lie in P_I and ||A-a0|| there is at most R,

```
|<chi,(A-a0)chi>|
 <= v^2 w0_upper/d0
    +R[2 sqrt(w0_upper)rho_err+rho_err^2].
```

Division by w>=w_lo bounds the conditional centroid's deviation from a0.
The actual neutral energy satisfies
8e-delta0<=Eg<=8e with delta0=v^2/(g_at-v). Thus the exact-ground
comparison source-excitation centroid, with reference cost4 explicitly
removed, lies in

```
[a0_lo-4-8e_hi-e_cent,
 a0_hi-4-8e_lo+delta0+e_cent].
```

Its width is below3.2e-9. This is an expectation of the actual full A on
the complete accepted bin, not a single energy eigenvalue or a claim
that the original maximum-residue E_star is in a smaller interval.

## Actual prepared-state error through the complete flagged instrument

Let S and T be the unnormalized admitted states obtained from the actual
rho and exact-ground comparison respectively. Their traces are s and w.
The pump followed by the bin flag is a complete instrument with success
K=P_I B and all pump/bin failures retained. Trace-norm contraction on its
full flagged output gives

```
||S-T||_1 + |s-w| <=2delta.
```

Indeed the full failure-block norm is at least the absolute difference
of its traces, which equals |s-w|. The normalized-state identity

```
S/s-T/w=(S-T)/w + S(1/s-1/w)
```

and ||S||_1=s therefore imply

```
D(S/s,T/w)<=delta/w<=delta/w_lo.
```

The effect of K has norm at most one, so |s-w|<=delta. In particular
s>=w_lo-delta>0 and the original-attempt admitted count has lower
probability q_min(w_lo-delta). The actual source/reference energy A has
spectral range **2r** on both normalized admitted states. Subtract its
lower spectral endpoint before applying the trace-distance inequality:

```
|tr(A S/s)-tr(A T/w)| <=2r delta/w_lo.
```

The absolute norm of the raw source Hamiltonian is not used here. The
preparation contribution is about3.012e-12. It accounts for the actual
conditional prepared state and is added to the exact-ground comparison
interval. The static reference cost and the same Eg subtraction are
constants and cancel in this comparison.

## Complete twenty-bit source-energy estimator

Let M=2^20 and b=c-2r. Use the explicit source-derived pointer phase b,
base source time t0=pi/(2r), and controlled active generator A. For an
actual energy E in I,

```
theta=(E-b)/(4r) in [1/4,3/4].
```

The ideal controlled unitary and finite pointer Kraus operators are

```
C=sum_j |j><j| tensor exp[-ij t0(A-b)],
K_k=(1/M)sum_j exp(2pi i jk/M) exp[-ij t0(A-b)],
Ehat_k=b+4r k/M,       0<=k<M.
```

Positive Fourier decoding is fixed by the negative source-evolution
sign. Every label is retained, including labels whose decoded energy
falls outside I because of finite-pointer tails. The b phase contains
the absolute atomic source energy and the reference4; it is not a free
phase equal to unknown Eg. The partner h and every reference continue
autonomously through the common exp(-iT Htotal). Only the extra controlled
interaction uses A-b. It is identity outside the admitted branch.

The general circle-grid tail bound is

```
Pr(d>=x)<=min(1,1/(2M^2 x^2)+1/(2M x)),
```

where d is the circular distance between theta and k/M. Integrating the
tail from1/M to1/2, and bounding the remaining interval by1/M, gives

```
E[d] <= [3/2+ln(M/2)/2]/M-1/M^2.
```

If the linear distance |k/M-theta| wraps around the circle, theta in
[1/4,3/4] implies d>=1/4. The probability of wrapping is therefore at
most2/M+8/M^2, while linear distance is at most d plus its wrap indicator.
For M=2^m and ln2<7/10,

```
E[|k/M-theta|]
 <= [7/2+7(m-1)/20]/M+7/M^2 = kappa_M.
```

The logarithm bound follows from the positive exponential series through
degree four at7/10, whose exact rational sum exceeds2. This bound is
uniform over every off-grid admitted source energy and extends to every
admitted density matrix and retained reference. The ideal estimator's
mean bias is at most4r kappa_M.

For an explicit finite rational realization take n=512 and

```
Uhat=R33(t0(A-b)/n)^n,
R33(X)=(a+i d)^-1(a-i d),
a=I-X^2/10, d=X/2-X^3/120.
```

The exact denominator is at least I. Since the admitted norm satisfies
||t0(A-b)||<=3pi/2<33/7, the full controlled-pointer operator error is at
most

```
gamma=(M-1)(33/7)^7/(100800 n^6)<3e-11.
```

All factors remain functions of the actual A and preserve its energy
eigenspaces. The finite record laws have total variation error at most
gamma. Because Ehat_k has range less than4r, its additional mean bias is
at most4r gamma. Define epsilon_clock=4r(kappa_M+gamma). This is below
1.04e-10. The actual prepared-state finite-pointer mean therefore lies
in the comparison centroid interval enlarged by

```
2r delta/w_lo + epsilon_clock.
```

The resulting source-excitation mean interval is approximately
[0.6916566792649156,0.6916566826118845], of width below3.4e-9.
Its exact rational endpoints are authoritative. The straightforward
binary-power realization uses n(M-1)=536870400 rational factors, with
largest coherent extra source time (M-1)pi/(2r). Finite coefficients,
functional inverses and native controlled actions are declared apparatus
services. No efficient gate synthesis or numerical execution of the
full source-plus-pointer matrix is claimed.

For a concrete source-energy display the decoder can subtract the
reference4 and the explicit rational baseline
Eg_report=4(e_lo+e_hi)-delta0/2 from each Ehat_k. This is a known source
bound midpoint, not the unknown exact Eg. The receipt separately gives
the absolute source-plus-reference mean and the mean of these decoded
source values. Interpreting the latter relative to the actual neutral
ground adds the explicit uncertainty
4(e_hi-e_lo)+delta0/2. That addition recovers the displayed excitation
interval exactly. The apparatus requires no hidden exact-ground phase.

## Complete records, work and clock

Let rho include the actual prepared active-source density and the fresh
static reference vacuum. The exact instrument law is

```
preparation failure:       1-q,
CAR pump failure:          q tr(K_fail rho K_fail^dagger),
outside energy bin:        q tr((I-P_I)B rho B^dagger(I-P_I)),
energy label k in bin:     q tr(K_k P_I B rho B^dagger P_I K_k^dagger),
                           every0<=k<M.
```

The actual rational-clock Kraus operators are used in this law. They
obey sum_k K_k^dagger K_k=I on the admitted sector by Fourier
orthogonality and exact unitarity, so the law sums to one. The post-bin
density in this expression is unnormalized. All failure flags, original
filter counters, source partners and unsupported zero-weight reference
occupations remain. Pointer mean bounds are conditioned on the explicit
bin record; the receipt separately reports its original-attempt count.

Every energy readout commutes with A. The native pump may not commute
with Htotal. Its complete Halmos unitary U_B nevertheless has a finite
Hermitian logarithm L_B. A one-unit interaction-picture implementation
uses V_B(t)=exp(-it Htotal)L_B exp(+it Htotal), producing
exp(-i Htotal)U_B. L_B is controlled by the preparation-success flag and
commutes with all Gauss constraints. All bin and energy operators commute
with Htotal, so this common free factor commutes through the later
instrument and leaves the stated probabilities and centroids unchanged.
No-op branches have zero extra interaction and keep their autonomous
clock. The pump can supply work; no zero-work preparation is claimed.

The source has dimension D^2, D=2^240 5^8 119^24. Retained charge and spin
references contribute2^30 and2^32, the original filter counter contributes
Mprep=1197721483502453176152084669, the preparation/pump/bin flags
contribute2^3 and the phase pointer contributes2^20. Controller and gate
synthesis workspace are outside this logical storage count.

This result is a source-resolved finite operational centroid with a
complete same-H detector law. It asserts neither one charged pole nor
stationarity of the prepared state, continuous physical spin, a
translation-invariant mass shell, a Standard Model identity or a GeV mass.
