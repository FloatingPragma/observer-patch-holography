# Dressed weak-current records on the actual native source

## Source, references and operator

Each source is the complete eight-cell Hamiltonian
`h=H0+2^-34 T` of `chiral_joint/interacting_ground_preparation`, with
`D=2^240 5^8 119^24`. The autonomous Hamiltonian is

```
Htotal=h_active+h_partner+4 N_R.
```

The second source and every spin reference and preparation counter are
retained. R is the complete thirty-mode static CAR reference of
`interacting_charged_response`; all `2^30` occupations are retained.
An eight-state position register X has zero Hamiltonian. Conditional on
`X=x`, the reference current is added to the source Gauss generator at x.
This direct-sum definition retains every Gauss constraint on the complete
carrier. The reference vacuum has zero current and energy. The unique
physical source ground g is unchanged, but the zero-energy position
register makes the full source-plus-apparatus ground degenerate.

The original physical two-copy preparation has success q bounded below
by its frozen exact certificate. Conditional on success, its reduced
two-source density rho12 satisfies

```
D(rho12, |g,g><g,g|) <= delta,
delta <= 2.35592248668e-7.
```

Here `D(rho,sigma)=||rho-sigma||_1/2`. This follows from the certified
ground-pair fidelity, before taking either marginal. The actual state
can be correlated with the retained preparation references. No additional
ground projection is applied and no stationarity of rho12 is assumed.
All initial preparation failures remain separate outcomes.

Use the original scalar/CAR polynomial exported by
`model/native_weak_current_interface`. At cell x it is

```
K_ij=phi_i phi_j^dagger,  K_a=tr(sigma_a K),
j_a=(1/2) sum_beta c_(beta,L)^dagger sigma_a c_(beta,L),
W_x=sum_a K_a j_a.
```

The exact identity `sum_beta J_beta^dagger J_beta=(tr K/2)N_L+W_x`,
with `J_beta=phi^dagger L_beta`, fixes the dressing and normalization.
All twelve Gauss commutators and all three native spin commutators vanish.
W_x is even and has raw L zero. Its full spectrum is
`{-1,-1/2,0,3/2,2}`. The eighty-state active scalar/lepton matrix has
projector ranks `(4,12,58,4,2)`, and the other twenty-six CAR modes are
identity spectators. Thus `||W_x||=2`, its spectral range is 3, and the
spectral range of W_x² is 4. These full-carrier bounds also hold after
tensoring all links and references.

Write `W_tilde=sum_x W_x tensor |x><x|`. The neutral setting prepares X
in its normalized uniform state. Actual antiperiodic translations send
W_x to W_(x+a): the two fermion phases cancel. The position register
transforms by the corresponding eight-label permutation. Hence the
neutral W_tilde apparatus is covariant under the actual native
translation group, including its parity-valued central relations.
Proper cubic operations permute x and rotate the summed native spin
index, under which W is scalar. The local spin statement does not assign
a continuous physical particle spin to an excited native eigenspace.
The charged setting uses the origin basis state of X and its static
compensator. It records this apparatus location; it makes no assertion
that the fixed-origin charged preparation conserves source momentum.

## Exact atomic contractions and their full-source transfer

The source provider supplies the complete sixteen-dimensional physical
C54 block for each local ground magnetic component. `local.py` maps W,
W², the composite density and its scalar term from the original CAR
polynomial using the actual particle-hole convention. For every operator
it verifies full image inclusion in the supplied basis and the metric
adjoint identity. No image is projected away. It also checks that the
mapped W² is the square of mapped W.

For each complete source Hamiltonian H with square-free annihilator p,
the simple least root e has projector

```
P_e=p(H)/((H-e)p'(e)),
[H^i] P_e = sum_(j=i+1)^degree(p) p_j e^(j-i-1)/p'(e).
```

This expression denotes the polynomial quotient, including at H=e.
Arithmetic is exact in `QQ[e]` modulo
`[e(e-1)]^3-313[e(e-1)]²+31708e(e-1)-1041468`.
The trace of P_e is exactly one in each magnetic block. Every listed
moment is `tr(P_e O)` and the two magnetic results agree exactly.
The complete five W-projectors are checked by polynomial multiplication.
Their probabilities sum to one and reproduce both moments exactly.

The local ground compression of W and W² is scalar on the doublet:
the diagonal values agree and magnetic selection kills off-diagonal
elements. Consequently the entangled eight-cell atomic singlet trial u
has the same one-cell moments. Exact root isolation gives

```
<W>_u  = 0.0501054229892607...,
<W²>_u = 0.1549753832212833...,
Var_u(W)=0.1524648298083505....
```

The frozen full-H bound is `||g-u||<=eta`, after phase alignment, with
`eta<=0.000169319250545`. Pure-state trace distance is at most eta.
For a Hermitian operator with spectral range R,
`|tr(O(rho-sigma))|<=R D(rho,sigma)`, by subtracting the spectral midpoint.
Thus the full interacting ground moments lie in

```
m_g in m_u + [-3eta,3eta],
s_g in s_u + [-4eta,4eta].
```

Subtracting the largest allowed m_g² from the lower s_g bounds the
actual variance. This lower bound is strictly positive. Every full W
projector probability differs from its atomic value by at most eta.
The actual prepared marginal contributes at most another delta to each
record, 3delta to the first moment and 4delta to the second. These bounds
hold at every source time because g is stationary and trace distance is
unitarily invariant. They do not replace the actual prepared state by g.

For unique source ground g, the spectral inelastic weight of the bounded
current pump W_tilde/2 is

```
||(I-P_g) W_tilde (g tensor uniform_X)/2||²
  = (<W²>_g-<W>_g²)/4.
```

Translation invariance of g makes all local means equal; projection
onto `P_g tensor I_X` leaves `m_g g tensor uniform_X`. The apparatus
ground degeneracy is therefore included in the elastic outcome.

## A complete finite inelastic experiment

The pump success operator is C_W=W_tilde/2, and its failure operator is
`F=sqrt(I-W_tilde²/4)`. Both are bounded physical operators on the full
carrier. After success, a fresh M-state pointer performs the Fourier
instrument of

```
C=(I-i h/B)(I+i h/B)^(-1), Z=C_active tensor C_partner^dagger,
K_k=(1/M) sum_(j=0)^(M-1) exp(2 pi i j k/M) Z^j.
```

Use the exact B, physical source gap Delta and M from the same physical
preparation certificate. `B>=||h||` includes the absolute energy origin.
For any excited physical source energy E, the Cayley chord from Eg obeys

```
|C(E) conjugate(C(Eg))-1|
 = (2|E-Eg|/B)/sqrt((1+E²/B²)(1+Eg²/B²)) >= Delta/B.
```

The geometric series therefore gives `|K_0(E,Eg)|<=2B/(M Delta)=epsilon`,
with `epsilon<=1/1024`. The ground-pair phase is exactly one, so K_0=1
there and every nonzero pointer label has zero amplitude there. No
unknown Eg is supplied as a free scalar phase. The input partner is an
actual retained source, whose preparation error is part of delta.

Let `e=m_g²/4` and `i=Var_g(W)/4`. For the exact ground-pair comparison,
the zero-label probability is in `[e,e+epsilon² i]`, and the sum of all
nonzero labels is in `[(1-epsilon²)i,i]`. Applying the complete flagged
instrument to the actual prepared pair changes each event probability
by at most delta. In particular the exact-rational lower bound on
the original-attempt nonzero-label probability is

```
q_lower * [(1-epsilon²) Var_lower/4 - delta] > 0.
```

This is an actual finite-current detector record. On the imperfect input,
a nonzero pointer label need not uniquely certify a one-source
excitation; the stated contamination allowance covers both sources.
The complete unnormalized law is

```
1-q,
q tr(F rho12 F^dagger),
q tr(K_k C_W rho12 C_W^dagger K_k^dagger), every k=0,...,M-1.
```

The reference and position factors are implicit identity extensions of
rho12. Fourier orthogonality gives `sum_k K_k^dagger K_k=I`; pump
completeness then proves normalization. Every initial and final pointer
label is retained. The finite Cayley inverse and state synthesis are
declared algebraic services; no efficient native gate synthesis or
execution of the full source-pointer matrix is claimed.

## An energy-selected charged insertion on the same Hamiltonian

The charged setting uses the native compensated e^c removal pump
`B_ec=sum_spin r_spin^dagger c_(origin,spin,ec)/sqrt(2)` on reference
vacuum, with the complete-carrier completion of
`interacting_charged_response`. Its success is 1/2 for exact g and its
reference energy is exactly 4. The full static reference has its
declared actual gauge/spin/raw-L representation. The selected bin I is
the frozen width-1/512 interval for `A=h_active+4N_R`. Its entire spectral
complement is retained. The readout uses the original native source;
the affine coupled architecture is a distinct Hamiltonian.

The exact atomic charged least root x satisfies

```
x^4-2x^3-212x²+213x+11114=0.
```

The complete charged target blocks at twice magnetic spin -2,0,2 have
dimensions 7,17,7. `charged.py` verifies all four original CAR removal
images from both magnetic source-ground blocks, complete W and W²
images, and their metric adjoints. The degree-nine square-free target
annihilator gives a rank-one projector at x in each block. Its W, W²
and five projector moments agree in all three blocks. Thus
`P_x W P_x=m_c P_x` in each entire target block. This identity applies
to every projected removal image regardless of its amplitude. Orthogonal
reference spin labels remove cross terms between removal components.
The source one-cell reduced state is I/2 on the ground doublet.
Therefore, for the complete selected atomic pump vector phi0,

```
||phi0||²=w0,  <phi0,W phi0>=w0 m_c,
m_c=0.0915665020745961... >0.
```

The raw w0 is about 0.419951524978232. It is the unnormalized bin probability,
including the pump's factor 1/2; there is no additional factor 1/2.
Every contraction is bounded by exact rational root isolation. The
selected current value is a density insertion, with no charged ground
or single-pole assertion.

For the exact native g, let `chi=P_I B_ec g`. The frozen complete-support
resolvent transfer bounds `||chi-phi0||<=nu`. It retains every other
atomic root and the full nonblank-link complement of h. Let a0 be the
selected atomic total energy, including the other seven cells and
reference energy 4. Since `||(A-a0)phi0||<=v||phi0||`, `v=360*2^-34`,

```
||exp(-itA)chi-exp(-ita0)phi0|| <= nu+|t|v sqrt(w0)=nu_t.
```

The declared pump and bin-admission services each take one source unit,
so the current record refers to t=2. This is a short-time bound; no
long finite-clock evolution is replaced by this comparison.
With `w=||chi||²` in its complete-source interval,

```
|<chi(t),W chi(t)>-w0 m_c| <= 2 nu_t (sqrt(w)+sqrt(w0)).
```

For each W-projector the same bound has coefficient 1. Dividing the
resulting joint intervals by the positive w interval yields conditional
current intervals. The resulting lower mean is positive.

For actual prepared rho, let S and T be the unnormalized admitted
densities from rho and |g><g|, with traces s and w. The complete flagged
instrument implies

```
||S-T||_1+|s-w| <= 2delta,
D(S/s,T/w) <= delta/w.
```

Indeed `S/s-T/w=(S-T)/w+S(1/s-1/w)` and `||S||_1=s`.
Thus each conditional projector probability contributes at most
`delta/w_lower`, and the mean contributes at most `3delta/w_lower`.
The unnormalized insertion contributes at most 4delta from the full
Hermitian norm bound. The actual bin count is in
`[w_lower-delta,w_upper+delta]`, multiplied by q on original attempts.

The complete ordered law is

```
1-q,
q tr(K_fail rho K_fail^dagger),
q tr((I-P_I) B_ec rho B_ec^dagger (I-P_I)),
q tr(P_w U(2) P_I B_ec rho B_ec^dagger P_I U(2)^dagger P_w),
  for every w in {-1,-1/2,0,3/2,2}.
```

W does not commute with A. The current measurement follows energy
admission and can change the subsequent energy distribution. No joint
commuting energy/current probability is asserted. The alternative
binary readout `E_plus/minus=(I plus/minus W/2)/2` has conditional plus
probability `1/2+<W>/4`; it is a separate square-root instrument, not
an identification with the five-outcome project's backaction.

## Finite-duration physical apparatus and scope

For the five-outcome current readout use
`U_read=sum_w P_w tensor S_5^index(w)` and retain all pointer outputs.
On the full carrier choose its finite Hermitian logarithm L. Over one
source unit set

```
V(s)=exp(-i Htotal s) L exp(i Htotal s).
```

The total evolution is `exp(-i Htotal) U_read`. It records the current
at the interval start. The same construction supplies the bounded
Halmos pump completion and the two-outcome bin flag. Their gates commute
with every combined Gauss constraint and parity. The time-dependent
interaction also does so, because Htotal preserves those constraints.
It adds a finite apparatus control, not a change of the autonomous h.
Every failure/no-operation branch undergoes the same autonomous Htotal.
After the neutral pump, the extra free factor commutes with every Cayley
filter. In the charged setting the two preceding free units are exactly
the t=2 included in nu_t. The free factor after the current pointer does
not alter its recorded probabilities. Apparatus work is not set to zero;
the source-energy change is bounded by `2||Htotal||` per service.

A common storage upper bound for these alternative settings is
`D² 2^30 2^32 Mprep² 8 5 2³`. It includes two complete sources, the full
static charge reference, preparation spin references, both finite
Cayley counters, position reference, current pointer and flags.
It excludes controller and synthesis workspace. The exact bin service
has its own declared spectral implementation premise. These results
give same-H dressed-current expectations, nonzero inelastic records
and energy-selected charged density insertions. They give no isolated
W-boson line, Lorentz-vector identification, weak cross section,
continuum limit, Standard Model matching or energy calibration in GeV.
