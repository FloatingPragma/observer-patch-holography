# Arrival followed by a finite energy record on the native charged source

## Source, admitted sector and retained preparation

The source is the **whole-bin** Hamiltonian of
`model/native_abelian_field_join/PROOF.md`, with eta=1/24 and no complement
penalty. Every original native intrablock Hamiltonian h, complete local
reference, 52-mode compensator, ternary longitudinal register and all 5456
transverse states are retained. The global compensator stays at vertex 1
and costs 4 in the one-reference sector. Opposite charges have the same
results by the proved photon-parity equivalence. Native G1 is an internal
finite-group representation. This result supplies no continuous spin,
Dirac dispersion, relativistic mass or full charged-sector ground label.

Write g and Eg for the actual native neutral ground and energy. The
one-bin sector has one native charged block in the full source spectral
bin J, all other active blocks in g, the compensating local port and
matched outer longitudinal registers. The fixed global reference has
one fermion. All native states in J remain present, including the
complement of the two normalized vectors

```
psi_s=P_J c_s g/sqrt(w_J),       s=1,2.
```

Here w_J=||P_J c_s g||^2 is the unnormalized per-spin weight. The native
source interval is [L,U], and c_s g has norm squared 1/2. The literal
compensated CAR pump succeeds with probability 1/2 and has bin weight
w_J, without an additional factor 1/2. The exact values and all inherited
source pins occur in the receipt.

The physical native preparation is performed thirteen times with its
partner copies, spin references, filter counters and complete failures
retained. Let Q denote the probability that all thirteen preparation
filters succeed. A separate exact spectral service P0 projects **all
26 native source copies**, including the partners, onto g^26. Its
conditional success probability is p. The parent proves explicit lower
bounds Q>=q_min^13 and p>=f_pair,min^13. P0 has its own retained binary
outcome and controlled-H apparatus cost; it is not an efficient or free
ground preparation. This use is unchanged from the parent source.

On success, the full neutral energy is 26Eg+Evac, with the thirteen
partner energies included. All extra preparation references have their
declared Hamiltonians, including the zero-energy spin references. Their
degeneracies do not imply uniqueness of a source-plus-apparatus ground.
Failed branches retain the same autonomous source Hamiltonian and clock.

## Actual ground in the whole admitted bin

Put Ec=(L+U)/2 and epsilon=(U-L)/2. For comparison only, replace
h-Eg by Ec I on J. The actual source h is unchanged. In the one-bin
sector this replacement has operator norm error at most epsilon. Call
the comparison Hamiltonian Hc and the actual one Ha, both reported above
the same neutral energy. On the psi subspace,

```
Hc = (Ec+2+H_affine,N=1) tensor I_G1.
```

The affine energy includes source cost 1 and reference cost 1. The shift
Ec+2 replaces these by the native excitation Ec and the retained reference
cost 4. Every other bin vector is immobile in Hc, has nonnegative electric
and photon energy, and lies at energy at least Ec+4. The affine first
excited energy exceeds 2.00219970691, while its actual ground EA lies in
the frozen full-70928 residual interval around 1.873809754966166. Thus Hc
has an isolated ground doublet at Ec+2+EA, and its complement is at least
Ec+4. Define

```
d0=2-EA_upper,       d=d0-2epsilon,
theta=epsilon/(d0-epsilon).
```

Minmax under the norm perturbation gives a rank-two lowest cluster for
Ha, ground energies in [L+2+EA_lower,U+2+EA_upper], and separation at
least d>0.1242370781 from every other state **within this admitted
one-bin sector**. All actual native h and the spin-contracted hopping
commute with the proved finite native cubic action. Along Hc+s(Ha-Hc),
0<=s<=1, the cluster remains isolated. Its finite-group character is
continuous and integer valued, hence equals the irreducible G1 character
of the comparison doublet. Schur's lemma makes the rank-two cluster a
single exactly degenerate source energy. A retained two-state compensator
spin factor multiplies its degeneracy. The singlet preparation specifies
a joint invariant vector, without removing that reference from the source.

Let Pc and P be the comparison and actual ground projectors. On the
actual orthogonal complement, Ha-(Ec+2+EA) is at least d0-epsilon.
The residual of Pc under Ha is bounded by epsilon. Functional calculus
therefore gives ||(I-P)Pc||<=theta. Equal finite ranks imply

```
||P-Pc|| <= theta <0.008.
```

No full-source gap is substituted for d, and no lower states outside J
are excluded by this argument. In particular, this is not a proof that
the admitted ground is the full outer-charge ground. No penalty is added.

## A coherent spatial CAR preparation

Let z be the normalized zero-photon part of the exact binary affine trial
stored in the frozen full-spectrum package. It has thirteen spatial
components. For fixed native spin s, let d_(y,s) be its literal native
charge operator with local-port vacuum projector and creation, appropriate
longitudinal shift and reference-vacuum guard. These are odd on distinct
native blocks. The port and longitudinal factors are even; different
vertices' longitudinal shifts commute. Hence distinct-site d operators
anticommute, and on the complete kinematic carrier

```
{d_(y,s),d_(y,s)^dagger} <= I,
A_s=sum_y z_y d_(y,s),       {A_s,A_s^dagger}<=I.
```

The inequality on one site follows directly from its two orthogonal
local-reference blocks: its two nonzero diagonal effects are c_s^dagger
c_s and c_s c_s^dagger, each at most I, with boundary-zero shift effects
bounded by I. The fixed global-reference creation is guarded by its
vacuum projector. Thus

```
B_z=(1/sqrt(2)) sum_s r_(-q,s)^dagger A_s
```

is a contraction on the full carrier, since its vacuum-reference effect
is (sum_s A_s^dagger A_s)/2<=I. Its Halmos unitary completion is an even,
Gauss-preserving, finite matrix-functional-calculus apparatus service.
Reference and native spin indices are contracted with the proved G1
singlet intertwiner. All complementary effects and pump failure remain.

On g^26, orthogonality of charge positions and the native scalar Gram
identity give total pump weight 1/2 and bin weight w_J. Conditioned on the
bin, the spatial/photon state is z tensor photon vacuum and the source/
compensator spins are the retained singlet. Let f0 be the certified affine
ground fidelity of this zero-photon spatial pump. The actual admitted
ground fidelity is at least

```
F=1-[sqrt(1-f0)+theta]^2 >0.9997343382.
```

This follows by the triangle inequality for (I-P) acting on the normalized
prepared state. The spatially coherent preparation is a distinct setting
from the localized pump used in the arrival experiment below.

## Direct finite clock and explicit raw-energy phase

The final energy detector acts directly on Ha, through the full source
Hraw. No long-clock Duhamel comparison with Hc is used. The native local
ground interval [e0lo,e0hi] and ||lambda T||<=45/2147483648 imply

```
Eg_center=4(e0lo+e0hi),
Eg_radius=4(e0hi-e0lo)+45/2147483648,
|Eg-Eg_center|<=Eg_radius,       -90<Eg<0.
```

The pointer has the explicit scalar phase
B=26Eg_center+Evac. Evac is the same exact finite transverse zero-point
energy used in Hraw. This is an apparatus phase; it does not subtract
source or reference energy from the autonomous clock. Relative to B the
actual admitted ground lies in radius

```
r=epsilon+10^-12+26Eg_radius
```

around c=Ec+2+lambda_affine, where lambda_affine is the frozen exact
binary trial energy center. All other admitted energies are at least
L+4-26Eg_radius. The full admitted shifted spectrum is in (0,15), from
EC<=6/31, Hph<9, and ||V||<3/10. In particular there is no unbounded
oscillator or discarded photon state. The corresponding raw sector and
the scalar B satisfy |Hraw|<2400 and |B|<2400.

Use M=2^20 pointer labels, t0=pi/8, and accept the integer labels with
|16k/M-c|<=0.03. The ideal controlled action and Kraus operators are

```
C=sum_j |j><j| tensor exp(-ij t0 Hraw) exp(+ij t0 B),
K_k=(1/M) sum_j exp(2pi i jk/M) exp[-ij t0(Hraw-B)].
```

The positive Fourier sign gives k for energy 16k/M above B.
Fourier orthogonality proves sum_k K_k^dagger K_k=I on the entire admitted
sector. Extend the extra controlled interaction by zero outside the
physical one-bin sector, so its controlled gate is identity there. Every
branch retains common autonomous exp(-i T Hraw), including preparation,
P0, pump, bin and position failures. The finite projectors selecting the
native admission are declared spectral services, as in the parent.

The inherited circle-grid proof gives

```
BM(delta)=min(1,1/(2 M^2 delta^2)+1/(2 M delta)),
delta_g=(0.03-r)/16,
delta_e=min((L+4-26Eg_radius-c-0.03)/16,
            1-(15-c+0.03)/16).
```

Both distances are positive. The second term controls wraparound.
Ideal ground rejection is at most BM(delta_g), and excited acceptance
at most BM(delta_e), uniformly on the actual admitted source eigenstates.

An explicit finite rational realization uses

```
a=I-X^2/10, b=X/2-X^3/120,
R33(X)=(a+ib)^-1(a-ib),
Uhat=R33(t0 (Hraw-B)/n)^n,       n=1024.
```

The exact denominator a^2+b^2=I+X^2/20+X^4/600+X^6/14400 is at least I.
The independently replayed integral identity gives a one-gate error
|t ||A|| |^7/(100800 n^6). Consequently the complete controlled pointer
has operator error at most

```
gamma=(M-1)*(165/28)^7/(100800 n^6)<3*10^-12.
```

Its functions of Hraw preserve every actual energy eigenspace and every
Gauss constraint. Both ideal and rational extra gates are identity off
the admitted sector; the sector norm is never applied to another source
sector. Define alpha=BM(delta_g)+gamma and beta=BM(delta_e)+gamma.
Applying the operator/trace-distance bound separately to each eigenstate
gives ground rejection at most alpha and excited acceptance at most beta.
For any unnormalized admitted input chi with ground weight g and excited
weight e, accepted ground weight is at least (1-alpha)g and accepted
excited weight at most beta e. The bound scales with that input's actual
weight, including rare arrival branches and any retained reference.

For the coherent spatial pump, g>=F and e<=1-F. Conditional on acceptance,
its ground infidelity is at most

```
(1-F) beta / [F(1-alpha)+(1-F) beta] <2.13021*10^-8.
```

This uses the bounded **extra control generator** Hraw-B. The scalar B
is the explicit source-derived pointer phase, and the common autonomous
Hraw evolution remains present. It does not replace the autonomous source
by a shifted generator. For comparison, the receipt also verifies the
separate raw/scalar implementation
R33(t0 Hraw/nr)^nr R33(-t0 B/nr)^nr with nr=2^19 and error
(M-1)*2*(6600/7)^7/(100800 nr^6)<10^-11. The direct shifted functional
calculus uses 1024 times fewer rational factors, with the same exact
matrix-service premise and explicit phase ancestry.

The exact finite coefficients, linear inverses and controlled source
actions are declared apparatus services. No efficient local circuit or
numerical execution of the full native-plus-pointer matrix is claimed.
The straightforward binary-power realization uses at most
n(M-1)=1073740800 rational factors. The parent retained storage is multiplied by M;
the arrival flag is the parent's existing arrival bit. Gate-synthesis
workspace is not included in that storage count. The inherited spectral
P0 and bin services have not been replaced by this final energy clock.

## Ordered arrival and energy, including measurement backaction

Use the localized native CAR pump at y=1, condition on its bin flag and
call the normalized state Phi. Evolve with Ha for t=1/4 or 1/64, then
measure arrival at x=0 (apex) or x=2 (adjacent boundary vertex). Its
projector Qx is source charge position, with all photon states retained.
It preserves the admitted sector, but need not commute with Ha or P.
The energy clock acts **after** the arrival record on

```
chi=Qx exp(-it Ha) Phi.
```

This is an ordered instrument with measurement backaction. It is not a
simultaneous position/energy projector, and the recorded short-time
position is not a claim about the position after the long energy clock.

The binary arrival readout has an explicit finite-duration realization.
On the admitted sector use

```
Uread=(I-Qx) tensor I + Qx tensor X,
L=(pi/2) Qx tensor (I-X),       exp(-iL)=Uread,
V(s)=exp(-is Hraw) L exp(+is Hraw),       0<=s<=1.
```

Extend L by zero on all complementary and no-op branches. The complete
Hamiltonian Hraw+V(s) has propagator exp(-i Hraw) Uread over this one
source unit. This follows directly in the interaction picture, and does
not assume [Qx,Hraw]=0. Every Gauss constraint and the admitted sector
are preserved. The recorded position refers to the start of the readout
interval. Its conditional arrived state is exp(-i Hraw) chi. Since the
subsequent finite energy Kraus operators and P commute with Hraw, this
extra common free factor leaves every probability and fidelity bound
below unchanged. The source clock runs on every branch. The controller
and its work are explicit apparatus resources: the readout can change
the source energy expectation, by at most 15 source units within the
admitted shifted spectrum. Synthesis and controller workspace are
separate from the retained binary logical pointer.

The comparison arrival vector has the same affine source matrix and
short-time Dyson expansion as the frozen mobile certificate. Write R
for its complete all-orders remainder norm and k_e for the exact finite
Peierls vacuum overlap. The first-order vacuum amplitude has lower bound

```
v=eta t k_e,lower [1-(Delta EC)^2 t^2/24],
Delta EC=26/217 for x=0, or 33/217 for x=2.
```

Indeed its time integral has magnitude eta t k_e |sinc(Delta EC t/2)|,
and sinc(u)>=1-u^2/6 in this interval. The first-order nonvacuum norm is
at most n1=eta t sqrt(1-k_e,lower^2), because H0 preserves photon number
and the Peierls matrix is exactly unitary on the full finite cutoff.
The total comparison arrival norm is at most A=eta t+R.

Let vtrial be the exact binary full affine trial normalized by its exact
rational norm, and s_aff its certified sine angle to the actual affine
ground. The actual ground's squared amplitude at |x,phvac> is bounded
below by |vtrial[x,0]|^2-s_aff. Let p_x,upper and q_ph,upper be its frozen
position and nonvacuum probability upper bounds. Cauchy-Schwarz on the
nonvacuum first-order term and on Qx times the remainder proves

```
a_aff=sqrt(|vtrial[x,0]|^2-s_aff) v
      -sqrt(q_ph,upper) n1 -sqrt(p_x,upper) R
 <= ||Pc Qx exp(-it Hc) Phi||.
```

This uses the full actual affine ground. No photon-vacuum compression
eigenvector is substituted. It is valid for each retained native spin
and hence for the entangled reference singlet.

Only the short propagation time uses Duhamel. The actual admitted
ground projector estimate gives

```
a=a_aff-theta A-epsilon t
 <= ||P Qx exp(-it Ha) Phi||,
||Qx exp(-it Ha) Phi|| <= A+epsilon t.
```

The theta factor multiplies the small comparison arrival norm A, not 1.
All four a are positive. Put G=a^2 and E=(A+epsilon t)^2-G. The actual
unnormalized arrived state's ground weight is at least G and its excited
weight at most E. Thus

```
unconditional arrival then accepted true-ground weight
 >= Q_min p_min w_J,lower G(1-alpha),
accepted ground infidelity <= E beta/[G(1-alpha)+E beta].
```

The true-ground component is a mathematically bounded part of the
accepted pointer event, not an additional experimentally measured ideal
projector. The observed arrival-plus-acceptance probability is at least
the displayed true-ground lower bound.

At t=1/4 the unconditional lower bounds are 4.5027825921*10^-6 for the
apex and 1.1346122060*10^-6 for the boundary site. The corresponding
accepted infidelity upper bounds are 0.001041 and 0.004354.
The exact receipt supplies the two t=1/64 settings as well.

## Complete law and source custody

Let a_arr=||Qx exp(-it Ha)Phi||^2, and use the implemented rational-clock
Kraus operators in the following exact finite record law:

```
native preparation failure:        1-Q,
P0 failure:                         Q(1-p),
native CAR pump failure:            Qp/2,
outside native bin:                 Qp(1/2-w_J),
bin nonarrival:                     Qp w_J(1-a_arr),
arrival with energy label k:        Qp w_J ||K_k Qx exp(-it Ha)Phi||^2,
                                    every 0<=k<M.
```

The six grouped entries sum exactly to one by clock completeness. The
preparation-failure entry groups all thirteen retained filter flags and
their counter records; none is erased. Other flags and zero-weight phase
labels remain part of the complete apparatus. For the coherent spatial
preparation setting, omit the position test, retain its unused flag at
zero, and use Qp w_J ||K_k Phi_z||^2 together with the same first four
failure groups. Reference charge and its cost 4 appear on every charged
branch. Outer flux, native gauge constraints, all complementary source
sectors and the common raw autonomous clock remain the parent's ones.

The code rebuilds all displayed rational inequalities from frozen source
receipts and the actual binary trial. Its finite independent controls
check projector perturbation, ordered backaction, Fourier completeness,
rare-state weighting and complete records. These controls do not replace
the source certificates or numerically execute the full native carrier.
