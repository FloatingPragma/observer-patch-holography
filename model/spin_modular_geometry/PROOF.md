# Same-source regional modular and geometric response

## 1. Actual marginal and modular generator

Use the eight vertices and twelve edges of the open cube in
`native_ordered_ground/source.py`, with J=1, B=1/4 and its staggered signs.
Write

\[
 H_B=\sum_{xy}h_{xy}-\tfrac18\sum_x\epsilon_x X_x,
 \quad h_{xy}=(X_xX_y+Y_xY_y+Z_xZ_y-I)/4,
 \quad e_x=\tfrac12\sum_{y\sim x}h_{xy}-\tfrac18\epsilon_xX_x.
\]

The complete 256-column integer eigensystem is accepted by the source's
polar-Weyl certificate. Undoing its Marshall signs gives a normalized
rational-amplitude proposal psi_hat within vector distance epsilon of the
unique actual selected ground Omega_B. The executable checks sum e_x=H_B
as an exact integer matrix equality. This finite reference is separate from
an infinite ordered state or a zero-field ground-state limit.

For A={0,1}, rho_hat=Tr_out |psi_hat><psi_hat| is rational. Four positive
leading principal minors prove rho_hat>I/40. Trace-distance contraction and
the certified vector error give

\[
 \|\rho_A-\widehat\rho_A\|\le\epsilon,
 \qquad \rho_A\ge I/41.
\]

No eigenvalue clipping or thermal reference enters this restriction. For a
faithful density rho, K=-log rho and delta S=Tr(delta rho K). The finite
remainder for halftrace distance eta<mu and rho>=mu I is

\[
 \Delta S=\operatorname{Tr}(\delta\rho K)-D(\rho+\delta\rho\|\rho),
 \quad 0\le D\le\eta^2/(\mu-\eta).
\]

The logarithm is evaluated by the convergent series sum_(n>=1)(I-rho)^n/n.
The executable rounds each entry to a grid 2^-b, b=96. Its matrix rounding
error is d=2/2^b. It recursively multiplies by Q=I-rho_hat_rounded, symmetrizes
and rounds. After n steps the power error is at most nd. Dividing each power
by n and rounding contributes a total at most 2Rd through order R=1024.
With q=1-1/40+d the tail is at most q^(R+1)/((R+1)(1-q)). The integral formula
for the logarithm bounds the reference change by
(epsilon+d)/(1/40-epsilon-d). Their sum kappa is below 6.05e-12.

All sixteen Pauli coefficients are c_P=Tr(K_hat P)/4. The dual norm of
(XX-YY)/4 is one, so the true coefficient difference differs from its
rational value by at most kappa. The certificate gives c_XX-c_YY<-0.13.
Every scalar multiple of the isotropic internal bond, plus arbitrary
one-site fields and a constant, has c_XX=c_YY. This particular physical
boost ansatz is therefore excluded. Nonlocal modular generators and other
continuum limits are outside this finite exclusion.

The native protected Pauli encoding is X->XX, Y->YX, Z->ZI. A symmetric
product of all nonzero K_hat terms, including its scalar term, implements
exp(-itK_hat). For S=sum_P |c_P| and n subdivisions, the product error is at
most S^3 |t|^3/(3n^2): the first two derivatives agree with the exponential,
and each third derivative is bounded by S^3. Literal-angle and full native
generator errors are added, along with |t|kappa. At t=1/64, n=8 the complete
error is below 2.35e-6. Four retained code-reference columns test the full
channel. The operator estimate permits any exterior reference. The modular
parameter t is not identified with the physical H_B clock or a Lorentz boost.

## 2. Source-current pulse and the static-area obstruction

Let V=i[H_B,e_0]. Its exact Pauli expansion is synthesized from the full
source action. Integer row norms give ||V||<=3/2. For
rho(lambda)=Tr_out exp(-i lambda V)|Omega_B><Omega_B|exp(i lambda V),

\[
 \dot S_A(0)=\langle i[V,K_A\otimes I]\rangle_{\Omega_B}.
\]

The parameter lambda specifies a driven native control. The experiment
retains H_B as its unperturbed energy reference, and supplies no identity
between this controlled pulse and autonomous H_B time evolution.

The rational tangent is evaluated before any restriction of the action.
The entropy slope error is at most 2||V||kappa+16||V||epsilon, since
||K_A||<=log41<4. Its certified value is below -0.053. Meanwhile
d< H_B >/d lambda at zero vanishes exactly in the actual H_B ground. This
last equality asserts first-order stationarity; a finite pulse can inject
energy at second order. The assigned density e_0 includes bonds across the
regional boundary and is not a one-region operator.

For lambda>0 with eta=lambda||V||<mu=1/41, the expectation Taylor estimate
and the relative-entropy bound give

\[
 |\Delta S_A-\lambda\dot S_A(0)|
 \le 8\lambda^2\|V\|^2+\eta^2/(\mu-\eta).
\]

At lambda=1/10000 the resulting upper bound for Delta S_A is below
-4.2938e-6. The executable uses the native source ground-preparation prefix
from blank, removes only its final X-pointer measurement, then executes the
symmetric Pauli current pulse on all 2^16 occupation rows. Every first
failure is retained. The preparation bound includes its original readout
error as a conservative overestimate; the ground-vector error is also
charged separately when comparing the rational input with Omega_B.

If eta_out bounds the complete output halftrace error, decode all failure
and off-code components into I_4/4. This fixed trace-preserving map retains
the protected block, decodes its spin factors and traces the exterior.
Contractivity gives
the same bound for the decoded density. Both endpoint densities and their
segment have eigenvalues above 1/100. Integration of the entropy derivative
then gives entropy error <=10 eta_out. The finite negative entropy witness
persists below -4e-6 for the actual decoded native output.

The same experiment is also measured by an actual native X_0 pointer.
Two complete histories start from blank, prepare the same reference, and
apply either the identity or the current pulse before the pointer coupling.
For a=Tr(dot rho_A X_0), the ideal plus-probability change differs from
lambda a/2 by at most lambda^2 ||V||^2. The rational a has error at most
4||V||epsilon. Subtracting both complete native output errors gives a
strict plus-probability increase above 17e-6. Both histories include all
2^17 work rows, both pointer outcomes, reserved pointer labels and all
first failures. This contrast uses unconditioned success probabilities.

Any area depending only on the unchanged cube graph and the fixed cut A
has delta area=0. Consequently Delta S_A+Delta area/(4G)<0 for this pulse,
for every fixed G>0. This disproves universal stationarity for that explicit
area assignment and pulse class. A geometrically restricted class of
variations, a responsive area, and an independently derived bulk/edge split
require further source results. Defining area=-4G S_A would assume the
desired relation and supplies none of those results.

## 3. Complete detector response metric

Four native fair coins select the sixteen two-spin Pauli strings uniformly
in the ideal instrument. Conditional parity measurements retain both
outcomes. For II the plus probability is one and the minus probability is
zero; both record labels are retained. The literal selector and all sixteen
conditional instruments are executed with all64 exterior columns of the
certified eight-source vector. The contraction kernel evaluates two groups
of32 columns without normalizing either group. Their outcome weights add.
The selector failures, all conditional first failures and the four reserved
six-level pointer labels form the complete output. No success division is
performed.

For an arbitrary trace-zero Hermitian regional tangent delta, write
r_P=Tr(rho P), a_P=Tr(delta P). The ideal joint probabilities and derivatives
are (1+-r_P)/32 and +-a_P/32. Their classical Fisher metric is

\[
 g_\rho(\delta,\delta)=\tfrac1{16}\sum_{P\ne I}
             \frac{a_P^2}{1-r_P^2}.
\]

If rho>=mu I on four dimensions, each nonidentity Pauli has two-dimensional
plus/minus spaces. Thus |r_P|<=1-4mu. Pauli orthogonality gives
sum_(P!=I) a_P^2=4Tr(delta^2), and hence

\[
 \tfrac14\operatorname{Tr}\delta^2\le g_\rho(\delta,\delta)
 \le\frac{\operatorname{Tr}\delta^2}{32\mu(1-2\mu)}.
\]

The metric is positive for every nonzero regional tangent. For the actual
source-current tangent, |r-r_hat|<=2epsilon and
|a-a_hat|<=4||V||epsilon. Exact interval arithmetic yields the actual source
Fisher interval [0.02126921343149,0.02126921343179]. The native finite
readout has a complete probability error bound below 1.63e-13 at the
reference input. A probability error at one input does not itself bound a
parameter derivative; the Fisher interval is certified from the exact
source tangent and ideal observables, with the native observables supplied
to arbitrary finite accuracy by their operator bounds.

This positive metric belongs to a space of observable density variations.
It supplies no causal cone, Lorentzian signature, spacetime area or Newton
normalization. Those identifications cannot follow from positivity alone.

## 4. Current filtering and finite-theorem attachment

The actual finite current is not conserved: an exact matrix entry of
i[H_B,V] equals 3/16. For source time T>0, define
V_bar(T)=T^-1 integral_0^T V(t)dt. Direct integration gives

\[
 V_{\rm bar}(T)=(e_0(T)-e_0)/T,\quad
 i[H_B,V_{\rm bar}(T)]=(V(T)-V)/T.
\]

The certified bounds ||e_0||<=13/8 and ||V||<=3/2 imply
||V_bar||<=13/(4T), ||i[H_B,V_bar]||<=3/T. This yields an exact approximate
Ward identity, but the current itself converges to zero in norm. Filtering
this finite derivative current cannot construct a nonzero momentum charge.
This statement concerns this specific current; an extensive momentum
construction or a controlled joint thermodynamic limit is not excluded.

`EinsteinBranch/Entropy.lean` separates the first law from the physical
2pi bulk-plus-center split and the identification of that center with an
edge variation. This construction supplies the first conjunct on the actual
faithful cube marginal. `SmallBall.lean` additionally requires an independent
G>0, a geometrical radius, stationarity, the continuum modular kernel and
the fixed-volume area formula. The static graph area fails stationarity on
the admitted pulse class. `Composition.lean` requires the common
four-dimensional chart, symmetric stress and geometry, null balance, Ward
and Bianchi identities, universality and calibration. No inhabitant of
those combined premises is constructed here. The finite observer can read
its source specification; this permission supplies the operators used
above without asserting these further physical identities.
