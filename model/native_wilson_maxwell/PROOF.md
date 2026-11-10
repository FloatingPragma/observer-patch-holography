# A finite spinor-current and transverse-field response

## Native source and declared spatial member

Each of eight vertices of the contractible open cube contains four copies
of the complete native eight-cell source h=H0+2^-34 T, with the full
old/dual/port reference of `native_abelian_field_join`. All240 fermion
modes, eight five-state scalars,24 PW119 links and the2^31+1-state local
reference of each copy are retained. Its actual neutral ground is
g tensor refvac at energy Eg. The source selector is the positive
elementary/Higgs cross-residue maximizer in the radius2v absolute bin of
`native_higgs_filtered_state`. Write m=E_H,sharp-Eg>0. Its bounds are
source quantities; no measured mass or coupling is inserted.

The two native charged G1 doublets at that energy and their opposite
charge partners supply the complete sixteen-occupation code at each
vertex. Two copies of each charge encode vacuum, both single occupations
and the pair singlet. The exact CAR, parity and cubic exterior-action
identities are those replayed by `native_block_dirac`. Let P be the global
16^8-dimensional code projector inside the full native physical carrier.
P commutes with the sum H_native of the32 constituent Hamiltonians.
The code annihilators a_xs,b_xs carry the actual outer port charges +1
and -1, respectively. Matrix units are selected at E_H,sharp, not at
the distinct frozen maximum-elementary-residue E*.

Assign spacing a=1/m, in the stated c=hbar=1 coordinate convention. This
is a declared length/time comparison. The open box has eight vertices,
twelve directed edges and six oriented faces. The executed cochains in
`field.py` have D_e,tail=-1, D_e,head=+1, CD=0, rank D=7 and
rank C=5. Their Green matrix satisfies
`D^T D G=I-11^T/8`, `G1=0`. The transverse stiffness has eigenvalues
4 (multiplicity3) and6 (multiplicity2); there are no harmonic one-forms.
All these matrices are consumed directly, rather than inferred from a
similar graph or continuum formula.

## Open Wilson gap, complete charge frame and free ground

Let beta and alpha_j be the four Hermitian Dirac matrices. On the full
32-dimensional site/spinor space, define

```
H_W=A+beta M,
A=sum_j alpha_j (T_j-T_j^dagger)/(2ia),
M=m I+L_graph/(2a),
B_e=-(beta+i alpha_axis)/(2a).
```

The diagonal Wilson coefficient at x is m+degree(x)/(2a); the tail/head
block is B_e and its adjoint is the head/tail block. Thus M>=mI,
M commutes with beta, A is Hermitian and anticommutes with beta. It is
unnecessary for M to commute with A. For any vector v,

```
Re <beta v,H_W v>=<v,Mv>>=m||v||^2,
||H_W v||>=m||v||.
```

Hence |H_W|>=m. The actual open-cube matrix additionally obeys the exact
charge-conjugation identity
`C H_W* C^dagger=-H_W`, `C=I_sites tensor i beta alpha_y`.
There are16 positive and16 negative modes. The producer verifies the
full characteristic polynomial and degree12 annihilator. The positive
energies in units of m are sqrt(3)+/-1/2, sqrt(7)+/-1/2 and
sqrt(13)+/-1/2, with multiplicities2,4 and2 for each sign.

Write P_+,P_- for the complete positive/negative one-body projectors and
E_+,E_- for the embeddings of the upper/lower beta subspaces. The frames

```
U=P_+ E_+ (E_+^dagger P_+ E_+)^(-1/2),
V=P_- E_- (E_-^dagger P_- E_-)^(-1/2)
```

are well defined. For example, a nonzero upper vector killed by P_+
would lie in the negative spectral space while having positive energy
expectation <v,Mv>, a contradiction. The analogous lower argument gives
the other positive Gram matrix. Thus [U,V] is unitary. Set
Omega_+=U^dagger H_W U>=m and Omega_-=-V^dagger H_W V>=m. Define all
four canonical site fields by

```
Psi=U a+V b^dagger,
H_F=a^dagger Omega_+ a+b^dagger Omega_-^T b,
rho_x=sum_(c=1..4) Psi_xc^dagger Psi_xc-2.
```

The complete CAR gives commuting integer densities rho_x in[-2,2] on
the code and `sum_x rho_x=N_a-N_b`. Normal ordering relative to the
a/b vacuum gives `:Psi^dagger H_W Psi:=H_F`. That vacuum is the filled
negative-band state in Psi coordinates. It is not the empty Psi Fock
state. The local rho_x differs from the original native onsite
N_a,x-N_b,x; its particle/antiparticle expansion retains pair terms.

The field change is implemented by a finite charge-preserving Fock
unitary U_FW on the code. Use site-major canonical variables
Psi0=(a_xup,a_xdown,b_xup^dagger,b_xdown^dagger)_x, with the conjugate G1
basis for the opposite-charge creators. The matrix[U,V] has grouped
input columns (all a, then all b^dagger), so its column permutation
must be retained. Define the conjugation convention explicitly by

```
M_site=U E_+^dagger+V E_-^dagger=[U,V][E_+,E_-]^dagger,
U_FW Psi0 U_FW^dagger=M_site Psi0=U a+V b^dagger=Psi.
```

This convention may use Gamma(M_site^dagger) when Gamma is defined
by its action on one-particle states. In particular,
`U_FW (Psi0^dagger P_x Psi0) U_FW^dagger
 =Psi0^dagger M_site^dagger P_x M_site Psi0`.
The inverse one-body conjugation is not substituted. An exact full-Fock
pair-rotation control checks both CAR and density identities; a separate
32-dimensional control checks the grouped/site permutation and rejects
its omission. Total
Psi0 number equals N_a-N_b+16, so this unitary commutes with total native
port charge. Extend it as identity on the complete native complement. With
`Q_native,x=sum_(four copies)(P_port_minus-P_port_plus)`, define

```
rho_x=U_FW Q_native,x U_FW^dagger
```

on the full carrier. Every rho_x commutes with P and with the other
densities. On the complement it is the original port charge, with
integer range[-4,4]. The total charge is unchanged on every source state.
All these operators preserve the original native Gauss kernel and
fermion parity. The extension is an explicitly defined finite spectral
and matrix-unit service; no efficient native circuit is asserted.

Add the nonnegative free kinetic operator

```
K_free=P [a^dagger(Omega_+-m)a
           +b^dagger(Omega_-^T-m)b] P.
```

It is zero on the complete complement. H_native+K_free therefore retains
the exact full native product ground and its positive physical gap.
This global spectral-code compression is part of the declared
architecture. It differs from both the uniform native lattice and the
frozen free Dirac member's complete complement extension.

## Finite Maxwell source, Gauss embedding and common Hamiltonian

Keep all56 states with total occupation at most three in the five
transverse modes. With nu=(2,2,2,sqrt6,sqrt6), the completed field
Hamiltonian is

```
H_ph=[Tr(nu)/2+dGamma(nu)]/a.
```

In electric/magnetic quadratures it includes the positive boundary
completion `P_top[dGamma(nu)+Tr(nu)]/(2a)`. The actual finite ladder
matrices, field equations, boundary correction and all occupations are
executed by the shared field provider. The link coordinate is
Q_e=a A_phys,e, and every edge satisfies

```
<Q_e^2>_0=1/16+1/(12sqrt6),
||Q_e||^2=13/48+5sqrt6/48<9/16.
```

Different finite Q_e need not commute at the top occupation boundary.
Each exp(i g Q_e) is its literal finite matrix exponential. No Weyl
identity or full-carrier commuting coordinate law is used.

At reference vertex r=0 retain R=-32,...,32 with Hamiltonian4|R|.
At each other vertex retain ell_x=-4,...,4. The normalized constraints
are ell_x=rho_x and R=-sum rho. They are independent of g, so g=0
does not enlarge the physical space. Let Y coherently attach these
joint eigenlabels, retaining their superpositions. The physical electric
flux has longitudinal part

```
P_L=g D G rho_ref,
rho_ref=rho+R delta_r,
H_C=g^2 rho_ref^T G rho_ref/(2a).
```

The actual electric divergence equals g times the full compensated
charge. Its transverse part has zero divergence on every finite photon
state. The normalized Gauss constraints are commuting operators on the
full finite carrier. All allowed labels and their complements remain.
Y is an isometry onto the outer integer-Gauss kernel and intertwines P, every
specified source observable and its full physical counterpart. This
label construction proves exact Gauss; it does not infer a canonical
continuum gauge-translation algebra at finite photon cutoff.
The native Gauss constraints are separate and are preserved by Y.
If their full kinematic sectors are included in the source carrier,
native-constraint-violating states stay in the native-code complement;
Y does not remove them. They have zero probability in the stated
native physical preparation. The additional readout failure label
below refers specifically to the outer integer constraints.

On the logical native/photon carrier declare

```
H(g)=H_native+K_free+H_ph+4|R|+H_C
       +P [H_W(Peierls)-H_W(0)]_Fock P.
```

The +phase exp(i g Q_e) multiplies
`Psi_head^dagger B_e^dagger Psi_tail`; its adjoint phase multiplies
`Psi_tail^dagger B_e Psi_head`. Full pairs of terms are retained. The
matrix derivative is J_phase,e=partial H/partial theta_e. The oriented
charge flow is I_flow,e=-J_phase,e, and
`dot rho=D^T I_flow` on the code. The Coulomb term commutes with each
density. In the physical program every operator is lifted through Y,
including the unchanged native and kinetic terms. Its clock is this
single H(g); it is not the ungauged free clock. On the complete kinematic
carrier choose exactly H_ext=Y H(g)Y^dagger direct_sum zero on the
outer-Gauss-orthogonal complement. Thus its unitaries extend by identity
there, retaining every constraint-failure record. The native complement
inside the physical kernel retains its specified native Hamiltonian;
it is distinct from that Gauss-orthogonal extension. The exact
label-shift controls include finite register boundaries but are not
assertions of an unguarded alternate Hamiltonian on invalid inputs.

Every term reduces P. On its complement the native h, port charges and
references are retained, the kinetic/current couplers are zero, and the
Coulomb energy is diagonal in the original port charges. No coupled
ground-state identification is made. The observed input below is
nonstationary under H(g). In particular, a free-vacuum response is not
radiation from the interacting ground.

## A preparation with no hidden ground projection

Prepare each of32 native sources by the deterministic unfiltered atomic
cube trial u and its retained spin reference J u. This is the original
native preparation service, with no finite agreement filter and no P0
measurement. All added R_loc registers are initialized in the shared
vacuum and are untouched by that preparation. The exact source ground
overlap obeys |<g,u>|^2>=p0, the frozen rational native bound.

Let Q_vac project all added R_loc factors onto their shared vacuum.
Every nonvacuum occupation-code word contains a charged port. Orthogonal
reference sectors therefore give the full-carrier support identity

```
P Q_vac=|g^32,refvac><g^32,refvac|.
```

The executed9-to4 doublet and81-to16 occupation embeddings check the
selected-space part exactly. On the full native complement it follows
from the same reference-sector orthogonality; no neutral excited state
can supply an additional code word. For the complete trial purification,

```
P |trial,all R_loc=vac,spinrefs> = |g^32,refvac> tensor chi,
p_code=||chi||^2>=p0^32.
```

The physical preparation is the lift by Y. Since Y intertwines P, the
same identity holds after Gauss dressing. The Bogoliubov unitary defines
rho and the readout, not an extra preparation applied to this trial;
it does not change the support identity. The initial empty native seed
has all port charges zero and lies outside P, so U_FW is identity there
and Y attaches zero labels. The lifted preparation can act entirely in
the physical kernel, with its full complement fixed.

All H(g) terms and final count effects reduce P. Thus code/complement
coherences cannot contribute to any count probability. The exact actual
law is p_code times the ideal code law, together with complement
probability1-p_code. On this prepared complement all port charges are
zero, so its photon vacuum remains unchanged. This is an unconditional
record decomposition; P is not measured or postselected initially.
The complement outcome can be more probable than the rare interaction
outcome without weakening its positive lower bound.

## Actual current coefficient on all fermion occupations

Differentiate the complete Peierls operator at g=0:

```
V1=sum_e Q_e J_phase,e.
```

In the ideal code vacuum Omega, J_phase has pair-creation and
pair-annihilation terms from the canonical Psi field. The target is one
a excitation, one b excitation and one transverse quantum. Its first
derivative norm satisfies the exact finite CAR identity

```
sigma^2/m^2=sum_(e,f) Cov(Q)_ef
      Tr[P_- J_phase,e P_+ J_phase,f],
Cov(Q)=P4/4+P6/(2sqrt6),
```

where the displayed one-body derivatives are in units of m. This is a
calculation on the entire32-mode fermion Fock space by its finite vacuum
contraction, not a deletion of many-body occupations. The free vacuum
current vanishes on all twelve edges by the executed time-reversal
identity. The field covariance is the actual finite-cutoff vacuum word;
two ladder operations never reach its top boundary.

One calculation reconstructs sign(H_W/m) with the exact binomial series
in I-(H_W/5m)^2. It sums indices0 through512, reduces the polynomial
modulo the complete executed degree12 annihilator, and bounds the tail
by `25(24/25)^513`. The open gap and ||H_W/m||<5 are the hypotheses.
The negative projector has half this error. Its Hermitian entries are
rounded to128 dyadic bits, with additional operator error at most
64*2^-128. If epsilon is the sum, the current-square error is bounded
by `144 epsilon(2+epsilon)`: each derivative has Hilbert-Schmidt norm2,
each covariance entry has magnitude at most1/4, and all144 pairs are
retained. Thus

```
1.0435128440992287 < sigma^2/m^2 < 1.0435186313865983.
```

The authoritative endpoints are the rational outward bounds in the
receipt. An independent exact Chinese-remainder calculation on the same
complete matrix gives

```
sigma^2/m^2 = 965/1456+5sqrt6/36-sqrt14/48+sqrt21/336
                +sqrt91/546+7sqrt546/1872.
```

The exact expression lies strictly inside the polynomial enclosure.
Neither calculation uses numerical eigenvalues or eigenvectors, and
neither uses an unrestricted photon coordinate identity.

## Entire-operator error and complete connected records

Fix g=1/64 and dimensionless time tau=mT=1/64. These are declared finite
action/observation choices. On the charge-zero code, every rho_x lies
in[-2,2] and sum rho=0, so ||rho||^2<=32 and ||G||=1/2. Consequently
||H_C||<=8m g^2. Each full-Fock Wilson bond has norm2m, and
||Q_e||<3/4 gives

```
||V1||<=18m,
||H(g)-H(0)||<=m nu,
nu=18|g|+(59/4)g^2.
```

The conservative quadratic coefficient includes both the Coulomb term
and the Peierls Taylor remainder. These are whole charge-zero-code
operator bounds, on every occupation and every retained photon state.
They do not restrict virtual intermediate states to the target mode.
The complete native complement is handled by the exact reducing
decomposition above.

The first-order target contains one pair and one photon, whose free
energy is less than11m. Hence its integrated leading amplitude is at
least `g tau (sigma/m)(1-11tau/2)`. Coulomb and the even Peierls terms
cannot take the photon vacuum into the one-photon target in a single
insertion. The remaining odd first-insertion error is bounded on the
full cutoff carrier by `(27/16) tau |g|^3`, using the finite functional
calculus remainder for sin(gQ)-gQ. Twice Duhamel bounds all further
insertions by `tau^2 nu^2/2`. Thus the ideal target amplitude is at least

```
a_*=g tau (sigma/m)_lower(1-11tau/2)
      -(27/16)tau |g|^3-tau^2 nu^2/2 >0.
```

The actual deterministic-preparation target probability is at least
`p0^32 a_*^2`. There is no trace-distance subtraction after squaring.
The first bound includes the complete nonlinear finite Hamiltonian,
all photon-cap boundary states and every multi-pair virtual state.

Read every (N_a,N_b,N_ph), with N_a,N_b=0,...,16 and N_ph=0,...,3,
and all four photon counts on the native-code complement. These
`(17^2+1)*4=1160` labels resolve the full physical source. Append one
Gauss-orthogonal constraint-failure label, giving a complete kinematic
pointer dimension1161. Its probability is zero for the stated physical
preparation. The target is(1,1,1); all other outcomes,
including(0,0,0), are retained. Charge-zero input forces N_a=N_b in
code but the unused count labels remain in the alphabet.

If E_alpha denotes this complete projector family, the bounded readout
unitary is `U_count=sum_alpha E_alpha tensor Shift^alpha` on the1161-state
cyclic pointer initialized at zero. It commutes with every Gauss
constraint and preserves the full source norm. The constraint-failure
branch shifts only the pointer; its autonomous source Hamiltonian is
the zero extension specified above.

Let X indicate any matter excitation in code and Z indicate nonzero
photon count in code; set both to zero on the explicit complement flag.
The ideal free input lies in each event's kernel. Duhamel gives both
ideal marginal probabilities at most tau^2 nu^2. The actual law therefore
has the strictly positive binary connected response

```
E[XZ]-E[X]E[Z] >= p0^32 [a_*^2-(tau nu)^4] >0.
```

It is a finite nonstationary current/field interaction. There is no claim
of scattering states, vacuum decay, an asymptotic photon emission rate
or a nontrivial four-dimensional continuum limit. A free-photon
refinement and its cutoff limit require their own field theorem; the
shared field result here is the exact finite cube and cap.

## Physical schedule, full controls and locality boundary

The complete source and references retain the same autonomous H(g).
The preparation is a specified finite gauge-preserving unitary S from
the physical blank seed to the full unfiltered trial, including spin
and Gauss references. To prepare that nonstationary state exactly at
the end of a duration-one pulse, declare the driven service

```
L=i Log(S),  V_control=L-H_ext,
exp[-i(H_ext+V_control)]=S,
```

with Hermitian branch phases in(-pi,pi]. This explicitly cancels the
autonomous generator during the prescribed preparation pulse. It is
not an implicit pause or a passive local controller. After the full
autonomous wait T, a duration-one counting pulse targeting U_count can
instead use `L=i Log(exp(-i H_ext)U_count)`. Its actual propagator is
exp(-i H_ext)U_count; the final common source evolution does not change
the retained pointer law. Elapsed source time is2+T on every branch.
Both constructions preserve the exact Gauss constraints. Native,
reference and vacuum energy constants belong to H_ext in these formulas.

For the32 main sources, an upper bound on the full raw physical carrier
is `32(B_native+120)+99m+18m|g|+(923/4)mg^2+128`, with
B_native=480+360*2^-34. It includes the complete reference cost, all
photon states and the worst full charge range[-4,4], not just the
neutral response sector. The pulse control norm is at most this bound
plus pi. Net bare-source energy work is bounded separately by twice
the raw norm; total absolute switching work for two pulses is at most
four times(raw norm+pi). Programmable matrix coefficients and their
additional workspace are declared resources. An implementation error
eta in the two complete pulse isometries changes every unconditioned
target amplitude by at most eta. The robust target lower bound is
`(sqrt(p0^32) a_*-eta)_+^2`, and a robust connected lower bound is this
quantity minus `(tau nu+eta)^4`. Such an error need not preserve the
trial's exact reference-vacuum support, so it cannot be inserted inside
a later factor p0^32 without a separate proof. Exact services here do
not constitute an efficient native gate compilation.

The displayed retained source/reference/pointer capacity in nats is

```
32[240ln2+8ln5+24ln119+ln(2^31+1)+16ln2]
 +7ln9+ln65+ln56+ln1161.
```

No agreement-filter pointers or discarded attempts are hidden in this
depth-zero count. The code Fock dimension is2^32; its full many-body
Hamiltonian is not numerically enumerated. Complete finite CAR traces,
matrix-function error bounds and whole-operator estimates cover it.

The existing finite causal observer theorem requires an actual local
interaction algebra and support bounds. A Wilson matter bond supplies
its full-Fock norm, but the global spectral-code compression, the
Coulomb-gauge charge frame and the total photon occupation projection
do not automatically meet that theorem's local tensor hypotheses.
No source-wide light cone or interacting Lorentz covariance is inferred.
The original native charge and parity supply the code's ancestry;
the regulator, common gauge coupling and reference are declared choices.
The result supplies no Standard Model species, weak chirality or GeV
identification.
