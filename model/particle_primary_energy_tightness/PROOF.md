# Sharper finite energy bands and the operator-norm scaling boundary

All results use the selected primary open Wilson source, its physical
Gauss polynomial density, the retained conjugate reference with idle
Hamiltonian zero, and its source clock. The native dimension-six block
has integer charge \(Q_6=1\). Its scalar beta-pair operator is gauge
neutral:
\[
 C_x={1\over\sqrt{d(d+1)}}\sum_{\alpha=1}^2\sum_{i=1}^{d}
             c^\dagger_{+,x\alpha i}c_{-,x\alpha i},\qquad d=6.
 \tag{1}
\]
The beta occupations and holes below are measured relative to a specified
Dirac matrix. They are not identified with protected topological defects
or measured particle species.

The [refinement response](../particle_primary_refinement_response/PROOF.md)
supplies the exact selected-state bound
\(\langle D_{\rm total}\rangle/V<9.657986<97/10\), uniformly in the finite
primary family. In particular the spatially averaged commutator of (1)
is at least \(1/20\). The finite scalar probe and its failure records
retain their original normalization. This proof tightens its second
moment using the literal source Clifford and hopping matrices.

## 1. Exact neighbour columns and the full Fock norm

Write the dimensionless source as \(\widehat H=aH\). For one native
block the forward hopping matrix is
\[
 S_j={-\beta_D+\beta_D\gamma_j\over2},\quad
 S_j^\dagger S_j={I-\gamma_j\over2},\quad
 S_jS_j^\dagger={I+\gamma_j\over2}.                         \tag{2}
\]
At fixed parent link coordinates, let \(h\) denote the Hermitian hopping
one-body operator. For any normalized beta-positive column \(p\) and
beta-negative column \(m\) at site \(x\), set \(u=hp\), \(v=hm\).
These neighbour vectors are orthogonal to all columns at \(x\), and
\[
             \|u\|^2=\|v\|^2={\deg(x)\over2}.             \tag{3}
\]
Link-coordinate unitarity cancels in each squared column norm, without
replacing the quantum links in the source. The cross inner product
\(\langle u,v\rangle\) need not vanish at a boundary.

For \(b=|p\rangle\langle m|\), the residual after the exact onsite
commutator is
\[
 [\widehat H,d\Gamma(b)]-6d\Gamma(b)
      =d\Gamma(|u\rangle\langle m|-|p\rangle\langle v|).
 \tag{4}
\]
The all-occupation four-mode calculation in the
[local response proof](../particle_primary_response_moments/PROOF.md)
gives squared norm
\(h_0+\sqrt{h_0^2-|\langle u,v\rangle|^2}\le\deg(x)\), where
\(h_0=\deg(x)/2\). It includes the two-particle sector; a one-body
norm alone would give the wrong constant. Summing the \(2d\) columns
in (1) therefore gives the boundary-safe full-Fock estimate
\[
 \|[H,C_x]-6C_x/a\|^2
               \le {4d\,\deg(x)\over(d+1)a^2}.             \tag{5}
\]

At a strictly interior site both neighbour directions occur in every
axis. Equation (2) then gives
\[
          h^\dagger h\big|_{\text{site }x}=3I_{4d}.        \tag{6}
\]
Thus all \(u_{\alpha i}/\sqrt3\) and \(v_{\alpha i}/\sqrt3\)
are mutually orthonormal. They are also orthogonal to the local
positive and negative columns. Each term in (4) is a sum of two
raising operators on disjoint fermionic mode pairs. After a harmless
phase change, the complete residual is
\[
 {[H,C_x]-6C_x/a}
    ={1\over a}\sqrt{{3\over d(d+1)}}\,J_+ ,
 \tag{7}
\]
where \(J_+\) raises \(4d\) independent CAR pairs. Each pair's four
states consist of a spin doublet and two singlets. In the complete
Fock decomposition the largest spin is \(2d\), and
\(\|J_+\|^2=2d(2d+1)\). Consequently
\[
       \|[H,C_x]-6C_x/a\|^2
          \le {6(2d+1)\over(d+1)a^2}.                     \tag{8}
\]
For \(d=6\), this is \(78/(7a^2)\). The spin algebra here computes
a finite operator norm; its auxiliary spin labels are not physical
particle spin.

The parent coordinate action is a direct integral, so these uniform
pointwise norm bounds apply to the quantum-coordinate carrier. The
source link projection commutes with \(C_x\). Its compressed commutator
is therefore the compression of (4), and the same upper bounds hold
at every selected word cutoff. They hold before, and hence after,
the Gauss restriction. No multiplication of independently compressed
link factors is substituted for the actual source.

## 2. Original selected-state moments and narrower windows

There are \(B=3L^2(L-1)\) bonds, \(V=L^3\) sites, and
\(I=(L-2)^3\) strictly interior sites. Summing (5) and replacing its
interior value by (8) gives the actual uniformly averaged original
response moment
\[
 \int(\omega-6/a)^2\,d\mu_d(\omega)
       \le {M_d(L)\over a^2},\qquad
 M_d(L)={8dB-6(2d-1)I\over(d+1)V}.                       \tag{9}
\]
This is a measure for the same selected stationary state, with no
conditioning on success. For dimension six,
\[
 M_6(L)={48B-66I\over7V}
 ={78+252/L-792/L^2+528/L^3\over7}.                       \tag{10}
\]
The primary sides are \(2,8,32,\ldots\). At \(L=2\) this equals
\(72/7\). For \(L\ge8\), the polynomial numerator is increasing in
\(x=1/L\le1/8\), since its derivative is at least \(54\).
Hence the whole primary family obeys
\[
                 M_6(L)\le {3141\over224}<15.             \tag{11}
\]

Reuse the exact [polynomial passivity window](../particle_primary_pair_channels/PROOF.md).
The selected populations satisfy
\(0\le-(\log p)'(E)\le\beta\), with their actual Taylor order and
dyadic power. A success measure with commutator lower bound \(\chi\)
and moment (9) obeys
\[
 \mu_d([\delta,R])
   \ge \chi-\beta\delta-{M_d(L)\over(aR-6)^2}.             \tag{12}
\]
For the dimension-six probe take \(\delta=1/(100\beta)\),
\(R=46/a\), and use \(\chi\ge1/20\). Equations (11)-(12) give
\[
 \boxed{\mu_6([1/(100\beta),46/a])\ge3/100.}              \tag{13}
\]
The receipt retains the stronger value from the actual member's defect
bound and moment. The full-source finite transfer-to-excitation coupling
then gives original final-excitation weight at least \(1/50\) in
\([1/(100\beta),46/a+100\epsilon]\). Complete two-energy clocks with
resolution \(1/(400\beta)\), one-clock tail and implementation
allowances each \(1/400\), retain reported transfer weight at least
\(1/50\) in \([1/(200\beta),46/a+1/(200\beta)]\). All failure and
initial-energy records and the idle reference are retained. The initial
energy instrument preserves the source marginal, not its purification.

For the dimension-one \(Q_6=6\) scalar,
\[
              M_1(L)=(4B-3I)/V.
 \tag{14}
\]
At level one \(L=8\), this is \(591/64\). The actual native-dependent
occupation estimate supplies \(\chi>0.0902814\). Equation (12) at
\(R=24/a\) consequently proves original scalar response weight greater
than \(1/20\). Its upper endpoint at this member is \(48\) source energy
units. No rest-momentum or fermion-species label is assigned.

## 3. The same retained-origin observation domain

Use the [retained-origin protocol](../particle_primary_origin_momentum/PROOF.md):
one classical \(X\) is uniform on the interior cube \(F_n\), and this
same \(X\) controls every intervention and source-time step in a history.
The source density and canonical reference are unchanged. The complete
local instrument is
\[
          C_X,\quad \sqrt{I-C_X^\dagger C_X},
 \tag{15}
\]
with \(X\) retained. Its failure map differs from a single aggregated
failure on the whole box, although the success measure has the indicated
average. No fresh origin or averaged observable replaces (15).

Nonnegativity of the native occupation defects gives the direct
commutator estimate
\[
 \chi_F={1\over|F|}\sum_{x\in F}
       \langle[C_x^\dagger,C_x]\rangle
 \ge {12-\langle D_{\rm total}\rangle/|F|\over42}.          \tag{16}
\]
This uses an occupation sum on the selected state rather than subtracting
the whole omitted volume from a probability floor. When the origin
collar is positive, all these sites are strictly interior, so (8) supplies
the moment coefficient \(78/7\). With the authenticated positive lower
bound \(c_F\), define
\[
 \delta_F={c_F\over4\beta},\qquad
 R_F={6+\lceil\sqrt{4(78/7)/c_F}\rceil\over a}.
 \tag{17}
\]
The original retained-origin transfer weight in this interval is at least
\(c_F/2\). If the bound in (16) is nonpositive, no positive weight is
inferred by this calculation. At level zero the origin set is the entire
box and the averaged coefficient (10) is used instead.

At level seven the retained fraction is \((15/16)^3\). Equations
(16)-(17) give a positive original interval on the very same relative
observation domain that has the canonical limiting spatial generators.
The receipt calculates its endpoints and original weight without
allocating the large source or origin register. The unfiltered local
probe is not silently identified with a regular limiting field. Fixed
spatial and time smears, their nonzero limiting response, and the
convergence of finite wavevector readouts to those generators are separate
requirements.

## 4. An actual physical obstruction to an operator-norm shortcut

Let \(\Omega\) fill all negative-beta rows of all five native blocks and
have constant links. This is a physical Gauss singlet, but it is not the
selected polynomial state. For one native internal column the Slater
quadratic identity and the literal Wilson matrices give
\[
 \|S_x^+\Omega\|^2=2,\quad
 \|([\widehat H,S_x^+]-6S_x^+)\Omega\|^2=\deg(x),\quad
 \langle S_x^+\Omega,([\widehat H,S_x^+]-6S_x^+)\Omega\rangle=0.
 \tag{18}
\]
The second relation sums distinct neighbouring particle-hole occupations.
The independent Haar link integrals retain their squared column norms.
The code checks the two-endpoint Slater traces for every boundary pattern.
Multiplying by \(d\) internal columns and normalizing as in (1) gives
\[
 \|([H,C_x]-6C_x/a)\Omega\|^2
          ={\deg(x)\over(d+1)a^2},\qquad
 \|[H,C_x]\Omega\|^2={72+\deg(x)\over(d+1)a^2}.            \tag{19}
\]
Every vector in (19) has at most degree one on each link and is retained
at every selected cutoff. Thus (19) is a lower bound on the actual
physical operator norm, without a parent-only or omitted-sector
interpretation. A bounded-energy operator-norm residual for these
unfiltered local probes is impossible as \(a\to0\). Improving the
constant in (9) cannot remove this inverse-spacing scale.

Equation (19) is not a lower bound on the selected-state weighted
residual. A source-derived dressing or spatial probe could have a
different weighted response. Such a candidate needs nonzero original
overlap in the selected preparation and a bounded weighted energy
moment, or a direct positive fixed-window record bound. The auxiliary
trial cannot supply that overlap.

## 5. Scope of the energy and particle conclusions

The shrinking lower endpoint and inverse-spacing upper endpoint in (13)
do not prove tightness on a fixed physical energy interval. The executed
two-level countercontrol uses a positive polynomial preparation, vanishing
mean ground excess, positive commutator and a populated response exactly
at \(6/a\). It satisfies the moment/passivity bounds while every fixed
finite-energy band eventually has zero response. This is a counterexample
to that inference, without replacing the primary source.

The primary common-algebra construction and the retained-origin
translation construction provide their scoped state and generator
interfaces. The finite ground-tail theorem controls downward energy;
it supplies neither fixed positive-energy overlap nor a particle pole.
The source clock may be assigned its declared SI/GeV conversion, but a
mass requires an identified excitation and its momentum/dispersion
attachment. No measured mass, topological-defect identification,
nontrivial interacting continuum field, or no-go theorem for such a
field is asserted here.

## 6. Selected hopping energy and a complete density sum rule

The [physical tiled trial](../particle_primary_refinement_response/PROOF.md)
supplies a full-source upper bound, without substituting its state for
the selected density. If \(u_0\) is its authenticated eight-cell
dimensionless energy upper bound and \(g=(u_0+720)/8\), then
\[
 a\langle H\rangle_\rho\le -90V+gV+a\epsilon .
 \tag{20}
\]
Write the actual source as
\[
 aH=-90V+3D_{\rm total}+\widehat H_{\rm electric}
                  +aH_{\rm hop}+aH_{\rm magnetic}.
\]
The first two operators after the constant are positive and
\(aH_{\rm magnetic}\ge-P/7\), where \(P=3L(L-1)^2\).
Consequently the selected source obeys
\[
 {-a\langle H_{\rm hop}\rangle_\rho\over V}
 \ge h_n:=-g-{P\over7V}-{a\epsilon\over V}
 \ge h_*:=-g-{3\over7}-{1\over128}>0.03459 .
 \tag{21}
\]
The uniform last bound uses \(a\epsilon/V\le1/128\).
This is a nonzero kinetic expectation in the actual selected state;
it is not an excitation residue in the trial.

Let \(N_x\) count all sixty native CAR occupations at site \(x\), and
put \(n_x=N_x-30I\). Then \(n_x\) is gauge invariant,
\(\|n_x\|\le30\), and the \(n_x\) commute. For real \(f_x\), let
\(A_f=\sum_x f_x n_x\). Every onsite, electric and magnetic term
commutes with \(A_f\). The CAR bilinear identity gives, bond by bond,
\[
 [A_f,[H,A_f]]
     =-\sum_{\langle x,y\rangle}(f_x-f_y)^2 T_{xy},
 \tag{22}
\]
where \(T_{xy}\) is the full Hermitian hopping bond including its
\(1/a\) factor. Both the link-word projection and Gauss projection
commute with \(A_f\); equation (22) therefore holds on the actual
compressed physical source. The finite double-commutator method used in
[the comparison electromagnetism calculation](../particle_comparison/electromagnetism.py)
is applied here to this source, without transferring its comparison
state or spectrum.

The [physical cube action](../particle_primary_channel_labels/PROOF.md)
first restricts every one-link ghost to its constant function.
Proper cube rotations then intertwine the interior graph, complete
plaquette products, source Clifford matrices and link cutoff.
They commute with the physical Hamiltonian and hence its selected
polynomial density. If \(H_{{\rm hop},j}\) is the sum of bonds in
direction \(j\), this proves
\[
 \langle H_{{\rm hop},1}\rangle_\rho
 =\langle H_{{\rm hop},2}\rangle_\rho
 =\langle H_{{\rm hop},3}\rangle_\rho
 ={1\over3}\langle H_{\rm hop}\rangle_\rho .
 \tag{23}
\]
There is no permutation of the asymmetric ambient dangling-link set
and no assumption of continuous spatial rotations or translations.

Fix an axis \(j\) and the coordinate phase increment \(\theta=\pi a\).
Define
\[
 A_c=\sum_x\cos(\theta x_j)n_x,\quad
 A_s=\sum_x\sin(\theta x_j)n_x,\quad \Lambda=30V .
 \tag{24}
\]
Their simultaneous spectral values satisfy
\(A_c^2+A_s^2\le\Lambda^2I\): the modulus of
\(\sum_x e^{i\theta x_j}n_x\) is at most
\(\sum_x|n_x|\le\Lambda\). Thus an actual complete instrument is
\[
 K_c=A_c/\Lambda,\quad K_s=A_s/\Lambda,\quad
 K_f=\sqrt{I-K_c^\dagger K_c-K_s^\dagger K_s}.
 \tag{25}
\]
All three labels and the idle reference are retained. This is one
three-outcome instrument. It is not an equally weighted choice between
two independently completed instruments. The phase coefficients are
finite algebraic numbers, since \(a=2^{-n}\).
They label a source-coordinate density apparatus, not a conserved
momentum outcome. The instrument is an intervention and need not
preserve its input density or purification.

In a physical energy eigenbasis define the combined original success
measure
\[
 \mu_n(d\omega)=\sum_{i,k}p_i
       \left(|(K_c)_{ki}|^2+|(K_s)_{ki}|^2\right)
       \delta_{E_k-E_i}(d\omega).
 \tag{26}
\]
Its total mass is at most one; no successful weight is divided out.
Stationarity and Hermiticity of \(K_c,K_s\) give the exact signed
first moment
\[
 \begin{split}
 M_{1,n}&=\int\omega\,d\mu_n(\omega)\\
 &=\frac12\left\langle[K_c,[H,K_c]]+[K_s,[H,K_s]]\right\rangle_\rho\\
 &=-{1-\cos\theta\over\Lambda^2}
            \langle H_{{\rm hop},j}\rangle_\rho
 =-{1-\cos\theta\over3\Lambda^2}\langle H_{\rm hop}\rangle_\rho>0 .
 \end{split} \tag{27}
\]
The cosine and sine squared bond differences add to
\(2(1-\cos\theta)\) on direction \(j\), and vanish on the other
directions. For a general complex operator the double commutator
instead gives the sum of its first moment and its adjoint's first
moment. Equation (27) avoids an unsupported factor of two.

The actual full-Fock bond norm is at most \(30/a\), including all
fifteen native internal rows. Thus
\[
 { (1-\cos\theta)Vh_n\over3a\Lambda^2}
 \le M_{1,n}
 \le{(1-\cos\theta)10B\over a\Lambda^2}
 \le {a\over6V}.
 \tag{28}
\]
The rational receipt uses the exact factors \(2\) at level zero and
\(1\) at level one. At every other level it uses
\(2a^2\le1-\cos(\pi a)\le\min(2,5a^2)\).
The lower inequality follows from the sine chord on
\([0,\pi/2]\); the upper uses \(1-\cos t\le t^2/2\) and
\(\pi^2<10\).

The positive finite moment does not establish survival in a fixed
positive-energy band. In fact the normalization in (25) supplies
an exact obstruction for this particular observation family.
Completeness and mean ground excess give
\[
 \int_{\omega<0}|\omega|\,d\mu_n
 \le\sum_i p_i(E_i-E_0)
       \langle i|K_c^\dagger K_c+K_s^\dagger K_s|i\rangle
 \le\epsilon .
 \tag{29}
\]
This is the first-moment version of the reused finite
[ground-record bound](../../Lean/QFT/GroundRecordTail.lean).
Combining (28)-(29), for every fixed \(\delta>0\),
\[
 \boxed{\ \mu_n((\delta,\infty))
       \le {a/(6V)+\epsilon\over\delta}\longrightarrow0.\ }
 \tag{30}
\]
The same estimate with the member-specific upper bound in (28) is
recorded. The threshold is fixed in source-clock energy units, or in
physical units under one fixed clock conversion. An independently
varying conversion changes the represented interval. Equation (30)
is an original-probability statement about the
volume-normalized density instrument. It does not concern the
collective pair instrument in (13), forbid other field normalizations
or dressings, or exclude a nontrivial field limit of the source.
Multiplying records by \(\Lambda^2\) changes their normalization and
cannot turn (27) into an unconditioned positive-band weight.
The exact sum rule fixes a weighted energy response; it does not locate
a spectral pole, supply rest momentum, or attach a mass in GeV.
