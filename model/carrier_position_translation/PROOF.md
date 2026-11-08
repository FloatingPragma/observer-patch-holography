# Native position covariance and its field boundary

## 1. The recorded coordinate

Use the exact twelve-port chart in `common_source/slow_geometry.py`.
Its unnormalized vectors are $u_a\in\mathbb Q(\sqrt5)^3$; the unit-vector
Gram is four times the rank-three slow projector. On the actual
78-dimensional nonnegative total-two occupation carrier define the
commuting self-adjoint matrices

$$X_i|n\rangle=x_i(n)|n\rangle,\qquad
x(n)=\sum_{a=0}^{11}n_a u_a.\tag{1}$$

The centroid is $X/2$; all equations below rescale accordingly. The chart
fixes dimensionless coordinates, without a metre or Planck-constant
calibration. Absolute $x(n)$ is a source-state observable. From a known
initial reference, its change is readable from the seam and current
records alone. On its declared binary-port source the complete local
E-record occupation decoder supplies a separate operational readout with
its coverage and erasure contracts. That decoder is not asserted to
measure every integer-load state in (1); this proof does not replace it
with access to a hidden load.

For a fixed seam $e=(a,b)$ and fair placement coin $c$, the actual native
rule emits current $j$ and changes
$n\mapsto n'=n+j(e_b-e_a)$. Let $D_r$ be the diagonal projector onto
the inputs emitting record $r=(c,j)$, and let $S_r$ be the partial
permutation mapping each such input to its actual output. Including $j$
makes the map injective, since $n=n'-j(e_b-e_a)$. Its Kraus operator is
$K_r=S_r/\sqrt2$. All six record slots are retained, with

$$K_r^*K_r=\tfrac12D_r,\qquad \sum_rK_r^*K_r=I.\tag{2}$$

## 2. Exact operator covariance on the branch support

Put $d_r=j(u_b-u_a)$. Applying (1) to each actual native column proves

$$X_iK_r-K_rX_i=d_{r,i}K_r.\tag{3}$$

Both sides vanish off $D_r$. Thus the result holds on every coherent
input, not only on occupation basis preparations. Functional calculus
for the finite diagonal operators gives the Weyl and position-projection
forms

$$e^{iq\cdot X}K_re^{-iq\cdot X}=e^{iq\cdot d_r}K_r,
\qquad E_X(B)K_r=K_rE_X(B-d_r).\tag{4}$$

For any joint source/reference density $\rho$, the unnormalized joint
probability of record $r$ and terminal position in $B$ is therefore

$$\operatorname{Tr}[E_X(B)K_r\rho K_r^*]
=\tfrac12\operatorname{Tr}[D_rE_X(B-d_r)\rho].\tag{5}$$

Identity operators on the reference are understood. Formula (5) keeps
the acceptance factor. Conditioning on $r$ filters the initial position
law by $D_r$ before shifting it; it does not translate the original law
without that filter.

Let $V=\sum_r|r\rangle\otimes K_r$ be the complete source isometry,
$X_{\mathrm{out}}=I_{\mathrm{record}}\otimes X$, and
$D_{\mathrm{record},i}=\sum_rd_{r,i}|r\rangle\langle r|\otimes I$.
Equation (3), before or after classical commitment of records, gives

$$(X_{\mathrm{out},i}-D_{\mathrm{record},i})V=VX_i.\tag{6}$$

For any finite adaptive history with the chosen seam retained, induction
gives (3)–(6) with $K_h$ and $d_h=\sum_{r\in h}d_r$.
Branch-dependent durations or stopping flags do not alter this algebra;
every stopped branch must remain in the complete instrument. Tensoring
these operator equalities with $I_R$ proves the reference extension. The
certificate additionally enumerates all 312 columns of the actual
two-event current-feedback observer, rather than substituting a different
scheduler.

## 3. Why this is not a unitary spatial translation

For each seam the six $D_r$ have ranks $1,11,66$, each occurring twice.
In particular, no nonzero-current branch has full support. A fixed
nonzero-current outcome can have probability zero or $1/2$ depending
on the source input. It has no input-independent inverse on the whole
carrier.

There is also a direct obstruction to replacing it by a finite unitary.
For any self-adjoint $D$-dimensional matrix $X_i$, unitary $U$, and real
$a_i$, set $R=U^*X_iU-X_i-a_iI$. Trace invariance and
$|\operatorname{Tr}R|\leq D\|R\|$ imply

$$\|U^*X_iU-X_i-a_iI\|\geq|a_i|.\tag{7}$$

Thus exact nonzero translation covariance cannot hold on the full
finite carrier. A restriction to a subspace can evade this trace
argument only through explicit boundary leakage or nonunitary
compression, whose contribution must be bounded.

The actual uniform-seam/fair-coin calibration gives a further operational
distinction. Its mean occupation map is $I-L/60$ and
$Lu_i=(5-\sqrt5)u_i$. Its complete channel consequently satisfies

$$\Phi^*(X_i)=\frac{55+\sqrt5}{60}X_i.\tag{8}$$

This is contraction toward the carrier's origin, not a fixed spatial
shift. Equation (8) concerns the declared uniform calibration; it is not
the mean law of the two-seam adaptive observer. Neither repair-attempt
time nor the slow-band eigenvalue is identified with the field's momentum.

## 4. The uncertainty statement supplied by these operators

For $R_r=(K_r+K_r^*)/2$ and $I_r=(K_r-K_r^*)/(2i)$, equation (3) gives

$$[X_i,R_r]=i d_{r,i}I_r,\qquad
[X_i,I_r]=-i d_{r,i}R_r.\tag{9}$$

The Robertson inequality on the same finite source Hilbert space gives

$$\Delta X_i\Delta R_r\geq\tfrac12|d_{r,i}\langle I_r\rangle|,
\qquad
\Delta X_i\Delta I_r\geq\tfrac12|d_{r,i}\langle R_r\rangle|.\tag{10}$$

These are state-dependent source identities. They supply neither an
implemented measurement of $R_r,I_r$ nor a particle momentum. A position
eigenstate has $\Delta X_i=0$ and is a legitimate input, excluding a
universal nonzero constant bound with any finite self-adjoint candidate
momentum on all source states. The complete committed record algebra
alone does not add the missing momentum measurement.

## 5. Total-load position is trivial in the actual protected spin sector

The spin construction uses a declared two-flavor field extension with
exactly one fermion per site in its protected subspace $P$. Write
$N_x=n_{x,\uparrow}+n_{x,\downarrow}$. Its actual full-domain repair
and inclusion obey

$$N_xP=P,\qquad M_{c,r}P=\delta_{r,0}P/\sqrt2.\tag{11}$$

The native occupations here belong to the field extension; the
twelve-port chart is not identified with the cubic site chart. Even if
a proposed spatial dictionary assigns arbitrary real weights $f_x$ to
these total occupations, its readout necessarily satisfies

$$X_fP=\Big(\sum_x f_x\Big)P,\qquad
X_f=\sum_xf_xN_x.\tag{12}$$

It is constant on **every** spin state and cannot locate a spin
excitation. Balanced coordinate weights give the zero operator. The
spin's discrete or limiting continuous translations do not change that
conclusion. An internal spin or defect detector requires a different
observable attachment.

This obstruction survives the specified source-to-spin limit. In
(13)–(14), $P$ denotes the canonical isometric inclusion of the spin
space into the source space. Let the same low-band isometry $W$ satisfy
$\|W-P\|\leq2\tau/U$, and let
$c_f=\sum_xf_x$. Since the self-adjoint $X_f-c_fI$ annihilates $P$,
both cross terms vanish exactly:

$$W^*X_fW-c_fI
=(W-P)^*(X_f-c_fI)(W-P).\tag{13}$$

As $N_x-1$ has spectrum in $\{-1,0,1\}$, this gives

$$\|W^*X_fW-c_fI\|
\leq4\Big(\sum_x|f_x|\Big)(\tau/U)^2
=\frac{4\sum_x|f_x|}{|E|^6}.\tag{14}$$

The last equality uses the actual source trajectory. For fixed local
weights it vanishes uniformly on all spin states. Rescaling the detector
by a divergent factor is a separate observable definition and requires
its own convergence, noise and operational implementation proof.

## 6. The required position–momentum join

An exact particle position relation requires an operational position
measure, its strongly continuous field translations, and covariance on
the same represented states. For a sharp unbounded position $\mathsf X$
and $U(a)=e^{-ia\cdot\mathsf P/\hbar}$, the desired covariance
$U(a)^*\mathsf X_iU(a)=\mathsf X_i+a_iI$ implies
$[\mathsf X_i,\mathsf P_j]=i\hbar\delta_{ij}I$ only on an appropriate
common differentiable domain. Finite approximants must supply the
position measure, the translated instruments, boundary and scaling
errors, and convergence in that same representation. Equation (3)
provides exact record displacement covariance. Equations (7), (12)
and (14) exclude identifying it directly with that particle relation
on the existing full finite or protected total-load carrier.
