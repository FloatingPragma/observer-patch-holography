# Source three-point response and Yukawa extraction

## The fixed source and the quantities calculated

The reference is the faithful distribution returned by
[`local_preparation.reference`](../native_modular_reference/local_preparation.py)
on all 66 binary configurations of total occupation two on twelve ports.
Write the occupations as a vector \(n\), the protected core as
`local_preparation.CORE`, and its occupation as \(N\). A configuration has
weight

\[
 w_0=\frac2{225},\qquad w_1=\frac1{36},\qquad w_2=\frac2{45}
 \tag{1}
\]

according to \(N=0,1,2\). These weights are checked against the native
preparation producer. They are not fitted to a fermion mass. Their
multiplicities are \(45,20,1\), so the complete core record probabilities
are \(2/5,5/9,2/45\). No outcome is discarded.

The actual source coordinate frame \(B\) is the twelve-by-three matrix from
[`slow_geometry.raw_frame`](../common_source/slow_geometry.py). Its exact
Gram matrix is

\[
 B^TB=(10+2\sqrt5)I_3.
 \tag{2}
\]

The recorded vector used here is \(q=B^Tn\). This is a finite occupation
observable in the source rank-three band. Physical fermion family labels
are not assumed.

The operator calculation uses the actual port currents and the declared
spin/exterior matter architecture in
[`native_matter_attachment`](../common_source/native_matter_attachment.py).
It recomputes that representation and solves its invariant tensor equations.
The architecture implements a chosen lift on an auxiliary encoded register;
unchanged native repair does not select that lift. This premise is retained
when calling its exact algebraic functions.

The finite charge-conjugation and spin-tensor identities reuse
[`WeylYukawaConventions`](../../Lean/Screen/WeylYukawaConventions.lean).
Its `right_charge_forced` theorem assumes both the all-left and Dirac
neutrality equations and concludes the singlet charge reversal. Its
`epsilon_SL2_invariant` theorem assumes determinant one over a commutative
coefficient ring and proves the two-component tensor identity. These
theorems constrain conventions and tensor transformation laws. Their
hypotheses supply neither a physical field assignment nor a coupling
coefficient. The native representation producer verifies the corresponding
matrix identities on the concrete source currents.

## Exact finite three-point amplitudes

For either \(H=N\) or \(H=N^2\), define

\[
 \mu=E[q],\quad \mu_H=E[H],\quad
 C=E[(q-\mu)(q-\mu)^T],\quad
 T_H=E[(q-\mu)(q-\mu)^T(H-\mu_H)].
 \tag{3}
\]

All expectations are complete 66-term traces with the exact weights (1).
The calculation gives

\[
 \mu=\left(0,\frac{14(1+\sqrt5)}{75},0\right)^T,
\]

\[
 C=\operatorname{diag}\left(
 \frac{353}{225}+\frac{19\sqrt5}{75},\quad
 \frac{936}{625}+\frac{1858\sqrt5}{5625},\quad
 \frac{19}{15}+\frac{19\sqrt5}{75}\right).
 \tag{4}
\]

Each diagonal entry is positive. For the linear core insertion,

\[
 E[N]=\frac{29}{45},\qquad \operatorname{Var}(N)=\frac{644}{2025},
\]

\[
 T_N=\operatorname{diag}\left(
 -\frac{112}{10125}-\frac{176\sqrt5}{3375},\quad
 -\frac{8864}{84375}-\frac{64\sqrt5}{253125},\quad
 -\frac{176}{675}-\frac{176\sqrt5}{3375}\right).
 \tag{5}
\]

For the squared core insertion,

\[
 E[N^2]=\frac{11}{15},\qquad \operatorname{Var}(N^2)=\frac{164}{225},
\]

\[
 T_{N^2}=\operatorname{diag}\left(
 -\frac{508}{3375}-\frac{28\sqrt5}{375},\quad
 \frac{3224}{28125}+\frac{7424\sqrt5}{84375},\quad
 -\frac{28}{75}-\frac{28\sqrt5}{375}\right).
 \tag{6}
\]

Since \(C\) is diagonal, its normalization is explicit:

\[
 A_H=\frac{C^{-1/2}T_HC^{-1/2}}{\sqrt{\operatorname{Var}(H)}}.
 \tag{7}
\]

The exact entries are retained in the receipt. Their numerical diagnostics
are

| Source insertion | First axis | Second axis | Third axis |
|---|---:|---:|---:|
| \(N\) | -0.10601870 | -0.08375413 | -0.36501956 |
| \(N^2\) | -0.17414553 | +0.16309769 | -0.34522621 |

These are dimensionless static connected response amplitudes. They are not
physical Yukawa values. Equations (3) through (7) use no measured mass,
Yukawa coefficient, Higgs expectation value or empirical clock.

### Conditional expectation and response checks

Conditional on \(N=c\), the configurations are uniform over a \(c\)-subset
of the two core ports and a \((2-c)\)-subset of the ten other ports. For
vectors \(b_1,\ldots,b_m\) and a uniform \(k\)-subset \(S\),

\[
 E\left[\sum_{i\in S}b_i\right]=\frac{k}{m}\sum_i b_i,
\]

\[
 E\left[\sum_{i,j\in S}b_ib_j^T\right]
 =\frac{k}{m}\sum_i b_ib_i^T
 +\frac{k(k-1)}{m(m-1)}\sum_{i\ne j}b_ib_j^T.
 \tag{8}
\]

Applying (8) separately to the core and its complement, and including their
cross term, reproduces (4) through (6). The test uses this sampling formula
independently of the 66-state moment loop. The receipt also reconstructs

\[
 T_H=\sum_c P(N=c)(c^r-E[N^r])
 E[(q-\mu)(q-\mu)^T\mid N=c],\qquad H=N^r.
 \tag{9}
\]

For the declared exponential intervention

\[
 \rho_h(x)=\frac{\rho(x)e^{hH(x)}}{\sum_z\rho(z)e^{hH(z)}},
 \tag{10}
\]

differentiation gives \(\partial_h E_h[f]|_0=\operatorname{Cov}(f,H)\).
Applying this identity to both terms of

\[
 \operatorname{Cov}_h(q_a,q_b)=E_h[q_aq_b]-E_h[q_a]E_h[q_b]
\]

proves

\[
 \left.\partial_h\operatorname{Cov}_h(q_a,q_b)\right|_0=(T_H)_{ab}.
 \tag{11}
\]

An independent numerical symmetric derivative of the complete finite
partition sum checks this equality. The derivative is an auxiliary
verification; equations (3) through (9) are exact.

### What normalization does and does not determine

Rescaling \(H\) by a positive constant rescales \(T_H\) and
\(\sqrt{\operatorname{Var}(H)}\) by the same factor. Positive rescaling of a
family observable similarly cancels in (7). A missing amplitude cannot be
inferred merely from an unnormalized statistical tensor.

The insertions \(N\) and \(N^2\) give a stronger witness: both are
well-defined records on the same source, with the same reference and
zero-field modular flow. Their normalized second-axis responses have
opposite signs. Thus the zero-field reference and modular flow alone do
not specify a unique normalized response to an unspecified intervention.
The two complete parameter-dependent laws (10) differ. No alternative
Yukawa is claimed at one fixed complete Higgs-dependent source law.

## Mixed fermion and composite-scalar three-point function

The M15 scalar-bilinear construction permits a more direct test than the
occupation statistic. Its physical domain is an independently declared
thirty-mode CAR module, Spin2 tensor M15, with its blank Fock state. This
domain is explicitly separate from both the 66-state occupation reference
and the auxiliary five-mode exterior grading used to construct M15.

For either component of a native scalar doublet, let \(C\) be its symmetric
internal coefficient and set \(\Delta=\epsilon\otimes C\), where
\(\epsilon_{12}=1=-\epsilon_{21}\). Let \(r\) be the directed internal
rank, equal to three for either color channel and one for the lepton
channel. Define the normalized composite

\[
 B^\dagger=\frac1{2\sqrt{2r}}\sum_{i,j=1}^{30}
 \Delta_{ij}c_i^\dagger c_j^\dagger.
 \tag{12a}
\]

Its vacuum norm is one because
\(\sum_{i<j}|\Delta_{ij}|^2=\operatorname{tr}(C^\dagger C)=2r\).
The exact mixed three-point function is

\[
 \langle0|c_jc_iB^\dagger|0\rangle
 =\frac{\Delta_{ij}}{\sqrt{2r}},\qquad
 \langle0|BB^\dagger|0\rangle=1,
 \qquad\langle0|c_ic_j^\dagger|0\rangle=\delta_{ij}.
 \tag{12b}
\]

The code checks all 435 pair contractions for all six native doublet
components using occupation-mask CAR signs. Every nonzero modulus is
\(1/\sqrt6\) in a color channel or \(1/\sqrt2\) in the lepton channel.
These normalized overlaps exist for a free Hamiltonian
\(H_0=\omega\sum_i c_i^\dagger c_i\), which contains no independent
Higgs field and no Higgs-dependent quadratic coefficient. Extending its
action by a spectator Higgs variable gives zero Higgs derivative.
Coupling an external source to \(B\) is a separate intervention. The
nonzero overlap (12b) expresses the definition of a composite interpolator;
it cannot determine the strength of a dynamical Higgs interaction.

The scalar-bilinear producer identifies one hypercharge-minus-half and two
hypercharge-plus-half doublet channels. It supplies no unique linear
combination, scalar binding kernel, pole residue or condensate. The
unit-normalized M15 tensor and its current intertwiner therefore do not
select the Higgs intervention omitted from the action. This mixed CAR
calculation is exact in its declared vacuum domain; no execution of a
native preparation or a transfer from the 66-state reference is asserted.

## Extracting a Yukawa matrix from a source Higgs variation

The actual native current equations produce three allowed charge origins

\[
 (-4,1),\qquad (1,2),\qquad (-3,6)
 \tag{12}
\]

in the \(6Y\) convention. Their invariant Higgs tensor spaces each have
complex dimension one. These are conditional internal representation
identities. Naming them up, down and charged lepton uses the existing
current intertwiner; it does not determine a numerical coefficient.

Let \(C_{s,k}\) be the symmetric fifteen-dimensional coefficient for
sector \(s\) and Higgs component \(k\). Both Higgs components determine
the complete weak-doublet and singlet supports. Let \(L_s,R_s\) be their
projectors and

\[
 D_{s,k}=L_s C_{s,k}R_s.
 \tag{13}
\]

Exact native arithmetic verifies

\[
 D_{s,k}^\dagger D_{s,k}=R_s,
 \qquad C_{s,k}=D_{s,k}+D_{s,k}^T,
 \qquad r_s=\operatorname{tr}R_s\in\{3,3,1\}.
 \tag{14}
\]

Consequently every directed nonzero singular value is one. Color supplies
multiplicity, with no extra sector-dependent tensor modulus. Choosing a
tensor to have unit singular values is a convention for its coefficient;
it supplies no physical coupling.

On a gauge-trivial three-dimensional family factor, the general allowed
real-component directional derivative of the all-left quadratic pairing
kernel, along H=v e_k with real v, is

\[
 V_k=\left.\frac{dK(v e_k)}{dv}\right|_{v=0}
 =\sum_s\left(D_{s,k}\otimes Y_s+
                   D_{s,k}^T\otimes Y_s^T\right).
 \tag{15}
\]

Here the directions are real. The stored conjugate-Higgs flag determines
whether the corresponding imaginary-component derivative carries +i or -i;
no Wirtinger derivative with respect to H alone is asserted.

Every \(Y_s\) is an arbitrary complex three-by-three matrix. Gauge
invariance supplies three tensor lines and hence 27 complex coefficient
coordinates. The transpose in (15) is essential; a nonsymmetric complex
family matrix is allowed. The API reports the all-left coefficient. The
separate native Dirac convention conjugates \(Y\) on entering this kernel.

Equation (14) gives an exact dual:

\[
 (Y_s)_{ij}=\frac1{r_s}\sum_{a,b}
 \overline{(D_{s,1})_{ab}}(V_1)_{3a+i,3b+j}.
 \tag{16}
\]

The oppositely directed term vanishes in this contraction because its
supports are reversed. Other sectors vanish because their oriented
supports are disjoint. The surviving term is
\(\operatorname{tr}(D_{s,1}^\dagger D_{s,1})Y_s=r_sY_s\), proving (16).
The executable substitutes 27 independent symbolic variables, recovers
each exactly, and reconstructs both \(V_0,V_1\). Any extra operator or
inconsistent coefficient between the two Higgs components is rejected.

This is a constructive completion of the coefficient readout: a supplied
same-source Higgs derivative removes all 27 coefficient freedoms. The
scalar tensors alone do not provide that derivative. The controls
\(Y_s=I_3\) and \(Y_s=2I_3\) satisfy the same current, canonical metric and
unit-Higgs premises, with different derivatives and spectra. They are
different interacting actions, not different outputs of one frozen action.

If the same unbroken irreducible native \(A_5\) triplet acts on both
family factors and the Higgs is family invariant, the exact existing
[`family_naturality`](../observer_particles/family_naturality.py) theorem
reduces \(Y_s\) to \(y_sI_3\) in common endomorphism coordinates. Its
bilinear is \(y_sG\) in a nonorthonormal family frame of Gram matrix \(G\).
Those extra hypotheses give three equal family singular values and retain
the three scalar coefficients \(y_s\). The protected-core reference used
in (1) is anisotropic; the unbroken-family theorem cannot be applied to it
without an additional invariance hypothesis.

### Coefficients of the existing declared Higgs action

The source-only trial in
[`chiral_candidate.action`](../common_source/chiral_candidate/action.py)
does supply a complete linear Higgs kernel. Its declared coefficient map is

\[
 Y_u=aF,\qquad Y_d=bF,\qquad Y_e=JF,
\]

\[
 a=\log\frac{w_1}{w_0}=\log\frac{25}{8},\quad
 b=\log\frac{w_2}{w_1}=\log\frac85,\quad
 J=a-b=\log\frac{125}{64}.
 \tag{18}
\]

Here \(F\) is the native one-body occupation compression onto the same
rank-three source band. Its diagonal entries in the native coordinate
order are

\[
 F=\operatorname{diag}\left(
 \frac{41}{225}-\frac{7\sqrt5}{750},
 \frac{41}{225}+\frac{7\sqrt5}{750},
 \frac{61}{450}\right).
 \tag{19}
\]

The producer chooses (18), a canonical fermion kinetic metric and a unit
Higgs coordinate. It supplies no dynamical Higgs kinetic residue. The reference and geometry determine the constants
appearing in that chosen map. No measured fermion mass is supplied.

The executable differentiates the actual linear dependence
\(K(v)=vK(1)\), carries both Higgs derivatives through the verified native
current intertwiner, and applies (16). It recovers all three matrices
in (18) exactly, including the sign convention of each native tensor.
Their ascending singular values in that declared Higgs coordinate are

| Declared action channel | First | Second | Third |
|---|---:|---:|---:|
| Up | 0.1544566473 | 0.1838502904 | 0.2314102039 |
| Down | 0.0637116031 | 0.0758361452 | 0.0954540664 |
| Charged lepton | 0.0907450442 | 0.1080141452 | 0.1359561376 |

These are dimensionless coefficients of the specified trial action. A fully
canonical Higgs field requires division by the square root of its kinetic
residue, which this trial does not supply. All
three channels share the eigenvalue ratios of \(F\), with maximum-to-minimum
ratio strictly below \(3/2\). They do not provide an observed charged-fermion
hierarchy or an independently established physical Yukawa prediction.
The declared map (18) is not forced by unchanged native repair. This
calculation establishes the exact conditional consequence of that map and
its actual Higgs derivative, rather than treating it as unspecified.

## Canonical kinetic normalization and physical attachment

Suppose the same source field construction supplies positive kinetic
metrics \(Z_L,Z_R\) and a positive Higgs kinetic coefficient \(Z_H\).
For an all-left bilinear, the field substitutions

\[
 \psi_L=Z_L^{-1/2}\psi_{L,c},\quad
 \psi_R=Z_R^{-1/2}\psi_{R,c},\quad
 H=H_c/\sqrt{Z_H}
\]

give

\[
 Y_c=Z_L^{-T/2}YZ_R^{-1/2}/\sqrt{Z_H}.
 \tag{20}
\]

The canonical tree masses at neutral component norm \(w\) are
\(w\sigma_i(Y_c)\). In the Standard Model convention \(w=v/\sqrt2\).
For real family Gram \(G\), inserting \(Y=yG\) and \(Z_L=Z_R=G\)
gives \(Y_c=yI_3/\sqrt{Z_H}\). Kinetic normalization cannot turn this
invariant bilinear into a split family spectrum. Tests also apply arbitrary
complex invertible field-coordinate changes and verify invariant canonical
singular values; this catches a transpose/adjoint substitution error.

The statistics in (7) use equal-time covariances. Equation (20) uses the
kinetic residues of propagating fermion and Higgs fields. Identifying the
two requires a same-source transfer theorem. The available occupation
observables are commuting records, and \(N,N^2\) are scalar source
insertions without a proved Higgs-doublet transformation law. Neither is
therefore substituted into (16).

For an interacting physical claim, the sufficient datum is the Higgs
variation of the inverse fermion propagator at stated kinematics, together
with the fermion and Higgs kinetic residues and their shared clock. A
complete connected three-point record can supply that datum after
amputation by the corresponding two-point functions and after establishing
the field dictionary. A static third cumulant alone is insufficient.
The finite inverse in
[`native_record_yukawa_field`](../common_source/native_record_yukawa_field.py)
recovers a supplied Hermitian \(Y\) from a declared static service
\(e^{-ihX\otimes Y}\); the service explicitly places \(Y\) in its
generator. It does not identify the occupation responses (5) or (6) with
that generator.

The distinct
[`record_two_point`](../common_source/chiral_candidate/record_two_point.py)
apparatus reconstructs a canonical two-mode transfer and its neutral-Higgs
multiplier derivative from complete detector probabilities for the supplied
quadratic action. That finite result assumes its apparatus, clock and alias
window, and supplies no interacting wavefunction renormalization or Higgs
propagator. The missing scalar residue and physical field attachment cannot
be supplied by reinterpreting its externally programmed Higgs multiplier.

The receipt therefore contains exact native tensor shape, exact normalized
finite-state responses and a complete coefficient extraction rule. It
contains no source-selected physical charged-fermion Yukawa matrix or
physical mass prediction. Target-mutation controls alter comparison
payloads while replaying this package's source calculations and verify
unchanged output and absence of reads from those target files. This tests
runtime input separation, not design history. Finite preparation errors,
physical pole matching and continuum refinement are outside this result.
