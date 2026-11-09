# Selected odd creation and its energy scale

The source, polynomial density, idle conjugate reference and clock are the
primary construction contract. At level \(n\), write \(a=2^{-n}\),
\(L=2/a^2\), \(V=L^3\) and \(\epsilon=a/16\). All expectations below use
the actual source marginal \(\rho\), with
\(\operatorname{Tr}\rho(H-E_0)\le\epsilon\). An anchor or Ritz vector is
used only to bound an operator or the source ground energy. Neither is the
preparation of this experiment.

## 1. A physical odd probe

The native weak doublet has integer charge \(-3\); the native singlet has
charge \(6\). In the positive-eigenvalue two-dimensional spin space of
\(\beta_D\), define

\[
 A_{x,\alpha b}=
 c^\dagger_{x,L_1,\alpha}c^\dagger_{x,L_2,\alpha}
 c^\dagger_{x,e,b},\qquad \alpha,b=1,2.                 \tag{1}
\]

The weak exterior determinant times the singlet is invariant under the
actual \(S(U(3)\times U(2))\) representation. The twelve native Gauss
generators annihilate this monomial, as checked in the local kernel of
`particle_primary_sector_threshold`. Its norm is one. It has odd CAR parity
and native occupation increment \((0,2,0,0,1)\), in charge order
\((-4,-3,1,2,6)\). These are source labels, with no Standard Model species
assignment. The four spin choices are an explicitly selected subset of the
six-dimensional spin-one weak pair times spin-one-half singlet. A pure
total-spin assignment is not made.

The native integer-mask basis and the coordinate-link exterior basis are
related by the row permutation
`(0,1,4,2,5,7,10,3,6,8,11,9,12,13,14)`. The code reconstructs this
permutation from the frozen subset lists and checks the weak determinant
times the singlet on a nonreal nondiagonal source link. The proof uses the
same exterior representation throughout.

A complete implementation uses a gauge-singlet idle apparatus memory with
the corresponding triple loaded. Exchange the source's empty and full
three-mode occupations against this memory, leaving partial occupations
unchanged. Measuring whether the memory emptied has source Kraus operators
\(A\) and \(I-A^\dagger A\), up to a known CAR sign corrected in the
apparatus. Since \(A^\dagger A\) is a projection, these are complete.
The exchange is even and preserves every combined source-memory native
occupation. The memory supplies the opposite increment to the source.
All memory and pointer records, including failure, are retained. This
declared apparatus input does not replace \(\rho\). An independent uniform
choice of \((x,\alpha,b)\) is retained as a classical record, so there are
\(8V\) complete probe outcomes.

## 2. Original selected-state success

Let \(N_q^+\) count positive-beta occupations and \(D_q\) count these
occupations plus negative-beta holes. Full-source charge conjugation is
linear, unitary, preserves \(H\) and \(\rho\), and sends the positive-beta
number to the corresponding negative-beta holes. The spin identity is
\(\beta_D M\overline v=-M\overline v\) for every positive-beta vector
\(v\). Therefore

\[
                  \langle N_q^+\rangle=\tfrac12\langle D_q\rangle.
                                                               \tag{2}
\]

For the three commuting occupations in (1),
\(A^\dagger A=\prod_{j=1}^3(1-n_j)\ge1-\sum_jn_j\).
Averaging the four spin choices and all cells gives the original success
probability

\[
 p_n=\frac1{4V}\sum_{x,\alpha,b}\langle A^\dagger A\rangle
 \ge 1-\frac{\langle D_L+D_e\rangle}{4V}.                \tag{3}
\]

This estimate requires neither site independence nor equality of individual
site probabilities. At level zero the existing full-source total-defect
bound is \(\langle D\rangle<4.992\). The exact bound used by the receipt
gives \(p_0>0.84402\).

For the whole family, retain the onsite fractions
\((0,21/50,0,0,21/50)\) in the electric Schur inequality. Allocate the
remaining onsite coefficient across at most six incident bonds. Put
\(m_q=(1-r_q)/2\) and

\[
 \delta=4\sum_qd_q(\sqrt{m_q^2+1/4}-m_q),\qquad
 w_q=(105C_q/2-\delta)^{-1}.                              \tag{4}
\]

The exact positive link margin is \(63/4-\delta>0\). Use the full-Fock
occupation Schur shifts for the four higher-dimensional native blocks and
the exact 256-state rational positive-pivot shift for the charge-six block.
If their sum is \(c\), the actual physical Hamiltonian obeys

\[
 aH\ge-90V-Bc-P/7+\frac{63}{50}(D_L+D_e).                \tag{5}
\]

Here \(B\) and \(P\) count actual interior bonds and plaquettes. Boundary
onsite slack is nonnegative. The tiled physical trial supplies
\(aE_0+90V\le Vg\), with the exact source-derived gain \(g\) from the
member-zero variational receipt. Since \(B/V\le3\), \(P/V\le3\) and
\(a\epsilon/V\le1/128\), (5) proves

\[
 \frac{\langle D_L+D_e\rangle}{V}
 \le\frac{g+3c+3/7+1/128}{63/50}
 <3.40203<\frac{341}{100}.
 \quad p_n>\frac{59}{400}.                              \tag{6}
\]

Every radical is bounded upward and the charge-six shift is certified by
exact rational pivots. No Haar expectation is substituted for the selected
state. The rational retained fractions select an operator inequality; they
do not alter the source action.

## 3. The full-CAR source moment

The onsite term satisfies \([H_{\rm on},A]=9A/a\). The electric and
magnetic terms commute with (1). For any link values, the one-particle
hopping image of each of its three orthonormal positive-beta modes is
orthogonal to the original cell and to the other two images. The squared
norm of each image is \(z_x/(2a^2)\), where \(z_x\) is the actual number
of neighboring cells. This follows from the literal Wilson matrices:

\[
 P_+S_i^\dagger S_iP_+=P_+S_iS_i^\dagger P_+=\tfrac12P_+,
 \qquad R(g)^\dagger R(g)=I.                             \tag{7}
\]

On the full link-function space, the residual consequently has the same
norm as \(\sqrt{z_x/(2a^2)}\,\mathcal R\), where

\[
 \mathcal R=b_1^\dagger a_2^\dagger a_3^\dagger
 +a_1^\dagger b_2^\dagger a_3^\dagger
 +a_1^\dagger a_2^\dagger b_3^\dagger.                   \tag{8}
\]

These six modes are orthonormal. On all 64 CAR states,
\(G=\mathcal R^\dagger\mathcal R\) satisfies
\(G(G-I)(G-2I)(G-3I)=0\), and eigenvalue three occurs. Thus
\(\|\mathcal R\|^2=3\). Multiplication by link coordinates is a direct
integral of these pointwise matrices. The physical cutoff and Gauss
projection commute with \(A\); compressing the residual cannot increase
its norm. This proves on the actual finite physical carrier

\[
 \|[H,A]-9A/a\|^2\le\frac{3z_x}{2a^2}.                 \tag{9}
\]

For the original positive measure of energy transfer
\(\mu_n(d\omega)=(4V)^{-1}\sum_{x,\alpha,b,i,f}
\rho_i|A_{fi}|^2\delta_{E_f-E_i}(d\omega)\), stationarity gives

\[
 \mu_n(\mathbb R)=p_n,\qquad
 \int(\omega-9/a)^2d\mu_n
 \le\begin{cases}9/2&n=0,\\9/a^2&n\ge1.\end{cases}     \tag{10}
\]

No success normalization occurs. Chebyshev's inequality and (6) imply

\[
 \mu_n([1/(2a),35/(2a)])
 \ge\frac{59}{400}-\frac{36}{289}>\frac1{50}.            \tag{11}
\]

At level zero, the sharper original band is

\[
                    \mu_0([3,15])>\frac{23}{32}.        \tag{12}
\]

The charge-conjugated hole operators have the opposite native increment
and exactly the same complete source energy laws. This follows from the
finite source unitary and the conjugated apparatus, without an active
operation on the retained reference.

## 4. Complete energy records and the refinement boundary

The existing two-clock finite spectroscopy uses the same \(H\), energy zero
and source time. Resolution \(1/(8a)\) and each retained tail and clock
budget \(1/4096\) leave original reported transfer weight at least \(1/64\)
inside \([1/(4a),71/(4a)]\). The receipt uses the full inherited error
bound. Initial energy recording preserves the source marginal \(\rho\);
its disturbance of the source-reference correlations is retained in the
complete quantum output. A classical record is copied, not an arbitrary
quantum state.

Since \(E_i-E_0\ge0\), every exact event of (11) has
\(E_f-E_0\ge1/(2a)\), with no additional probability loss. At level zero,
cooling also gives original final-excitation weight greater than \(21/32\)
in \([3,16]\), by discarding at most \(\epsilon/1=1/16\) original weight
with initial excitation greater than one. This inference does not require
computing \(E_0\).

Equation (11) proves that at least \(1/50\) original odd-probe energy
weight escapes every fixed upper source-energy band as \(a\to0\). The
unscaled family of these local composite measures is therefore not
uniformly energy tight. A fixed finite-energy component can coexist with
this escaping weight; (11) proves neither its absence nor its survival.
It does not decide the critical-strength density Weyl response. Those
Weyls commute with every native number and CAR parity, so they cannot
supply the odd occupation change of (1) from a fixed occupation input.

The literal source scales as \(H_n=a_n^{-1}\widehat H_{L_n,\ell_n}\)
with fixed dimensionless action coefficients. Fixed physical source-energy
bands require dimensionless differences of order \(a_n\). The populated
band (11) has order-one dimensionless energy. Identifying a relativistic
particle requires a controlled low-energy branch with the corresponding
dispersion, interactions and clock conversion. No rest mass in GeV, pure
Lorentz spin or continuum particle identity follows from this populated
finite gauge-singlet odd response.
