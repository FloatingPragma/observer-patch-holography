# Coherent charge readout on the same source

The [source charge construction](../protected_defect_catalog/PROOF.md)
supplies a strongly continuous circle action on the transverse configuration
time-history space. It fixes the selected ground and commutes with the
actual positive energy generator and spatial translations. Write
`R(theta)=exp(i theta C)`. Since `R(2pi)=I`, the spectral values of C are
integers. They are internal spin charges. No physical electric charge or
relativistic particle label is assigned.

## 1. A complete detector that preserves a charge class

For a power of two N at least four, prepare a blank log2(N)-bit pointer in
the uniform superposition. Apply R(2pi j/N) controlled by its integer value
j, then the inverse finite Fourier transform, and record r. The complete
Kraus operators are

    K_r=(1/N) sum_(j=0)^(N-1) exp(-2pi i rj/N) R(2pi j/N)
       = 1_(C congruent to r modulo N),       0<=r<N.       (1)

The finite geometric sum proves (1) on each integer charge. Thus the
K_r are pairwise orthogonal projections summing to I, and the coherent
map V psi=sum_r K_r psi tensor |r> is an isometry. It commutes with the
source energy in each branch. These identities survive tensoring any
untouched reference. Every residue is retained, including zero-probability
residues. Measuring every local X and grouping the resulting strings has
the same residue probabilities but generally dephases within a residue;
that different instrument is not used here.

The finite eight-spin source has
C=sum_x epsilon_x X_x/2 with four positive and four negative source signs.
Its charge range is -4 through4. The source Walsh matrix W satisfies
W*W=256 I and diagonalizes this exact C. Therefore

    K_r=(1/256) W diag(1_(q congruent r)) W*.                (2)

The executable computes all16 matrices for N=16, their ranks,
orthogonality, completeness and commutation with the full selected H_B.
It computes the original probabilities in Z0 Omega_B, with the inherited
ground-vector error. No fine X strings are measured in this construction.
For this particular finite cube N=16 resolves every actual charge without
aliasing. The source's complete raw repair is the identity on the protected
spin code; its off-code outcomes are part of the implementation errors.

## 2. A literal source program and its error

For q protected spins with balanced signs and N=2^m, decode each spin to
its logical bit. Apply Hadamard to those bits. Charge then equals
`-sum_x epsilon_x b_x` in the logical computational basis. Prepare the m
pointer Hadamards. For pointer bit j and each spin x apply the controlled
phase with angle `-epsilon_x 2pi 2^j/N`. Apply the inverse Fourier circuit,
undo the logical Hadamards, and re-encode. No arbitrary rotation oracle
or new measurement service appears.

Every controlled phase exp(i theta |11><11|) consists of the actual native
J(-theta) and scalar I(theta) pulses. A Hadamard consists of I(-pi/2) and
H(pi/2). The inverse Fourier circuit retains its bit-order swaps, each
three of the existing phase-exact CNOT words. Counting decode/encode,
basis changes, all controlled phases and the inverse Fourier transform
gives the explicit source-pulse count

    L=14q+4m+2qm+m(m-1)+9 floor(m/2).                      (3)

It is55 for two spins with N=4, and222 for the actual eight-spin cube with
N=16. The two-spin execution streams every literal six-level source pulse
on all64 work rows and four coherent reference columns. It retains all
196 first-failure slots and all36 physical pointer labels, including32
reserved labels. Its complete ideal charge probabilities are1/2,1/4,0,1/4.
The record map preserves within-class amplitudes and all reference columns.
The eight-spin statement uses the exact matrix certificate and emitted
program; its full native tape is not executed in this package.

The source's existing dyadic integrated-port algorithm gives per-pulse
complete-isometry error at most29*2^(-b). Telescoping the actual whole
isometries gives operator error29 L 2^(-b) and unhalved trace error
58 L 2^(-b), uniformly in every input and reference. Literal binary64
commands instead receive their own exact angle and port-coefficient
bounds; their host floating execution is checked with its explicit
numerical admission. The complete coefficient bound is independent of
that diagnostic tolerance. Neither fixed binary64 precision nor a finite
register ceiling is used to assert arbitrary accuracy.

## 3. Quantified aliasing for the continuum preparation

Let F be a real, unit-L2 Schwartz test whose smooth Fourier support lies
in `1/2<=|k|<=1`, with the inherited normalized Fourier measure. The same
source continuum configuration fields Q3(F) and Q2(F) strongly commute
and obey `||Q2(F)Omega||^2<=6`. Put A_s=sin(s Q3(F)). Source rotations give

    R(theta) A_s Omega
       = sin(s[cos(theta)Q3+sin(theta)Q2]) Omega.

Use their joint spectral measure. Differentiation at zero in Hilbert norm
is justified by the finite second moments and dominated convergence.
Consequently A_s Omega belongs to Dom C and

    ||C A_s Omega||^2
       =s^2 ||Q2 cos(s Q3) Omega||^2 <=6s^2.               (4)

This is an unnormalized moment of the actual bounded preparation, not
an assumption of finite charge support. If f(L) is a contraction and
commutes with C, the same bound holds for f(L)A_s Omega. It applies in
particular to the complete resolvent-success branch lambda/(L+lambda).
The original sine preparation and filter probabilities are retained.

In the residues r=1 and r=N-1, any unwanted charge differs from +1 or -1
respectively by a nonzero multiple of N. Its magnitude is at least N-1.
The spectral Markov inequality and (4) therefore give, for the union of
both wanted residue outcomes,

    Pr(original preparation, unwanted aliased charge)
       <= min(1,6s^2/(N-1)^2).                            (5)

This is an absolute original-weight error. If the corresponding total
preparation probability is certified at least b>0, the conditional alias
error is at most min(1,6s^2/[b(N-1)^2]). A complete detector trace error
delta contributes at most delta/2 to an event probability, before any
division by a preparation probability. The same rules apply to energy
scores jointly measured with this charge instrument.

Equation (4) has a same-source prelimit estimate. For the actual lattice
Q2_p,Q3_p and its invariant selected state, differentiation of the bounded
rotation and a unitary Duhamel identity give

    ||C sin(sQ3_p)Omega||
       <= |s| ||Q2_p Omega||+(s^2/2)||[Q2_p,Q3_p]||.        (6)

For the Fourier test above and0<p<=1, the inherited transverse structure
bound `g(k+Q)<=sqrt(e(k+Q)/(2e(k)))<=3/|k|` gives
`||Q2_p Omega||^2<=6`. Poisson summation gives
`p^3 sum_x F(px)^2=1`: the Fourier support of F squared lies within radius2,
whereas nonzero reciprocal sampling vectors have length at least2pi/p.
Thus the exact source commutator bound gives `||[Q2_p,Q3_p]||<=p/2`.
The right side of (6) is at most `|s|sqrt(6)+s^2 p/4`.
For nonzero s and `p<=min(1,1/(2|s|))` its square is below7s^2; s=0 is
trivial. This supplies a concrete
prelimit charge-tail estimate; it does not supply a modulus for preparing
the infinite selected ground or for the separate continuum spectral limit.

Finally, sin(sQ3)Omega/s converges strongly to Q3 Omega. Projection onto
charge±1 therefore converges to Phi_±Omega/2. The inherited positive-energy
charged Stone vectors prove nonzero positive-energy charge±1 weight for
some rational small s. Charge conjugation gives equal±1 energy measures
for every such bounded sine preparation and every commuting real energy
filter. No effective small-s threshold, sharp dispersion or physical mass
is inferred solely from these second moments.
