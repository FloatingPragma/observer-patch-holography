# Speed concentration from one source susceptibility

The actual source is

    H_B=J sum_<xy>(S_x.S_y-I/4)-B sum_x epsilon_x S_x^1,
    K_B=H_B-E_B>=0,  J>0, B>0.

Use its selected positive-field ground and its own time evolution. The
joint prescription is p=n^-2, B=J n^-5. The concentration statistic below
uses an exact source first moment and one weighted inverse susceptibility.
Its smallness is not proved for the growing cubic ground. It supplies a
finite certificate with original probabilities, rather than a physical
particle assignment or a dispersion inferred from a broad energy band.

## The static inputs and the positive ratio

Let psi=A_f Omega, A_f=sum_x f_x S_x^a for a transverse component a=2,3.
Take a real packet with spatial Fourier support in a small annulus away
from zero, in the even-translation chart used in
[the source shell diagnostic](../spin_shell_diagnostic/README.md).
Write

    e(k)=sum_j(1-cos k_j), s(k)=sqrt(2e(k)),
    0<s_min<=s(k)<=s_max on that support,
    X=K_B/[J s(P_space)],
    w=||psi||^2, H=<psi,X psi>, chi=<psi,X^-1 psi>.       (1)

P_space denotes spatial translations; it is not the canonical internal-spin
field P. K_B and these translations commute. Their joint spectral calculus
defines X on the packet support. Assume w>0 and finite H,chi. The inverse
takes value infinity at zero, so finite chi excludes a zero-energy part of
psi. For each positive B this also follows from the source gap and the
vanishing transverse ground mean. The same ratio can be used on a genuine
finite periodic translation grid; open-cube reflection characters supply no
small-momentum identification.

The exact positive-field ground identity in
[ARCTAN_CHARGED_BANDS.md](ARCTAN_CHARGED_BANDS.md), Eq. (1), applied to the
profile with Fourier transform fhat/sqrt(s), gives

    H=(kappa_B/2) integral s(k)|fhat(k)|^2
        +(B m_B/(2J)) integral |fhat(k)|^2/s(k).         (2)

Here kappa_B=-omega_B(S_x^1 S_y^1+S_x^2 S_y^2) on a nearest bond and
m_B=omega_B(epsilon_x S_x^1). Integrals use d^3k/(2pi)^3, or the matching
normalized torus sums. Indeed the graph Dirichlet form is integral
2e(k)|fhat|^2, and e/s=s/2. Equal-time two-point correlations determine w.
Thus H and w are static source quantities; chi is the sole dynamical
statistic in (1). The static correlations themselves need certified
enclosures; their values are not supplied by symmetry alone.

Equation (2) is a microscopic source identity. The inequalities for the
projected continuum Stone vectors do not justify inserting it as an
equality after a reducing projection or exchanging limits. Fourier support
and the first-moment identity extend from finite profiles through their
stated source form bounds.

## An uncut original-weight inequality

Let nu be the positive joint-ratio measure of psi, with total mass w.
For any specified c>0 define

    D_c=H-2cw+c^2 chi
       =integral (x-c)^2/x dnu(x)>=0.                  (3)

For every epsilon>0,

    nu({|X-c|>=epsilon}) <=(c+epsilon) D_c/epsilon^2.  (4)

For x>=c+epsilon the function (x-c)^2/x is increasing and has minimum
epsilon^2/(c+epsilon). If 0<x<=c-epsilon, its minimum is
epsilon^2/(c-epsilon), which is larger. This second set is empty when
epsilon>=c. These observations prove (4) without an upper energy cutoff.
For R>=2c, the same integrand is at least x/4 on x>=R, hence

    integral_(x>=R) x dnu(x)<=4D_c.                   (5)

Taking c=H/w gives D_c=H d with d=H chi/w^2-1. Thus d tending to zero,
bounded positive c, and a common limiting c across the admitted packet
profiles and directions give speed concentration and first-moment tail
control. For a fixed common c the sufficient condition is D_c/w tending
to zero. On k=pq in a bounded annulus, s(pq)/p tends uniformly to |q|,
so the speed tube becomes the rescaled cone tube. A positive D_c for the
whole response does not exclude a sharp subbranch in that response.

If an admitted complete preparation has successful vector psi/Lambda,
divide (4) by Lambda^2. Its original inside probability is at least

    max(0,w-(c+epsilon)D_c/epsilon^2)/Lambda^2.         (6)

Every preparation and detector complement is retained. There is no
division by observed success. Given intervals w in [w_l,w_u], H in
[H_l,H_u], chi in [chi_l,chi_u], the executable uses

    D_+=H_u-2c w_l+c^2 chi_u,
    outside<=min(w_u,(c+epsilon)D_+/epsilon^2).        (7)

It rejects reversed, negative or inexact intervals, invalid preparation
normalization, and data contradicting positive Cauchy-Schwarz or (3).
These are necessary admission checks, not a full moment-realizability
test. A declared speed is an explicit hypothesis being tested, not a fit
to an empirical particle speed.

## A complete finite-field inverse score

For zeta>0, put chi_zeta=<psi,(X+zeta)^-1 psi>. A source gap g>0 gives

    X>=g_X=g/(J s_max),
    chi_zeta<=chi<=(1+zeta/g_X)chi_zeta.               (8)

After the complete preparation, the joint energy/translation detector
assigns score

    q_zeta(E,k)=zeta J s(k)/(E+zeta J s(k))            (9)

to successful records in the packet's momentum support, and zero to
all preparation failures. It retains every energy bin, spatial label,
and native failure. The ideal original mean is zeta chi_zeta/Lambda^2.
This is a bounded score, not a successful filter branch reweighted to one.

The complete periodic translation/energy instrument is given in
[native_joint_spin_spectroscopy](../native_joint_spin_spectroscopy/README.md).
Its raw generators are g_j=T_j Z_all. A physical S^3 Fourier packet has
the ordinary raw label, whereas an S^2 packet has that label shifted by
(pi,pi,pi). Recenter the latter labels before evaluating s, or use S^3
and the internal-circle equality. The projection includes both charge
representatives in the stated even-translation chart. Raw reflection
labels or unrecentered S^2 labels cannot be substituted in (9).
The bounded resolvent score, finite residual enclosure and explicit
all-record error conventions are in
[spin_relaxed_response](../spin_relaxed_response/README.md), Eqs. (5),
(7) and (10). They supply the finite ingredients; this module executes
no large joint native tape. A bound on the norm of A_f alone does not
supply a controlled native preparation. Its Pauli-sum or Weyl compilation,
ground preparation and label certificate remain required inputs.

Use exact finite spatial labels and clamp negative numerical energy labels
to zero. If the true-label energy error is at most d_E except with original
probability ell, an absolute score error bound is

    xi=d_E/(zeta J s_min)+ell+d_out+d_in+d_sample,      (10)

capped at one. Here d_out is the unhalved complete-output L1 error and
d_in is the unhalved input trace error; extra compilation errors must be
included. The energy-origin error belongs in d_E, and clock aliases belong
in ell. The derivative in E of (9) has absolute value at most
1/(zeta J s_min). Approximate spatial labels need a separate certified
error; the module does not silently assign them exact values.

If a measured interval for the mean score is [q_l,q_u], enlarge it by xi
and intersect with [0,1]. Multiplication by Lambda^2/zeta encloses
chi_zeta. Apply (8) to obtain chi_u in (7). Only declared instrument
normalizations and regularizers are divided out.

The actual uniform source bound is g=B/8, including even periodic boxes
and the positive-field ground representation. If |k|<=p r with the joint
prescription, then

    g_X>=1/(8r n^3).

Choosing zeta=rho/(8r n^3), 0<rho<=1, makes the multiplicative inflation
in (8) at most 1+rho. At the outer symbol bound the corresponding source
energy regularizer is J rho/(8n^5). This explains the growing time and
precision costs; it is not a uniform susceptibility or state-convergence
modulus. The finite gap cannot be inserted after the B=0 limit.

## Connection to bounded charged preparations

Apply (4) to the normalized canonical P packet used in the arctan
construction, with its actual source scaling included in w,H,chi. Let Pi
be the speed-tube projection, completed with the full internal charge
orbit and invariant under charge conjugation. The two transverse
components have the same leakage bound by the source internal circle.
If eta^2 bounds ||(I-Pi)P Omega||^2, Eq. (13) of
[ARCTAN_CHARGED_BANDS.md](ARCTAN_CHARGED_BANDS.md) yields

    ||Pi B_+ Omega||^2=||Pi B_- Omega||^2
       >=max(0,d_comm-eta)^2/(2||P Omega||^2).         (11)

Thus a proved source limit with eta tending to zero would supply a
nonzero original-weight charged component on the cone, with all other
records retained. It would not require the nonlinear bounded creator's
entire response to be a single pole. Passage to the actual bounded-history
representation needs its established limit order and domain/projection
control; (2) is not an asserted equality for a projected field. No source
estimate forcing (3) to vanish is supplied here.

## Exact source control and the remaining dynamical datum

`finite_source_control()` reconstructs the full four-dimensional dimer
from `native_ordered_ground.source.matrix`, J=1,B=2/3, and verifies its
physical frame. Its exact ground has E0=-4/3, m=kappa=2/5. With
A=S_0^x-S_1^x and Lambda=1 the success norm is 4/5; I-A^2 is the other
Kraus map. The complete two preparations times three distinct energies
give six records, including zero weights and total probability one.

For this finite control the ratio divisor is fixed to one. It is not a
spatial momentum label. The exact moments are

    H=16/15, chi=3/5, c=4/3, D_c=0.

The source Dirichlet and selection terms give the same H. At zeta=1/6
the original score is 4/45; using the exact gap 4/3 in (8) returns the
exact upper chi=3/5. Using the valid coarser B/8 gives a larger interval,
as it must. This finite equality is not a cubic or scaling-limit result.

The source first moments, CCR, reflection bounds and force-slack lower
bound also admit the sharp and two-speed implication controls in
[RESOLVENT_CLOSURE.md](../spin_dynamical_memory/RESOLVENT_CLOSURE.md),
Eq. (12). With kinetic coefficient 12kappa/m and gradient coefficients
kappa/(2m),3kappa/(2m), the canonical P ratio has

    H chi/w^2-1=7-4sqrt(3)>0.

One gradient coefficient kappa/m gives zero defect with the same H and
chi, but a different static w. Both satisfy the stated source bounds;
neither construction is substituted for the cubic ground. Inverse
susceptibility saturation together with the source structure factor is
the discriminating datum. The source current Ward rule fixes a total
inverse-current moment, without fixing its frequency distribution.
No additional concentration identity follows from those constraints.
