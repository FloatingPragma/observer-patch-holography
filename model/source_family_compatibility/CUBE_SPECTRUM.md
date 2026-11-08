# Complete source-energy records for cube wavepackets

The eight-site full occupation source at c=1/8 and U=1 has a certified
wavepacket-dependent addition/removal response on its unique centered
ground. The calculation retains the original probability of each charge
branch and every record of a finite clock instrument. Its energy windows
are bounds on the actual spectral measure. The open-cube wavepacket
labels are prepared coefficients, without an identification as momentum
eigenvalues or a relativistic particle mass.

## One source, ground and positive excitation clock

Use the exact signed cube action K=H-N/2, the interleaved sixteen-mode CAR
basis, and the full N=8 singlet ground g from
`model/centered_source_ground/gap.py`. This source is unitarily expressed
in its hidden-spin frame; no spin-occupation projection or aligned
reference replaces its ground. The shared `cube_response.certified_ground`
replays that full-gap certificate and verifies a rounded 4900-component
trial p by exact Rayleigh, residual and spectral-separation arithmetic.
Let p_hat=p/||p||. It supplies

    D(|p_hat><p_hat|,|g><g|) <= delta,
    E_lower <= E0 <= E_upper,

where D is half-trace distance. The fixed reported allowance is
delta=1/1920; a trial that exceeds it is rejected. The energy enclosure is
[-206203/50000,-82481/20000]. Each reported number includes these reference
errors. Ground preparation is a declared service, without a
claim that the rounded vector was prepared by an executed native tape.

For k in {0,1}^3 define

    f_k(x)=(-1)^(k.x)/sqrt(8),
    a_(k,s)=sum_x f_k(x) a_(x,s),   L=K-E0.

Each f has norm one, so a a*+a* a=I and the two Kraus operators a*,a
form a complete addition/removal instrument on every source input. The
two successful charge labels have no common conditional normalization.
The odd transfer/reference and coherent clock services are supplied.
Physical service-failure records and their error allowance are additional
to the ideal table; no such outcome can be discarded in its implementation.

The exact source particle-hole involution obeys C K C*=K and
C a_(x,s) C*=(-1)^(x1+x2+x3) a_(x,s)*. Uniqueness fixes its action on g
up to phase. Thus the removal measure for k equals the addition measure
for k xor (1,1,1), at every excitation-energy set. Their total weights
sum to one. Spin-singlet invariance makes the two spin-component tables
equal. The producer checks the full centered clock under the literal
involution and all sixteen CAR transformations. Its inherited gap proof
checks the full-source spin algebra. This does not identify a physical
charge-conjugation or Lorentz-spin representation.

The calculation acts on the complete N=9, Sz=1/2 block of dimension 3920.
It verifies that this is a reducing block of the original integer matrix
16K. The other charge and spin responses follow from the exact source
symmetries, rather than a comparison of fitted spectra.

## Original energy moments without an extensive norm loss

Let v=a* g and B=[K,a*]. Then

    mu0=||v||^2,
    mu1=<a* g,Bg>,
    mu2=||Bg||^2.                                        (1)

These are original moments of the positive excitation measure, since
L a* g=B g. In the trial expectations the two terms of the commutator
are evaluated separately with the same source matrix; p_hat is never
assumed to be an exact ground eigenvector.

The source-generated composite CAR identity in `COMPOSITE_CAR.md` gives
B=a*(h f)+(U/2)b*(f). On this cube the literal one-body numerator obeys
(16h)^2=3I. Consequently

    ||B|| <= sqrt(3)/16+1/2 <39/64 = b.

The full-source b family is canonical on the same occupation space.
It adds no independent fermionic modes. The trial errors in (1) are
bounded by delta, 2 delta b and delta b^2, respectively: the first and
third are expectations of effects with ranges [0,1] and [0,b^2], and the
Hermitian part of a B has norm at most b. Interval division by the
strictly positive original mu0 encloses the conditional mean and variance

    mean=mu1/mu0,    variance=mu2/mu0-(mu1/mu0)^2.          (2)

The conditional quantities are reported alongside their original mu0;
they do not replace its preparation probability. The code also encloses
the first four original moments of the shifted source clock by direct
integer actions and the separate extensive operator norm bound.

## A complete finite-resolution clock measurement

Use the exact rational clock reference

    X=K+33/8,    Uclock=exp(-i X/2),    n=256.

On both charged response blocks, 1/8<X<=45/8. The lower charge bound is
the inherited full-source theorem. The upper bound follows from the
twelve bond norms c and the nonpositive centered onsite term. The
producer additionally checks a row bound of 90 on the actual integer 16X.
The excitation energy is X-(E0+33/8); the certified ground-energy
interval is retained when converting its records to excitation windows.

Prepare a uniform n-address clock register, apply its controlled powers
of Uclock, and Fourier analyze that register. Every record j=0,...,n-1
has the source Kraus operator

    F_j(X)=(1/n) sum_(t=0)^(n-1) exp(2 pi i jt/n) Uclock^t.

Fourier orthogonality gives sum_j F_j*F_j=I exactly. Thus the full
operators F_j a* and F_j a retain all 512 charge/clock outcomes for each
profile and spin setting. A uniform choice of the sixteen profile/spin
settings adds its explicit factor 1/16. Coarsened tables retain every
underlying j; their reported groups partition all 256 records.

For the original addition measure dmu(X), the actual record probability is

    q_j=integral |F_j(X)|^2 dmu(X)
       =[n c0+2 Re sum_(l=1)^(n-1)(n-l)c_l exp(2 pi i jl/n)]/n^2,
    c_l=<a* g,Uclock^l a* g>.                             (3)

The Fejer-shaped record is not a sharp spectral projection. In particular,
the source can produce several clock labels on a single eigenvector;
the table does not interpret those labels as distinct particles or poles.
The supplied finite instrument compiler can realize this finite menu to
a specified complete-output tolerance. No large source clock tape is
reported as executed by this calculation.

## Integer propagation and numerical enclosure

Only the trial eigensolver is used as a proposal. Admission of the source
reference, subsequent clock propagation, moments and Fourier coefficients
uses integer or rational inequalities. With p's inherited denominator 2^23,
form each unnormalized charged vector w=sum_x (-1)^(k.x) a_x* p/2^23.
Its original probability uses the denominator 8||p/2^23||^2. The code
checks ||w||<=3, 3/4<=||p/2^23||^2<=2 and stores clock vectors with
denominator 2^46.

One step uses the degree 32 Taylor polynomial in -i(16X)/32, rounding
each complex term back to that dyadic grid. Every sparse product and sum
has an int64 accumulation prebound; norms and correlation dots are summed
over Python integers. For dimension d=3920 put

    q=ceil(sqrt(2d))/2^47,
    epsilon_step=4*3^33/33!+32*21*q.                      (4)

The clock-step norm is at most 3. For an input of norm at most 4, the
unitary Taylor integral remainder is at most 4*3^33/33!. The accumulated
rounded-term error is at most 32 exp(3)q<32*21q. The norm bound 4 is checked
on every stored step. Unitarity then bounds the error at time l by
l epsilon_step, without exponential growth in l.

Machin's formula pi=16 atan(1/5)-4 atan(1/239), alternating rational
remainders and a rational Taylor expansion enclose each Fourier root.
Its stored 42-bit coefficient has complex error below 2^-41. No sampled
trigonometric error is treated as an outward certificate. If
r=||p/2^23||^2, (3) has per-record numerical error at most

    epsilon_num=epsilon_step*(n^2-1)/(8nr)+2^-40.          (5)

Indeed |c_l-computed_c_l|<=3l epsilon_step/(8r), and
sum_l (n-l)l=n(n^2-1)/6. The computed correlation magnitude is at most
3*4/(8r)<=2, so Fourier rounding contributes the last allowance.
Admission checks epsilon_num<14/10^9 and reports this fixed allowance.
An arbitrary group J of clock outcomes therefore has total error
|J| epsilon_num+delta. The ground error occurs once for the group effect;
it is not multiplied by the number of retained records. Every reported
rational interval is rounded outward. Fixed rational output envelopes are
admitted only when they contain the freshly computed exact interval; the
envelopes themselves supply no spectral premise. Individual records, complementary
records and their grouping are all included in the machine result.

## True spectral windows from the finite clock records

The clock period is 4pi. For any energy and any tolerance d>0, the
probability that a clock record's wrapped energy differs by at least d
is at most

    tail(n,d)=min(1,1/[2(nd/(4pi)-1)]),                   (6)

when the denominator is positive; otherwise use 1. To prove (6), the
record probability at wrapped grid distance u<=1/2 is at most
1/(4n^2 u^2), using |sin(pi u)|>=2u. On each side, grid distances
have spacing1/n. The sum beyond a=nd/(4pi)>1 is bounded by
(1/2)(1/a+1/a^2)<=1/[2(a-1)]. The code uses pi<22/7, so its rational
denominator replaces a by 7nd/88 and enlarges the tail.

For a contiguous selected set J with extreme source-energy centres A,B,
every source energy outside [A-d,B+d] has selected-record probability at
most tail. In the displayed windows there is no wrapping across the
source support [1/8,45/8]. Thus, before any conditioning,

    mu([A-d,B+d]) >= q(J)-tail(n,d)*mu0.                  (7)

Equation (7) consumes the original q(J), its complete-output error and
the upper bound on mu0. Ground-energy uncertainty enlarges the final
excitation interval. Clipping its lower endpoint to 1/8 is justified by
the full-source charged gap. These positive window weights belong to
the actual source measure, without an eigenvalue or Ritz-node assumption.

## Certified source energy table

The following rational decimal intervals are admitted for every indicated
profile and either spin component. Units are those of the supplied source
K at c=1/8,U=1. The weight and last column are original probabilities for
the chosen setting. Only the two middle energy statistics condition on
the reported nonzero addition probability.

| Addition profiles k | Original weight | Conditional excitation mean | Conditional energy variance | Original weight in [0.125,0.93539] |
|---|---:|---:|---:|---:|
| (0,0,0), (0,1,0), (1,0,0), (1,1,0) | [0.58337,0.58444] | [0.47976,0.48282] | [0.00986,0.01390] | at least 0.46633 |
| (0,0,1), (0,1,1), (1,0,1), (1,1,1) | [0.41556,0.41663] | [0.52298,0.52738] | [0.00830,0.01456] | at least 0.33525 |

For removal at a given k, use the opposite row, by k xor(1,1,1).
The complementary charge and all 256 phase records are retained. The
displayed true-energy window is obtained from phase records 6 through 15
and tolerance 1/5, with tail allowance 55/338. It is not the uncorrected
histogram bin interval. A uniform selection among the sixteen profile/spin
settings multiplies each table weight by 1/16.

The source energy table and every phase-record enclosure are returned by
`cube_spectrum.certificate()`. The finite variation across prepared
profiles demonstrates a wavepacket-dependent response. A physical
dispersion requires operational spatial momentum on a growing source,
an appropriate state and clock limit, and an observable-sector attachment.
Positive energy variance of a wavepacket does not exclude a sharp
spectral subbranch or a fixed mass with varying momentum.
No continuum pole, scattering cross section or physical mass follows
from the finite table.
