# Finite source strength in complete collective geometry

Let N be a positive integer, q=1/8 and omega=1. In each retained native
source/probe energy sector the declared operator is

```
Q=(q/N) sum_i Z_i,
H=(omega/2) sum_i(I-X_i)+Es exp(Q)+Ep exp(F Q),  0<=F<=1.
```

All eight native energies at each matter site, their references, and all
2^N geometry states belong to the source. The spectral statements below
concern a fixed matter-energy sector. They do not replace the full matter
spectrum by the geometry gap.

## Two spin blocks determine the full geometry gap

Total spin commutes with H. Its spin-j multiplicity is
binomial(N,N/2-j)-binomial(N,N/2-j-1), with the second binomial zero at j=N/2.
In the m basis the diagonal is omega*N/2+Es exp(2qm/N)+Ep exp(2qFm/N),
and the adjacent off-diagonal entry is
`-(omega/2)sqrt((j-m)(j+m+1))`.

The single-bit flip graph is connected and all its off-diagonal entries
are strictly negative. Perron positivity gives a unique positive ground
vector. Permutation invariance places it in j=N/2. Each spin block also
has a positive ground vector. Extend the spin-j vector by zero into the
spin-(j+1) m grid. Every internal negative hopping increases in magnitude,
so the Rayleigh energy decreases. For j=0, adding a small positive boundary
component lowers the energy to first order. Thus sector ground energies
strictly decrease with j. No ordering of their higher eigenvalues is used.

For N>=2 the full geometry gap is exactly
`min(E1(j=N/2),E0(j=N/2-1))-E0(j=N/2)`. For N=1 only its spin-half block
exists. All other spin sectors and their multiplicities are covered by the
ordering proof. Outward diagonal bounds enclose the two required spectra
by matrix ordering; exact Sturm sign variations use rational off-diagonal
squares. No square-root rounding of those off-diagonal entries enters the
inertia certificate. Independent complete bit matrices check the result.

## Uniform finite-energy comparison

Put B=Es+F Ep and M>=Es+Ep. The comparison
`Hlin=Hg+Es+Ep+B Q` has a product ground with
`Omega=sqrt(omega^2+(2Bq/N)^2)` and
`m=<Z>=-2Bq/(N Omega)`. Convexity gives a positive remainder

```
R=Es(exp(Q)-1-Q)+Ep(exp(FQ)-1-FQ),
0<=R<=M exp(q) Q^2/2.
```

For a product spin state with mean m,

```
mu2=q^2[1/N+(1-1/N)m^2],
mu4=q^4/N^4 [N(3N-2)+N(N-1)(6N-8)m^2
             +N(N-1)(N-2)(N-3)m^4].
```

Use their values at m²=x²/(1+x²), x=2Mq/(N omega), as uniform upper
bounds. Min-max gives gap at least `d=omega-M exp(q)mu2/2`. If d>0, the
centered product-ground residual bounds its angle from the actual ground:
`sin(theta)<=epsilon=M exp(q)sqrt(mu4)/(2d)`.

For A=exp(FQ), the independent-copy variance identity and its scalar
Lipschitz constant give `sigma(A)<=F exp(q)q/sqrt(N)`. Write the actual
ground as cos(theta) times the product ground plus an orthogonal vector.
Centering A bounds the cross term by 2epsilon*sigma(A). The quadratic
term is at most `2||A-I||epsilon²<=2Fq exp(q)epsilon²`. This centering
retains the small observable scale; a full-operator norm would discard it.

The product expectation divided by the zero-matter vacuum expectation is
`[1+m tanh(Fq/N)]^N`. Put k0=2q²/(N omega) and
`kappa(F)=2q tanh(Fq/N)/omega`. Taylor's theorem on [-1,0] gives
`0<=(1+z)^N-1-Nz<=N(N-1)z²/2`. Also
`1-1/sqrt(1+x²)<=x²/2`. These yield, uniformly in Es+Ep<=M,

```
|<A>/cosh(Fq/N)^N - [1-kappa(F)B]| <= F C,
C=2 exp(q)q[epsilon/sqrt(N)+epsilon²]
  +k0 M x²/2+N(N-1)x²(q/N)²/2.
```

Let E(Es,Ep,F) be the unique geometry ground energy. Hellmann--Feynman
gives its Ep derivative <A>. Integrate the difference between source-on
and source-off ground expectations from probe energy zero to Ep. The
additive source and probe self-energies cancel. The mutual energy differs
from `-k0 F Es Ep` by the relative bound encoded in `uniform_certificate`:

```
(b-1)+b q²/(3N²)+2b C/(k0 Emin),
b=1/[1-q²/(2N)],  Es>=Emin>0.
```

Here `cosh(Fq/N)^N<=b` and
`F k0(1-q²/(3N²))<=kappa(F)<=F k0`. All 49 positive native pairs lie
in the certified M=7 domain. Differences between the inherited first and
third profile values bound their central forces, including the two error
budgets. The calibrated huge N is evaluated analytically, without building
its carrier. The M=10^6 result is a declared energy-domain extension, with
no identification of additional physical particles.

## Complete autonomous quench and two observation contexts

Prepare the geometry vacuum and retained source/probe energy labels.
Under Hlin each spin has
`m(t)=-2h omega/(omega²+4h²)[1-cos(sqrt(omega²+4h²)t)]`, h=qB/N.
For every time up to t, |m|<=4Mq/(N omega). The same fourth-moment bound
and Duhamel's identity give
`sin(angle(psi_H(t),psi_lin(t)))<=min(1,t M exp(q)sqrt(mu4)/2)`.
The centered observable bound applies with this angle. All coordinate
labels q(2k-N)/N are retained. The product comparison has exact binomial
probabilities, and the actual complete distribution lies within the stated
total variation bound.

Use two independently prepared runs, at delays zero and tau=355/113, with
equal weights. Measuring the first run does not alter the second run.
The linear response factor differs from one by `(1+cos(omega*tau))/2`.
For finite B its additional error is at most
`x²(1+omega*tau/4)`. The binomial quadratic remainder and full-H Duhamel
error give the complete source-on/off normalized-lapse bound in the receipt.
This constructs a finite-time coarse record without assuming equilibrium
preparation. The two delays, context bit and blank output memory are
supplied apparatus resources; no autonomous acquisition of those resources
or nonlinear spacetime limit is inferred.

The initial energy is
`Es cosh(q/N)^N+Ep cosh(Fq/N)^N`. A complete Q measurement removes every
X expectation, so its work is exactly `omega*N/2-<Hg>`. Positivity and
energy conservation enclose this work by `N omega/2-Mb` and `N omega/2`.
It is extensive even though the gravitational response is small. The
finite-control records retain this work and every coordinate outcome.

The inherited N was chosen using observed G. Finite-strength universality
does not turn that calibration into an independent prediction, identify the
strong-sector spatial carrier, or establish an Einstein continuum.
