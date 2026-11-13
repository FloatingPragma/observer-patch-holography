# Complete geometry response, fluctuations and information

Use the declared geometry Hamiltonian `Hg=(omega/2) sum(I-X)`, its unique
vacuum `|+>^N`, and `Q=(q/N) sum Z`, with q=1/8, omega=1. The native matter
energy labels remain controls of `H(Es)=Hg+Es exp(Q)`; the information
statements concern this geometry state family at Es=0 and Ep=0. They do not
bound a different observer architecture, extra parameter-bearing references,
feedback, or an interacting continuum.

## An exact normalization shortcut

All N+1 coordinate outcomes occur with probabilities binomial(N,k)/2^N.
Their moment generating function is `M(f)=cosh(fq/N)^N`. Thus

```
<Q>=0, v=Var(Q)=q^2/N,
Cov(Q,exp(Q))=M'(1)=q cosh(q/N)^N tanh(q/N),
Cov(exp(fQ),exp(gQ))/(M(f)M(g))
    =[1+tanh(fq/N)tanh(gq/N)]^N-1.
```

The inherited complete spectral-response coefficient therefore satisfies
`g_N=S Cov(Q,exp(Q))/(2 pi omega)` exactly. This is a relation between
fluctuation and weak response in the same declared source. An independent
coordinate-variance measurement, together with q, fixes N; an independently
known covariance fixes g_N. Inserting a covariance computed from an
observed-G-calibrated N is not an independent derivation of G.

For a=q/N and f,g in [0,1], use `tanh(t)>=t-t^3/3`, `tanh(t)<=t`, the
binomial theorem, and `exp(z)<=1/(1-z)` for z<1. These give the rational
all-N covariance bounds in the executable. Their leading covariance is
`fg v`: all profiles share the same Q. They are not independent local
geometry fluctuations. The leading covariance matrix has rank at most one;
the full nonlinear covariance need not have rank one.

A complete coordinate measurement removes every X expectation. Its vacuum
work is W=N omega u/2 GeV, with the unchanged electron anchor u. Hence
`W Var(Q)=q^2 omega u/2` exactly. When the model light speed is matched to c,
`G=g_N/u^2`, and `G W` equals `S q^2/(4 pi u)` times
`cosh(a)^N tanh(a)/a`. This factor lies between `1-a^2/3` and
`1/(1-v/2)`. Without speed matching the inherited invariant is
`G/(c_model/c)=g_N/u^2`; the strong-sector length is not silently equated
to the geometry length.

## Exact information in the entire carrier

In the X basis an r-excitation state's matrix element of exp(Q) from the
vacuum is `cosh(a)^N tanh(a)^r`; its excitation energy is r omega, and
its multiplicity is binomial(N,r). Nondegenerate ground perturbation theory
and normalization give

```
|ground'>=-sum_{r>0,labels} |r,labels><r,labels|exp(Q)|0>/(r omega),
I_ground=4 cosh(a)^(2N)
             sum_{r=1}^N binomial(N,r)tanh(a)^(2r)/(r omega)^2.
```

Here I is the pure-state quantum Fisher information,
`4(<psi'|psi'>-|<psi|psi'>|^2)`. For any complete measurement, applying
Cauchy--Schwarz to each outcome derivative bounds its classical Fisher
information by I. This is the standard quantum information bound of
[Braunstein and Caves](https://doi.org/10.1103/PhysRevLett.72.3439).
All coordinate-basis ground amplitudes are real and strictly positive near
zero source by Perron positivity. Their derivatives are constant on each
coordinate-count class by permutation symmetry. Summing all N+1 coordinate
outcomes therefore retains exactly I_ground. No exponentially many
unobserved bit labels are needed to achieve this local information.

Retaining r=1 gives the lower bound
`4v(1-a^2/3)^2/omega^2`. For the upper bound discard 1/r^2, then use
`cosh(a)^(2N)<=1/(1-v)` and
`(1+tanh(a)^2)^N-1<=v/(1-v)`. Consequently

```
4v(1-a^2/3)^2/omega^2 <= I_ground <= 4v/[omega^2(1-v)^2].
```

For R independently prepared copies, quantum information adds. The
Cramer--Rao inequality for a locally unbiased source-energy estimator gives
`Var(Es_hat)>=1/(R I_ground)`. Thus unit source-energy variance requires at
least `ceil(omega^2(1-v)^2/(4v))` such copies. This is an inference about
this parameter-encoding family, not a general lower bound on resources
available to all OPH observers. It does not infer a finite-difference
distinguishability bound from a local derivative.

For autonomous evolution from the same zero-source vacuum, integration of
each r-excitation phase gives

```
I_quench(t)=16 cosh(a)^(2N)
  sum_{r=1}^N binomial(N,r)tanh(a)^(2r)sin^2(r omega t/2)/(r omega)^2
  <=4 I_ground.
```

The diagonal vacuum phase cancels from the pure-state information. Waiting
arbitrarily long in this isolated geometry family cannot evade this bound.
No source-dependent phase reference was supplied to that calculation.

## Scope and verification

Independent complete bit matrices verify the spectral sum, all coordinate
probability derivatives, the equality of quantum and coordinate information,
finite-difference ground derivatives, cross-profile covariance, and work.
The huge calibrated N is evaluated by exact formulas and rational bounds;
its full carrier is not enumerated. The same positive complete architecture
at N, 2N and 4N keeps its native matter source and free gap, but changes G,
noise and measurement work. Declaring one hardware count is a valid premise;
no uniqueness requirement is imposed. To obtain an unfitted G, the count
or covariance requires independent calibration ancestry. No physical
metric-noise identification, common strong/geometry Hamiltonian, or Einstein
continuum follows from these finite-source identities.
