# Detecting a surviving positive-energy response

The [local spectral-access theorem](../source_family_compatibility/SPECTRAL_WINDOW_ACCESS.md)
retains complete original spectral weights in a window that includes zero.
The [selected cooling theorem](../particle_joint_reconstruction/POSITIVE_ENERGY.md)
controls negative transfers and supplies a continuous time-filtered response
subsequence. Neither statement excludes a limit concentrated at zero.
The following bound tests that distinction using an unnormalized correlation
contrast, with ultraviolet and negative-transfer allowances explicit.

Let a positive finite spectral measure nu have mass at most M0 and first
absolute moment at most M1. Set

```
F(t)=integral exp(it omega) dnu(omega),
D(t)=F(0)-Re F(t)=integral (1-cos(t omega)) dnu(omega).
```

Suppose a complete experiment bounds D(t) below by d, and |t|<=t0.
For 0<delta<R suppose also `nu((-infinity,-delta])<=eta`. Partition the
line into `|omega|<delta`, `omega<=-delta`, `delta<=omega<=R` and `omega>R`.
On the first set use `1-cos x<=x²/2`; on the other sets use `1-cos x<=2`.
The positive ultraviolet tail has mass at most M1/R. Hence

```
nu([delta,R]) >= max(0,
  [d-M0*(delta*t0)²/2-2*eta-2*M1/R]/2).                 (1)
```

Every weight is original. Conditioning on a successful filter or deleting
failure records changes the measure and is outside (1). If estimated
F(0) and Re F(t) have absolute errors e0 and et, then their difference
has error at most e0+et. Sampling, state, quadrature and coherent-control
errors must be included in those errors before using (1). A trace-distance
bound e for a bounded signed score in [-1,1] contributes at most 2e.
The two source clocks, preparations and filter definitions must agree.

## The actual primary allowances

For the selected member n, a physical source probe of squared norm at most
C² and triangular time filter of half-width T have

```
M0=C²,  M1=2*C²/T,
eta=min(C²,C²*2^(-n-4)/delta).                          (2)
```

These follow from the existing finite downward-transfer and filter bounds
on the same source. `selected_requirement` authenticates that member and
computes (2). The contrast itself is explicitly supplied; no primary
positive contrast has been computed by this package. Its selected example
uses zero contrast and derives zero certified positive-band weight.

A uniform strictly positive right side of (1), for a fixed compact interval
[delta,R], gives a nonzero positive-energy component of any weak spectral
limit. Indeed this interval is closed, and the portmanteau inequality gives
`nu_limit([delta,R])>=limsup nu_n([delta,R])`. Such a component prevents both
complete zero-energy collapse and loss of all response to infinite energy.
It does not isolate a pole, fix a particle mass or establish an interaction.
Extending it to a relativistic species needs the actual momentum, spin,
current and scattering attachments.

## Exact controls and failure cases

Finite diagonal Hamiltonians with integer energy differences, evaluated at
time pi/2, give exact cosine values 0 or ±1. The rational bound |t|<11/7
keeps (1) entirely rational in the executable controls. Original measures
with a positive excitation, a zero atom, an ultraviolet excitation, a
small negative tail and a mixed response are evaluated without removing
any row. The first and mixed controls give positive lower bounds; the
zero and ultraviolet controls give none. These finite controls are
separate from the spatial source and do not supply its missing witness.

The central next measurement is therefore a certified nonconstant
time-filtered correlation on the selected growing family, with enough
original weight to exceed the allowances in (1). The bound is a diagnostic
for that physical premise, not an assumption that it succeeds.
