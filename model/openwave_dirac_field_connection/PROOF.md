# The electric field in a signed-band transformation

A dynamical compact gauge field cannot in general retain its bare electric
operator after a field-dependent matter diagonalization. Let E=-i d/dtheta
on periodic functions and U(theta) be a smooth periodic unitary. On smooth
vectors the product rule gives U* E U=E+A with A=-i U* U'. Consequently the
complete electric quadratic form transforms as alpha(E+A)^2. The derivative
of A and its square are part of this expression. Smooth U and its bounded
first two derivatives preserve the periodic Sobolev domains, so this identity
extends to the self-adjoint operators. This is a transfer identity, independent
of a choice of particle interpretation.

The executed diagnostic uses all integer fluxes and two signed internal bands:

```
D(theta)=[[m,t exp(-iq theta)],[t exp(iq theta),-m]],
epsilon=sqrt(m^2+t^2), H=alpha E^2+D+epsilon I.
```

Its alpha=1/32 and the three native charged mass parameters are inherited
inputs. The hopping t=1/16 is a declared diagnostic choice. This rotor is not
asserted to be a reducing sector of the complete native fermion source.
No particle-antiparticle dictionary or physical pole claim is transferred.

Set c^2=(1+m/epsilon)/2, s=t/(2 epsilon c), and
U=[[c,-s exp(-iq theta)],[s exp(iq theta),c]]. Then U*DU=epsilon sigma_z and

```
A=q [[s^2,cs exp(-iq theta)],[cs exp(iq theta),-s^2]].
```

For every integer n, the span of upper flux n and lower flux n+q is reducing.
The complete source is the orthogonal sum of these two-dimensional blocks.
Each has exact eigenvalues

```
epsilon+alpha[n^2+(n+q)^2]/2
  +/- sqrt([m+alpha(n^2-(n+q)^2)/2]^2+t^2).
```

All labels, signs and both charge choices are retained. Since D+epsilon I is
positive, omitted block energies are bounded below by
alpha min(n^2,(n+q)^2). This bounds the full integer-flux ground without a
cutoff-defined source. Rational square-root intervals certify the finite
minimum and the complementary floor.

In the free-band basis the field off-diagonal entry is
alpha c s [(n+q)^2-n^2]. It is nonzero for the executed controls. Omitting it
changes the spectrum and suppresses a transition present in the exact source.
Complete +/- readouts at the specified times include every outcome; all other
flux-block outcomes have zero probability for these reducing preparations.
Their floating probabilities and terminal work are numerical diagnostics;
the block formulas and completeness identities are exact. Preparation work
uses the full-source ground and a rational interval for its initial energy.
Band dephasing has its actual energy cost. It is not an energy-free record.

The calculation identifies the field connection needed by a signed-Dirac
completion of a positive-root matter model. It proves neither that the
positive-root member is inconsistent as a declared architecture nor that the
corrected diagnostic is a local Lorentz-covariant quantum field theory.
