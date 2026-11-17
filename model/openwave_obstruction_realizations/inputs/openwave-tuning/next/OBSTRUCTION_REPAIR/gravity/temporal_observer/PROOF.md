# The proper spatial regulator has a preferred-time ultraviolet sector

The exact continuum source at zero transverse momentum is

    H(p)=alpha p+beta R³p⁴,
    E(p)=sqrt(p²+R⁶p⁸),
    P_s(p)=[I+sH(p)/E(p)]/2,  s=±1.

Each band has multiplicity2. All calculations below retain both bands.

## Weighted-energy boost: an exact failed algebra transfer

On a common smooth core in momentum representation, let x=i partial_p and
K={x,H}/2. Then

    [H,K]=−i{H,partial_p H}/2
         =−i partial_p(H²)/2
         =−i[p+4R⁶p⁷]I.

The physical translation generator remains p. The relative boost-algebra
defect is4(Rp)⁶, which vanishes in the infrared but equals4 at Rp=1. The
source therefore does not have the Lorentz boost algebra at the momenta
that control its finite ultraviolet area coefficient. No nonlinear momentum
relabeling repairs the SAME source/position/clock map.

The two-band projector nevertheless has the formal property [K,P_minus]=0,
because {H,partial_p P_minus}=0. Vacuum preservation alone therefore does
not establish a geometric Lorentz boost. No globally self-adjoint boost
generator or complete boost flow is asserted here; the displayed identity
is a core commutator, sufficient to show failure of the desired algebra.

## A physical tilted field restriction retains extra real branches

Sample the actual evolving field on the plane t=u x, with0<u<1. Its spatial
Fourier variable is

    kappa=k−u s E(k).

The original Fourier Hilbert space uses real k. Each real root contributes
the positive CAR Gram weight

    P_s(k)/|1−u s E'(k)|.

The occupied covariance before Gram normalization contains the same weights
only for s=−1. Squaring the equation gives the degree8 polynomial

    u²R⁶k⁸+(u²−1)k²+2 kappa k−kappa²=0.

`run.py` retains all its roots, checks the unsquared equation for each real
root, and marks complex roots as outside the real-momentum source Hilbert
space. No positive-norm source branch is discarded.

As kappa tends to zero from above, besides the two infrared roots there are
two additional real ultraviolet roots

    k=±R^−1(u^−2−1)^(1/6).

At either ultraviolet root the Jacobian magnitude is3(1−u²). Their
projectors sum to I+u alpha. The two infrared roots have combined Gram

    G_IR=(I+u alpha)/(1−u²),

which is the local Dirac value. The complete Gram is instead

    G_total=(4/3)G_IR.

After local Dirac normalization, the limit is4I/3, not I. The discrepancy
persists as u tends to zero: the ultraviolet roots move to infinity but
their positive Gram weight does not vanish. This is a nonuniform ultraviolet
limit, not a numerical zero-mode effect.

If only the four sampled field components are Gram-whitened, their occupied
covariance tends to

    C_sample=I/2−3 alpha/8−beta/8,
    eigenvalues=1/2±sqrt(10)/8 (each twice).

It is not a projection. The four sampled components therefore do not form a
complete pure-vacuum Cauchy field algebra. They combine eight positive-norm
real-root band modes and lose branch information. A local invertible spin
change followed by correct Gram normalization cannot change these covariance
eigenvalues.

One can retain a complete branch register instead: decompose the original
real-momentum Hilbert space into monotonic branches of kappa(k), with the
usual square-root Jacobian change of variables. All weights are positive;
the band-occupation covariance on that enlarged direct-integral space is a
projection. This restores complete records, but it has variable branch
multiplicity and no demonstrated interpretation as a local four-spinor
field on the observer's geometric slice. It is not a derivation of the
Einstein observer family.

## What the infrared limit does and does not supply

The dispersion and projector have expansions

    E/|p|=1+(Rp)⁶/2+O((Rp)^12),
    P_minus=(I−alpha sign(p))/2−beta(R|p|)³/2+O((Rp)⁶).

The corresponding flat preferred-time action contains the explicitly
declared dimension7 spatial operator

    −R³ psi_bar(−spatial_Laplacian)² psi.

Writing a generally coordinate-covariant completion requires a preferred
timelike normal/foliation in addition to the metric, so that the spatial
Dirac operator is specified. Its normal-dependent stress and dynamics are
additional obligations. This packet does not derive their gravitational
coefficients or show that they disappear from an induced low-energy gravity
law. Discarding the ultraviolet branches to define a restricted infrared
service would also require remeasuring entropy/area and records; the old
eta cannot simply be carried over.

These are failures of a full temporal/observer transfer for THIS declared
regulator. They do not prove that OPH has no other spacetime source. Spatial
finite entropy, the responding coframe and the measured infrared modular
convergence remain valid at their stated scope. A selected lapse, attraction
of resting massive probes, all-observer stationarity and physical G remain
unresolved.
