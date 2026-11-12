# Full spatial domain and controlled trap release

## Declared continuum-coordinate source

Replace the finite oscillator shells of the charged quantum block by
L²(R³;C²), with p=-i grad, r the coordinate and
N_m=(p²+r²-3)/2. Retain the full three-mode photon Fock space and the same
mass16, coupling1/8 and Gaussian profile exp(-r²/5). The Weyl operator is

    B=sigma.p+sum_i M_i(r)b_i,
    M_i(r)=exp(-r²/5) sigma.(e_i cross r)/2,
    b_i=(a_i+a_i*)/8.

In photon coordinates each B(x) is the free self-adjoint Weyl operator plus
a bounded spatial potential. Their direct integral defines a self-adjoint B.
The nonnegative closed form sum

    H_Omega=sqrt(16²+B²)+Omega N_m+N_photon

defines the source for Omega>=0. The vector potential is bounded in the
spatial coordinate and at most linear in the field coordinate. With
N=N_m+N_photon, its bounded Gaussian coefficients give
||B psi||<=C||(N+1)^(1/2)psi|| on Hermite/Schwartz vectors. Their bounded
derivatives also give the Nelson commutator form bound
|<Bpsi,(N+1)psi>-<(N+1)psi,Bpsi>|<=C<psi,(N+1)psi>.
The oscillator reference operator therefore supplies a Hermite core for B;
this closes the core premise in the complete-shell argument below.

For Omega>0, Cauchy--Schwarz and Young's inequality make the square-root
energy form infinitesimally bounded relative to N:
<sqrt(16^2+B^2)><=16||psi||^2+||psi||||Bpsi||
<=epsilon<N>+C_epsilon||psi||^2. The full form norm is thus equivalent to
the oscillator form norm, so finite Hermite/Schwartz vectors are a common
form core. At Omega=0 the closed positive form sum is used directly; the
threshold argument below needs only its Gaussian trial vectors and the
zero-energy kernel conditions. For Omega>0 the lower bound
H_Omega>=16+Omega N_m+N_photon gives compact resolvent. The ordered positive
CAR lift retains all charged fermion sectors, whose N>=2 floor is32. This
argument concerns the charged block; the finite native neutral floor1/2 is
not transferred into a different spatial domain without its own proof.

## Uniform confined spin-half band

For both spin states take spatial and photon oscillator vacua. Their entire
two-dimensional Rayleigh space obeys

    <B²>=3/2+(3/256)(5/7)^(5/2),
    E<=16+v0,
    v0=3/64+(3/8192)(5/7)^(5/2)<1/16.

The equality uses sum_i M_i²=r² exp(-2r²/5)/2 and the exact Gaussian radial
integral. The scalar inequality sqrt(16²+t)<=16+t/32 supplies the energy
bound. Every vector orthogonal to the joint vacuum doublet has
Omega N_m+N_photon>=min(Omega,1). Thus if min(Omega,1)>v0, exactly two
levels lie below16+min(Omega,1). Projection to the vacuum doublet is an
equivariant bijection. Schur irreducibility gives exact J²=3/4 and exactly
degenerate energy, including the full matter and photon angular action.
For every state in this ground band,

    P(N_m>L)<=v0/[Omega(L+1)].

This bound covers arbitrary spin mixtures. It is a sufficient confined-band
domain, not a claim that the same gap persists at all positive Omega.

## Source-fixed complete-shell refinement

Let P_N project onto every spatial oscillator state of total degree at most
N, and use B_N=P_N B P_N with the exact unchanged Gaussian matrix elements.
For a fixed complete photon cap, finite Hermite vectors are a core for B;
B_N psi converges to B psi on that core. The self-adjoint extensions that
are zero outside P_N therefore converge in strong resolvent sense.
Continuity of functional calculus gives the same statement for
sqrt(16²+B_N²). Positivity converts this to Mosco convergence of its forms.
On the common Hermite core the energy forms converge: the uniform second
spectral moment from ||B_N psi||² controls the linear growth of the square
root. The common oscillator form supplies recovery vectors and the liminf
bound after adding Omega N_m+N_photon. Its compact sublevel embedding gives
convergence of each fixed low eigenvalue for Omega>0. The same core and
coercivity argument applies to simultaneous complete photon refinements.

The finite source is extended by mass16 on the spatial complement before
adding the unchanged oscillator operator. The complement therefore moves to
high energy as N grows and cannot hide a low level. No truncated square root
is identified with the compression of the full square root. Finite source
energies need not be monotone. The displayed shell/cap rows are full-mode
diagnostics; convergence alone does not provide a numerical continuum error
bar for the final row.

## Explicit trap-release bound and threshold

For a normalized spatial Gaussian of width L and a photon vacuum,

    E_Omega <=16+3/(64L²)
       +3L²/[8192(1+2L²/5)^(5/2)]
       +(3Omega/4)(L²+L^-2-2).

Every term is an exact expectation bound on the same source. Set
Omega=1/(16L^4) and let L grow. The upper bound tends16, while positivity
gives the lower bound16. At Omega=0 the same Gaussian sequence proves
infimum H_0=16.

No normalized vector attains this threshold. Equality in the nonnegative
sum H_0-16 requires N_photon=0 and B psi=0. Projecting the latter equation
onto the photon vacuum gives sigma.p psi=0. Fourier transformation implies
support only at p=0, which is impossible for a nonzero L² vector. Thus there
is no normalizable self-bound ground at16. Free particles ordinarily also
lack a normalizable zero-momentum eigenstate; this result does not exclude a
free-particle rest-mass threshold or all electron identifications. It excludes
using the confined dressed ground as an established self-localized electron.

The unchanged finite-cap-two electron unit assigns the threshold about
0.000509529040252GeV. This is a conditional conversion of a changed spectral
observable, not a second calibration or physical electron prediction. The
Gaussian field profile has a specified spatial centre, so a translation-
invariant physical dispersion and its source/observer dictionary require
additional work. No neutral, Gauss-tree or magnetic theorem is silently
transferred to this larger spatial source.
