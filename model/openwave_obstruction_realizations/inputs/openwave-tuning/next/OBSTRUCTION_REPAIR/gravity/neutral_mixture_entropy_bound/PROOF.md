# Actual non-Gaussian neutral-mixture entropy from a low-rank bound

Let rho0 be the actual regional Gaussian vacuum, rho1 the actual regional Gaussian restriction of the neutral four-particle/four-hole Slater source, and rho_e=(1-e)rho0+e rho1. The source preparation has two explicit global branches, both with exactly the same conserved charge. Its external classical preparation flag remains in the complete source record; it is not part of the specified regional observer algebra. Including that accessible flag would add its Shannon entropy and would be a different observable.

The mixture rho_e is generally not Gaussian. We do not replace it by the Gaussian state with its covariance. If K0=-log rho0, define J=Tr[(rho1-rho0)K0]. The exact relative-entropy identity is

S(rho_e)-S(rho0)=e J-D(rho_e||rho0).

For faithful finite rho0, the relative modular operator and scalar Jensen inequality give

D(rho||rho0) <= log Tr(rho^2 rho0^-1).

Indeed D is the expectation of the logarithm of the inverse relative modular operator in the vector rho^(1/2), and log is concave. Expanding the actual mixture yields Tr(rho_e^2 rho0^-1)=1+e^2 chi2, where chi2=Tr(rho1^2 rho0^-1)-1. Thus

e J-log(1+e^2 chi2) <= S(rho_e)-S(rho0) <= e J.

This controls the actual non-Gaussian finite-e entropy, not only its formal derivative.

For number-conserving Gaussian branch states with faithful covariance C0 and C1=C0+Delta, use rho(C)=det(I-C) Gamma[C(I-C)^-1] and Tr Gamma(A)=det(I+A). Algebra gives

Z=Tr(rho1^2 rho0^-1)
 = det[I-2C1+C1^2 C0^-1]/det(I-C0)
 = det[I+([C0,Delta]+Delta^2)(C0(I-C0))^-1].

The commutator is required when the two covariances do not commute. All unaffected regional modes cancel exactly in this identity. For Delta=X Jsig X†, with r columns and Jsig diagonal signs, let

M=X†X,
G0=X†[C0(I-C0)]^-1 X,
G1=X†(I-C0)^-1 X,
G2=X† C0(I-C0)^-1 X.

The determinant lemma reduces Z to the following 2r-dimensional determinant:

det { I + [[Jsig G1, Jsig G0],
            [-Jsig G2+Jsig M Jsig G1, -Jsig G1+Jsig M Jsig G0]] }.

The neutral four-pair branch has r<=8, so this is at most16 by16. This cancellation does not discard physical bands or approximate the mixture as Gaussian. The matrices G are computed using every actual regional covariance sector. For displaced balls all irreducibles and representation components contribute.

The analytic bound is exact. Floating-point determinants and the numerical approximation of the continuous observer covariance must still be controlled in any quantitative use. Explicit small Fock-space comparisons test the algebra and the finite-e entropy bracket before production use. Geometry at finite amplitude remains an independent obligation; a controlled matter entropy bound alone does not establish joint entropy equilibrium or gravity.
