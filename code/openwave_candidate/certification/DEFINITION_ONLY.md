# Twelve-direction scalar propagation: independent computation packet

This packet supplies the definition, requested outputs and verification
rules. It contains no expected expansion coefficients. Mathematical
recomputation is distinct from empirical prediction or OpenWave certification.

Let phi=(1+sqrt(5))/2 and take the twelve unit vectors obtained by normalizing
the vectors (0,+/-1,+/-phi), (+/-1,+/-phi,0), (+/-phi,0,+/-1), with independent
signs. For a positive length a, a unit direction n and scalar momentum k, define

    omega_squared(k,n,a) = sum_u [1-cos(a*k*(u dot n))]/(2*a*a).

Use every direction once and equal weights. Define P6 as the degree-six
Legendre polynomial. Define I6(n) as sum_u P6(u dot n), divided by its value
when n is any one of the twelve unit vectors. These are declared mathematical
inputs; no identification with a physical field or physical lattice scale
is supplied.

Compute from these definitions:

1. The exact homogeneous quadratic, quartic and sextic polynomials.
2. Their rotationally invariant parts and the first nonconstant angular part.
3. The coefficients in the stated I6 normalization and their scale-free ratios.
4. The change under a -> a/2.
5. An independent direct-symbol numerical check over several directions and
   decreasing k, with a bound on the Taylor remainder rather than cancellation-
   sensitive subtraction alone.

Before reading any target values or OPH result files, save the implementation,
outputs and a statement of which files were read. Check identities by exact
algebra where practical. For the direct symbol use at least 60 decimal digits.
Deliberately remove one port or change one weight and demonstrate that the
isotropy test rejects the altered operator. Do not choose weights or
normalizations by comparison with an experimental target.

Acceptance is internal mathematical consistency of these outputs, exact
normalization, the direct-symbol remainder check, and successful detection of
the broken control. Agreement with another implementation is evaluated only
after both outputs are saved. A later external referee controls its own
certification and empirical-data gates.
