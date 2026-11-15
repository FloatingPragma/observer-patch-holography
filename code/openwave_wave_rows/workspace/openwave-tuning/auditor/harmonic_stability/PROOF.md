# Uniform stability of the carrier field's harmonic response

The face symbols have36edge rows and26face columns per cell. Let A(k) be
their exact Bloch matrix and K(k)=A(k)A(k)*. Reconstruct its finite integer
convolution coefficients K_(ij,d) from the actual face columns. Define the
nonnegative36by36matrix M by M_ij=sum_d |K_(ij,d)|. Then
|K_ij(k)|<=M_ij at every momentum. The same convolution estimate holds
before periodic folding, so it covers every admitted periodic refinement.

Choose the positive integer vector w=M1. Exact arithmetic gives

M w <= (1504/17) w.

For an eigenvector v of K(k), select an index maximizing |v_i|/w_i.
The componentwise triangle inequality yields
|lambda| |v_i|<=sum_j M_ij |v_j|<= (1504/17)|v_i|.
Thus every curvature eigenvalue is at most1504/17. The receipt contains
M, w and every source pin, so this bound needs no sampled momentum grid.
The unweighted maximum row sum158 is a valid weaker estimate.

On a curvature mode lambda, the electric/magnetic symplectic step has
cos(omega)=1-alpha beta lambda/2. For a nonzero mode,
0<alpha beta lambda<4 places both eigenvalues on the unit circle. The
scalar matched coefficients alpha1/8,beta1/4 obey

alpha beta lambda <=47/17<4.

This supplies a uniform harmonic stability certificate for that setting.
Zero-curvature modes retain their separate Gauss, free and global-sector
interpretation; the estimate does not assign them an oscillator gap.

For the proposed Dirac normalization rho=2/3, coefficients
alpha=1/(2N), beta=3/N give c²=N²rho alpha beta=1 and the same canonical
electric and magnetic coefficients1/3and3 for every N. At N=6 the
curvature product is at most188/51<4. At N=8 it is at most141/68<4.
These are different Floquet schedules, despite sharing their leading
speed and coupling. Their accumulated time-splitting errors and native
word/readout costs require their own comparisons.

The pi corner gives an exact integer matrix. A rational trial vector in
the receipt has Rayleigh quotient33989643/500087>200/3. Consequently
the N=5schedule has alpha beta lambda>4 at that momentum. That mode
exists on even periodic refinements, so a uniform N=5stability assertion
fails. Six is the smallest integer N compatible with harmonic stability
throughout this same-speed coefficient family, as certified by these
upper and lower bounds.

Harmonic stability does not bound cosine-to-quadratic quantum-state error.
It does not prepare a photon vacuum or prove an interacting Maxwell limit.
Norm conservation of the finite quantum circuit is a separate property.
