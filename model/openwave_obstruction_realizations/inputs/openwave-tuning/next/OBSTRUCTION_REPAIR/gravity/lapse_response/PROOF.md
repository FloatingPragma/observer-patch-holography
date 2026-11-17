# Actual lapse response on the named two-band regulator

The declared lapse coupling is H_N={N,H}/2 with N=exp(n). It is the linear
coupling of the actual source clock, not a field equation selecting n.
At each momentum, H(k) has the two energies±E(k), each twice.

For a lapse mode q, the perturbation vertex is V=(H(k+q)+H(k))/2. Its
filled/empty divided-difference response simplifies EXACTLY to

    delta P_N(k+q,k)
      =[E(k+q)−E(k)]/[2(E(k+q)+E(k))] [P(k+q)−P(k)].

To derive it, project V between the opposite bands. The two terms are
proportional to(E_q−E)(P_plus,q P_minus−P_minus,q P_plus)/2, and the bracket
equals P_plus,q−P_plus. The negative spectral divided difference then gives
the formula above. Each factor is O(q), so delta P_N starts at q². The
positive ground-energy susceptibility has another square of the interband
vertex and starts at q⁴. Uniform N merely rescales H and leaves P unchanged.

The script computes the regional modular trace of this ACTUAL covariance
derivative, retaining all source momenta. A generic eigenvector
divided-difference check agrees to3.15e-16. Independently, a complete500-mode
Hamiltonian with finite N=exp(epsilon cos(qx)) gives the entropy derivative
to2.67e-10; it does not use the simplified vertex in that calculation.

On the same improved R/h1.5,L60 source as the responding-coframe packet,

| Radius | q1 lapse entropy derivative | q2 derivative |
|---|---:|---:|
|4|−.05461069733|−.20129739171|
|4.5|−.07265219262|−.26283904458|

The uniform response is exactly zero. chi_N/q⁴ is.0388220 at q1 and.0379102
at q2. A nonzero regional q² entropy response must not be confused with a
q² Newtonian ground-energy kernel.

At radius4.5 the independently measured fixed-volume spatial kernel is
K_sigma(q1)=−.2683805164. The leading joint scalar equation would therefore
read K_sigma sigma+K_N n+J=0. This packet supplies its actual lapse component,
but no independent temporal/observer equation selecting n. Setting
n=−sigma, solving a target Einstein equation for it, or identifying this
free clock with the electron clock would be additional input.

The full source's temporal observer transfer fails independently in
../temporal_observer/. A low-energy limit with preferred-frame corrections
is a possible research route; no lapse profile, attractive field or physical
Newton constant has been derived by this computation.
