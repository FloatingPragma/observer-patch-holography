# Regional modular energy of the routed scalar

## What the compiler supplies

After heralded preparation, one stationary marker and one mobile excitation occupy the protected register. The data word acts on the mobile amplitudes by \(U=S(2ee^*-I)\), where \(e=(1,\ldots,1)/\sqrt8\) and \(S\) is the flip-flop shift. At a fixed period the pure state on any station cut has Schmidt decomposition
\[
 \sqrt{p}\,|\psi_R\rangle|0\rangle
 +\sqrt{1-p}\,|0\rangle|\psi_{\bar R}\rangle .
\]
The marker is fixed on one side. The formula holds for a union of routed slot sets at restored period boundaries. The entropy of that state on \(R\) is exactly the binary entropy \(H_2(p)\) in nats. Resolving or averaging classical herald and clock outcomes specifies different instruments. [Compiler](../../../../model/experiments/seam_program/compiler.py) · [Preparation and return of work rails](../../../../model/experiments/seam_program/PROOF.md)

Use the three-dimensional \(2^3\) periodic chart with private size one, and prepare an equal-slot wave of momentum \((\pi,0,0)\) or \((0,\pi,0)\) in lattice units. Both have scalar eigenphase \(\pi/3\). In each case \(a=e^*Se=1/2\) and \(U^2e=2aSe-e=Se-e\). The resulting normalized coin lives entirely on the two slots of the chosen axis. The two resulting coin rays on one station are orthogonal. The mobile occupation probability is \(p=1/8\) in each case. Mixing the two inputs with weights \(1-u,u\) gives station entropy \(H_2(p)+pH_2(u)\). The right derivative at zero is infinite because a previously empty regional state gains weight \(pu\). Under the registered shift-only program \(S^2=I\), both local coin rays equal \(e\) and the same mixing leaves the station entropy fixed. The test compiles the pulses on the native data carrier and checks those ray overlaps to \(10^{-13}\). This identifies what is missing from a source transfer based solely on a match of one-particle eigenphases. [Native calculation](compiled.py)

The frequency-based Gaussian construction assigns an oscillator to every nonzero Fourier mode, including unrestricted oscillator occupations, and makes the state their joint vacuum. Native one-mobile detection does not implement these field or momentum quadratures. The calculations below use the canonical completion explicitly specified by the source. [Vacuum prescription](../../../../model/experiments/seam_program/scalar_vacuum.py) · [Scope of the measurement](../../../../model/experiments/seam_program/CLOCKED_FIELD.md)

## Finite regional algebra with the spatial constant removed

Let the periodic side be \(L\), the mesh \(h\), and the number of sites \(N=(L/h)^3\). On the scalar branch with private angle zero, \(\omega_h(k)=E(hk,0)/h\) with low-frequency speed \(c=1/2\). Site coordinates \(q_x=h^{3/2}\phi_x,\ p_x=h^{3/2}\pi_x\) have vacuum covariances given by the discrete inverse Fourier transforms of \(1/(2\omega_h)\) and \(\omega_h/2\), respectively. Both Fourier arrays have zero in the omitted mode. The symmetrized \(q,p\) covariance vanishes. The commutator is \(i(\delta_{xy}-1/N)\). A proper \(n\)-site region has the invertible matrix \(D=I-\mathbf1\mathbf1^T/N\); its canonical momentum is \(D^{-1}p_R\). The raw site generators in a region and its set-theoretic complement need not commute after the global mode is deleted. The finite read is the restriction to the indicated canonical subalgebra.

Let \(X\) be the regional position covariance and \(P\) the covariance of the regional *canonical* momentum. If \(X^{1/2}PX^{1/2}=O\,\mathrm{diag}(\nu_j^2)O^T\), the independent oscillators in coordinates \(O^TX^{-1/2}q,\ O^TX^{1/2}D^{-1}p\) have variances \(1,\nu_j^2\). Their modular frequencies are \(\epsilon_j=\log[(\nu_j+1/2)/(\nu_j-1/2)]\). A coherent mean shift whose coordinates in this basis are \(a_j,b_j\) therefore gives relative entropy, equal to its modular-energy increase,
\[
 C_R(a,b)=\tfrac12\sum_j \epsilon_j(\nu_j a_j^2+b_j^2/\nu_j).
\]
The covariance and entropy of a displaced finite Gaussian state stay fixed. The calculation rejects unresolved gaps between \(\nu_j\) and \(1/2\). A constant *local* mean can be continued outside the region with total spatial average zero. No global zero-frequency oscillator has been restored.

## What converges to the conformal ball weight

Put \(w_R(x)=(R^2-|x|^2)/(2R)\). Write the coherent classical means in chart units as \(\phi,\pi=\partial_t\phi\), with \(T_{\rm can}=(\pi^2+c^2|\nabla\phi|^2)/2\). Under relativistic time \(x^0=ct\) the normalized field is \(\Phi=\sqrt c\,\phi\), so the relativistic energy density is \(T_{\rm can}/c\). The free massless continuum scalar on its double cone has the coherent relative entropy
\[
 C_R=\frac{2\pi}{c}\int_{B_R}w_RT_{\rm can}\,d^3x
            +\frac{\pi c}{R}\int_{B_R}\phi^2\,d^3x .
\]
This is the three-dimensional spatial case of Corollaries 5.5 and 5.6 of Longo and Morsella. For compact smooth waves the integrals are finite; smooth rapidly decreasing data give the same local form by approximating the data on a neighbourhood of the ball. [The massless modular Hamiltonian](https://arxiv.org/pdf/2012.00565)

For momentum-only data \(\phi=0\) and a smooth \(\pi\), this gives \(C_R=8\pi^2 R^4 T_{\rm can}(0)/(15c)+o(R^4)\). With \(\phi=1,\pi=0\) on the ball, both the canonical and improved local energy densities vanish there, whereas \(C_R=4\pi^2 cR^2/3\). In relativistic units \(T^{\rm imp}_{00}=T^{\rm can}_{00}-(1/6)\Delta\Phi^2\). Integration by parts, using \(\partial_n w_R=-1\) and \(\Delta w_R=-3/R\), writes the same modular expression as \(2\pi\int_B w_R T^{\rm imp}_{00}+(\pi/3)\int_{\partial B}\Phi^2\). The edge term depends on the scalar algebra. On this algebra a premise with only smooth bulk stress omits a term of lower order in \(R\). Changing the algebra or transferring this term to an area response requires an explicit readout. [Boundary term for the conformally coupled scalar](https://arxiv.org/abs/1610.02261)

**Refinement theorem for finite observation menus.** At fixed physical radius, choose finitely many linearly independent real functions in \(C_c^\infty(B_R)\); smear both the canonical field and momentum with these functions. Sample the smearings on each lattice and discard its global zero mode as above. For a smooth coherent source sampled with the same limiting local means, each menu's relative entropy converges to that of the corresponding continuum menu as \(h\to0,\ L\to\infty\).

To see this, orthonormalize the fixed continuum test functions in spatial \(L^2\). The sampled commutator Gram matrix, including the zero-mode subtraction, tends to the identity. Its inverse therefore supplies canonical momenta on every sufficiently fine lattice. For test functions \(f,g\), the covariance sums are Fourier Riemann sums of \(\overline{\hat f_h}\hat g_h/(2\omega_h)\) and \(\overline{\hat f_h}\hat g_h\,\omega_h/2\). On the Brillouin zone \(|k|/\pi\le\omega_h(k)\le C|k|\), by the sine and arcsine bounds for the compiled eigenphase. The Fourier transform of a sampled compact smooth function is a sum over translates \(2\pi n/h\) of its rapidly decreasing transform, hence decays faster than any fixed inverse power uniformly over the zone for bounded \(h\). Away from zero the frequencies converge locally uniformly to \(c|k|\). The normalized low-frequency sum of \(1/|k|\) below \(\delta\) is \(O(\delta^2+L^{-2})\) by counting lattice shells. Splitting off small and large momenta proves covariance convergence. Smooth means converge by direct quadrature. The free continuum vacuum is separating for a ball algebra, hence faithful on its finite oscillator subalgebras. Each fixed menu has minimum symplectic eigenvalue strictly greater than \(1/2\). Its finite-dimensional Gaussian formula is consequently continuous in these covariances and means.

Exhaust smooth compact smearings by increasing finite menus dense in the test-function topology. Their Weyl operators generate the continuum double-cone algebra by the time-zero construction. Relative entropy on an increasing weakly generating family of subalgebras tends to its value on the full algebra. For each fixed menu, monotonicity on the finite lattice and the convergence just proved give
\[
 \liminf_{h\to0,\ L\to\infty} C^{\rm full\ lattice\ ball}_{R,h,L}
 \ \ge\ C^{\rm continuum}_R .
\]
Choosing a sufficiently fine regulator for each successively larger menu gives a diagonal sequence of *those menu readings* converging to \(C_R\). One regulator per menu can serve the first \(j\) members of a prescribed countable test family to accuracy \(1/j\), by applying finite-menu convergence to each member. This is an existence statement for increasing coarse observations in the Gaussian completion. A uniform upper bound on all unresolved modes of the complete ball, or an effective physical mesh-rate for the signed gravity residual, requires further information. [Bostelmann, Cadamuro and Del Vecchio, Lemma A.1](https://arxiv.org/pdf/2012.14401)

**A two-sided source derivative.** For a finite Gaussian regional read let \(K_R\) be its modular generator and let \(u\) be one normalized global oscillator wave. Denote the coherent modular increase at amplitude 1 in that wave by \(C[u]\). The curve \((|0\rangle+s|2_u\rangle)/\sqrt{1+s^2}\) is a state for both signs of \(s\). Expanding the coherent state at amplitudes \(a\) and \(ia\) to order \(a^2\) gives its modular derivative exactly as \((C[u]-C[iu])/\sqrt2\). In fact the diagonal one-particle terms cancel, the off-diagonal vacuum/two-particle terms add, and coherent costs scale quadratically. A normalized real even Fourier wave supported in a smooth annulus can be chosen with position mean nonzero at the centre and zero momentum at the initial time. Its phase rotation by \(i\) has zero position mean. Both waves have smooth bounded energy densities. The *same* diagonal sequence of finite readouts for these two phases thus gives a signed modular derivative with leading term \(4\pi^2c\,\phi(0)^2R^2/(3\sqrt2)\) in the continuum small-ball expansion. The bounded stress contributions start at \(R^4\). This argument uses polarization on finite registers, without defining a trace-class continuum regional density. The two lower bounds for the unrestricted complete lattice balls cannot be subtracted to infer a signed limit there.

## A single quadrature detects the extra term

On the field algebra measure the average \(\bar\phi\) of the positions in a lattice ball. A local mean shift by one leaves its Gaussian variance \(V_h\) fixed. Its classical relative entropy against the vacuum reading is \(1/(2V_h)\), a lower bound on the full regional relative entropy. If \(W_h(k)\) is the Fourier transform of the site weights summing to one, then
\[
 V_h=L^{-3}\sum_{k\ne0}\frac{|W_h(k)|^2}{2\omega_h(k)} .
\]
This read has a direct continuum limit even for the sharp indicator of a ball. On bounded momentum annuli both the normalized ball indicator transform and the dispersion converge uniformly by Riemann summation. At low momentum \(|W_h|\le1\) gives the \(O(\delta^2+L^{-2})\) estimate above. At \(|k|>K\), Parseval and the frequency lower bound bound the tail by \(\pi/(2K h^3 n)\), with \(n\) the number of retained sites. For fixed \(R\), \(h^3n\to |B_R|\). This proves joint decompactification and refinement, and by rescaling proves the analogous dimensionless asymptotics with \(h/R\to0,\ L/R\to\infty\).

At equal time the continuum field covariance is \(1/(4\pi^2 c |x-y|^2)\). The overlap volume of two radius-\(R\) balls at distance \(d\le2R\), divided by the ball volume, is \(1-3d/(4R)+d^3/(16R^3)\), by the elementary spherical-cap formula. Integrating against \(4\pi\,{\rm d}d/|B_R|\) gives the average of \(|x-y|^{-2}\), \(9/(4R^2)\). Hence
\[
 V_h\sim\frac{9}{16\pi^2cR^2},\qquad
 D_{\bar\phi,h}\sim\frac{8\pi^2c}{9}R^2 .
\]
One scalar observable recovers two thirds of the full continuum cost for the flat field shift. The omitted Fourier zero mode does not remove locally constant means or this lower bound.

For \(R=1\), \(h=1/8\), \(L=8,12,16\), the measured variance divided by \(9/(16\pi^2c)\) is \(0.95450,0.98831,1.00025\). For the full finite canonical restriction at fixed 48 sites per side the computed ratios to the analytic field and momentum moments are as follows. These small balls show the mesh-shape dependence and do not certify a uniform error rate.

| \(R/h\) | Retained sites | Flat field / continuum | Flat momentum / continuum |
| ---: | ---: | ---: | ---: |
| 1 | 7 | 1.4090 | 2.3779 |
| 1.5 | 19 | 1.2048 | 1.5614 |
| 2 | 33 | 1.0053 | 1.0525 |
| 2.5 | 81 | 1.1700 | 1.3928 |
| 3 | 123 | 1.0899 | 1.1332 |

These results establish the momentum-only normalization for increasing coarse readings of the supplied vacuum. The extra field observable has to be addressed when a gravity premise admits general sources. Native implementation of the vacuum algebra and its measurements, proper-ball geometry, and the signed area response are separate inputs. [Numerical method](measure.py) · [Consistency tests](test_kernel.py)
