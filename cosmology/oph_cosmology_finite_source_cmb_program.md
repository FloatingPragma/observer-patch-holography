**Paper release:** `r1549`
**Released:** July 18, 2026

# Claim Boundary

A physical cosmological prediction requires one construction that maps finite OPH records to a covariant source, propagates that source through a declared background, and evaluates the resulting observables with a complete likelihood. No such construction is supplied. Simulator outputs may test numerical interfaces and expose incompatible assumptions. They carry no cosmological verdict.

The source ledger is shared with the dark-gravity, cosmological-vacuum, and compact OPH papers.

# Finite-Source Interface

Let $`X_r`$ denote a finite observer-screen state at refinement level $`r`$. A cosmological source map must emit
``` math
\begin{equation}
\mathsf S_r(X_r)
=\bigl(g_{\mu\nu},T_{\mu\nu},\mathcal I,\mathcal B\bigr)_r,
\label{eq:source-map}
\end{equation}
```
where $`T_{\mu\nu}`$ is conserved, $`\mathcal I`$ contains gauge-invariant initial data, and $`\mathcal B`$ states the background and boundary conditions. The map must be defined without using the observed sky quantities that its outputs are compared with.

The complete one-shell covariance determines the three-dimensional correlation only for chord lengths $`0\le s\le2R_\star`$. The resulting radial kernel is infinite-dimensional and persists under positivity. A source-derived lift therefore requires either a scale-natural embedding that transports refinement to physical log-wavenumber dilation or complete radial cross-covariances that support spherical-Hankel tomography. The physical source map must also specify the radial measure, null-space treatment, clock, normalization, and non-fitting forward residual. Construction of one finite source DAG that passes this packet is work in progress.

## Conditional source-spectrum packet

On a source-stress-complete uniform-density release cut, the screen field is
``` math
\begin{equation}
q_r=\Pi_{\ge2,r}\left[\frac13\log\frac{J_{X,r}}{\bar J_{X,r}}\right].
\end{equation}
```
The primitive collar ledger selects a relative-MaxEnt release law before the screen energy is evaluated. It excludes spectra, residuals, likelihoods, fitted amplitudes, and measurement-calibrated proxies. The release energy then fixes
``` math
\begin{equation}
A_{q,r}=\frac1{d_r}\mathbb E_{\nu_r^{\rm rel}}
[q_r^{\mathsf T}M_rK_{\theta,r}q_r]
=\frac{2E_{q,r}^{\rm src}}{d_r}.
\end{equation}
```
The conformal precision and angular spectrum are
``` math
\begin{align}
K_\theta
&=\frac{\Gamma(\sqrt{-\Delta_{S^2}+1/4}+3/2+\theta/2)}
{\Gamma(\sqrt{-\Delta_{S^2}+1/4}-1/2-\theta/2)},\\
C_\ell^q
&=A_q\frac{\Gamma(\ell-\theta/2)}
{\Gamma(\ell+2+\theta/2)},\qquad \ell\ge2.
\end{align}
```
The infinitesimal reserve-generator receipt supplies a full-collar density $`P_\star/24`$. Orientation reversal gives its source-facing half, so
``` math
\begin{equation}
\theta=\frac{P_\star}{48},\qquad
n_s=1-\frac{P_\star}{48},\qquad
\kappa_{\rm rep}^{\rm edge}=\frac{P_\star}{48(P_\star-\varphi)}.
\end{equation}
```
A finite one-step survival probability has exponent $`-\log u_q(\log b)/\log b`$ and does not supply this infinitesimal density.

On the source-dilation branch, the finite source embedding must satisfy
``` math
\begin{equation}
D_{s,r}U_r=U_rR_{s,r},\qquad
C_{\zeta,r}=U_rC_{{\rm src},r}U_r^*.
\end{equation}
```
Strong finite-to-continuum convergence, a uniform covariance norm bound, and a vanishing operator residual then give
``` math
\begin{equation}
D_s^{-1}C_\zeta D_s=e^{-\theta s}C_\zeta.
\end{equation}
```
This identity forces
``` math
\begin{equation}
\Delta_\zeta^2(k)=A_\zeta(k/k_\star)^{-\theta},\qquad
A_\zeta=\frac{A_q}{\pi^{3/2}Z_q^2(k_\star R_\star)^\theta}
\frac{\Gamma(3/2+\theta/2)}{\Gamma(1+\theta/2)}.
\end{equation}
```
The tomography branch obtains uniqueness from complete radial cross-covariances. A finite prior selects a conditional radial continuation. Every finite branch publishes its singular spectrum, null basis, finite-window error, positivity check, and forward residual. Temperature and polarization spectra belong to the transfer interface.

# Dark-Sector Interface

The dark-matter paper proposes a repair-charge rotor with occupation $`n`$ and phase $`\theta`$. Conditional on that action, a source-free dilute branch has
``` math
\begin{equation}
n\propto a^{-3},
\qquad
\rho_R\propto a^{-3},
\qquad
w_R\simeq0.
\label{eq:dust-branch}
\end{equation}
```
Its cubic condensed branch gives the deep radial-acceleration law and the baryonic Tully–Fisher scaling under static spherical premises.

Equation <a href="#eq:dust-branch" data-reference-type="eqref" data-reference="eq:dust-branch">[eq:dust-branch]</a> does not determine the abundance, initial perturbations, relativistic stress, gravitational slip, sound speed, nonlinear evolution, or coupling to the visible sector. A cosmological use requires a covariant completion that emits these quantities from the same action and obeys the full energy-momentum ledger. That completion is work in progress.

# Primordial Interface

A physical primordial branch must supply all of the following:

1.  a finite source ensemble and a clock;

2.  a conserved covariant stress tensor;

3.  gauge-invariant scalar, vector, and tensor initial data;

4.  a screen-to-volume lift with a controlled null space;

5.  phase, amplitude, and isocurvature conditions;

6.  refinement stability and an error model.

The conditional source theorem selects $`n_s=1-P_\star/48`$ from the infinitesimal edge-center generator receipt. Physical primordial status requires one finite source DAG that supplies this receipt together with the release law, stress, clock, freezeout, mode, scale, and radial receipts. Numerical agreement with a measured summary does not supply the missing source map.

# Transfer Interface

Given a completed source map, the transfer calculation must declare the background species, equations of state, collision terms, recombination model, initial modes, gauge convention, numerical tolerances, and observable projection. Standard Einstein–Boltzmann software can test this interface with imported inputs. Such a run tests the solver path rather than an OPH source theorem.

The minimum transfer outputs are
``` math
\begin{equation}
\bigl(C_\ell^{TT},C_\ell^{TE},C_\ell^{EE},C_L^{\phi\phi},
P_m(k,z),H(z),D_A(z),f\sigma_8(z)\bigr).
\end{equation}
```
All outputs must arise from one compatible background and perturbation system.

# Likelihood Interface

A data comparison must state the datasets, masks, covariances, nuisance parameters, priors, and combination rule. Overlapping samples and shared calibrations must be represented in the covariance model. Visual overlays, diagonal residual sums, and comparisons made after selecting a formula are diagnostics.

The current cosmology materials contain diagnostic CMB, galaxy, and background comparisons. They do not define cosmological falsifiers or forward targets.

# Status

<div class="center">

| Object | Status |
|:---|:---|
| Finite observer-screen carriers | Conditional finite construction |
| Screen scalar covariance | Phase III conditional theorem |
| Source-functional screen amplitude and edge-center tilt | Phase III conditional theorem; finite receipt work in progress |
| Radial-lift mathematics | Phase III conditional theorem packet; one-shell unrestricted inversion is impossible |
| Physical source-dilation or radial-tomography receipt | Work in progress |
| Primordial covariant source | Work in progress |
| Repair-charge cosmological abundance | Work in progress |
| Relativistic repair stress and slip | Work in progress |
| Boltzmann initial-state bridge | Work in progress |
| Joint CMB, BAO, lensing, growth, and cluster likelihood | Work in progress |
| Cosmological falsification target | Excluded from the program |

</div>

# Conclusion

Finite OPH screen data become cosmology only through a covariant source map, a three-dimensional lift, a conserved stress tensor, compatible transfer equations, and a complete likelihood. These constructions are work in progress. Existing cosmological calculations are diagnostics and carry no falsification verdict.
