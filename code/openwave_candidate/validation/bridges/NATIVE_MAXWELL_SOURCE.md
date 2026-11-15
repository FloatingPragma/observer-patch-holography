# Native coherent source into the existing Maxwell executor

This bridge imports `code/electromagnetism/seam_maxwell_continuum.py` from RER and calls its existing `blocks`, `propagator`, and `constant_forcing` functions. It does not supply a replacement Maxwell derivation. Source hashes retain that provider alongside the native construction. The byte-identical upstream module is packaged at `providers/seam_maxwell_continuum.py`, with Apache-2.0 licence text, upstream licensing notice and a provenance/hash record. Both development and public invocations import this same local provider; the upstream receipt-building CLI is not invoked. All executable receipt pins are bundle-root-relative.

The source is the actual zero-commit coherent history of the modular local-coarse triangle in its explicit uniform-Gauss-fiber preparation. The exact one-load intertwiner permits execution on its six reduced port/transit coordinates. The supplied initial amplitudes are `[1,2i,0,1−i,0,0]/sqrt7`, mass .23, modulus6. There is no energy potential or source/field backreaction.

For seam e=(a,b), let j_a and j_b be the Hamiltonian probability currents from its two ports into its transit. They are evaluated from the actual H matrix elements and density. Define the whole-seam current `J_e=(j_a−j_b)/2`. The corresponding observer gives each transit half to each endpoint: `q_a=n_a+sum_e incident(a) n_transit(e)/2`. Exactly `dot q_a=−sum outgoing J+sum incoming J`. This effective-port observer differs explicitly from the whole-transit regional-charge instrument; its fractional transit readout does not establish sharply fractional native charge.

The declared embedding places the three ports on a unit equilateral triangle in a three-dimensional coordinate chart. For Fourier convention `qhat=sum q_a exp(−ik.x_a)`, map each seam current to the straight-edge line current

`Jhat=sum_e J_e (x_b−x_a) exp(−ik.(x_a+x_b)/2) sinc(k.(x_b−x_a)/2)`.

Here sinc(z)=sin(z)/z. The endpoint exponential identity gives exactly `dot qhat+ik.Jhat=0`; a finite-difference test verifies its sign, and omitting transit fractions is a rejecting control. The embedding and its units are supplied, not derived spacetime.

At the declared momentum (.7,.4,.3), the existing transverse Maxwell propagator receives forcing `(−Jhat_T,0)`. It uses supplied continuum dispersion omega=|k|. A zero initial field and zero-current control stays zero. Each interval holds the current at its exact native left endpoint. The longitudinal field starts at `E_L=−ik qhat/|k|²` and advances with `dot E_L=−Jhat_L`, preserving Gauss for the held-source charge. That held charge is compared with the actual native endpoint charge; 20-to-40-step refinement reduces their difference. No exact native-history field solution is claimed by a finite hold discretization.

The receipt retains every interval source and transverse mode energy. The work integral `Re integral F† forcing dt` is evaluated independently by quadrature and compared with the transverse energy change. Longitudinal energy change is compared with the exact affine-update work. This source is externally prescribed, so field energy gain is supplied source work; no conserved joint native/field Hamiltonian has been constructed.

One Fourier mode's norm is an algebraic energy observable. It does not resolve a spatial shell, near/far separation, distant flux, quantum radiation quanta, outgoing scattering, calibrated electromagnetic coupling or physical clock. Those require a further explicitly matched source/field construction. The bridge closes a reusable source-to-existing-propagator connection, without promoting a physical OpenWave criterion.
