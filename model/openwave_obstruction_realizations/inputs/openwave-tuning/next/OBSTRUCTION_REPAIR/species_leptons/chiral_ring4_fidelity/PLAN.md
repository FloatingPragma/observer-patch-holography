# Complete heat-trace eigen-count test

The source and every preparation are the frozen L4,S2,16-CAR member.
Use the exact split `H=(E+qbar)^2+H_matter+B`. Golden–Thompson and the
shifted-rotor diagonal heat kernel give

`Tr exp(−beta H) <= Z_rotor(beta) average_theta Tr_Fock exp[−beta(H_matter+B)]`.

Retain all12870 N8 occupations at every one of1024 holonomies. Positivity
of B bounds every excited pointwise level by its complete free Fock
counterpart; only the ground contribution uses the frozen Coulomb
two-block lower bound. Degenerate and high-energy levels remain in the
trace. The global holonomy derivative bound controls the quadrature.

Predeclare beta in {1/4,1/2,3/4,1,5/4,3/2,2,3}, and cut offsets in
{1/16,1/8,1/4,1/2,1} above each preparation's Rayleigh upper bound.
These are verification filters, not source parameters. A count upper
bound below2, together with a variational state below the cut, proves
exactly one neutral eigenvalue there. If this gate passes, use the full
stored infinite-flux spectral residual for Temple/fidelity bounds and
the bounded charged currents to transfer prepared responses to the
actual vacuum. If it fails, retain the failure as insufficiency of this
counting bound. It is not evidence for several vacua.

Transfer the electron-input case only with the declared4|m| operator
bound. Retain every beta/offset test, every Fock occupation and the full
per-grid band and trace data. A finite-volume vacuum result does not
certify a physical charged-family mass or continuum chiral protection.
