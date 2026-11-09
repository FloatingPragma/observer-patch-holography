# A complete three-input charged-lepton controller

The [complete attenuator](BARE_EFFECTIVE_CONTROLLER.md) can reproduce the
three quoted charged-lepton central masses as conditional action coefficients.
All three masses are calibration inputs. This construction supplies no
held-out mass prediction, source selection or interacting mass matching.

Let the input energies be `m_e < m_mu < m_tau`, in GeV, and put
`x=(sqrt(m_tau),sqrt(m_e),sqrt(m_mu))`. Define their discrete Fourier
coefficients

    A=(x0+x1+x2)/3,
    B=(x0-(x1+x2)/2)/3 + i sqrt(3)(x2-x1)/6.

The inverse Fourier identity gives
`x_j=A+2 Re[B exp(2 pi i j/3)]`. For the existing native constants
`a=sqrt(10/72)`, `b=sqrt(5/72)=a/sqrt(2)`, choose

    r=sqrt(2)|B|/A,  z=B/|B|,  s_GeV=3 A^2/a^2.

Then `C_r=aI+rb(zR+conjugate(z)R^2)` has eigenvalues `(a/A)x_j`.
Consequently the existing family action `Y_r=C_r^2/3` satisfies exactly

    s_GeV * spec(Y_r) = {m_e,m_mu,m_tau}.

Parseval gives `r^2=3Q-1`, with `Q=sum(m)/(sum(sqrt(m)))^2`.
Thus this attenuation-only apparatus admits a nondegenerate target when
its derived amplitude lies in `[0,1]`. This is a substantive realizability
condition; a target with `Q>2/3` cannot be fitted by this attenuator while
retaining the stated positive roots. Other source actions are outside that
exclusion.

For the supplied CODATA2022 electron and muon central values and PDG2026
tau central value, the exact radical enclosures certify

    r       ≈ 0.9999966950984938674,
    theta   ≈ 0.222224761894592598 radians,
    s_GeV   ≈ 6.779180152508157998 GeV,
    Q       ≈ 0.666664463402636703.

The resulting conditional action energies are
`0.00051099895069`, `0.1056583755`, `1.77693` GeV. The overall factor
`s_GeV` combines a supplied action normalization and energy calibration;
it is not an identification of the physical Higgs vacuum expectation value.
The original sixteen outcomes, including loss flags, are retained. There
is no postselection or division by a successful event probability.

The code supplies rational approximations to `r`, the phase half-angle
tangent and `s_GeV` at forty decimal places. Its phase lies exactly on
the unit circle. Outward rational radical enclosures certify every output
within `10^-36 GeV` of its supplied central value. That is a numerical
realization bound, not experimental or physical theory accuracy. The
complete-controller theorem applies; a native tape at this particular
fitted point is not replayed here.

This test differs from the [one-ratio Koide comparison](../charged_lepton_response/CALIBRATED_SANITY.md):
there the muon/electron ratio selects the phase, exact balance is retained,
and tau is held out. Here tau selects attenuation and also changes phase
and scale. The reference mass central values and the independently published
mass-ratio entries need not have identical rounded quotients. No joint
confidence region or covariance fit is claimed.

Connecting this action to measured electrons, muons and taus requires the
same interacting chiral state and physical clock, canonical kinetic terms,
and infrared-threshold or resonance matching. A source derivation must
calculate the controller or the effective coefficient without reading these
three target masses. Measuring complete same-state correlators and deriving
their kinetic and self-energy kernels is one concrete route; allowing every
symmetry-compatible matrix correction alone does not select a value.
