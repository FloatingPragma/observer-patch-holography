# A finite charged energy cluster

The declared Wilson cube has an isolated pair of charged energy levels at
coupling \(1/64\), photon capacity three and completion cost \(21/2\).
Exact matrix identities and rational bounds give a gap above the pair
greater than \(0.17175185m\), and a reference-subtracted charged transition
between \(1.11134603m\) and \(1.12278560m\). These are finite source energies;
their interpretation as particle masses requires a physical clock and a
particle correspondence. The [proof](PROOF.md) states the construction and
the role of the reference energy.

Run `python -m model.native_charged_cluster_minmax.audit` to replay the
exact result. Run `python -m pytest -q model/native_charged_cluster_minmax`
to check the retained result and arithmetic rejection controls.

The [self-energy calculation](self_energy/PROOF.md) preserves exact
coefficient arithmetic, a residual-controlled one-body eigenframe and
virtual-photon contractions. Its proposed many-body transfer allowance is
unproved; the binary64 second-order energies are diagnostics. The modules
`self_energy.diagnostic` and `self_energy.contractions` run from the
repository root without a separate research checkout.
