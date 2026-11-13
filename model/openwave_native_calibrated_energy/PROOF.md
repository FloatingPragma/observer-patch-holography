# Same-source spectator interaction and calibrated transition energies

Use the complete corrected C/E source and preparation/readout domain of
[the coherent clock](../openwave_native_coherent_clock/PROOF.md). Its exact
compression can be written

    h = I_n tensor E0(s) + n tensor Delta(s),
    Delta = (1/256,-1/1024,-1/3072).

The difference `Delta(1)-Delta(0)=-5/1024` is nonzero. This is a genuine
interaction between the binary source coordinate and its spectator. It is
calculated from the same H2, matching code and native recovery, without
installing a second interaction Hamiltonian. H2's two-hop phase remains a
declared architectural service.

Prepare binary X-plus and a coherent superposition of spectator values zero
and one, with the matching syndrome included. The ideal evolved state has
binary-versus-spectator negativity

    N_ideal = |sin(t (Delta(1)-Delta(0))/2)|/2.

At t=256 the sine argument is 5/8. Let epsilon be the inherited complete
source/reference trace-norm error at 65,536 intervals. Partial transposition
on a two-dimensional factor has induced trace norm at most two, hence
`|N(rho)-N(sigma)|<=||rho-sigma||_1`. The actual complete finite source has
negativity at least `N_ideal-epsilon`. Directed rational sine bounds prove
this lower bound positive for C and E. Every spectator and syndrome row is
included. A zero-phase control has no entanglement. The stored floating
partial-transpose spectrum is a diagnostic; the positive lower bound is
analytic and applies to the exact declared source.

## One exposed energy scale

The smallest absolute transition gap is 1/3072. Equating this gap to the supplied
electron anchor `0.00051099895069 GeV` fixes the conversion
`scale = anchor/(1/3072)` GeV per source unit, giving
absolute transition energies in the ratio `1:3:12`. Exact sine-record
intervals from the clock producer propagate through that same conversion.
They are calibrated source transition energies. They are neither stable
rest poles nor three charged leptons. The literal OpenWave lepton-spectrum
comparison `1:207:3477` fails for this assignment. The electron and lepton
hierarchy were exposed before the comparison, so it is retrospective.

The weighted transfer-ground protocol in
[`native_transfer_control`](../native_transfer_control/PROOF.md) provides a
separate positive source-derived logarithmic energy and complete native
readout, with interval-checked coherent implementation. Its weighted source,
ground and time assignment do not identify it with the corrected phase
program used here. Neither construction proves that no coherent energy
bridge exists. A common spatial particle, active charge, energy clock and
readout requires an explicit interface on the same declared program.
