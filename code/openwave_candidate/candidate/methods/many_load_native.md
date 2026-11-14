# Shared-occupation native extension

`ManyLoadNative` is an explicitly supplied nonnegative Bose extension of the existing one-load generator and local coarse instrument. It reuses the already executed Fable L04 `Occ(bose=True)` half-hop convention, independently implemented here without a Fable runtime import. No occupation or register measurement is inserted during evolution.

A basis state is `(u,F)`, with integer nonnegative occupations on every port and transit, fixed total L, and one Z_N register per half-seam. Start from the supplied occupations and definite zero registers. Breadth-first search retains the complete coherently reachable sector, failing rather than truncating when the declared state budget is exceeded.

For half-seam h connecting port x to transit e, the forward hop is

`|u,F> -> sign(h) sqrt(u_x (u_e+1)) |u−1_x+1_e,F+1_h>`.

The reverse is its adjoint; first/second endpoint signs are +/−. The diagonal is `m sum Gamma_z u_z`, with Gamma +1 on ports and −1 on transits. Thus L=1 recovers exactly the existing one-load H, including all record labels and holonomies. This is a Bose statistics choice, not a derivation of physical statistics. There is no signed vacuum, pair creation, potential, or hard-core interaction.

For drawn seam (a,b) and fair coin c, compute `d=u_a−u_b`, `j=(d−(2c−1)(abs(d) mod2))/2`. Transfer j from a to b and update `F_a+=j,F_b−=j`; transit occupations are unchanged. C retains (coin,j), E retains (coin,d), and MIX additionally draws C/E fairly. Each indistinguishable local outcome is one coherent partial permutation, scaled by the coin/law probability. For fixed j, the final occupation and unchanged total reconstruct the input; for fixed d, j is already fixed. Hence every branch is injective, the Kraus completeness sum is I, and the complete instrument is CP including arbitrary idle quantum references. There is no rejection in this complete modular sector. `commit` retains all outcome probabilities; `step` includes the no-commit complement.

The sector preserves `G_v=u_v+sum incident_port F`, and `G_e=u_e−F_ea−F_eb`, modulo N. Summing any region gives exactly `Q_R+cut(F)=G_R mod N`, including transit vertices. Individual representative registers do not define exact integer Gauss charges. `decode_bounded` requires an independently valid absolute bound B and N>2B. Nonnegative fixed L bounds regional charge changes by L; net signed total alone would not supply such a bound for a signed extension. Three units alias zero modulo3; modulo7 uniquely decodes signed changes in [−3,3]. Historical winding remains unobservable from residues alone.

Tests exercise the exact cyclic one-load generator and MIX instrument, all C/E/MIX completeness sums on a 35-state three-unit star, trace/positivity/event accounting, every tested cut, and N3 rejection/N7 decoding. The two-unit event (2,0)->(1,1) is deterministic under both coins, rejecting a substitution by independent one-particle commits. Individual-record/definite preparation is kept distinct from the one-load uniform-fiber preparation. No many-load uniform-fiber reduction or physical species identification is established.
