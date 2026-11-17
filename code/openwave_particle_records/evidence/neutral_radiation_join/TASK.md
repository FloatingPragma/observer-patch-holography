# Frozen finite field experiment

Read DESIGN.md. All source code is hash bound in SOURCE_PINS.json. No measured
values of nature, fitted constants, tail trimming or selected successes.

1. Run every contact branch (30 seams ×2 orientations ×2 coins) on each even
ring L=8,16,32 through first antipodal detection. Execute unchanged R-work
charge attempts, both actual local swap layers, and the complete inverse
detector menu. Preserve public references, every native current, swap and
detector outcome. All field readouts use only those records, not a source state
snapshot. Require zero reference/current/guard/state-reconstruction defects.
2. Before and after every transport cycle verify exactly, for each i,
`u'_i=u_(i-1)`, `v'_i=v_(i+1)` and
`ε'_i-ε_i = J_(i-1/2)-J_(i+1/2)`, where
`J_(i+1/2)=u_i²-v_(i+1)²`. Verify the equivalent first-order E/B update with
integer arithmetic after multiplying by2. On all three-time windows without
a source/detector operation verify `F_i(t+1)+F_i(t-1)=F_(i-1)(t)+F_(i+1)(t)`
for both F=E,B. Explicitly separate source, transport and absorption phases.
3. Every native emission must transfer load-square -2 to field norm +2 and
every absorption reverse it. Require initial/flight/absorbed combined budget2,
charge0, first antipodal response L/2 and zero field after absorption. These
are the old declared dimensionless budget and clock, not physical units.
4. Derive the full 2L×2L one-cycle field transfer from executed single-packet
impulses at every lane and site for each L. Retain those planted references
and their native attempts/swap/detector-disabled outcomes. Require exact
permutation orthogonality, expected one-edge opposite routing, and integer
wave-matrix identity `T²-S T+I=0`, S being spatial nearest-neighbor sum.
For every discrete Fourier k=2πm/L, require observed transfer eigenphase
residual <1e-10 against exp(∓ik). Record all modes, including zero/Nyquist;
do not fit a speed. The same reader must recover all nine planted local
amplitude pairs {-1,0,1}² and known E/B values exactly.
5. On every L run no-pair, transport-off, detector-off, and one-current-tag
scramble controls (seam0, orientation+1, coin0). Keep all their outcomes.
No-pair gives zero field; transport-off must fail the unit-speed field
transfer and finite-arrival test; detector-off and scramble preserve field
norm but fail original-dipole response. A scramble is an explicit intervention,
not a native event. No control failure is relabeled a source pass.
6. At L=8 only, for each contact branch, read both packet slots at source
post-emission, site1 after first transport, antipodal site just before
detection, and antipodal site just after detection: 960 symbol services.
For service index h and symbol s, initialize carrier r=0,1 to full78 indices
`(s+13h+17r)%78`. Use independent reproducible coin streams with prefix
`characteristic-field-v1:<h>:<s>:<r>:<stage>` and SHA256 counter bits, low bit
first in each byte. Each reset/read uses four seed scans and32 retries;
each write edge allows32 retries. Require every conclusive symbol/E/B read
correct, and report all inconclusive services in the denominator. Recompute
the analytic per-service bound and use a union bound for the entire set.
7. At each of the 120 post-emission checkpoints run the same symbol services
with writing disabled (240 further services, separately denominated); every
conclusive result must be blank. Force one write timeout with32 holding coins.
Force a read census seed timeout by supplying the holding coin for every
requested edge through four scans on a known binary control state. These
deliberately changed-coin controls use actual native calls and must return
no guessed symbol. All complete records and failure work are retained.
8. Independently replay every symbolic event and every memory fine record,
checking native floor/ceiling arithmetic and current continuity. The headline
must read the retained files, verify source/manifest hashes and all checks;
it must not regenerate success claims from stored summary Booleans alone.

Exact checks use integers; only Fourier phases use the fixed1e-10 bound.
The execution saves every case, every control, all checkpoint memories and
their outcomes. Any implementation error is retained in a numbered attempt
before correction. No new evidence may be executed before this packet's
DESIGN/TASK/source/code hash freeze is sent to the parent for registration.
