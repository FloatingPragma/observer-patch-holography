# Charge, rest energy and frame on one signed source

One construction carries all three results: the twelve-port carrier, a signed
source on it, and the native repair law read by a complete current pointer.
Declared architecture: the positive-amplitude memory lift of the repair branches,
the finite carrier graph, the preparation and pointer services, and the physical
dictionary that reads the primitive boundary source as electric charge, the
transfer spectrum as energy and the first-current mean and covariance as frame
and field.

## Carrier, source and instrument

On a seam (a,b) native repair records

    d = q_a - q_b,
    j = (d - (2c-1)(|d| mod 2))/2,
    q' = q - j e_a + j e_b,

with coin c. Given q', the seam, c and j, the input q is unique. Every branch of
the positive-amplitude lift is therefore a partial injection, and the branch
effects together with the invalid sector sum to the identity: the source carries
a complete quantum instrument. The one-particle chart has 256 addresses and the
pair chart 65,536. The same work register is prepared, evolved and read; the
complete current pointer has preparation-plus-readout error below 10^-6 over the
whole experiment, reset failures included. The sealed run makes 435,080 literal
native calls; a standalone copy of the published tree replays it from its own
root.

## 1. Charge comes in whole units

Integrated current F obeys BF = q_initial - q_terminal, with B the port-seam
incidence. On the twelve-port cell (30 seams, 20 triangular faces) rank B = 11
and the face matrix D has Smith factors 19 units and one zero. Hence
ker B = im D over the integers and

    Z^30 / im D = H_1(S^2, V; Z) = reduced H_0(V; Z).

The class of one signed pair e_p - e_r is primitive. Two sharp half-valued
daughter charges lie outside the measured source image. A state that gives
probability 1/2 to each of two integer outcomes is a probabilistic alternative
between whole charges; it is not a pair of fractional charges.

The response depends on the law. Native expected seam current is d/2, so the
stopped mean flow is

    E[F_tau] = B^T L^+ q,

an exact rational Bellman equation on every source, with almost-sure absorption
on the finite transient domain. All 264 native trials absorb and synchronize; all
48 changed-law controls stay nonabsorbing; a planted half-current is read by the
same incidence decoder as two sharp halves and rejected. On the full 240-port
graph the integer image C_1/ker B = im B is the zero-total source lattice; the
sphere-face certificate is stated for the cell.

## 2. Rest energy of a compact charged state

The full graph has 240 ports and 630 seams; its one-particle transfer is
T = I - L/1260. The four-port charged ray f = (|4> - |5> - |6> + |7>)/2 satisfies
Lf = 6f. With H = -1260 log T and the constant mode at zero energy,

    E_6  = -1260 log(209/210) = 6.014331228348673,
    E_12 = -1260 log(104/105) = 12.057508280349891

for f and for the disjoint pair f (x) g (g on ports 16 to 19). The value 511 keV
is the permitted calibration that fixes the unit; it is not a predicted number.
With that one unit the pair comes out at 1024.4508487 keV, 2.45 keV above twice
511 keV; the excess E_12 - 2E_6 = 1260 log(43681/43680) belongs to the finite
clock. The spectrum also has lower nonconstant modes (the difference of the
first two cell indicators has Rayleigh quotient 1/3), and the energy pointer
|f><f| certifies the eigenray while every complementary outcome is kept as an
unresolved-energy outcome. The energy uses only the full 630-seam kernel;
addressing the 30 cell seams for the charge and frame assays leaves it unchanged.

## 3. Field and frame under a half-turn

For the antipodal preparation e_p - e_antipode(p), the first current is
Delta = j(v_b - v_a). Summing all 60 alternatives for each of the 12 sources
(720 branches) gives, with golden ratio phi,

    m = ((phi - 3)/15) v_p,
    C_parallel = (3 - phi)/9,   C_transverse = (phi + 2)/15.

Mean and covariance are measured independently from the current records. The
executed proper half-turn R = diag(-1,-1,1) returns C and reverses m; R^2 returns
both. The covariance differs between non-antipodal axes, the fixed-seam native law
fails the source-axis and shape gates, and an off-axis plant is detected, so the
return is a property of the records and not of the reader. The executed rotation
is a finite symmetry of the addressed cell; it is kept distinct from any
continuous interpolation and carries no exchange-statistics claim.

## Reproduction

The experiment records, receipts, manifests and the standalone replay are
published in `code/openwave_static_signed_source/`.
`python3 -B code/openwave_static_signed_source/replay.py` reruns the complete
experiment (435,080 native calls) in a disposable copy and compares every
regenerated receipt and record. `python3 code/openwave_signed_source/verify.py`
checks the published evidence against its manifest and recomputes all three
results from the carrier graph and the native rule in a few seconds
(`recompute.py`, which imports none of the experiment scripts).
