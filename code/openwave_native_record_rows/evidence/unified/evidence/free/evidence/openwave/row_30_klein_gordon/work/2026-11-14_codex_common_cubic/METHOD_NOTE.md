# Method and headline

The record-derived identity is P=I-L/(2W), W=6318. The full native census
checks L=L_macro tensor I12+I81 tensor L_internal. The internal bands
are0,5-sqrt5,6,5+sqrt5 with multiplicities1,3,5,3. The finite occupation
covariance and the period3 first-hit matrix are written in `PROOF.md` and
computed by `probe.py` from the branch counts.

Actual source work is169,364 calls:37,908 endpoint/spectator branches,
384 full-vector checks, and131,072 actual events across four components.
All four tapes independently reconstruct every occupation, scalar-band
readout, curl value and clock tick. Their completed clock counts are
4,5,3,0. All four final open intervals remain in the data; another-group
style same-group replay found lengths9442,5276,1717,32768. In particular,
the component with zero completed ticks was not removed.

The native period3 first clock has mean12636 attempts and tail0.07477077
at32768 attempts. These are population statements, not inferred from four
samples. The finite-block population comparison at n=10^8 has maximum
contact-accounted covariance relative error6.7631e-5(gate0.001), source
Fourier residual2.418e-14(gate1e-10), and first-hit agreement4.180e-14
(gate1e-10). The massless electric/magnetic relative difference, with
contact retained, is at most0.00046081(gate0.001). Its small longitudinal
electric eigenvalue is negative because of the finite contact; it is
retained, not called a positive finite-volume polarization.

Independent-location resampling has only contact covariance and zero
nonconstant first-hit transform under the same readers. The finite
positive Green ratio serves as an analytic reader reference. A full
rerun reproduced receipt and raw compressed records byte for byte;
`REPLAY.json` binds the before/after hashes.

**Scope:** actual occupation-only cubic common-source member with shared
field/band/clock records and exact compressed population checks; its
free-wave continuum interpretation is supplied by the bound source
theorem, with pin, in-platform filing and another-group review pending.
