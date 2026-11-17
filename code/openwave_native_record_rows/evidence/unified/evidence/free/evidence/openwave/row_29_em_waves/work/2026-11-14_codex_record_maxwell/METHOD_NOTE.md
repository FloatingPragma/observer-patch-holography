# Method and results

For the branch-derived macro transfer P=I-L/(2W), the normalized n-event
occupation covariance on an eigenmode z of P is

    c_n(z) = (1+z)/(4W(1-z))
             - z(1-z^n)/(2Wn(1-z)^2).

`experiment/probe.py` constructs every nonidentity macro count from the
retained native current branches. It verifies equality of the twelve
port-conditioned macro rows. It evaluates this finite geometric-sum identity
using the eigensystem of that216-state count matrix. `PROOF.md` maps its
internal constant-band readout and nearest-neighbor differences to the target
tensor covariance; the shared occupation proof supplies the identity above.

The three spatial indices are(1,0,0),(1,1,0),(1,1,1), on a3x3x3x8 periodic
array. With k_j=2sin(p_j/2), Q²=sum k_j², the target slice energy is
omega=2asinh(Q/2). The measured population reader is
acosh((g(1)+g(3))/(2g(2))). This calculation uses no fitted coefficient.

At n=10^9, the largest separated covariance relative error is3.0695e-5
(gate1%), the largest recurrence relative residual is6.7394e-6(gate1%),
and the largest energy relative error is4.6932e-6(gate1%). The reflected
electric longitudinal residue is at most3.2092e-7; the two transverse
eigenvalues agree with their magnetic partners within the registered1% gate.
The n=10^7 comparison is retained in full and shows the finite-block bias.

Actual source work is167,008 native endpoint calls:101,088 compressed
endpoint/spectator branches,384 full-vector equivalence checks and65,536
events on four disclosed preparations. The latter contain29 nonzero currents.
Every event, including every zero current, is retained. A second observer
reconstructs the occupation counts from the tape and initial reference; all
counts match exactly. Six actual field-strength readouts are retained.

The full rerun reproduced receipt and compressed raw records byte for byte;
the before/after SHA256 values are in `experiment/REPLAY.json`. The numerical
population result is not an estimate from65,536 events. Blocks of10^7 and10^9
events were summed through the exact transfer formula; they were not executed
as long tapes. No finite-block Gaussianity, finite-sample precision, physical
phase experiment, clock or charged interaction was measured.

**Scope:** conditional native-record Maxwell covariance and Gaussian-limit
reconstruction on an unpinned N+ product protocol; all physical row gates open.
