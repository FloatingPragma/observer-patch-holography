# Method and receipt

The first-step native record equation is

    (A+2w)h_z-w h_(z-1)-w h_(z+1)=0,

with the reflecting-bottom equation(A+w)h_-R-w h_(-R+1)=0 and top value1.
Its transfer is cosh(omega(R+1/2))/cosh(omega(R+3/2)), with
omega=2asinh(sqrt(A/w)/2). `PROOF.md` derives this directly from the native
current probabilities and proves convergence to the occupation slice transfer
exp(-omega). `probe.py` constructs its pre-hit matrices from the retained
endpoint census, then solves their absorption and capped-time equations.

All405,755 pinned native calls are accounted:48,600 census branches,
316,803 actual all-incident attempts and40,352 vertical-only control attempts.
The two schedules each have4096 actual trials. No run was censored. The
actual first-hit means at the three spatial indices are
0.22143555,0.14709473,0.09582520, against capped population targets
0.20879121,0.12702079,0.09167368. Each discrepancy is below the registered
simultaneous99% halfwidth0.04412877. These intervals use bounded independent
trials and a union bound over all six comparisons; no completed-only
conditioning or postselection is used.

The exact native expected attempt cost is76, with observed mean77.34448.
The cap4096 has exact population censor tail3.526e-15. The changed schedule
gives Fourier means1,1,1, exceeding each source mean by more than0.5.
Every first-hit record replays from its starting position and retained
current letters. A complete rerun reproduces receipt and compressed raw
records byte for byte, bound in `REPLAY.json`.

Fixed finite bottom needs care: a subsequent1→2 hit above the same bottom
has h_(R+1), not h_R. The identical repeated transfer in `PROOF.md` section3
is stated only for the infinite-depth cylinder. No moving boundary or hidden
reset is used in this experiment. The infinite-depth expected attempt cost
diverges; finite strips and caps retain their resource bounds.

**Scope:** actual native first-hit record transfer and its Euclidean clock
dictionary on a registered unpinned N+ protocol; Lorentzian waves are the
reconstructed effective description, not primitive endpoint updates.
