# One periodic source supplies the fields and the first-hit clock

This corollary joins the finite first-hit experiment to the periodic native
product family, without declaring the two full tapes identical.

Use a vertical cycle with H>=3 slices, starting at slice0 and stopping on
the first native hit of slice1. These are relative to the disclosed initial
slice, so a uniform initial preparation is preserved. For an internal/spatial eigenmode A, cut
the cycle at that target. Its unfolded path has boundary values1 at heights
1-H and1, and the same native harmonic recurrence

    (A+2w)h_z-w h_(z-1)-w h_(z+1)=0.

For omega=2asinh(sqrt(A/w)/2), its unique solution gives

    h_cycle(A)=cosh(omega(H/2-1))/cosh(omegaH/2).              (1)

For A=0 take h_cycle=1. This first-hit transfer is a function of the same
internal/spatial generator as the occupation field. In particular, H=2R+3
makes(1) exactly the reflecting-strip transfer
cosh(omega(R+1/2))/cosh(omega(R+3/2)). The first-hit packet's R=2 native
experiment therefore checks the stopped endpoint transform of the periodic
H=7 family, using a reflection fold and deletion of holds. Attempt times and
full histories remain different; this is an equality of the specified
stopped observable laws, with an explicit mathematical map.

More generally,

    0 <= h_cycle/exp(-omega)-1 <= exp(-omega(H-2)).            (2)

Thus the cycle transfer converges to the occupation field's infinite-slice
transfer with a controlled finite-period correction. At finite H the
one-hit transfer is not substituted for the exact thermal slice generator.

A finite cycle also gives an exact repeated record clock: after first hitting
plane j modH, stop next on the first hit of plane(j+1) modH, and retain j in
observer memory. By strong Markov and cyclic translation, its endpoint
transfer is h_cycle^n after n such clock ticks. This is an unwrapped stopping
counter on a finite source, not a claim that every raw unit hop is future
causal. It is different from repeated planes above a fixed reflecting bottom,
which have successively h_R,h_(R+1),... and were not claimed equal.

For global IID seam draws the exact mean attempt cost per cycle tick is
W(H-1)/w; for the incident scheduler on the regular graph it is
d(H-1)/w, d=5+8w. The nearest-neighbor vertical walk on an H-cycle has mean
H-1 vertical jumps to its adjacent target, and the corresponding vertical
jump probabilities are w/W and w/d. These finite costs diverge in the
infinite-depth limit. Cutoffs and missed-clock outcomes must be charged in
an execution; no free infinite-resource clock is inferred.

The common proposed field candidate can therefore use only periodic routed
graphs and the globally IID native event law. Its observer's retained tape
supplies both the occupation statistics and the stopped endpoint clock.
Extra clock copies, a separate clock Hamiltonian, a moving boundary and
primitive real-time unitaries are not needed for this source-time dictionary.
The four record channels needed by the Maxwell readout remain part of the
declared common apparatus.
