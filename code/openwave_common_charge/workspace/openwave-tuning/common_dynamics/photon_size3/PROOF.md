# Size-three native free-field member

This is an additive member of the frozen noncompact quadratic v4 source.
It consumes the unchanged opposing photon audit and its transitive pins.
The actual periodic LOCAL12 complex has 324 ports, 972 edges and 702 faces.
The parent spatial-domain audit proves rank B=323 and rank A=646, with
three independent harmonic periods. `GEOMETRY.json` stores its integer
incidence, face and winding data. The code rechecks BA=0, period A=0,
B winding=0 and period winding=3I. Positive and negative axis neighbours
are spatially distinct at this size.

`FRAME.json` freezes one rational proposed frame Q and positive curvature
list Lambda. Floating eigendecomposition selected these integers once;
it is not the certificate and is not rerun when consuming the member.
The exact certificate rechecks Q^TQ-I and KQ-QLambda, K=AA^T, using signed
18-bit limbs. Every integer product and absolute partial dot-product sum
is below 2^53, so binary64 GEMM here performs exact integer arithmetic;
Python integers recombine the limbs and form squared Frobenius sums.
Integer square roots give outward rational norm bounds epsilon and r.

Put delta=r/min Lambda and

    gamma=epsilon+2(1+epsilon)delta+delta^2,
    eta=r+(1504/17+max Lambda)(delta+gamma).

Y=KQ Lambda^-1 lies exactly in the coexact range and differs from Q by
at most delta. Its Gram error is at most gamma<1/2, so U=polar(Y) exists,
spans all 646 positive modes and obeys ||U-Q||<=delta+gamma. The exact
curvature C=U^T K U then obeys ||C-Lambda||<=eta. The integer source has
no missing positive mode because the frame rank meets the exact upper
bound. With a=1/2,b=3, the square-root operator estimate bounds

    ||sqrt(ab C)-sqrt(ab Lambda)||
       <= (5/4)eta/sqrt(min Lambda).

The certified cap-one Hamiltonian includes the zero-point scalar and all
646 one-photon states. Its operator error is at most (646/2+1) times the
frequency bound, namely about 1.15270048e-6. The cap is a compression of
the full free Hamiltonian, not a literal square of capped field matrices.
The three harmonic Gaussian factors stay in the exact source; this
terminal radiative instrument traces them, by the source-contract proof
in ROW_29, rather than assigning them an oscillator ground state.

Physical labels are vacuum 0 and photons 1..646; labels 647..1023 are
retained padding with declared zero generator. `walsh.py` implements the
exact complete Walsh transform through ten bits. The emitted 25606-pulse
diagonal keeps all scalar and zero coefficients. Its certificate consumes
every primitive against all 1024 outward phase targets and includes the
native controls. Integer direct-character sums and transform-square
identities test the compiler; unsupported domains reject.

`native_member.py` executes native response matrices with the independent
parent audit executor, including selected physical and padding columns.
Native preparation resets the ten-bit spectral code and applies the
actual Hadamard to prepare vacuum plus photon1. Both field X/Y analyses
and both common-clock X/Y analyses are emitted. Disjoint field/clock
factorization is used exactly, with no post hoc phase alignment. The
result retains the visible full zero-point scalar and every physical and
padding output. Both native pointer laws supply all 2048 coarse outcomes
and the one-bit fine effects whose tensor contractions specify the full
fine law. Eleven blank and eleven pointer-ready reset failures, first
export complements and pointer errors retain their own outcomes. No
success renormalization occurs.

The emitted words have 25615..25617 pulses before pointer reads. The
maximum complete-record TV upper bound is 1.26596789e-6. The recovered
first excitation frequency is about 0.604305766664; the actual source
certificate supplies its unaliased band. Both X/Y field settings are used.
A native counterphase erasure executes the deleted-frequency control;
the Y contrast changes by about 0.56819. No absolute vacuum-energy
measurement, interacting Weyl implementation or local edge-register
encoder is inferred. This member and the controlled physical-field limit
in ROW_29 share the actual face source and spectral-code architecture.

The interrupted integer-GEMM attempt remains recorded in RUN_ATTEMPTS.
Its frozen frame was retained and the exact arithmetic was optimized;
it was not replaced by a favourable spectral sample. The manifest is a
producer freeze for independent review, not a self-issued row validation.
