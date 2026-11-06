# Native preparation followed by the complete phase/work instrument

Let S be the 78-dimensional total-two native carrier and B its 66-dimensional
binary subspace, with inclusion J. Let R be any finite untouched reference.
The input rho_SR is an arbitrary density on S tensor R. The total-two domain,
the fair placement coin and the quantum C/E instrument semantics are premises.
No diagonal, binary-support, product-state or stationary-input premise is used.

Every preparation sweep visits the thirty source seams in the fixed order of
`model/native_equilibration/sweeps.py`. The preparation tape and any permitted
law-only controller are ignored at its boundary; subsequent controls cannot
access them. There is no conditioning on their values. The all-law theorem in
that module permits the stated native current-Gram family and its specified
law-only feedback; the numerical executor here uses fixed C or E.

## Preparation bound

The exact sweep theorem proves that one covering sweep removes collision
support. Let P be the exact 66-state diagonal kernel of one sweep, with output
rows and input columns, and let U be the uniform kernel. After m further
sweeps set Q_m=P^m-U. For the entire admitted current-Gram law family the
off-diagonal Hilbert--Schmidt contraction is at most c^m, c=1/32. Fixed C and E
have stronger bounds, which are unnecessary for this composition. The theorem
gives

```
1/2 ||Phi_(m+1)(rho_SR) - J (I_B/66) J* tensor rho_R||_1 <= epsilon_m,
epsilon_m = min(1, (66/2) max(||Q_m||_2,c^m), ||Q_m||_1/2+(66/2)c^m).
```

Here matrix norm 1 is maximum absolute column sum. Exact integer matrix
multiplication computes Q_m; row/column bounds and a rational upper square
root bound its spectral norm. This proof does not infer a singular norm from
the eigenvalues of a nonreversible scan. It preserves the ORIGINAL reference
marginal. One may therefore trace out the preparation tape even when the
retained field is entangled with R, without treating its conditional state as
an independent input.

The first all-law full-carrier certificate below 2^-24 uses 19 sweeps, or 570
native attempts. Its exact epsilon is

```
1009634291956847693669001532064517518244609244697758873135
/27150825195799899695556991512196283332215141782954067745252048896
```

This is approximately 3.71861365051e-8. This finite error replaces exact
stationarity; the certificate does not assert exact mixing in finite time.

## Complete instrument composition

Let Gamma_N be the exact ideal complete instrument of
`model/clocked_quantum_readout`, including its independent work preparation,
native herald, stopping slots, phase detector, all sixteen work records,
feedback, failed herald and two unmeasured timeouts. The balanced coherent
operation and the work dynamics are the parent's supplied services. Finite
Kraus composition (`Lean/Dynamics/ChoiCPTP.lean`) makes this a trace-preserving
map to a classical direct sum of terminal quantum states. Tensoring with R
does not change this assertion. Contractivity gives

```
1/2 || (Gamma_N tensor id_R)(Phi_(m+1)(rho_SR))
       - (Gamma_N tensor id_R)((I_B/66) tensor rho_R) ||_1 <= epsilon_m.
```

The inclusion J and its inverse on binary support are implicit after the
covering sweep. The same argument allows any common CPTP continuation with
its declared initial resources, including adaptive interventions on suffix
records. No factor depends on cap N, record count, herald yield or reference
dimension. The number of joint suffix records is 128N+3. Comparing only the
classical records or one event cannot increase the distance. Since the ideal
outputs are normalized, the error of an event probability is at most epsilon.
Appending the field work register therefore needs no additional preparation
estimate, and it does not require executing its large tensor matrix again.

This statement concerns the full SUFFIX instrument. A retained warm-up tape
can remain correlated with S and is excluded from the comparison. Accessing
that tape during continuation changes the map and invalidates the use of its
unread marginal. Physical erasure is unnecessary: ignoring a stored tape and
erasing it are different operations. Conditioning on a successful suffix
herald also changes the error contract; epsilon is unconditional.

## Binary source mixer and numerical reference

The parent compiles the balanced operation on its full 512-address carrier,
including padding leakage, and gives an all-input two-call half trace-norm
allowance delta. The native binary experiment is embedded into that carrier
only after the absorption sweep. The parent's CPTP padding extension applies
between mixer calls. Its estimate already charges the possible norm gain of
the binary mixer matrices. Apply that estimate to the ACTUAL prepared input,
then use ideal contractivity above. Thus the half trace-norm error from the
stationary ideal reference is at most epsilon+delta. There is no division by
the probability 5/33 of a useful herald.

Binary source matrices need not be exactly trace preserving. Their individual
event weights therefore have the bound epsilon+2delta relative to exact
stationary ideal event probabilities: the epsilon term compares normalized
ideal states, while 2delta uses the full norm for the binary comparison.
For C-minus-E event weights, add both candidate allowances, 2epsilon+4delta.
The receipt uses exact rational phase-match weights from the parent's oracle.

If b_N bounds the parent's reported numerical complete output against its
stationary exact ideal reference, the total half trace-norm comparison to the
prepared source-mixer instrument is b_N+epsilon+delta. The parent's b_N is
consumed with its source-pinned receipt, not recomputed from decimal displays.
Finite cap is essential. No unlimited binary-source terminal-state limit is
asserted, and the parent's exact stationary failure and timeout probabilities
become bounded probabilities after finite preparation.

The executor composes actual source per-seam C/E Kraus maps with the actual
native herald, current and pointer maps and the ideal balanced service. It
accepts arbitrary input/reference matrices and retains subnormalization.
Its floating residuals are diagnostics. The analytic composition uses the
exact sweep theorem and exact CPTP instrument; floating execution is not an
outward error certificate for its own returned entries.
For a subnormalized input of mass a, linearity multiplies the analytic
preparation and source-comparison allowances by a. Even a requested tolerance
of one executes the covering sweep required to enter the binary readout domain.

## Physical scope

This construction removes independently supplied exact I/66 preparation and
initial binary support from this phase/work experiment. It does not derive
the total-two sector, fair coin, quantum state/effect probability rule,
balanced target design, coherent compiler services, effective work action,
work vacuum preparation, occupation cutoff or model clock. The uniform
reference has trivial modular flow on the full binary algebra. Its derivation
does not identify a physical thermal energy, temperature, pressure or equation
of state. Scalar record consensus alone does not establish Born affinity or
a physical instrument attachment.
