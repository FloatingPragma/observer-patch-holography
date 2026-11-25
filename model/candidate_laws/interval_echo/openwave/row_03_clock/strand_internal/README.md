# Can the closed strand tick at its measured rest frequency?

The two commissioned Interval Echo rings have reconstructed rest frequencies of 15.000082 and 14.999985 radians per dimer-clock unit. No free internal measurement on either whole ring has a spectral line between 14.99 and 15.01. The distinction can be checked on observer records. Pair spectroscopy recovers the rest parameter by combining three measured sum frequencies. A single dimer, prepared and read in its own exchange basis with no drive during flight, produces a band near twice that parameter. [Rest measurements](../../submission/row_30_cautions/README.md#the-ring-source) · [Exposure history and controls](PROTOCOL.md) · [Count banks](records.json)

## The complete internal spectrum

The native blank dimer and the coherent unread connector reduce the guarded ring to the physical orientation Hamiltonian
$$
H=J\sum_j Z_jZ_{j+1}-h\sum_j X_j,\qquad J=8,\quad h=1/2 .
$$
An additive constant is common to all the words and both flip sectors. The earlier coupling meter independently gives a 95% counting interval of [14.996789, 15.004576] for the zero-wave-number constituent frequency. The fresh replay recommissions the same two 80-cell boards and verifies both source circuit digests. The bound below uses the declared microscopic couplings, and the measured interval tests the identification. [Native restriction](../../row_30_collective/PROTOCOL.md) · [Independent coupling meter](../../row_30_collective/coupling_meter.json)

In a sector of the physical flip $F=\prod_j X_j$, the cyclic wall transform has one-body block $2JI-2hS$. The signed shift $S$ has either the periodic or the antiperiodic momentum grid, depending on $F$ and ring length. The quadratic fermion Hamiltonian gives physical energies $\sum_k (b_k-1/2)\epsilon_k$ up to the common constant, with occupations $b_k=0,1$ and
$$
\epsilon_k=\sqrt{257-32\cos k}.
$$
In either sector $\sum b_k=n\pmod 2$. To see the sign, turn the transverse flip $h$ down to zero while retaining $J$. At zero the excited bonds are aligned neighbours. There are $n$ bonds, and the number of opposite neighbours on a closed binary ring is even. Complementing every physical spin leaves these bond records unchanged, giving each flip sector the same allowed occupation parity. Fermion parity is preserved by the quadratic Hamiltonian, and none of the positive singular frequencies passes through zero on $0\le h\le 1/2$, since each is at least $2(J-h)$. The occupied parity therefore persists on both grids. The replay independently compares **every physical energy** in both flip sectors on 5-, 6-, 7- and 8-dimer rings with the corresponding parity-restricted occupation sums, including their zero-point energies. [Ring reduction](../../row_30_collective/rest/README.md) · [Physical orientation check](measure.py)

Since $\epsilon_k\ge15$, rationalizing gives $|\epsilon_k-16|\le |1-32\cos k|/31$. Both momentum grids have $\sum\cos k=0$ and $\sum\cos^2 k=n/2$. Cauchy–Schwarz gives $\sum|\epsilon_k-16|\le n\sqrt{513}/31$. Each energy hence differs from $16(\sum b_k-n/2)$ by at most half this bound. In a difference of two energies the centres are separated by an integer multiple of 32. **This holds even for a coherent superposition of the two flip sectors** because the bound includes the zero-point shift on each grid. If $\omega$ is any physical Bohr frequency, then
$$
\operatorname{dist}(\omega,32\mathbb Z)\le n\sqrt{513}/31 .
$$
The right-hand sides are 12.42070 and 14.61259 for $n=17,20$. Even the wider enclosure excludes the entire rest test interval. The replay checks the restrictive endpoint by integer arithmetic. For any initial density matrix and fixed terminal effect, the probability is a sum of phases at these Bohr differences. This includes a named joint terminal outcome on several dimers and a preparation herald completed before release. A later ratio of postselected probabilities, phase relative to a distinct flight arm, or an algebraic combination of detected spectral lines is a different observable. The bound concerns these two closed, commissioned rings with their coherent connector memories and free evolution in the guarded orientation space. Larger rings need their own bound.

## The local record

Choose a random physical X-product word at the registered total flip parity, retaining the X sign of one dimer. Let the whole ring evolve freely until one terminal X read on that dimer. On independent copies a local dimer hand supplies two timing quadratures. Each register bank retains all four source/terminal sign outcomes over 1,024 waits, including both source signs and every opposite-sign return. The linear statistic is the mean product of these signs. In a maximally mixed flip sector it equals the two-time trace of the physical $X_j$. In the dual representation $X_j$ is one Majorana bilinear across a neighbouring bond. The trace is the two-by-two minor of the propagated Majorana basis on that bond, since distinct quadratic products are trace-orthogonal in these parity sectors. The replay checks the trace directly against physical spin matrices on each small ring. [Generator and decoder](measure.py)

On the 17- and 20-dimer rings, the measured native pair-band maxima are 0.03787 at frequency 31.878 and 0.04724 at 32.046. The largest amplitudes in the [14.9,15.1] band are 0.000534 and 0.000295. These are finite-window measurements with shot noise and leakage. The spectral exclusion follows from the full-ring bound. Exchange-only flight leaves each X sign fixed. A planted diagonal contact chosen to give a 15-unit transition shows that the very same read resolves that frequency, at amplitudes 0.50015 and 0.50005. Frequencies come from the timing counts and all settings use 65,536 launches per appointment. The register words are simulated from the exact free-sector probability law. The one-dimer detector has frequencies at most 34 on the native source, below the sampling Nyquist frequency. The all-effect proof refers to unsampled physical frequencies.

The source does produce observable internal oscillations. On these rings they cannot supply the constituent rest frequency required by OW-03. The prepared medium, instrument and calibrated rest definition retain the scope of the original OW-30 study. Interval Echo on other preparations and a full charged-particle rest-energy clock remain open. [OpenWave criterion](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L103) · [Upstream M10 clock entry](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L437)

From the repository root, reacquire and check every retained count and the spectral gates with

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_03_clock.strand_internal.measure
```
