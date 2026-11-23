# Preparation and recovery of three charged levels

The experiment uses the Gauss-projected field extension, two identical charged-cell Hamiltonians related by the physical port rotation, and the material bridge from port 3 to port 48. Its inputs are the lowest three isolated-cell levels of each sign, ordered by native energy. The energy meter, phase reference and delay controls are supplied apparatus. Simulated energies and waiting times use the extension's units. The pulse schedule has no calibration to an observer clock or to keV.

At native time π, each of the six prepared inputs can carry charge through the bridge. For either sign the receiving cell tends to record its lowest charged energy, even when the source prepared an excited level. The table reports charge arrival per selected input and correct energy recovery per charge arrival. Herald-and-recover counts both the source-selection probability and the correct destination outcome in every attempted preparation.

| Signed input, level | Native energy | Source herald | Charge arrival at π | Correct level given arrival | Herald and recover per attempt |
| --- | ---: | ---: | ---: | ---: | ---: |
| −, 0 | 4.343760 | 5.8234% | 3.9452% | 76.3832% | 0.17548% |
| −, 1 | 4.620847 | 3.3437% | 1.9247% | 0.00977% | 0.00000628% |
| −, 2 | 4.720239 | 2.7209% | 2.6055% | 2.3006% | 0.001631% |
| +, 0 | 5.450006 | 6.1932% | 76.4282% | 85.5277% | 4.0483% |
| +, 1 | 5.730242 | 3.8218% | 70.1289% | 1.2870% | 0.03449% |
| +, 2 | 5.919776 | 2.8306% | 68.8515% | 1.1377% | 0.02217% |

A reference test detects entanglement across the fully active bridge. At π the normalized three-level output has partial-transpose negativity $1.3605\times10^{-5}$ for the negative charge and $3.4580\times10^{-5}$ for the positive charge. The moving seam-only control has zero negativity. For equal source priors, even the most successful reassignment of the three receiving energy labels identifies the source correctly in only 48.28% and 37.09% of the accepted three-level records, respectively. These rates refer to this energy measurement. An additional calibrated sender and receiver filter can transmit the unknown three-level state exactly on a heralded branch, with probabilities $6.7650\times10^{-10}$ and $6.1585\times10^{-6}$ at π after ideal source selection.

The negative-sector level energies divided by their lowest value are 1, 1.063790 and 1.086671. The positive-sector ratios are 1, 1.051419 and 1.086196. These are finite one-cell energies. This two-cell measurement supplies neither a momentum-resolved rest energy nor a repeated-cell trajectory or evidence of three lepton species.

## Preparation

The sender starts with one signed load at the charged anchor and a uniform superposition of accessible field words. The source may imprint the parity of the edges in one dual path between two faces. Each target uses the largest spectral weight among the 190 face pairs and the unmodified source. The source record includes this weight, the weights at the other two levels and the probability discarded by an energy herald. Different filters can extract all three targets from the same underlying law.

The ideal herald projects onto level $i$, with probability $w_i$. A finite alternative prepares a uniform superposition of $N$ controlled delays separated by $s=\pi/B$, where $B$ bounds the charged spectrum from above. The sender applies a phase $e^{iE_i k s}$ to delay branch $k$ and projects the delay register onto its initial uniform state. The amplitude response is

$$
F_i(E)=N^{-1}\sum_{k=0}^{N-1}\exp[-i(E-E_i)ks].
$$

At $E_i$ the response is one. With nearest level separation $\delta_i>0$, the off-target amplitude is at most $1/[N\sin(s\delta_i/2)]$ over the rest of the finite spectrum. The receipt gives branch counts sufficient for 99% fidelity after a successful herald and the corresponding interval for the herald probability. These filters require coherent delays and energy references. The two-cell calculation propagates the ideal selected inputs. Source-inclusive transport yields refer to that limiting instrument.

## Receiver and reference

The two-cell state consists of a left-charge block and a right-charge block. In the latter write $R_i(q,a)$ for charged configuration $q$ and departed neutral configuration $a$, given a prepared source level $i$. The receiving charge counter reads the squared norm of this whole block. On a right-charge click, the energy instrument tests the same three isolated-cell eigenstates $u_j$. Its joint weight is $\sum_a |K_a(j,i)|^2$, where $K_a(j,i)=\sum_q u_j(q)^* R_i(q,a)$. The remaining outcomes are another energy in the receiving cell or the charge staying in the sending cell. Reading the sending neutral ground state as well gives $|\sum_a n(a)^* K_a(j,i)|^2$.

Each $K_a$ acts on an arbitrary coherent input in the three selected levels. Tracing out $a$ gives the entire receiving three-level operation, including its input-dependent acceptance. Applying it to half a maximally entangled three-level state supplies a reference test. Negative eigenvalues of the partial transpose certify that the accepted output can retain entanglement with the sender's reference. Their absolute values sum to the negativity. It is zero for a separable output and one for a maximally entangled normalized three-level state.

For the more restrictive departed-ground outcome the amplitude matrix is $K_g=\sum_a n(a)^*K_a$. If $K_g$ is invertible, the receiver can apply the contraction $\sigma_{\min}(K_g)K_g^{-1}$ as a heralded filter. Its product with $K_g$ is a scalar times the identity, including for entangled inputs. After ideal source selection, the probability of both the departed-ground and reversal heralds is $\sigma_{\min}(K_g)^2$. This is the largest uniform probability of exact reversal for that specific branch. The extra instrument is calibrated to the whole transfer map.

Leaving the bridge off gives exactly zero charge arrival. Opening only the bridge for $\pi/8$ swaps the gateway field words. The arriving charge sees the same uniform gateway field for every source, while the departed cell retains the source-dependent field. This moving-charge control gives a separable receiver-reference state. The first derivative of full neutral-ground transfer with the two complete cells also has rank one, since its matrix is $-4 a a^T$ with $a_j=\langle g,n|u_j\rangle$. The calculation at finite native times includes every internal charged and neutral transition on both sides. It compares against the archived two-input records at the two overlapping times.

## Reproduce

From the root of the reverse-engineering-reality checkout use Python with the repository requirements installed. The output path must be unused.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 \
python -B -W error model/candidate_laws/interval_echo/openwave/electron_field/lepton_family/transfer/measure.py --output /tmp/echo-three-level-transfer.json
```

The executable verifies the relevant archived constructors and two-input receipt against the electron-field SHA-256 manifest. All transport uses the complete joint state. The [protocol](DESIGN.md), [source](measure.py) and [numerical records](results.json) specify the sampling, operators and outcomes.
