# Low interacting energy windows for a signed Wilson spinor transfer

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

On the unchanged interacting Wilson cube, the lowest exact energy with nonzero short-string sender-to-receiver amplitude lies in **[−727.5, −648.375]** in the native clock units, for either sign. A separate clock filter rejects the interacting ground and still permits a full-rank two-spin transfer with the signed charge, holder and written field checked at the next cell. The filter localizes energy to a numbered window. No individual spectral line is resolved. [Exact receipt](results.json) · [Source and prior receiver measurements](../README.md) · [Earlier common-energy proof](../ENERGY_BRIDGE.md)

## Complete receiver records

Use the admitted beta stock, signed native CAR pair and original idle stock reference of the preceding study. The moving excitation starts at $(1,0,0)$, its conjugate holder is $(1,0,1)$, and the next-cell receiver is $(0,0,0)$. The source spin may be entangled with an independent idle Bell reference. The clock acceptance requires the success port and return of every clock and arithmetic register to blank. After that read, hold under the original $H$ for a separately fixed native wait. The receiver checks the moved charge, unchanged conjugate holder, expected link words and spectator occupations, and the Bell spin match without an earlier spin-component read. No term of $H$ is shuttered.

The ground bracket comes from the archived full physical occupation bound and variational upper bound on this same member, without using the beta stock energy as the ground. Rounded outward, it gives $E_0\in[-727.5,-723.75]$ at the archived energy zero. The narrow existence bound and the runnable clock plan serve different purposes. The tapered clock window is wider to leave a finite guard between the energy-dependent acceptance and rejection. [Ground enclosure][ground] · [Native clock and carrier][contract]

| Exact window or held record, either charge sign | Short string | Detour |
| --- | ---: | ---: |
| Lowest energy with a nonzero spin-preserving sender-to-receiver amplitude | [−727.5, −648.375] | [−727.5, −632.625] |
| Tapered low clock window | [−727.5, −640.375] | [−727.5, −624.625] |
| Low-window clock and full receiver probability at hold 2⁻¹⁸, lower bound | 3.00832 × 10⁻⁷ | 4.11151 × 10⁻⁷ |
| Low-window transfer minimum singular value, lower bound | 0.00054848 | 0.00064121 |
| Exact energy outside low window and accepted clock, upper bound | 3.64232 × 10⁻¹⁰ | 2.96115 × 10⁻¹⁰ |
| Tapered excited-only clock window | [−723.5, −462.375] | [−723.5, −446.625] |
| Excited-window clock and full receiver probability at hold 2⁻²¹, lower bound | 1.02047 × 10⁻⁷ | 9.81145 × 10⁻⁸ |
| Excited-window transfer minimum singular value, lower bound | 0.00031944 | 0.00031323 |
| Exact energy outside excited window and accepted clock, upper bound | 5.01436 × 10⁻⁹ | 5.01130 × 10⁻⁹ |

Subtracting the certified ground bracket puts the excited-only windows at **[0.25, 265.125]** and **[0.25, 280.875]** in ground-relative source units. Both ideal filters exclude the ground. The leakage bounds refer to a sharp energy check after clock acceptance and before the terminal position-and-spin check, on separate copies. The transfer probabilities already propagate the amplitude error from the approximate clock. Each event probability uses the accepted-stock denominator, with all clock and receiver refusals included. Multiply by the original stock admission $s\ge 2^{-1{,}073{,}742{,}520}$ to express it per original primary-polynomial start. The admitted full-spin channel and the bounds apply to either sign. The archived charge-conjugation unitary also gives exact equality of the signed comparison at each energy when the holder and the two spin columns are conjugated as specified in the earlier bridge. [Conjugate spin and holder comparison](../ENERGY_BRIDGE.md#same-energies-on-the-exact-conjugate-preparation)

## Window witness

Let $V_i,V_f$ be the two isometries from an input spin qubit into the sender and expected receiver code. The signed spin transport fixes their relative CAR phases. Their charge addresses make $V_f^\dagger V_i=0$. For $c=(\mu_i+\mu_f)/2$, $d=|\mu_i-\mu_f|/2$ and $X=H-cI$, the exact mixed moments are

$$
V_f^\dagger X V_i=-\tfrac12 I_2,\qquad V_f^\dagger X^2 V_i=0.
$$

The second identity uses the two-hop parity and gauge-word selection of the native receiver calculation, whose diagonal contribution is $(\mu_i+\mu_f-2c)V_f^\dagger H V_i=0$. This argument also holds between different spin columns. The first moments at each end are the scalar $\mu_j I_2$. Each two-by-two beta block of the archived three Wilson hop matrices has at most one nonzero entry per row and column. The replay verifies these matrix patterns. The two spin basis columns at a given address differ in one occupied positive-beta mode or one negative-beta hole. A common Fock output from one hop on each column could occur only by removing the distinct occupied modes or filling the distinct vacant modes through the same other orbital. That would require two nonzero entries in one row or column of a single beta block. The magnetic term leaves the occupations unchanged and is orthogonal to the one-hop outputs. Thus each centered second-moment Gram matrix is diagonal, and its largest sharp-spin variance bounds its norm. [Exact full-law moment calculation](../README.md#why-the-finite-bands-follow)

For a real quadratic $r(x)$ with linear coefficient one, $V_f^\dagger r(X)V_i=-I_2/2$. Write $\theta(x)$ for a taper between zero and one. On its excluded tails assume

$$
|(1-\theta(x))r(x)|\le k_*(x^2-d^2),
\qquad |x|>d .
$$

Since $(E-\mu_i)(E-\mu_f)=x^2-d^2$ at $x=E-c$, factoring the two centered Hamiltonians on the left and right bounds the tail matrix in operator norm by $k_*\sqrt{v_i v_f}$, where $v_i,v_f$ are the maximum variances across the sharp-spin columns. If $\theta r$ has a fixed sign, normalize its absolute value by a rational bound $R$ to make a positive contraction $F$. Its transfer matrix differs in norm from a signed $I_2/(2R)$ by at most $k_*\sqrt{v_i v_f}/R$. This certifies both spin directions and a positive Bell overlap in a coherent window. It neither selects a rank-two individual eigenspace nor counts how many levels lie inside.

For the low window use $r(x)=x-b-k(x^2-d^2)$ at $b=48$ and $k=b/[4(2b^2-d^2)]$. On $x\ge b$, the maximum of $(x-b)/(x^2-d^2)$ is $1/[2(b+\sqrt{b^2-d^2})]$, at most $2k$. This makes $k_*=k$ and the polynomial is negative throughout the pass region. Taper linearly to zero across the two source energy units after $b$. Taking $b=42$ with a discontinuous cut already forces positive exact common energy support at or below $c+42$. This is the first table row.

For the excited-only window use $r(x)=(x-a)(1-k(x+a))$ with $a$ halfway between $-727.5-c$ and $-723-c$, and $k=(b-a)/(2b^2-a^2-d^2)$ at $b=224$. On $x\ge b$ the quotient $(x-a+k(a^2-d^2))/(x^2-d^2)$ decreases from $2k$ to zero, giving the high-tail bound $k$. On the interval from the ground floor through $E=-723$, the zero at $a$ and the positive minimum of $x^2-d^2$ give the separate low-tail bound. The clock starts its linear turn-on at $E=-723.5$ and completes it at $-723$. The polynomial stays positive inside the selected band, and the upper taper ends four units after $b$. Both tails are accounted for by the larger weighted bound. The observed threshold is still very broad. This construction proves spectral access strictly above the ground, without asserting where the lowest such excited line lies. [Rational intervals and taper parameters](results.json)

## Finite clock and hold allowance

Use the original absolute clock $U(t)=\exp(-itH)$. The source norm radius bounds the base tick by $\tau=1/(4\lVert H\rVert_{\mathrm{cap}})$, so there is no aliasing of real source energies. Reuse the archived Fourier-clock tail with energy tolerance $\eta=1/4096$. For scalar energy $E$, each controlled exact clock register produces a distribution of energy labels. Use nine coherent copies of its register, compute their median, and rotate a success ancilla by amplitude $F$ at that median. Uncompute the registers and accept success only if every work register returns blank. On an energy eigenvector the accepting amplitude is the expectation of $F$ in the median distribution. This stock-compatible operation is scalar also on degenerate energy eigenspaces.

The bound on the probability of five or more bad estimates is the binomial upper tail. The probability estimates are independent conditional on exact energy. Clipping the energy label below the ground floor to the ground floor gives a globally Lipschitz amplitude as a function of the displayed energy. The uniform Kraus error bound is the median tail plus $\eta$ times the checked Lipschitz bound plus a complete implementation allowance of $2^{-17}$. The plan needs nine 36-bit Fourier registers and 648 controlled native clock powers including inverses. The longest controlled native wait is about $8.43\times10^5$ source units. It also needs the stated arithmetic, rotation and counterclock precision. The replay evaluates the exact input moments, tail and probability bounds for this finite protocol. It does not execute the full interacting time evolution or the complete large clock circuit. [Fourier clock and alias guard][fourier] · [Clock precision budget](results.json)

After clock success, the ideal filtered amplitude is supported in $[l,u]$. At hold $h$ the difference from the zero-hold transfer, up to a calibrated common phase, is at most $h(u-l)/2$ in operator norm. Adding this error and the complete clock allowance leaves the positive minimum singular values and Bell probabilities in the table. The hold verifies a transferred component at a fixed later time with the same terminal gauge and spin checks. The clock itself uses $H$ and can contribute to the transfer. The bound assigns neither a speed nor the part of the amplitude acquired during the hold alone.

The detector comparison would fail if the calibrated spin-transfer matrix lost rank, the accepted receiver probability fell below the certified band, or leakage exceeded its allowance. The current finite proof establishes the shared interacting support and the finite receiver implication on this exact held compensator construction. Level-by-level energy resolution, an isolation gap, a freely propagating endpoint with calibrated momentum and a physical rest mass all remain open.

From the repository root:

~~~sh
.venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.theorem_interfaces.endpoint_recoil.low_branch.measure
~~~

[contract]: https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/particle_common_contract/PROOF.md
[ground]: https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/particle_primary_occupation_bound/PROOF.md
[fourier]: https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/particle_primary_joint_spectroscopy/PROOF.md
