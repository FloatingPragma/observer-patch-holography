# One flux quantum connects two measured spin-active lines on four sites

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

For either signed charge on the existing interacting four-site contact ring, the lower bright level's translation-character-3 component continues in a fixed spin-½ block to the second measured bright level. The level is simple in that block before one complete supplied twist. The terminal energy has an orbital doublet at fixed spin projection. The continuous approach picks one direction in that doublet, with local same-ground charge-and-spin residue **0.08430170036**. The translation-character-1 copy of the initial energy reaches a different higher level under the same twist. [Exact receipt](results.json) · [Original spectra and records](../results.json)

## Instruments and endpoints

Use the archived contact $U=1$, original hopping $1/16$, original onsite spin frame and half-filled singlet. Spin-up creation and spin-down removal give opposite signed excess $s$ but both have $S_z=+1/2$. A Fourier read tags the translation eigenvalue $i^k$ on the four-site ring. The half-filled ground has $k=2$, so the detector's relative wave index is $k-2$ modulo 4. The same spin-independent schedule acts on the other spin component and preserves an idle spin reference.

Multiply each forward hop by $z=\exp(i\phi)$ and each reverse hop by $\overline z$, with $\phi$ increasing from zero to $\pi/2$. The contact is unchanged throughout. The two ends of this schedule are gauge equivalent. On a word with site occupations $n_x$, the exact gauge relabeling is $G=i^{\sum_x x n_x}$. With $N=4+s$, the calculation verifies $H(\pi/2)=GH(0)G^\dagger$ and $TG^\dagger=i^N G^\dagger T$. Returning with $G^\dagger$ shifts the recorded translation index by $s$ modulo 4. The ground's number four gives zero shift on its own reference. [Exact source and gauge replay](check.py)

| Signed excess | Initial $k$ | Returned $k$ | Terminal raw energy | Second pump-accessible line |
| --- | ---: | ---: | --- | --- |
| +1 | 3 | 0 | $3/2-\sqrt{21}/8$ | Yes |
| −1 | 3 | 2 | $1/2-\sqrt{21}/8$ | Yes |
| +1 | 1 | 2 | $3/2-\sqrt{13}/8$ | No |
| −1 | 1 | 0 | $1/2-\sqrt{13}/8$ | No |

The first two rows connect the archived reservoir-corrected energies 0.4175457854 and 0.4723570583. Each starting momentum component has local CAR residue 0.08363250802 on the normalized singlet. These are unscaled residues, before the earlier pump, ground-admission and any finite clock refusals. [Original admission convention](../PROTOCOL.md)

## No crossings before the terminal orbital doublet

Write $u=\tan(\phi/2)$ and $X=16H$. The exact signed six-dimensional momentum matrices have a one-dimensional image under spin raising. Dividing out that spin-$3/2$ factor leaves the degree-five characteristic polynomial of the complete spin-$1/2$ block in each tracked momentum. The replay verifies the intertwining coefficient by coefficient in $z$, including the contact and both orientations separately. For the monic polynomial the discriminant is

$$
\frac{2^{44}(1-u^2)^4}{(1+u^2)^{16}} B(u).
$$

The exact degree-24 polynomial $B$ has positive Bernstein coefficients on $[0,1]$. The verified coefficient floors are 24,879 for $k=3$ and 25,857 for $k=1$. Thus no lines within a tracked spin-and-translation block can cross before the terminal doublets. Every compact subinterval $0\le u\le u_*<1$ has a quantitative gap. Since the spectral width of $X$ is at most 36, each pairwise energy gap for $H$ is at least the square root of this discriminant divided by $16\cdot36^9$. The bound permits crossings between different spin or translation sectors and goes to zero at the endpoint.

Projecting the derivative of $X$ with respect to $\phi$ onto the final orbital doublet of the $k=3$ route gives slopes $\pm\sqrt{(40+8\sqrt{21})/21}$. The lower eigenvalue before the endpoint selects the positive derivative as $\phi$ approaches $\pi/2$ from below. After the gauge return the compressed derivative is purely imaginary in a real two-dimensional basis. A real local CAR column then divides its exact weight equally between the two complex-conjugate slope directions. The archived second-line residue is 0.16860340071 on the normalized singlet, so the selected direction carries half. The full endpoint still contains two orbital states per spin sign. The terminal slope reports a continuation of one particular component. [Exact slope and residue check](check.py)

These are static spectral-continuation and probe statements on one supplied control path. The replay does not certify a finite-rate ramp, a uniform endpoint gap, a physical electromagnetic operation or an autonomous relativistic shell.

From the RER root:

~~~sh
.venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.theorem_interfaces.interacting_chain.flux_continuation.check
~~~
