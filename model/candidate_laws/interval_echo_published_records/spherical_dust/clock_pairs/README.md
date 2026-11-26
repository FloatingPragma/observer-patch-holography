# Two local clock readings on a common material pair

At each endpoint an observer reads one fair bit at the local clock appointment. The entropy of one bit is \(\log 2\) on every stock. Counting both bits measures a change in their correlations. This declared material channel uses shared classical randomness deposited before the clock measurements. Its endpoint law accepts the local conformal-time stamp, comoving pair separation and the local scale reading. It can therefore serve different observers passing the same material. [Finite source](channel.py) · [Independent word counts](experiment.py)

## One finite experiment specified before the measurements

Write \(d>0\) for the pair separation, \(b=2^{-k}\) for the lag spacing and \(S\) for the number of subphases. The endpoint with stamp \(\eta_i\) records at fine tick \(t_i=\lfloor S\eta_i/(bd)\rfloor\). A shared uniform integer \(z\in\{0,\ldots,S-1\}\) fixes its table index \(j_i=\lfloor(t_i+z)/S\rfloor\). The material chooses its class before either appointment.

* Half the copies have one fixed, common fair sign. One quarter have independent fair signs at every distinct local clock ticket and endpoint.
* In the remaining quarter, most copies use independent sign tables. A copy in lag class \(m\) uses entries \(H_{j_A}\) and \(-\operatorname{sgn}(w) H_{j_B-m}\) from the same table of independent fair signs. Each end multiplies its table sign by an independent local sign of mean \(\sqrt{a_i/a_{\max}}\), rounded to the finite dyadic grid. The independent-table class uses the same local gates. Here \(a_i\) is the clock scale, \(a_{\max}\) its registered bound, and \(w\) the material's signed entropy weight.

Each class contributes its unconditional probability. For \(|m|\le J<1/b\) the lag-class probability is \(|w|a_{\max}F(mb)/\kappa\) up to dyadic rounding, with \(F(x)=1/(1-x^2)\) and \(\kappa=\log(3)/2\); all other offsets have zero weight. The code checks that the lag slots fit in the available quarter and increases the integer ticket precision for smaller weights. At zero weight the two-bit vacuum has correlation \(c_0=1/2\). Given any finite clock slab, every required sign table and gate coin can be allocated in advance. The implementation assigns independent entries on first access and retains each value for repeated queries. Distinct pairs draw independently. A local record is diagonal in its two sign states. A local replacement may discard arbitrary initial endpoint memories and use these independent material tapes while preserving the marginal state of an outside reference.

One endpoint's complete sign history, for the same local clocks and scales, has the same law for every weight. Within its table class an index shift or sign reversal preserves the law of an independent fair table. The weight only moves probability among these classes. Each possible history at \(N\) distinct local addresses (endpoint, fine tick and scale ticket) has probability at least \(2^{-N}/4\) from the independent fine-tick class. In particular the vacuum and each sourced density matrix have full support on the current two-bit register. A change of the second appointment requires no operation at the first endpoint.

## Signed entropy, error and the boost coefficient

The two-bit probabilities are \(P(a,b\mid c)=(1+cab)/4\) for \(a,b=\pm1\). Given a lag class, its table entries correlate precisely if \(j_B-j_A=m\). If \(t_B-t_A=S q+s\), \(0\le s<S\), this selects offsets \(q\) and \(q+1\) with probabilities \(1-s/S\) and \(s/S\). This gives the exact joint law for the discrete source. On any compact lag interval \(|x|\le\beta_*<1\) with one spare bin and \(x=(\eta_B-\eta_A)/d\), each sampled bin differs from \(x\) by at most \(b(1+1/S)\). The derivative of \(F\) is bounded on the padded interval. With \(D\) tickets, the class rounding error is at most \((\sup F+2)/D\); the two gate means each have error at most \(1/D\). The code chooses \(D^{-1}=O(|w|b^2)\) at fixed positive scale envelope. Consequently

\[
\delta c=-w\sqrt{a_Aa_B} F(x)/\kappa+O(|w|b).
\]

This estimate is uniform in absolute appointment and lag on that interval. An uncertainty \(e\) in either conformal stamp adds \(O(|w|e/d)\). The proof applies to any fixed finite boost; approaching the light cone requires increasing \(Jb\) toward one and decreasing \(|w|\) enough to fit the lag slots.

For the same regional register and its unpulsed reference,

\[
 \Delta S=-\kappa\,\delta c-D_{\rm KL}(P_{c_0+\delta c}\Vert P_{c_0}),\qquad
 0\le D_{\rm KL}\le \delta c^2/(1-c_0^2).
\]

The upper bound follows from the \(\chi^2\) divergence of the two distributions. Only a region holding both ends receives the entropy signal. Entropy and divergence add over independent pairs. The leading response per whole pair is \(w\sqrt{a_Aa_B}F(x)+O(|w|b+w^2)\). This compares the same public slot algebra on both arms of the experiment. Comparing raw entropies on two different sets of slots introduces changes from the reference entropy of those slots. The public register here contains one current bit per included endpoint, conditional on the clock and stock preparation. A retained record of a longer history has a different entropy.

In the rest coordinates of a tangent dust chart with scale one, a ball moving at speed \(\beta\) projects to the ellipsoid with principal radii \((\ell,\ell,\gamma\ell)\). For separation in direction \(n\) making angle cosine \(\mu\) with the velocity, integration over positive separation \(d\) of the two-endpoint overlap volume gives \(\pi\ell^4\gamma/\sqrt{1-\beta^2\mu^2}\). The endpoint time difference on the tilted slice is \(d\beta\mu\). Isotropic averaging of the channel's factor therefore gives

\[
 \gamma\int_0^1\frac{d\mu}{(1-\beta^2\mu^2)^{3/2}}
 =\gamma\left[\frac{\mu}{\sqrt{1-\beta^2\mu^2}}\right]_0^1
 =\gamma^2.
\]

This calculation uses the specified source in a tangent chart. Curved proper balls, cofinal spatial coverage, the common reference mapping and independent area and dust-stress readouts enter the regional theorem separately.

## Readings

The independent counter enumerates the issued ticket, subphase, shared or private sign entries and both local gate coins. Each of four synthetic count banks contains \(2^{25}\) complete reads. At the enrolled lag the vacuum, full material, locally attenuated material and independent-table control give correlations \(0.5,\ 0.48997132,\ 0.49498566,\ 0.5\). The simultaneous 99% half-width on each terminal word frequency is \(0.000347\). For the fastest held tangent survey the 4-bit and 8-bit readings give 2.92219 and 2.90948 against \(\gamma^2=2.90951\). A constant lag response gives 1.98837. These surveys use finite angular quadrature. The integral above fixes the limiting coefficient analytically. [Counts and boost results](results.json) · [Checks including both weight signs and faster lags](test_channel.py)

~~~sh
python -m model.candidate_laws.interval_echo_published_records.spherical_dust.clock_pairs.experiment
python -m pytest -q model/candidate_laws/interval_echo_published_records/spherical_dust/clock_pairs
~~~
