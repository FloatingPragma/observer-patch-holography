# Can the observed three-label decisions be read as a lossless flavour transition matrix?

A complete three-state unitary transition matrix sends equally populated independent inputs to equally populated outputs. The three material winding recipes and their classified local/distant coincidences violate that necessary condition on **all six** enrolled receiver appointments. The largest supplied output population at owner 12, time 256 was 2.67506 times the equal-population value, with simultaneous 99% limits [2.64148, 2.70665]. This identifies a property of the classifier and its selected material clicks. The production and detector response needed to infer a latent mixing matrix from them is a separate measurement. [Three-recipe experiment](../joint_receivers/THREE_SOURCES.md) · [Original counts](../joint_receivers/blind_counts.json) · [All 18 intervals](columns.json)

## Finite identity

Let each source prepare one member of an orthonormal three-state basis, let every source use the same trace-preserving unitary flight V, and let a receiver resolve a complete orthonormal triple in the output space. Counting each source trial with the same source- and channel-independent acceptance gives

$$
P_{ab}=|\langle b|V|a\rangle|^2,\qquad
\sum_a P_{ab}
=\langle b|V\Big(\sum_a|a\rangle\langle a|\Big)V^\dagger|b\rangle=1.
$$

Row sums also equal one. A common mixture of such flights or reference settings preserves both identities by linearity. These hypotheses cover an ideal three-flavour vacuum oscillation probability at a common baseline and energy, with unitary PMNS mixing. Unequal spectra, detector efficiencies, tagging errors or source-dependent selection can change the column sums of the reported decisions. A column-sum test supplies only a necessary condition even for ideal lossless records.

## Independent reconstruction

This is a retrospective audit of the simulated trial banks at source commit `d810f85ac`. It uses the registered three-recipe classifier, including its equal source priors conditional on a received pair. For each owner, wait and phase setting, the independent decoder chooses for each material exit the source with the largest training frequency among received pairs. Exact rational arithmetic resolves the ordering and the smallest source index wins ties. Validation counts supply the three conditional probabilities per actual source after labelling. Average the two phase settings with equal weight within each source, then sum over sources. Reporting the same label for several distinct material exits is permitted by the instrument.

For a fixed trained decision and a nonempty validation pair bank, the number assigned to a chosen label is binomial conditional on the pair count. Exact two-sided binomial intervals with 4,096 Bonferroni allocations cover the 108 relevant source/setting/label probabilities jointly with at least 99% confidence in this declared family. A panel assigning either none or all of its exits to one label has a known conditional probability of zero or one. Summing interval endpoints bounds each column of the validation decision matrix. The audit uses every receiver appointment in the source acquisition. It has no prospective holdout for selecting this particular diagnostic. [Independent decoder](column_sums.py)

| Receiving owner | Free wait | Output 0, 99% interval | Output +, 99% interval | Output −, 99% interval |
| ---: | ---: | ---: | ---: | ---: |
| 12 | 256 | [2.64148, 2.70665] | [0.14670, 0.17928] | [0.14665, 0.17924] |
| 12 | 512 | [1.36679, 1.38891] | [1.47496, 1.51927] | [0.11394, 0.13613] |
| 12 | 2048 | [2.66694, 2.69638] | [0.15146, 0.16616] | [0.15215, 0.16690] |
| 17 | 256 | [1.41482, 1.43030] | [1.41684, 1.43218] | [0.13752, 0.16834] |
| 17 | 512 | [1.25684, 1.28913] | [0.42387, 0.48458] | [1.25801, 1.28839] |
| 17 | 2048 | [1.24091, 1.27034] | [0.22948, 0.25893] | [1.47375, 1.52782] |

For column totals sᵦ, the total variation distance of the output label distribution from uniform at equal received-pair source weights is Σᵦ|sᵦ−1|/6. Taking the smallest absolute deviation in each certified interval gives lower bounds from **0.17171 to 0.55565** across the six settings. Source labels, conditional receipt and the trained decision all enter this quantity. The distant neutral receiver alone has identical probabilities on the opposite winding pair. [Receiver theorem](../theory/PROOF.md)

From the repository root, run with a single numerical thread.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_10_oscillations.kernel_check
python -Bm "$p.column_sums"
python -Bm pytest -q --pyargs "$p.test_columns"
~~~
