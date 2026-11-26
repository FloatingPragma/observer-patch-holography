# Which of three launch positions does a distant observer distinguish?

Owner 13 can place the same sharp vacancy–excess pair on any of three consecutive local bonds. On the clean 14-port guide a distant signed-contrast detector distinguishes all three launch positions in repeated trials. One common spatial response predicts every measured probability when the launch or the detector moves. This control measures how a choice of source position can change the apparent composition of a neutral signal.

At bond 8 in owner 17, the registered observer distinguished the central launch from either neighbour. On the native law the two neighbouring launches differed by at most 0.000017 in the observed signed frequencies across the four appointments. At half contrast this difference was at most 0.000177. Their two possible flight separations are +6 and −6 around the ring. Reflecting the ring about the launch bond and exchanging the two constituent signs fixes the source and preserves the two signed outcome labels. The integer native Hamiltonian commutes with that permutation at both tested strengths. The two distributions there are equal at all times. [First protocol](PROTOCOL.md) · [Complete simulated records](counts.json)

A separately enrolled detector one bond away, also in owner 17, distinguished each pair of source positions. The largest absolute signed-bin differences across four appointments for source pairs (0,1), (0,2) and (1,2) were respectively 0.32901, 0.27559 and 0.05342 on the native law. The simultaneous 99% counting half-width for any two-bank difference is 0.001242. Independent calibrations made with source 0 predicted all 36 held probabilities at this fixed detector and all 36 when it moved together with the source. The comparison at strength 4 also passed all 72 predictions and distinguished all three sources. No opposite-charge coincidence occurred in the blocked and filled-host shams. The sealed comparisons use both ordered-pair outcomes and their complement. They retain misses. There are 120 banks of 2²⁴ simulated starts in each experiment. [Second protocol](SECOND_RECEIVER.md) · [Independent repeat](second_counts.json)

Run from the repository root with one BLAS thread. The final receiver decodes the integer records; the optional acquisition takes a fresh output path.

~~~sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_10_oscillations.translated_pairs
python -Bm "$p.measure" --receiver 8
python -Bm "$p.measure" --receiver 7
python -Bm "$p.measure" --receiver 7 --acquire /tmp/second-neutral-position.json
python -Bm pytest -q --pyargs "$p.test_positions"
~~~
