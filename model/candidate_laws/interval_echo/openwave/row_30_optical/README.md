# A matter oscillator on the blank sphere

One load, a stationary uniform reference and one gateway detector suffice to measure a rest beat on the intact board. The estimated angular frequencies on 20, 80 and 320 cells are 2.992, 2.999 and 3.010 per local clock radian. Doubling the internal exchange gives 6.000, 5.997 and 6.010. The local optical profile has both signs of amplitude within a cell. Every occupied port carries the same type of load. [Preparation and apparatus](PROTOCOL.md) · [Complete detector records](results.json) · [Counts-only frequency decoder](readout.py)

| Cells | Rest frequency, counting-only 95% bootstrap interval | Dimension of constructed stationary space | Stationary weight at the first positive wave number | Localized loads two or more cells away at exposure 32 |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 2.992 [2.984, 3.001] | 70 | 93.20% | 16.30% |
| 80 | 2.999 [2.983, 3.014] | 280 | 97.65% | 22.40% |
| 320 | 3.010 [2.978, 3.041] | 1,120 | 99.35% | 23.30% |

Four relative source phases measure the phase of the local beat. Their two detector quadratures give its unwrapped angle. On the short calibration appointments the dimer hand stays within its first half cycle, so its terminal port frequency gives elapsed time unambiguously. The slope of the beat against these inferred times gives the table. The longer flights and return experiments use specified apparatus exposures. The clock and the signal use independent material copies. Absolute energy and a particle mass in laboratory units are outside this measurement.

## Why most of the longest wave stays put

In the one-load space the uniform reference subtracts the common energy 48. An isolated cell acts by half its port Laplacian. On even antipodal amplitudes of zero sum the operator is exactly 3 times the identity. Each port has one neighbor on each of the five other axes. Give every axis that contains a gateway amplitude +1 and balance the private axes with a common negative amplitude. Opposite sides of every shared link then agree on this wave at zero spatial wave number.

The board generator is the sum of the cell generators and a projector P onto the antisymmetric port combination of each shared link. These links are disjoint. Let a be a spatial wave of the cell graph with its vertex weight equal to the squared norm of that cell's internal profile. For the corresponding normalized load wave, the graph Rayleigh quotient is k². The expected value of P is k²/2, measurable by the antisymmetric port read on a randomly chosen intercell seam. The expected energy of the wave is 3 + k²/2. Its exact energy variance is (k²/2)(1−k²/2). The strongest response at a given spatial wave number can therefore include several temporal frequencies.

The subspace that agrees across each gateway and sums to zero in each cell has frequency 3. Its dimension follows by equating amplitudes along the gateway links on the six-axis graph and imposing one independent zero-sum constraint per cell. Every cell has private axes, which makes those constraints independent. Orthogonal projection onto this subspace gives a certified stationary weight. An inverse source preparation followed by a home-port read records departures from the initial wave at every positive wave number in the test. On the 320-cell board the five spatial waves retain at least 97.88% stationary weight. A single localized source excites a broader set of modes and reaches distant detectors. This transport alone cannot attach the rest frequency to a dispersing branch. [Native source and free run](experiment.py)

## Directly preparing the complementary wave

A follow-up removes the stationary component using the same spatial constraints, normalizes the remaining source, and prepares it with local splitters. This changes the source and its wave-number distribution. Its squared overlap with the original spatial wave and an independent bond count are recorded. The time meter recombines 128 controlled release times on a Fourier register against the uniform reference. Its period 8 covers the complete relative spectrum in [0,6]. Every outcome is counted. The frequency resolution and test were chosen after exploring the component's spectrum. [Complementary preparations and timing instrument](complement.py) · [All fifteen timing histograms](complement.json)

Any single frequency places at least 8/π² of its probability in the two nearest bins of this meter. Across all fifteen sources the largest observed pair share is 71.7%. Even a simultaneous 99% Hoeffding allowance over 15 × 128 pairs adds less than 0.1 percentage point. Removing the stationary wave therefore gives another resolved mixture. Narrower preparations within that mixture require their own spatial and temporal measurements.

From the checkout root run the following commands using the model environment.

```sh
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_30_optical.experiment
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_30_optical.readout
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_30_optical.complement
```
