# A distant dimer removes ordinary energy from the already connected material

Prepare the ground of the four-dimer Interval Echo material with its contacts already open. A pulse on the two source dimers exchanges their separately calibrated highest and lowest levels. It raises the mean energy by 33.1281 native units. The last dimer has no active contact with either source dimer. Its two orientation words have the same calibrated electric charge. The ordinary free generator, the source connection and the cold material all stay fixed during the flight. The pulse deliberately prepares an excited parent. [Frozen experiment](PROTOCOL.md) · [Native material, calibrated levels and charge](../passive_decay/README.md)

| Four million simulated launches per entry | Excited source | Idle source | Cut between core and receiver |
| --- | ---: | ---: | ---: |
| Local far exchange-minus at launch | 446,480 | 446,684 | 447,418 |
| Local far exchange-minus at wait 16 | 934,619 | 446,911 | 446,089 |
| Far probe records positive energy loss, aperture starts at 0 | 291,182 | 355 | 301 |
| Far probe records positive energy loss, aperture starts at 16 | 291,568 | 349 | 320 |

Each bank contains exactly 4,194,304 starts. The source pulse acts only on the core, so the reduced far-dimer preparation at launch is exactly that of the corresponding ground. Its subsequent local exchange-parity response at wait 16 is 22.28%, compared with 10.66% on the intact cold control. In the cut material the source and receiver generators factor. A pulse on the source then has zero effect on any far record. The source/far joint histograms and both signs of the frequency observation are retained. [All records](records.json)

**How the far probe reads the sign of an energy transfer.** The probe flips the local orientation phase at the last dimer at one of 16,384 coherent delay labels. Every branch has the same total duration. Write δ for initial minus final ordinary energy and A for that single local operation. Its transition amplitude from level i to level j in Fourier bin k is the matrix element of A times the mean of exp(2πin(k/16384 − δ/128)) over all labels n, apart from a common final propagation. A positive bin labels energy removed from the material. The complete band from 0.25 to 2 lies below the independently fixed light threshold of 3.36060 energy units. The full level span is smaller than half the Fourier period of 128, so the sign is unambiguous up to the calculated finite-window leakage. The long absorption trial spans about 804 native time units. The separate local trials resolve when the far signal is present. [Acquisition and response matrix](measure.py)

Delete from the local probe every matrix element that lowers the ordinary energy by any amount between zero and 3.36060. Stack the remaining, windowed matrices of all accepted bins and take the largest eigenvalue of their squared norm operator. For every normalized incident amplitude the resulting light-free response is at most 0.00031117 on the intact material. The measured loss fraction is 0.06952 at the second start, including every original launch. With a simultaneous 99% Hoeffding allowance of 0.00113 per probability, both starts exclude the light-free response and the cold and cut controls. The positive read therefore requires accessible low-energy transitions at a disjoint receiver. The opposite energy-gain band had 2,424,762 counts at the second start. The loss result describes that resolved outcome. Its value is an inclusive absorption probability, which cannot assign one exact energy to an individual time-resolved arrival. The material geometry, cold preparation, level pulse and coherent time reader are supplied. The short material permits reversible energy exchange. [Direct time-domain reader test](test_joined.py)

Run the experiment with one numerical thread from the repository root.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.candidate_laws.interval_echo.openwave.row_09_neutral.joined_decay.measure
~~~
