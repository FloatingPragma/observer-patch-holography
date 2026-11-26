# Interval Echo with partition repair

A change to the local echo law makes the entropy of overlapping observer records stationary at fixed recorded volume. Each exposed endpoint has either one or two possible partners. That integer determines the entropy response to an open aperture, including for perturbations that correlate different faces. An ordinary contrast write can open the aperture, so this record is accessible through a native echo experiment. Interval Echo is the main candidate. Partition repair is an experimental alternative.

On one four-cell cycle, summing all 10,000 global words gives one bit of entropy loss per added open aperture at fixed volume. The same calibration predicts a loss of 8/5 bits per unit of a held perturbation involving two faces. A preparation that biases the two open endpoint symbols relative to one another fails the prediction even at fixed area. [Exact transition, regional restrictions and controls](REGIONAL.md)

The geometric comparisons use independent measurements of the boundary. On a cubical survey, rods and word counts satisfy a fixed-volume first law on overlapping regions. On a supplied smooth survey with ten moving or resting frames, a further rule ties the number of possible endpoint partners to proper area. It enforces an error of order the fifth power of radius for the first variation, smaller than the fourth-order curvature signal. Both the smooth geometry and its effect on partner counts enter this correction explicitly. Recovering either from unmodified Interval Echo and independently identifying regional stress remain open. [Law and equivalence criterion](SPEC.md) · [Cubical rods](surveyed_cells/README.md) · [Volume-matched moving observers](proper_balls/README.md)

Run the checks from the repository root.

~~~sh
python -m model run --law interval_echo_partition_repair --check
python -m model run --law interval_echo_partition_repair --full
python -m pytest -q model/candidate_laws/interval_echo_partition_repair
~~~
