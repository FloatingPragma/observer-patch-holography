# A signed sum of local analyzer readings

On the four-fan guide use the fresh local launch with the preparation clock bypassed. At terminal hub 1045 read the signed load and the face-0 pair with positive-chart words 334 and 365. The only other observer reads the untouched source reference. The same local procedure applies to a vacancy after complementing the words and calibrating the two local axes. This one selected pair can be measured by a dedicated detector on each trial. It does not require the multiplexed port and face detector.

The exact model forecasts in the [local face bank](word_reads.json) were known when choosing this hub, pair, wait and witness. They were selected to give an explicit test in terms of detector counts. Predictions at different waits and rates use the same witness. A future independent acquisition can freeze this procedure before observing any data.

Order each observer's axes as identity, X, Y, Z. For simultaneous accepted readings let each sign be +1 or −1. Use $t_{ij}$ for their mean product on the receiving axis $i$ and the reference axis $j$. An identity means ignoring that observer's sign. Evaluate $w=\sum_{ij} C_{ij}t_{ij}/32$ with

~~~text
       I   X   Y   Z
I      8   2   2   2
X      2  -2   4   6
Y     -2  -4  -6   2
Z      2   6  -2   4
~~~

For every separable pair this statistic is nonnegative. Indeed take $v=(2,1,1,-1-i)$ in the positive-chart receiving/reference convention. The matrix $C$ is $v^\dagger(\sigma_i\otimes\sigma_j^T)v$, and $v^\dagger v=8$. Thus the statistic is the expectation of the partial transpose on $v/\sqrt8$. For a product input the numerator is an absolute square, and convex mixtures preserve its sign. The same argument includes a classical message that tells the receiver which local state to prepare. A negative observed value with a resolved margin certifies the transfer of quantum information written at the source.

For a concrete count replay use nine Pauli settings, $2^{23}$ original preparations per setting and generator seed 981045. Record the four local sign pairs and refusal. Read one-body means from the diagonal settings. Use separate acquisitions for native flight at waits 27.0625 and 30.8125, for half-rate flight at the later wait and for a Z-dephased reference at that same native wait. No launch uses an intermediate station's report. Give a family-wide 99% two-sided Hoeffding enclosure conditional on the observed numbers of accepted trials, and include the numerical readout enclosure. These are prospective draws from the computed instrument laws for the selected witness. [Counting code](witness.py) · [Simulated records](witness.json)
