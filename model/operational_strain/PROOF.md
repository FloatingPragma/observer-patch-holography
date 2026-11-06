# Finite work identification and native sector weights

## Operational contract

An energy record supplies a positive volume setting, an energy reading, an
absolute error bound and a held-state identifier. The identifier expresses
a protocol requirement: all settings concern the same canonical state. It
does not certify that requirement physically. Volume settings are exact in
the mathematical interface; volume calibration error, imperfect holding,
and meter error require bounds before applying the result. Energy error
may include an independently bounded state-holding contribution.

The intervention, energy meter, chosen response class and work convention
`p=-partial_V H` are supplies. The certificate never labels generated strain
examples as native observations. The scalar source supplies its verified
fixed-volume histories and a value of `Q`; it supplies no strain apparatus.
The physical scale relating `Q` to an energy amplitude is separate.

## Three settings identify the power response

Assume, over the volume settings and near the reference volume,

\[
 H(V)=C+A(V/V_0)^{-\alpha},\qquad A>0,\quad\alpha\ge0,\quad V_0>0.
\]

The constant `C` is independent of volume and the same across all records.
The family is positive on all positive volumes if `C>=0`; local pressure
identification itself requires no sign for `C`. A positive density or its
ratio additionally requires positive total reference energy.

For a chosen `r>1`, take three energies `E_j=H(V_0 r^j)`, `j=0,1,2`.
Writing `t=r^-alpha` gives `E_j=C+A t^j`. If `alpha>0`, then

\[
 D_0=E_0-E_1=A(1-t)>0,\quad D_1=E_1-E_2=tD_0,
\]
\[
 t=D_1/D_0,\quad \alpha=-\log(t)/\log(r),\quad
 A=D_0/(1-t),\quad C=E_0-A.
\]

Consequently `p(V_0)=alpha A/V_0`. For `E_0>0`,
`p/rho=alpha A/E_0`. An additive energy offset changes this ratio even
though it changes neither pressure nor the differences `D_j`. The recovered
`C` is the asymptote of the stipulated family; this does not select an
absolute physical vacuum energy outside that family.

Three settings are generically minimal when both amplitude and offset are
unknown. For any two energies `E_0>E_1` and any `0<t<1`, the assignments
`A=(E_0-E_1)/(1-t)` and `C=E_0-A` reproduce both readings. For
`E_0=144, E_1=90, r=64, V_0=3`, contractions `1/4` and `1/8` give,
respectively, `(A,C)=(72,72)` and `(432/7,576/7)`. Both have positive
energy at every positive volume. Their reference pressures are `8` and
`72/7`, while their third energies are `153/2` and `333/4`. If either
the amplitude or the offset is independently known, fewer settings can
suffice; minimality here concerns three unknown parameters.

When `alpha=0`, all three energies coincide. Given the assumption `A>0`,
exact equality identifies zero pressure and `alpha=0`, but does not
separate amplitude from offset. Nonzero reading errors around constant
data receive an unresolved result, without an inferred zero derivative.
If amplitude zero were permitted, constant data would not identify alpha.

The implementation performs outward rational interval operations on the
three energy intervals. It divides only after proving both decrements
positive and the contraction strictly between zero and one. Monotonicity
of the logarithm then encloses alpha, amplitude, offset and pressure.
These are guaranteed enclosures **if the family and record error bounds
hold**; interval dependency can make them conservative. They are not a
test that three data points obey that family at other volumes.

For `1<=x<=2`, let `z=(x-1)/(x+1)`. The exact lower sum
`2 sum_{k=0}^{N-1} z^(2k+1)/(2k+1)` encloses `log x` with remainder at
most `2 z^(2N+1)/((2N+1)(1-z^2))`. This follows by bounding every omitted
denominator below by `2N+1` and summing the geometric series. Range
reduction by powers of two and `log(1/x)=-log x` cover all positive
rational arguments. No floating transcendental evaluation certifies the
intervals.

The source example independently replays all 240 retained rational scalar
updates through `code/equation_of_state/native/verify_native_eos.py`.
It uses the initial `Q=72` as a declared amplitude in two added
calibration examples, with `C=72`, `r=64`, and exponents `1/3,1/2`.
Each reading has error at most `2^-20`. Their exponent and pressure
intervals are disjoint. The scalar histories are native; the deformation
laws and energy readings in these examples are supplied mathematical
interventions. No simulator source checkout is needed for this independent
arithmetic replay, and its external source bytes are not checked here.

## General response with a regularity bound

Let `H` have two derivatives on `[V_0-h,V_0+h]`, with `h>0` and the
interval contained in positive volumes. Suppose `|H''|<=L` there. By the
fundamental theorem of calculus,

\[
 \left|\frac{H(V_0+h)-H(V_0-h)}{2h}-H'(V_0)\right|
 \le\frac1{2h}\int_{-h}^{h}L|s|\,ds=\frac{Lh}{2}.
\]

If instead a third-derivative bound `|H'''|<=M` is supplied, the symmetric
Taylor remainder improves the bound to `Mh^2/6`. Endpoint reading errors
`epsilon_-` and `epsilon_+` add
`(epsilon_-+epsilon_+)/(2h)`. Negating this interval encloses pressure.
The offset cancels exactly, so no energy zero is needed for pressure.
An energy density or pressure/density ratio requires a separate reference
energy. The derivative bound applies throughout the interval; checking
the sampled points does not establish it.

Without a response class or quantitative regularity, even positive smooth
energy functions defeat every finite stencil. For arbitrary finite sampled
positive volumes `V_i`, set

\[
 u(V)=(V-V_0)\prod_{V_i\ne V_0}\frac{V-V_i}{V_0-V_i}.
\]

Duplicate settings are used once. Then `u(V_0)=u(V_i)=0` and
`u'(V_0)=1`. For any real `d` and `B>0`,
`H_d(V)=B exp(d u(V)/B)` is positive and smooth on positive volumes,
has the same energy `B` at every sampled setting, and satisfies
`H_d'(V_0)=d`. Thus finite exact energy samples alone give no finite
bound on pressure. The certificate checks the rational polynomial
identities; the displayed exponential proves positivity analytically.

## Native readout of calibrated sector forces

The binary total-two carrier has 66 states and core charge
`q=n_0+n_1`, with multiplicities `(45,20,1)`. Reuse the existing E
instrument census on seams `(0,1)` and `(0,5)`, retaining both fair-coin
values, contrasts and disturbed endpoints. Its decoder reads only the
two-event tape and reconstructs the initial core charge. This is the
zero-prefix specialization of
`model/native_equilibration/history_entropy.py`.

For every complete tape `y`, the Kraus branch is one half times an
injective partial permutation. The exact execution checks injectivity
and all 264 input branches, as well as the decoder for every branch.
Its effect is diagonal. Summing the effects with decoded charge `q`
gives exactly the initial charge projector `P_q`. Thus, for any source
state, including one entangled with an untouched reference, this native
measurement has outcome weights `Tr(P_q rho)`. All outcomes are retained;
no postselection or terminal-state restoration is asserted. The E
selective-instrument convention and binary total-two carrier are premises.

Suppose independently calibrated sector energies specify a diagonal
strain Hamiltonian `H(V)=sum_q E_q(V) P_q`, with the projectors fixed
under the declared variation. If `-E'_q(V_0)` belongs to an interval
`I_q`, then

\[
 -\operatorname{Tr}(\rho H'(V_0))\in
 \sum_{q=0}^2\operatorname{Tr}(P_q\rho)I_q.
\]

The native charge outcome is therefore a bounded per-run estimator of
the calibrated generalized force, with the displayed expectation for
arbitrary coherent inputs. The certificate's probabilities are exact
model probabilities, not empirical frequencies from a finite sample.
Experimental frequency estimates require their own sampling error.

The two native moments in `common_variation.py` give
`p=(1-m+c,m-2c,c)`. Its global reference has
`p=(17/25,23/75,1/75)`. The calibration example uses three supplied
power responses, requiring nine energy settings, and combines their
pressure intervals with those exact native weights. This supplies an
operational readout of a specified charge-diagonal work operator. It
does not select its strain coefficients, make the conserved moments
physical energy, or turn changes of discrete cap into volume variations.

The Whitney Maxwell action has a supplied metric and canonical variation;
its conformal pressure relation is separate. This package derives no
Maxwell pressure formula, thermal equilibrium, total apparatus energy,
continuum limit or physical equation of state. The positive result
replaces a supplied derivative by finite calibrated energy data within
the stated response or regularity class.
