# Baryon and lepton labels at two material stations

Two Interval Echo stations distinguish the declared signed baryon and lepton labels of ten one-generation fields by interference. They agree on every signed charge and on the eight oriented field multisets that pass the dimension-six selection. A second interferometer measures the sign for each of the 24 orders of four distinct exterior generators. With these signs the four classes and their conjugates each support one independent contraction. Omitting the signs erases the contraction with fields Q, Q, u_R and e_R, and its conjugate. [Fixed protocol](PROTOCOL.md) · [All material readings and coefficient vectors](results.json)

## What is assigned to a field

The five-dimer register has three declared colour modes and two weak modes. For colour and weak occupations c and w the integer hypercharge weight is q = 6Y = −2c + 3w. The empty and all-full words are neutral reference lines. Excluding them leaves 30 internal words in ten components, each with a distinct q. The two stations choose the even nontrivial exterior degrees for left-Weyl labels. In the order Q, L, u_R, d_R, e_R and their Hermitian conjugates, the q values are 1, −3, 4, −2, −6, −1, 3, −4, 2, 6. The right-singlet names designate the conjugates of the corresponding left-Weyl singlets in the even sector. Field names, Weyl conventions and dimensional assignments are inputs to this construction. [Native exterior register](../README.md) · [Finite Weyl table and theorem](../../../../../Lean/Screen/BaryonDimensionSix.lean)

Three times baryon number is +1 on colour fundamentals, −1 on antifundamentals and 0 on singlets. This chooses the signed representative of colour triality and the conventional unit. L and e_R carry lepton number +1. Their conjugates carry −1. The other fields carry zero. In the native word chart these assignments read

$$
3B=(0,1,-1,0)_c,\qquad
L=\mathbf{1}_{c=3,\ w<2}-\mathbf{1}_{c=0,\ w>0}.
$$

Baryon and lepton labels are declared accidental charges. The gauge charges and exterior census alone permit a different sign or normalization of either assignment. The neutral phase origin and the three candidate integers −1, 0, +1 are part of the meter calibration.

## A local and a remote number reading

For each nonempty word x, cell 13 prepares the coherent sum of the empty word and x. A Hadamard on one marked vacancy and controlled flips onto the other selected dimers suffice. A reference dwell t on one marked rail gives relative phase exp(8it) to its vacancy branch. Any two marked rails can visit an adjacent pair of private-reference ports. The refocused native controlled phase there and its two inverse routes restore the temporary registers. Dimer Hadamards convert that phase to a controlled flip. The check resolves every ordered gate on the 32-dimensional code and records physical paths for all ten unordered pairs in both stations. [Local pulse and routing checks](apparatus.py)

For either f = 3B or f = L the sender programs the diagonal exp(iπf/2). If pₛ(x) is the parity of the occupied bits selected by a nonempty mask s, its coefficient is

$$
c_s=-\tfrac{1}{16}\sum_x(-1)^{p_s(x)}f(x),\qquad f(x)=\sum_{s\ne0}c_s p_s(x).
$$

The second identity follows from the orthogonality of the 32 binary characters and f(0) = 0. A sequence of controlled flips computes the parity onto an enrolled marked rail. One positive reference dwell applies its phase and the reverse sequence restores the inputs. The 3B program uses four nonzero masks. Each lepton orientation uses sixteen. Every computational-basis column of the compiled operator is tested. Phase alone identifies a label modulo four relative to the chosen origin.

On separate launches the first station reads the result locally or transports the five coherent dimers to cell 0 using the 438-contact material permutation. Each analyzing station reverses the coherent fanout and measures X or Y on the remaining dimer. The bright rates for label f are (1 + cos(πf/2))/2 and (1 + sin(πf/2))/2. A raw reading comprises twenty echo clicks from both rails of all five dimers. In **372 independent simulated settings of 4,096 launches each**, both stations recover all assigned labels, including the reversed lepton program, for every nonempty word. The joint 99% Hoeffding radius is 0.037004 in bright probability, including the sampling-grid allowance. The two fractional charge signs give opposite Y readings. The neutral label gives a deterministic bright X reading. All raw count bins of nonzero count are retained and their totals include every launch. [Two-station acquisition](measure.py)

## The 715 four-field ballots

A second, separately initialized sharp frame supplies a binary spinor-component index. A component record is its internal occupation word followed by this index frame. The two-frame code uses the product of the two complete deterministic transmissions. All 64 frame pairs have their sender, receiver and blocked readings in the receipt. Four contain neutral lines outside the field table. The other 60 distinguish all one-generation field components and their two spin tags. Suppressing exchange on the external links gives a refusal on every frame at the receiver. The spin tags are two formal basis indices with declared antisymmetric pairing. A physical spin measurement would require an independent spin and spacetime attachment.

Each reader identifies a fundamental colour by its occupied slot and the conjugate by the missing slot. The conjugate weak basis is indexed by the complementary positive word. Component phases and epsilon orientations are chosen consistently in these bases. The formal chirality indicator is one on even nontrivial exterior degrees. Giving each field the declared canonical dimension 3/2, four fields exhaust dimension six. The grammar uses one generation of these ten field labels, with repetition allowed and sorted names. It has no right-handed neutrino, extra field insertion or flavour index. [Eligibility and contraction reader](grammar.py) · [Source grammar and independent checks](../../../../../code/a5_closure/baryon_dimension_six_census.py)

Both stations enumerate the 715 sorted four-field multisets from their own decoded field records. The successive counts are 554 with nonzero baryon number, 12 also neutral in hypercharge, 12 with triality zero modulo three, 12 with even weak parity, and eight with even left-Weyl count. The surviving positive rows are QQQL, QQu_Re_R, QLu_Rd_R and u_Ru_Rd_Re_R. Every one has B = L = +1. The four complementary rows have B = L = −1. Reversing only the lepton phase program retains the eight eligible rows, with B = −L on every one. This comparison measures the convention on which the signed charge equation depends.

## Antisymmetric contraction signs

For a four-field monomial with four distinct component indices, assign its sorted indices in order to four vacancy dimers. A fifth dimer is the coherent reference. Start the data empty and the reference in its plus state. Before each insertion both arms pass the required local absence read. Controlled phases from the reference to preceding modes give the ordered exterior sign. A controlled flip then inserts the chosen index only on the marked reference arm. Four ordinary controlled flips return the occupied arm's data to the empty word without changing its sign. Transport the register and measure the fifth dimer in the X basis. The sign equals the sorting parity of the four component indices. Earlier occupied smaller indices contribute the sign on insertion, and the sum of their count and the inversion count is six. Reversing all four factors also has even parity.

The receiver distinguishes twelve even and twelve odd orders. With all prefix phases suppressed it reads the positive sign on all 24. Four separate attempts to insert an occupied index yield the local presence word 0001 and refuse. These are measurements of the finite four-generator rule. Relabeling its four slots for a different monomial uses the same gates and read directions.

The two readers combine their own field catalogues with this shared sign calibration. For each eligible field multiset they enumerate every weak and left- or right-spin epsilon pairing against the three-index colour epsilon. Repeated component indices contribute zero. For distinct indices they collect equal sorted monomials with the measured sorting sign. Exact rational elimination gives rank one for each of the eight oriented spaces. Counting Hermitian conjugates together leaves four independent classes in disjoint field-content sectors. The control with ordinary unsigned sorting gives zero on both QQu_Re_R representatives. The full signed coefficients, including exact zero polynomials for two repeated-right-quark pairings, are compared to the source polynomial in its own ordered basis. [Independent coefficient comparison](test_baryon_labels.py) · [Original contraction certificate](../../../../../code/a5_closure/receipts/baryon_dimension_six_census.receipt.json)

The stations implement the declared labels of PR-55 and the finite eligibility and contraction grammar of PR-56. This is an encoded operator census. It assigns no free baryon-number dynamics, physical proton, decay occurrence or lifetime. Those applications need the separate matter, gauge and effective-action attachments, including Wilson coefficients and QCD matrix elements.

From the repository root, with one numerical thread:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
nice -n 15 .venv/bin/python -m model.experiments.interval_echo_gauge.exterior_register.baryon_labels.measure
nice -n 15 .venv/bin/python -m pytest -q model/experiments/interval_echo_gauge/exterior_register/baryon_labels
(cd Lean && nice -n 15 lake env lean -j1 Screen/BaryonDimensionSix.lean)
```
