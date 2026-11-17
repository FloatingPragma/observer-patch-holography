# Conditional polar-analyzer map and balanced-path refinement

This is an additional preregistered diagnostic after the balanced96-event
fringe missed its -0.9 target. Original failures stay on file. Test matched
weight histograms at dwell4 and8 against twice-around references at dwell2
and4, with visibility<-0.9 as the unchanged useful-fringe threshold. These
are new comparisons of retained kernels, not new primitive-event experiments.

For a conditional apparatus interpretation, declare an equilateral face and
an independently recorded polar analyzer n(theta)=(cos(theta),sin(theta)).
Its three outward seam axes have angles0,2pi/3,4pi/3. Give seam e
64(1+cos(theta-phi_e)/2)/3 tickets before rounding by largest remainder.
All64 tickets are retained and weights remain positive. Rotate the declared
analyzer through96 evenly spaced samples per full turn, starting after0 and
ending at2pi. The angle records are supplied apparatus coordinates, not a
native physical-clock or spatial-rotation derivation.

Control: replace polar response by apolar response cos(2theta-phi_e), retaining
the same recorded physical-analyzer angle theta and ticket arithmetic. This
tests whether the proposed bridge depends on choosing a polar coupling. It
does not claim that the native process selects that coupling.

Execute all3 inputs×64tickets×2coins at all96 points for each coupling.
Retain all source outcomes and eraser instrument probabilities. Compute true
native transfer products at dwell1 and2. Desired diagnostic: polar one-turn
source overlap<-0.99 and apolar one-turn>0.99. All failures retained; exact
rational independent products check the signs. No source phase is inserted.

The signed geometric axes and external polar/analyzer coupling are supplies.
Changing their calibration is not a native physical rotation. The test can
produce a concrete conditional apparatus dictionary and expose its premises;
it cannot by itself fill the original row's physical-source gap.
