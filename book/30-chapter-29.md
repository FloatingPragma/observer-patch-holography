# 29. Why Do the Masses Land Where They Do?

Yoshio Koide was an associate professor at Shizuoka Women's University at the end of 1981, working on composite models of quarks and leptons. Models of that kind put something smaller underneath the particles on the roster, so that the shape of the roster becomes a consequence of how few pieces there are to build from. The field was busy. In 1979 Haim Harari and Michael Shupe published the same proposal four pages apart in the same issue of *Physics Letters B*, without either knowing what the other was doing: every quark and every lepton built from triples of two objects, one carrying a third of the electron's charge and one carrying none. Harari called them rishons. Nobody has ever seen one. One line of arithmetic came out of Koide's model and stayed.

Take the masses of the three charged leptons, the electron, the muon and the tau. Add them up. Then take the square root of each of the three masses, add the three square roots together, and square that total. Divide the first number by the second.

The answer is two thirds.

The three masses have no business producing a number like that. The muon is two hundred and seven times the electron and the tau is three and a half thousand times the electron, and a quantity assembled out of three separately measured masses is supposed to come out as a decimal with nothing in it. Koide was in a position to notice, because his model had given him a reason to compute exactly that combination.

Two of the three masses were known to many digits by 1981. The tau was the trouble. It had been found at Stanford in the mid-seventies, it dies in less than a picosecond, and its mass had to be reconstructed from what came out of the decay rather than weighed. The world average sat at 1784.2 million electron-volts, give or take 3.2. Koide's relation, given the electron and the muon and asked for the third, returned 1776.97.

The gap is 7.2 million electron-volts, and physics reports a gap in units of the experiment's own quoted uncertainty. A measurement is a claim about a range, which is how chapter four read two different Higgs announcements as one result. The half-width of that range has a symbol, the Greek letter sigma, and a distance between two numbers gets quoted as a multiple of it. Koide's tau mass sat 2.25 sigma from the measured one. Gaps that size open and close as the measurements improve. The word discovery is reserved for five, which is the threshold ATLAS and CMS cleared in July 2012.

The unit in that comparison belongs to somebody. Sigma is whatever the experiment declared its own uncertainty to be, so a distance in sigma is a distance measured with a ruler that a particular group of people manufactured in a particular year, and rulers get replaced.

In the spring of 1992 a group running the Beijing Spectrometer, a detector on China's electron-positron collider, went at the tau mass from a different direction. Below a certain beam energy you cannot make a pair of taus at all and the production rate is zero. Above it the rate climbs. The corner between the two is the mass, twice over. They set the beam using the standing world average of 1784, took twelve measurements across the threshold, and found the corner at 1776.9 million electron-volts, 7.2 below where the average had put it, with an uncertainty smaller by a factor of seven.

Koide's 1776.97 had not moved. The world average had walked 7.2 million electron-volts and arrived on it.

Something holds the three masses at two thirds to three parts in a million. It is a triangle.

## Three corners, two numbers

Chapter twenty-six put the three families of matter on the three corners of an oriented face of the icosahedron, which is where the count of three came from. The electron, the muon and the tau sit one to a corner.

An oriented face has one rotation left in it. Rotate the corners round the face: one goes to two, two goes to three, three goes to one. Do that three times and everything is back where it started. Of the sixty rotations of the icosahedron, that turn, its repeat and leaving the solid where it stands are the only three that keep both the face and its sense of circulation.

Take three identical rooms arranged in a ring, each joined to its two neighbors by identical doorways, and put a heater in one of them. Wait, then measure the temperature in all three. Doing that for each room in turn gives you nine numbers: how warm room one gets when room one is heated, how warm room two gets when room one is heated, and so on across and down. Nine numbers, and seven of them are redundant, because the ring has no room that differs from any other room. There is a number for how much heat stays where it was put and a number for how much crosses a doorway. The whole table is those two numbers arranged in a pattern.

The face works the same way. Push at one corner and read what comes back at all three, which is again a three-by-three table.

A response that respects the rotation cannot tell corner one from corner two, because the rotation carries the first onto the second and the response has no way of objecting. So the second row of the table is the first row shifted one place along, and the third row is the second shifted again. Nine numbers collapse into two: one for a corner's response to itself, written a, and one for its response to a neighbor. The neighbor entry carries a size and an angle, because a response can arrive turned as well as scaled, and chapter nineteen's amplitudes are that kind of number. Call the size b and call the angle the twist.

A table whose every row is the row above it shifted one step round the cycle is called a **circulant**.

Chapter fifteen took a twelve-by-twelve grid of repair weights and hunted for the grid's own directions, the ones it merely scales without swinging them round, and called the scale factors eigenvalues. Finding those for a general grid means solving an equation whose degree is the width of the grid. A circulant hands them over for nothing, because its own directions are fixed by the cycle before anybody says what a and b are: they are the three ways of stepping a phase round the triangle in equal parts.

The three scale factors come out as a, plus twice the size b, multiplied by a cosine. The three cosines are the cosine of the twist, the cosine of the twist plus a hundred and twenty degrees, and the cosine of the twist plus two hundred and forty.

Three cosines, a hundred and twenty degrees apart. The face has fixed their spacing and left the starting angle free.

## Two lines about cosines

Draw three arrows of length one from a single point, a hundred and twenty degrees apart. They are the spokes of an equilateral triangle drawn from its center. They add up to nothing, because there is no direction for the sum to point. The cosine of an angle is the shadow one of those arrows casts on a chosen axis, and the shadow of a sum is the sum of the shadows.

Three cosines a hundred and twenty degrees apart sum to zero.

For the squares, use the fact that the square of a cosine is half of one plus half the cosine of the doubled angle. Double the three angles and you get the twist doubled, that plus two hundred and forty degrees, and that plus four hundred and eighty. Four hundred and eighty degrees is a hundred and twenty with a full turn discarded, so the doubled angles are three directions a hundred and twenty degrees apart, and by the line above their cosines sum to zero. Three halves of one, plus half of nothing.

Three cosines a hundred and twenty degrees apart have squares summing to three halves.

Add the three eigenvalues. The three copies of a give three a. The cosine terms cancel by the first line. Total three a, whatever the twist happens to be.

Add their squares. Each square is a squared, plus twice a times twice b times a cosine, plus four b squared times the square of a cosine. The middle terms cancel by the first line. The last terms come to four b squared times three halves, which is six b squared. Total three a squared plus six b squared.

What a corner supplies is an amplitude, and chapter nineteen's rule for turning an amplitude into a measurable quantity is to square it. Each eigenvalue is therefore the square root of its corner's mass, up to one overall scale that fixes the units and cancels out of any ratio. The square root of a mass is a positive number, so the only responses describing anything that exists are the ones whose three scale factors all come out positive.

Koide's ratio, written Q, is the sum of the eigenvalues' squares divided by the square of their sum. The two lines above have evaluated both of those.

$$Q = \frac{m_e + m_\mu + m_\tau}{\left(\sqrt{m_e} + \sqrt{m_\mu} + \sqrt{m_\tau}\right)^2} = \frac{1}{3} + \frac{2}{3}\left(\frac{b}{a}\right)^2$$

On the left are the three charged lepton masses, measured in laboratories. On the right is one number, the size of a corner's response to its neighbor divided by the size of its response to itself. Everything else has cancelled, the twist included, which appears nowhere on the right and therefore cannot be tuned to help.

Q equals two thirds exactly when b over a is one over the square root of two.

The left side of that equation is trapped, for any three positive masses anybody cares to write down. Set all three equal and it gives one third, which the right side agrees with by sending b to zero. Let one mass run away from the other two and it climbs toward one, without arriving. So every possible triple of masses in the universe lands somewhere in the range from one third to one. Two thirds is the exact midpoint of that range. The three charged leptons, whose masses cover a factor of three and a half thousand, sit on the midpoint.

Two pieces make up the response. There is the direction in which all three corners move together, which carries weight three a squared, and there is the two-dimensional plane in which they move against one another, which carries weight six b squared. Setting b over a to one over the square root of two sets those two weights equal. One condition on the response of a single face, and Koide's two thirds is that condition written out in masses.

The twist gets a small allowance rather than a value. At balance, the three eigenvalues stay positive as long as the twist is within fifteen degrees of one of the three spokes, and anywhere inside that band the ratio sits at two thirds without moving.

Take the three square roots of the masses as the three coordinates of a point in space. There is a line through that space along which all three coordinates are equal, the line where the three leptons would weigh the same. The balance says the point stands at exactly forty-five degrees from that line, and forty-five degrees is why the relation looks so peculiar written in masses. A clean statement about square roots turns into an agreement to five decimal places once everything in it has been squared.

## Seventy-two electron-volts

Feed in the two masses known best. The electron is 0.51099895 million electron-volts and the muon is 105.6583755 million electron-volts. Requiring the balance to hold turns the relation into a quadratic in the square root of the third mass. A quadratic has two roots. One of them is 3.3174 million electron-volts, which would put the third charged lepton between the electron and the muon rather than above both of them. The tau is the heaviest of the three, which is why it can decay into a muon and a muon cannot decay into it, so the upper root is the one.

Chapter twenty-eight turned a computation into a proof by carrying intervals rather than points. The same machinery applies here: give the electron and the muon their measured widths, do the arithmetic on the widths, and what comes back is an interval guaranteed to contain the answer.

It runs from 1776.968991 to 1776.969063 million electron-volts, centered on 1776.969027.

That window is 72 electron-volts wide. The measured tau mass is 1776.93 million electron-volts, give or take 0.09, and 0.09 million electron-volts is ninety thousand electron-volts, so the uncertainty on the measurement is twelve hundred and fifty times the width of the window it is being compared against. The center of the window is 0.039 million electron-volts from the measured value, which is 0.43 sigma.

Fed the measured tau mass, Koide's ratio comes out at 0.66666446.

The sharpest tau measurements are threshold scans of the kind the Beijing Spectrometer ran, reaching ninety thousand electron-volts. What limits them is how precisely a beam energy can be calibrated at the point where a pair of taus becomes possible, and that is where the next three decimal places are.

A measured central value more than three standard uncertainties from 1776.969027 destroys the balance condition. There is no free parameter left to absorb the move: the twist cancelled, the scale cancelled, and the one remaining quantity was set to one over the square root of two by the requirement that the two pieces of the response weigh the same. At an uncertainty of 0.09 million electron-volts, three of them is 0.27, so a central value below 1776.70 or above 1777.24 ends it.

## Thirty-one axes

In June 1963 Nicola Cabibbo published three pages in *Physical Review Letters* under the title "Unitary Symmetry and Leptonic Decays". Decays that changed strangeness were running slower than decays that did not, by a consistent factor, and the theory of the day had no room for a consistent factor. Cabibbo's move was to let the weak interaction couple to a rotated combination of the down quark and the strange quark, turned by a single angle. One number, fitted from decay rates, and the discrepancy went away.

The strength of the transition from a down quark to a strange one is measured at 0.2250. The angle whose sine that is comes to 13.0029 degrees. Rates go as the square of that strength, so a decay that changes strangeness runs at about a twentieth of the rate of one that does not. That factor of twenty is the suppression sitting in the data he was accounting for.

An icosahedron has axes of three kinds. Six of them run through pairs of opposite corners. A fifth of a turn about one of those leaves the solid looking untouched. Ten run through pairs of opposite faces, at a third of a turn. Fifteen run through pairs of opposite edge midpoints, at half a turn. Six and ten and fifteen is thirty-one.

Take the angle between every pair of those axes. Four hundred and sixty-five pairs, one angle each, and the smallest nonzero angle anywhere in the list is 20.9052 degrees.

13.0029 is not in the list. The reason it cannot be got at is that 20.9052 is the floor. The smallest angle the geometry has to offer is sixty percent larger than the one the quarks require, and relabeling which axis is which does nothing, because the floor is a property of the solid rather than of the labeling. The angle that sets quark mixing is not an axis angle of this geometry.

The list has an end, which is what makes an empty search a result: four hundred and sixty-five angles, written out in full and read down to the last one, at an arithmetic cost below what a phone spends drawing its own screen once.

Koide's own model reached past the leptons to the quarks. The quark side did not hold. Two structures built out of entirely different material, a set of subparticles and a twenty-sided solid, keep the charged leptons and stop at the quarks.

The wiring fixes the electric charges to every decimal anyone has measured. It fixes the number of colors at three and the number of families at three. It fixes the balance at the corners of one oriented face to within 72 electron-volts, and prints the value at which the balance dies. It does not hand over the angle by which the weak interaction turns a down quark into a strange quark: 20.9052 against 13.0029 proves that the angle is no axis angle of this solid. Had 13.0029 turned up anywhere on that list, the charges and the family count would have been worth less, because a solid that yields whatever angle is asked of it has yielded nothing.

The settled entries and the untouched ones end up on the same line. Working physicists carry it in four terms, holding every charge, every family count, every coupling and every mass on the roster, the forced ones and the typed-in ones side by side, with nothing marking which is which. Nobody has ever derived that line.