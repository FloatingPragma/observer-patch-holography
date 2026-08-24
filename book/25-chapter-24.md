# 24. Why Does Gravity Look Like Geometry?

In the autumn of 1907 Albert Einstein was examining patents in Bern and owed Johannes Stark a review article on relativity for a yearbook of radioactivity and electronics. He was assembling it at his desk in the office when the thought arrived that he described fifteen years afterward, in a lecture in Kyoto on 14 December 1922, as the happiest of his life. If a person falls freely, he will not feel his own weight.

Every other influence in physics discriminates. A magnet lifts a steel paperclip and does nothing at all to an aluminum one of the same mass, because the response depends on what the object is made of. Gravity pulls harder on a heavier object by exactly the factor that makes a heavier object harder to shift, the two cancel, and everything falls at one rate, which is why David Scott could stand on the Moon in August 1971 with a hammer in one hand and a falcon feather in the other and drop them together for the camera. Write the description in a frame that falls along with everything else and gravity drops out of it, leaving nothing anywhere in the equations for it to do.

The sharpest test of that cancellation flew between 2016 and 2018 in a French satellite called MICROSCOPE, which carried two hollow concentric cylinders, one of a titanium alloy and one of a platinum alloy, around the Earth in an orbit continuously corrected so that nothing but gravity acted on the pair. Electrostatic forces held both cylinders in place. The size of those forces reports any difference in how the two of them fall. Five months of usable free-fall data came out of two and a half years in orbit. The difference in acceleration between titanium and platinum came out at one and a half parts in ten to the fifteenth, with a statistical uncertainty of two and three tenths in the same units, which is to say a null result at the level of one part in a thousand million million.

So gravity belongs to the place rather than to the thing. Which leaves two questions. Why does an influence that treats everything identically turn out to be a shape, and what decides how much shape a given lump of matter makes?

## Two walkers on the equator

Two people stand on the equator, a kilometer apart. Both face due north. Both are told to walk straight ahead and never turn, and both are conscientious about it: no course corrections, no navigating, nothing but one foot in front of the other in the direction they set off in.

They collide at the North Pole, about ten thousand kilometers later.

Both of them walked straight the whole way. The gap between them closed from a kilometer to nothing, at an accelerating rate over the last stretch, and no force acted on either walker at any point in the trip. A cord tied between them would have gone slack rather than taut. The convergence is a fact about the surface they crossed, and they can establish it without ever seeing the Earth from outside, because they have their own two pairs of boots and a tape measure. Send a third walker off from a point between them and the same thing happens to each of the three gaps at once, which is what makes the effect a property of the ground rather than an agreement between two people.

The path you get by going straight ahead is called a **geodesic**. On a flat floor it is a straight line, on a sphere it is a great circle, which is why the flight from London to Tokyo goes over the Arctic and why that route looks bent on a wall map. Two geodesics that set off parallel and end up crossing is what curvature means, and the walkers measure it from the inside.

Free fall is walking straight ahead. Two ball bearings released side by side above the Earth converge, and the amount they converge by is a fact about where they are rather than about what they are made of, which is what MICROSCOPE measured to one part in a thousand million million. Nothing in the account calls for a sheet with a weight on it, and nothing bends into a further dimension for the bending to happen in.

Release a small ball of particles, all at rest with respect to each other, and let them fall. The ball's volume changes. John Baez and Emory Bunn, in the *American Journal of Physics* in 2005, wrote Einstein's equation as a sentence about that ball: the rate at which it begins to shrink is proportional to its volume times the energy density at its center, plus the pressure in each of the three directions at that point.

## The equation, read rather than solved

The whole of that, in the notation the textbooks use, takes one line.

$$G_{ab} + \Lambda\, g_{ab} = 8\pi G\, T_{ab}$$

That says the shape of the world at a point is fixed by what is sitting at that point.

The letters $a$ and $b$ each stand for one of the four directions, one of time and three of space. Four directions taken two at a time gives sixteen pairs, and swapping the two directions of a pair changes nothing, which leaves ten of the sixteen independent, so that single line is ten equations. The object $g_{ab}$ is the **metric**: ten numbers at every point of the world that turn a pair of directions into their contribution to chapter seventeen's interval, which makes the metric the ruler, the thing that says how much time and how much distance separate two nearby events. The object $G_{ab}$ is a particular combination of the rates at which those ten numbers change from point to point, assembled so that its own bookkeeping closes automatically, with nothing leaking in or out. On the right, $T_{ab}$ is the stress: energy density, momentum flow and pressure, everything that is there. And $8\pi G$ carries Newton's gravitational constant, which fixes how much bending a given amount of energy buys.

The one term that spoils the arrangement is $\Lambda$, the cosmological constant, which multiplies the metric itself and has nothing on the other side of the equals sign to account for it.

Einstein reached the equation without that term in November 1915, by requiring that it reduce to the old law of gravity for weak fields, that it not depend on the coordinates anybody chose, and that both sides conserve what has to be conserved. Those requirements narrow the possibilities to one. They do not say why geometry should respond to energy in the first place.

## Four pages in 1995

Ted Jacobson's paper occupies pages 1260 to 1263 of volume 75 of *Physical Review Letters*, titled "Thermodynamics of Spacetime: The Einstein Equation of State".

He needed three things. An observer who accelerates steadily has a horizon behind it, a boundary past which signals never arrive. A horizon carries an entropy equal to a quarter of its area. And an accelerating observer finds empty space warm, at a temperature set by its acceleration. Then Clausius's relation, which is the first thing anybody learns about heat: the heat flowing into a system, divided by the temperature it arrives at, equals the entropy that heat buys.

Jacobson's move was to demand that this hold on every one of those horizons. A black hole has one. So does the air above your kitchen table, and so does every other point of the world, one horizon for every direction an observer passing through that point might accelerate in.

Matter crossing such a horizon carries energy, so it delivers heat. The horizon's entropy is a quarter of its area, so the entropy change is an area change. The temperature is the one the accelerating observer reads. Insist that heat over temperature equals the area change at every point and in every direction at once, and the geometry has no freedom left: it has to distort in exactly the way that keeps the accounts straight everywhere, and that distortion is Einstein's equation, coupling constant included.

Jacobson had to bring the area law and the temperature in from outside. The area law came out of chapter twenty-three's collar, where what a region can hold is a count of the seams a cut along its boundary severs. The temperature came out of chapter twenty-one, where the flow a state determines runs an accelerating observer along and reads warm on a thermometer riding with it.

The repair rule drives every observer to the state carrying the most entropy compatible with the constraints it has to satisfy. Take the total entropy of a small region, meaning the entropy of its contents plus a quarter of the area of its boundary. At the state the repair rule settles on, that total is at a maximum, and the first-order change of anything at a maximum is zero. Stationarity of the total entropy on every small region, for every small variation of the state or the shape, is Clausius's relation with nothing thermodynamic assumed.

## Where the eight and the pi come from

Take a small ball around a point, in the rest frame of one observer, small enough that the geometry inside it is nearly flat, and run the stationarity condition on it.

The contents side first. Chapter twenty-one's generator, restricted to this ball, weights the energy density by a factor that falls from a maximum at the center to zero at the edge, namely the difference between the squared radius of the ball and the squared distance from the center, divided by twice the radius. That shape is what a boost centered on the ball does. Integrate that weight against the change in energy density over the whole ball and out comes four pi times the fourth power of the radius, divided by fifteen. That is an integral a first-year student can do in a line, and the fifteen comes out of it.

One factor is owed, and chapter twenty-one supplied it. Entropy is heat divided by temperature, and the temperature an accelerated observer reads is its acceleration divided by two pi. Dividing by that temperature multiplies the heat by two pi, which takes the matter side to eight pi squared times the fourth power of the radius, divided by fifteen.

The geometry side. Alfred Gray and Lieven Vanhecke, in *Acta Mathematica* in 1979, worked out how the size of a small geodesic sphere depends on the curvature around it. The piece of that needed here says that if you hold the volume of a small ball fixed and raise the curvature, the area of its boundary falls by four pi times the fourth power of the radius, divided by fifteen, times the change in exactly the curvature combination that stands on the left of Einstein's equation. Curving the surroundings costs boundary area.

Set the sum of the two to zero, remembering that the area enters the total entropy divided by four times Newton's constant. Every appearance of the radius cancels. Both fifteens cancel. What is left is eight pi squared on one side against pi over Newton's constant on the other, and the ratio of those two is eight pi times Newton's constant, which is the coupling that stands in front of the stress in the field equation. Nobody chose it and nobody fitted it. The eight pi that every textbook writes in front of Newton's constant without comment is the ratio of two elementary integrals over a small ball.

That gives one of the ten equations, in one observer's rest frame. The other nine come from the overlap rule. Chapter seventeen gave every observer a three-dimensional space of velocities, and observers passing through the same point at every velocity in an open range each supply the same relation in their own frame. Written out in terms of the velocity, each relation is a polynomial of degree two, and a polynomial that vanishes for every velocity in an open range has every coefficient equal to zero. Ten components, from a ball of observers who have to agree with each other.

## The gas law does not know about molecules

Jacobson ended his abstract with a sentence that has been irritating people for thirty years: it may be no more appropriate to canonically quantize the Einstein equation than it would be to quantize the wave equation for sound in air.

Pressure times volume equals a constant times temperature. Émile Clapeyron put it in that form in 1834, and it held up in every laboratory in Europe through the rest of the century without containing a single molecule. In 1908 Jean Perrin put grains of gamboge resin in water under a microscope, counted how the population thinned out with height, and got Avogadro's number out of the count, which finished an argument about whether matter came in pieces at all. Wilhelm Ostwald, who had spent a career arguing that atoms were a bookkeeping convenience, gave that position up in the preface to the fourth edition of his own textbook in 1909, citing Perrin by name. The gas law did not change a symbol, because it was never a claim about molecules. It is a claim about a gas that has settled down, and it stays true whatever the gas turns out to be made of.

Einstein's equation sits on that shelf. It is the condition a network of observers at its settled state has to satisfy, written in the language of the shape they collectively report. Quantizing the metric in the hope of finding the atom of space is quantizing the pressure of a gas in the hope of finding the molecule. The pressure is real, it is measurable, it will lift a piston, and it belongs to a settled state rather than to a constituent.

The 1907 thought comes out as a theorem. Every observer builds its local frame the same way, out of its own region's flow, so the frame of a falling observer and the frame of an observer floating in empty space are the same construction, and gravity can always be made to vanish locally by choosing the right one. And the walkers come out too. For matter with no pressure, the conservation that the field equation forces on the stress says that each speck goes straight ahead. Free fall along geodesics is a consequence of the equation rather than a separate law bolted to it.

## Nine light rays

Everything in the argument crossed a horizon, and horizons are made of light rays.

An observer's access to anything is bounded by chapter seventeen's light cone, so the boundary of what can reach it is traced by rays, and any quantity read off at that boundary is read by evaluating the stress along a ray direction. The stress takes two directions and returns a number. Feed it the same light-ray direction twice and you get the one number a horizon can report.

That restriction comes from the world rather than from a choice of coordinates. A record held by one observer turns into a fact for a second observer when something crosses between them, the edges of every region anybody can name are set by what can cross, and the boundary case travels at the speed light does. Slower probes report on a region from inside it, after the fact, and faster ones are unavailable at any price. So the geometry of what can be established is the geometry of rays, and a quantity the rays cannot report is a quantity no local measurement will ever hand anybody.

Ten numbers describe the stress at a point. Light rays report one number each. How much of the stress do the rays determine?

A symmetric object that reads zero along every light ray has to be a multiple of the metric, and nothing else. So the rays pin down nine of the ten numbers, and what they leave undetermined is a multiple of the metric.

Nine readings do it. Take the six rays that travel along one axis at a time, forward and backward along each of the three. Then take three more that share their motion equally among all three axes, each space component equal in size to one over the square root of three, with the signs plus plus plus, plus plus minus, and plus minus plus. The six axis rays never pick up the cross terms, the entries pairing one space direction with a different one, because a ray running along a single axis has no second direction to pair it with. The three tilted rays are there for those, and their three sign patterns read the three cross terms off one another. Nine light rays, nine numbers, one system of linear equations to invert. Its determinant is exactly eight thousand one hundred and ninety-two over twenty-seven, which is two to the thirteenth over three cubed, and a determinant that misses zero is a system that inverts. The inversion is stable into the bargain: an error in the nine readings grows by a factor of at most two plus the square root of three, which is under three and three quarters, in the nine numbers recovered from them.

## The vacuum and the tenth number

The vacuum has a stress. Whatever is in empty space when nothing is in it, the same everywhere, looking identical to every observer no matter how fast they are moving, has to be built out of the only object with that property, and the only object with that property is the metric. So the stress of the vacuum is a number times the metric.

Contract it with a light ray. The metric evaluated on a light-ray direction twice is the interval of a light ray with itself, and that is zero: chapter seventeen defined a light ray as a separation whose interval is exactly zero. So the vacuum's contribution to what a light ray reports is a number times zero.

Exactly zero, for any vacuum energy density whatever, however enormous, and zero by the definition of the probe rather than by a cancellation between large terms.

The single direction the light rays cannot see is exactly the direction the vacuum's energy points in. That is why the cosmological constant survives in the field equation as a term nobody derived: the small-ball calculation that produced everything else on the left-hand side is built out of horizons, horizons are made of rays, and rays are blind to it.

Add up the zero-point energies of the quantum fields of empty space, cut the sum off at the shortest length anybody takes seriously, and you get a vacuum energy density around a hundred and twenty orders of magnitude larger than the one the sky reports, the sky's value of the cosmological constant being about 1.1 times ten to the minus fifty-second per square meter. That is the worst prediction in the history of physics, and being wrong by that factor takes a kind of commitment a merely careless calculation cannot reach. The arithmetic in it is sound. The instrument it was carried out with reads zero on the quantity it was trying to compute. Every local calculation in physics is built out of light cones, so no local calculation was ever going to fix that number.

What fixes it is a count that no small region contains: the total number of distinguishable records the horizon around an observer can hold. The constant comes out small because that count is enormous.

## Newton's constant is an area

The area law says the entropy of a boundary is a quarter of its area divided by Newton's constant. Read that from right to left instead. Entropy is a count of records, area is an area, so Newton's constant is whatever converts one into the other. It is the number of square meters per record. One square meter of boundary holds about ten to the sixty-ninth of them.

Chapter thirteen said the size of a screen is a count of readings and not a width in meters. This is where the two get joined. Chapter twenty-three's collar says what the readings are: the seams a cut along a region's boundary severs. Take one cell of that screen. Whatever area the cell covers, that area is a whole number of one small unit, the same unit everywhere, and on the wiring this world has the cell carries a quarter of that same whole number in severed seams. Divide the cell's area by four times the entropy it carries and the whole number cancels, top and bottom, exactly, leaving the unit standing on its own.

Newton's constant is that unit of area, one square of the screen, and its being an area is the reason gravity has a length in it at all.

The area is 2.61228 times ten to the minus seventy square meters, a square 1.6163 times ten to the minus thirty-five meters on a side. Multiply it by the cube of the speed of light and divide by Planck's constant over two pi, which is the conversion any area has to go through to be quoted in the units gravity is measured in, and out comes 6.67430 times ten to the minus eleven cubic meters per kilogram per second squared. The laboratory value is 6.67430 times ten to the minus eleven, with an uncertainty of 22 parts per million.

That uncertainty is where the comparison stops being informative, and worth being plain about which side it sits on. Newton's constant is the least precisely known of the fundamental constants by a wide margin, measured to five figures where the electron's magnetic moment is known to twelve, because gravity is feeble and everything in the laboratory has mass. A derived value landing inside a window that wide is a weak test passed, not a decimal place matched. MICROSCOPE measured gravity's indifference to what falls at one part in a thousand million million, and the strength of the thing being indifferent is known to twenty-two parts per million, ten orders of magnitude worse, because a difference between two cylinders nulls out everything the two have in common and a magnitude has to be assembled out of a kilogram, a meter and a second.

So the geometry is complete. Ten numbers at every point, a rule fixing how they answer to what is there, a coupling that fell out of an integral over a small ball, and a constant of proportionality that turned out to be a unit of area. Curvature is what a settled network of observers looks like from a distance, free fall is going straight ahead, and the two walkers on the equator were doing physics the whole way to the pole.

And nothing is standing in it. A geometry can be told how much energy and momentum is present at a point, which is all the field equation ever asks. It cannot be told why what is present comes in exactly three kinds of interaction, one of them with eight varieties, or why the twelve ports of chapter fourteen have room for those and for nothing else.