# 15. Why Does Space Have Three Directions?

In September 1879 a man left his house in Essex carrying a rack of brass resonators and walked to the parish church to hit a bell.

The house was Terling Place, where John William Strutt, third Baron Rayleigh, had a wing fitted out as a laboratory. The bell was the second of the five in the village tower, counting down from the highest, cast by a founder named Gardiner in 1723. The brass was a set of Helmholtz resonators made in Paris by Rudolph Koenig, who had been apprenticed to the violin maker Vuillaume in 1851, set up on his own account in 1858, and spent the next forty years turning Helmholtz's ideas into apparatus a physicist could order by post. Each resonator is tuned to a single frequency. Hold one to your ear inside a noise and that frequency comes up out of the mess loud enough to name, while everything else drops away. Rayleigh's verdict on them, printed a decade later, was that without the security the resonators afford, deciding which octave a bell is sounding in is very uncertain.

He published in January 1890, in the *Philosophical Magazine*, under the title "On Bells". Much of the paper is a list of what he heard, bell by bell, with the founder and the casting date beside each one, going back to a Terling bell cast by Graye in 1623.

One entry is a hemispherical bell of about three hundredweight, a little over 150 kilograms of bronze, lent to him by the London founders Mears and Stainbank. Four tones came out of it plainly and he named all four. Then this: "The gravest tone has a long duration. When the bell is struck by a hard body, the higher tones are at first predominant, but after a time they die away, and leave e-flat in possession of the field."

For a fraction of a second the bell has four voices. The loudest of them are the high ones. After that it has one. Which of those two descriptions is the bell depends on how fast the listener is. Give the job to an instrument that cannot register anything shorter than a few seconds, and the four voices never happen. It reports one note, the lowest, and it reports that note as the whole bell.

You can walk forward, walk sideways, and climb a ladder, and any motion you are capable of is some combination of those three. Three numbers pin a body down and two will not: an air traffic controller is given latitude, longitude and altitude, and losing any one of the three leaves a whole line of places the aircraft could be sitting. A fourth independent number has never been available to anybody. General relativity generalizes the shape of space in every respect except the count of directions, which it takes as a starting datum. Nothing in its field equations prefers three to five or to ten. The same equations take it for granted that the three run smooth, with a place between any two places and no smallest step anywhere in them.

It comes out of a bell.

## Twelve readings and a grid of weights

Chapter fourteen left the carrier with twelve ports, thirty seams, five seams meeting at each port, six axes of paired ports, and sixty rotations carrying the whole arrangement onto itself. Chapter nine left a rule running on it: each port's reading becomes itself, less a sixtieth of the total gap between it and the five readings around it. Sixty because there are thirty seams and a midpoint move closes half of a gap rather than all of it.

Write the twelve readings as a column of twelve numbers, top to bottom, one line per port. A column like that is a **vector**, and the word carries no more than that here: an ordered list, kept in order so that the third entry always means the third port.

The rule turns one such column into another. Look at what a single output line is made of. Port three's new reading is eleven twelfths of port three's old reading, plus a sixtieth of each of its five neighbors' old readings. Every output line has that shape, with different neighbors in it. Each new number is a weighted sum of the old numbers, and the weights can be laid out in a square grid: twelve rows for the twelve outputs, twelve columns for the twelve inputs, and in row three, the number eleven twelfths in column three, a sixtieth in the five columns belonging to its neighbors, and zero in the other six. A square grid of weights used this way is a **matrix**. This one has 144 entries in it, of which 72 are zero, and eleven twelfths plus five sixtieths is one, so every row of it adds to one and so does every column. That is the arithmetic statement of something chapter nine established by a different route: repair moves the total around the twelve observers without creating or destroying any of it.

Run it once on something. Put 12 units at one port and nothing anywhere else. That port keeps eleven twelfths of what it had, or 11, and each of its five neighbors picks up a sixtieth of 12, which is 0.2. Eleven plus five lots of 0.2 is 12, so the total sat where it was, and the gap between that port and each of its neighbors shrank from 12 to 10.8. That is one pass.

The question the bell asks is what the grid does after many passes. That means applying it a hundred times, which by hand means a hundred rounds of multiplying and adding across 144 entries, and after all of it a column of twelve numbers that has to be interpreted.

## Two observers and a number raised to a power

Alice holds 10 and Bob holds 0. The only rule in operation is that each of them moves a tenth of the way toward the other.

One pass: Alice holds 9, Bob holds 1. Second pass: 8.2 and 1.8. Third: 7.56 and 2.44. Do the fourth on paper. Nine tenths of 7.56 plus a tenth of 2.44 gives Alice 7.048, and the same sum with the tenths swapped gives Bob 2.952. The numbers are unpleasant and getting worse, and their sum is 10 at every stage, and their difference goes 10, then 8, then 6.4, then 5.12, then the 4.096 you just produced.

The arithmetic behind that is two lines. Let Alice hold A and Bob hold B. After one pass Alice holds nine tenths of A plus a tenth of B, and Bob holds a tenth of A plus nine tenths of B. Add those two expressions together and the tenths cancel, leaving A plus B, which is what the pair started with. Subtract the second from the first and the tenths reinforce, leaving eight tenths of A minus B.

Look at what the rule did to those two combinations. It multiplied the sum by 1 and the difference by 0.8, exactly, and it will do the same on the next pass and on every pass after that, from any starting values whatever, because the arithmetic that does it never refers to what the numbers are. The whole two-by-two grid, read in those two combinations, is two ordinary numbers.

A direction that a grid of weights merely scales, without swinging it round to point somewhere else, is called an **eigenvector**, and the number it gets scaled by is its **eigenvalue**. The prefix is the German for "own", and these are the grid's own directions, the ones it was built out of whether anybody noticed or not.

What is left of the gap between Alice and Bob after a hundred passes is 0.8 raised to the hundredth power of what it was, which is one button on a calculator rather than a hundred rounds of arithmetic. Ten passes leave 0.8 to the tenth, which is 0.107, so a gap of 10 has closed to a little over 1. The same substitution works with twelve observers instead of two, on one condition.

The weight carrying port three into port seven has to equal the weight carrying port seven into port three, and a grid where that holds everywhere is **symmetric**. When it holds, a full set of such directions always exists: as many of them as the grid is wide, sitting at right angles to one another, and every possible column of readings is a sum of pieces along them. That statement is the **spectral theorem**. It holds for every symmetric grid of numbers of any size, and the repair grid qualifies, because a seam pushes both of its ends toward the middle by the same amount.

So the twelve-port rule is twelve numbers wearing a costume. Getting them out of the 144 by force means solving an equation of the twelfth degree.

## Sixty rotations and four numbers

The wiring hands them over instead.

The rule cannot tell one port from another. It is written entirely in terms of which port shares a seam with which, and a rotation of the arrangement carries seams to seams, so rotating the readings and then repairing gives the same twelve numbers as repairing and then rotating. Either order, same answer, all sixty rotations.

That is a heavy constraint, and it turns into a count. The sixty rotations shuffle the twelve readings around among themselves, and under that shuffling the twelve directions of reading-space fall into pieces: collections of directions that the rotations mix internally and never mix with anything outside. Some such pieces can be cut into smaller pieces of the same kind. Some cannot, and the ones that cannot are where the argument lands. For twelve readings on this wiring the uncuttable pieces have sizes one, three, three and five, and one plus three plus three plus five is twelve.

A rule that commutes with all sixty rotations acts on each uncuttable piece as multiplication by a single number. Suppose it did something else, treating one direction inside a piece differently from another. The rotations carry the first direction to the second, so the rule would have to agree with itself about two directions it was distinguishing, and it does not get to do both. That argument is **Schur's lemma**, which Issai Schur published in Berlin in 1905 to put the character theory of finite groups on a new footing. It is chapter five's symmetry doing arithmetic on 144 unknown weights.

So a rule respecting this wiring, however elaborate its 144 entries look, has at most four distinct answers in it. One for each uncuttable piece.

The two pieces of size three are not copies of each other. The sixty rotations move the directions inside one of them along a different pattern from the other, which is why the two get separate numbers and why the count is four rather than three. Two pieces of size three went into the wiring, and the world is built on one of them.

## The three rates

The piece of size one is the direction in which all twelve ports shift together. Nobody disagrees with anybody, there is nothing for a repair to close, and the rule leaves it exactly where it is: eigenvalue one, forever. It also carries no comparison, which means no observer inside the arrangement can read it, and it drops out. Eleven directions are left.

The other three numbers are the fraction of a disagreement that survives one pass of repair in each of the three remaining pieces. Two of the three are irrational.

$$\frac{55+\sqrt{5}}{60} \approx 0.9539 \qquad \frac{9}{10} = 0.9000 \qquad \frac{55-\sqrt{5}}{60} \approx 0.8794$$

Three directions keep 0.9539 of whatever they were holding, five keep exactly nine tenths, three keep 0.8794, and those three numbers, with the 1 that no observer can read, are the entire content of a twelve-by-twelve grid of weights. The 55 over 60 in the outer two is the eleven twelfths a port keeps of its own reading, and the square root of five over sixty is what its five neighbors put back or take away. Slow band, middle band and fast band, named for the rates rather than the sizes.

The square root of five is the same square root of five a reader met in chapter two, in one plus the square root of five all over two, which a calculator kept landing on from any starting number. It is here because the solid carrying the twelve ports is built out of it, and because a fraction cannot be its own reciprocal plus one. The distance between the top rate and the bottom one is the square root of five divided by thirty, which comes to 0.0745 per pass, and thirty is the number of seams. One of the two pieces of size three fades slower than the other by that margin, which is a small number to have a body in.

Both of them can be written down. Set the arrangement on a table, choose any flat surface to measure heights above, and give each port a reading equal to its own height above that surface. Those twelve numbers fade at 0.9539 a pass and at no other rate. There are three independent versions of the pattern because there are three independent directions to measure a height along. The other three-dimensional piece is the same recipe run on a different solid. Chapter two's calculator settled on 1.618034 and the quadratic behind it had a second solution, minus 0.618034, a fixed point that nothing ever arrives at. Build the twelve corners using minus 0.618034 everywhere the golden ratio stood, read the heights off those, and the patterns that come back are the ones that fade at 0.8794.

## Slow ears

Strike the arrangement, meaning disturb it away from agreement in no particular direction, and the disagreement lands across all eleven directions at once. Repair works on all three bands simultaneously, at the three rates above.

Follow it. A disturbance in the fastest band is down to half in 5.4 passes. The middle band takes 6.6 passes to halve, and the slow band takes 14.7. Spread the initial disagreement evenly across the eleven directions and three elevenths of it, 27 percent, sits in the slow band at the start, which is a minority position. After ten passes the slow band holds 42 percent of what is left. After thirty passes, 73 percent. After sixty, 94 percent. After a hundred, 99.5 percent.

Those are shares of a total that is itself collapsing. After a hundred passes the slow band holds nearly nine parts in a thousand of what it was handed, the fast band holds under three parts in a million, and the ratio between them is 3,414. Nothing was tuned to produce that. Change the rates and you would have to change the wiring, and the wiring is what chapter fourteen showed there is only one of.

Which is Rayleigh's hemispherical bell, with the founder's name filed off. Four voices at the strike, one voice a moment later, and e-flat left in possession of the field.

Observers here are built out of repair records, and a repair record is not made in one pass. Chapter nine's limit on how fast influence travels says one seam per repair, and the two most distant ports on this wiring sit three seams apart, so the quickest possible confirmation between that pair is three repairs out and three back. Six repairs is longer than the fast band's half-life of 5.4 passes, so more than half the fast content is gone before the quickest question this arrangement can ask about itself has its answer back. Any observer inside this arrangement is therefore reading it late, after the fast content has gone, which makes observers built out of repair records slow ears. Over those same six passes the slow band loses a quarter of what it holds and no more.

That settles which content can carry a position. A position is a quantity two separated observers can compare and both get the same answer for. The comparison takes rounds to run. Content that has faded to nothing before the comparison completes cannot anchor anything: by the time the second observer reports, the quantity it was reporting on has relaxed away. Only the slowest band lasts long enough, and its size is three.

The five-dimensional piece is the largest one on the wiring and it loses anyway. Its content halves in 6.6 passes against the slow band's 14.7, and after sixty passes the five directions are holding about a twentieth of what the three are holding. Whatever an observer writes into those five directions has drained most of the way out before a second observer can be told about it. A quantity like that describes nothing that has a place.

That is where three comes from. It is the number of independent directions in the one band that outlives every other band, on a wiring that has no alternatives.

## The table of angles

The slow band is three directions wide. Look at what the twelve ports are doing inside it.

Each port leaves a shadow in the slow band, a direction with a length. To compare two directions, take the number that says how much they point the same way: one when they coincide, zero when they are at right angles, minus one when they point exactly opposite. It is the cosine of the angle between them, the only thing about two directions that survives forgetting where they are.

Do that for all 144 pairs of port shadows and the table has exactly four different numbers in it. A port against itself gives 1. Two ports sharing a seam give one divided by the square root of five, which is 0.4472, an angle of 63.43 degrees. Two ports that are two seams apart give minus that. A port and the port directly across the arrangement from it give minus 1.

Those are the angles the twelve ports of chapter fourteen's arrangement stand at, arriving here out of a rule about disagreement rather than out of anybody having measured anything. What went into the calculation was a list of which port shares a seam with which. What came back was a solid.

The count of independent directions in a table like that is called its **rank**: how many of the twelve shadows you have to be given before the other shadows are forced. Two will not do it. Four is more than necessary. The rank of that table is three. Every entry in it is 1, or minus 1, or one over the square root of five with a sign in front, so the whole calculation runs in arithmetic built out of the square root of five and nothing anywhere in it is rounded.

That three is the slow band's size arriving a second time, since shadows lying in a three-directional band cannot supply more independent rows than the band has directions. A world with four directions in it would need a slowest band of size four. The sizes available on this wiring are one, three, three and five.

## Six counters

A three-dimensional space in the abstract has nowhere in it to stand, and getting into it takes whole numbers.

What an observer holds about where it is, at any moment, is a tally. Every comparison it completes through a port adds one to a count. The twelve ports come in six pairs facing opposite ways, and a comparison recorded out through one port is the same event as a comparison recorded inward through the port facing the other way, so the two counts are one count with a sign on it. Twelve ports, six pairs, six running totals, and every one of them a whole number, because a comparison either completed or it did not. Nobody ever made three and a half comparisons. An observer that has completed forty comparisons out through one port and thirty-seven inward through the port facing it holds a 3 on that axis, and five more numbers like it, and that list is everything it has ever known about where it is.

Six whole numbers, dropped into a space with three directions in it. That looks like it should lose information, and on ordinary coordinates it does: six real coordinates fed through the table of angles collapse onto three. On whole numbers it loses nothing. For two different tallies to name the same point, some collection of whole numbers other than all zeros would have to be flattened by the table, and the flattening needs the square root of five to be a ratio of two whole numbers. It is not, and no two histories of comparison ever name the same place.

What those tallies do instead is fill the space. Take any point of the three-dimensional band and any tolerance you care to name, and some finite history of comparisons lands inside that tolerance of it. The landing places can be numbered off one after another the way the fractions can. They are packed everywhere. No history of comparisons lands on the point itself.

## Nobody ever holds pi

Write down 3. Then 3.1. Then 3.14, then 3.141, then 3.1416.

Every one of those is a ratio of two whole numbers: the fifth is 31,416 over 10,000. The list crowds in on something, and the something it crowds in on is not on the list and never will be. Nobody has ever written pi down. Every circle any engineer has ever cut was cut to a finite number of digits, and the number the digits are converging on is not a thing anybody has held.

The argument that a list like that is a problem is older than the number line. Zeno of Elea made it in the fifth century BC, in a book that is lost. It survives because Aristotle set it out in the *Physics* in order to answer it: the first argument "asserts the non-existence of motion on the ground that that which is in locomotion must arrive at the half-way stage before it arrives at the goal." Cross a room and you cross half of it first. Then half of what is left, then half of what is left after that. The halves do not run out, so the crossing is a list of tasks with no last entry on it, and Zeno took that to mean the crossing never happens.

Add the list up by hand. A half. A half and a quarter is three quarters. Put in an eighth and you have seven eighths, and the gap left over is an eighth, the very piece you just put in. That happens at every stage: each total falls short of 1 by exactly the piece last added, which is also exactly what all the pieces after it come to. Writing n for the number of pieces added, the pattern is one line.

$$\frac{1}{2} + \frac{1}{4} + \cdots + \frac{1}{2^n} = 1 - \frac{1}{2^n}$$

At n equal to three the left side is the three pieces you just added, seven eighths, and the right side is 1 minus an eighth, which is the same number. Name a tolerance and the line says which stage gets you inside it. A thousandth needs ten pieces, because 1 over 2 to the tenth is 1 over 1024, and every total after the tenth is closer than that one.

That is the whole of what "the sum is 1" says, and the surprise is in what it does not say. The limit is not the last term. There is no last term. None of the totals is 1, no amount of adding produces one, and the statement is about where the totals pile up rather than about any total on the list. Calling it a trick is fair, provided the trick gets named accurately: it replaces an act nobody can perform with a test on finite quantities that anybody can run. For any tolerance somebody names, some finite stage gets inside it and every stage after that stays inside. That is what an infinite sum is, all of it. Cauchy had pushed the argument into inequalities in the *Cours d'analyse* of 1821 and stopped short of putting it that way; Weierstrass put it that way to a differential calculus class in Berlin in the summer of 1861, twenty-three centuries after Zeno.

The fix for pi is a single operation and it is the operation that built the number line. Take the ratios of whole numbers, which have gaps in them wherever a sequence crowds together without arriving, and add in the limit of every such sequence. What comes out has no gaps left. It is the real numbers: the fractions together with every destination the fractions were pointing at. The name for the operation is **completion**.

Run that operation on the record points inside the three-dimensional band, using the distances the table of angles supplies, and what comes out is continuous three-dimensional space. Every point of it is the destination of some sequence of ever longer comparison histories, in the sense that 3.1416 and its successors are a destination.

The sixty rotations come along. Each of them permutes the twelve ports, so each leaves every entry of the table of angles precisely where it was, so each preserves every distance in the completed space. A transformation that preserves all distances is a **rigid motion**: what a machinist does picking up a part and setting it down the other way round, with nothing stretched. The group chapter fourteen counted at sixty elements acts on the space it built as the turns of a rigid body.

## What holds at a finite resolution

Both constructions have a limit in them, and nothing in the world is at either limit.

Run the repair grid any finite number of passes and count the independent directions in what it does. The answer is eleven. Every time, at ten passes and at ten million, because all three bands are present, none of them has reached zero, and a number smaller than another number is not zero. Three appears when the overall size is divided out and the number of passes is allowed to grow: the fast and middle bands vanish against the slow one, and the rank of what survives is three. At any actual number of passes it is eleven.

The tally has the same shape. Cap the comparisons at any figure you like, a thousand or a thousand billion, and the places an observer can be sitting are finite in number, a scatter with gaps between them, and raising the cap never closes a gap. They close in the completion, which is a statement about where the sequence is heading rather than a place any tally has reached. An observer a billion comparisons old holds six integers, sits at a finite number of passes, and is separated from the continuum by exactly as many steps as it was at the beginning, which is all of them.

So the three-dimensional continuum is what both constructions converge on, and nothing inside the world is sitting at it. Every observer sits at some finite number of passes and holds some finite number of counts. Every rate, rank and angle quoted here holds at that resolution and at no other. Space is smooth in exactly the way 3.1416 is pi: in a limit, exactly, and at no stage anybody occupies.

Every quantity in both constructions is one somebody could write down: twelve ports, thirty seams, six counts on six axes, a whole number of passes. Smooth three-dimensional space is where descriptions like those accumulate, in Weierstrass's sense of the word, and nobody occupies it, the way no total on the list of halves is ever 1. Zeno's runner crosses the room at walking pace, because crossing a room was never the same act as finishing a list with no last entry on it.

Which leaves one word doing more work than it has been paid for. Every rate here is per pass. Nine tenths per pass, half in 5.4 passes, 94 percent after sixty. A pass is an event that has finished happening, and the shares only mean anything if the events fall into an order. What a pass is, why passes fall into an order at all, and why that order has a preferred way round are three questions a wiring diagram does not answer. Rayleigh's bell rang down and did not ring up.