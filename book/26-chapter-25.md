# 25. Why Exactly These Forces?

In February 1764 the French army sent a twenty-seven-year-old engineer named Charles-Augustin de Coulomb to Martinique to build a fort. Eight years on Fort Bourbon got him a good deal of what is taught as soil mechanics and a succession of tropical illnesses that damaged his health for the rest of his life. He came home a captain in June 1772.

The Academy of Sciences in Paris wanted a better mariner's compass, and Coulomb's entry, which took a share of the grand prize in 1777, hung the needle from a fine thread instead of resting it on a pivot, because a pivot has friction in it and a thread has almost none. That left him needing to know exactly how hard a twisted thread pulls back, so he measured it, at length, and published the numbers in 1784: the pull rises with the fourth power of the wire's diameter, falls off with its length, and holds exact proportion to the angle of twist, which makes a long thin one as weak a spring as an experimenter could want and an exact one.

In 1785 he hung a light horizontal arm from one, mounted a small ball on its end, charged it, and brought a second charged ball up beside it. The two pushed apart, the arm swung, and the twist in the thread read off how hard they pushed.

Three readings from that first memoir on electricity carry the whole argument. Left alone, the thread let the balls settle thirty-six degrees apart and carried thirty-six degrees of twist. He turned its top through a hundred and twenty-six degrees, forcing them closer, and they came to rest eighteen degrees apart, held by a hundred and twenty-six plus eighteen, which is a hundred and forty-four degrees of twist. Half the separation, four times the push. He wound the top through five hundred and sixty-seven degrees, better than a full revolution, and the balls held out at eight and a half degrees, with five hundred and seventy-five and a half degrees of twist in the thread. The separation had come down by a factor of a little over four, and sixteen times thirty-six is five hundred and seventy-six.

Which is a fact about distance: every version of the law written down since 1785 says only how far apart two things are. What does a force law look like on a carrier that has no distances in it, and how many forces will that carrier hold? The twelve observers of chapters six through fourteen have a wiring diagram, thirty pairings saying which port is joined to which, and not one number anywhere saying how far apart anything is.

Two centuries of better instruments have left the exponent where Coulomb put it. Williams, Faller and Hill, writing in *Physical Review Letters* in 1971, wrote it as two plus an unknown small amount and then measured the amount. It came to two point seven parts in ten million billion, give or take three point one, which pins the exponent at two out to the fifteenth decimal place. Fifteen decimal places of pure distance, and the machine has no ruler in it.

## A value at every place

Between 1698 and 1700 Edmond Halley took a small Royal Navy vessel called the Paramore twice across the Atlantic, on two voyages fitted out for the single purpose of measuring the Earth's magnetism, and recorded the angle between where his compass pointed and where north actually was. The angle changes as you sail. In 1701 he published the results as a chart, and rather than print a table of positions and readings he drew curves through the places where the deviation came out the same. Sailors called them Halleyan lines for the next hundred years.

What that chart holds is a number attached to every place on it. A captain anywhere in the Atlantic could read his local deviation off the paper. The reading changes from place to place, the rule for finding it does not, and a ship crossing the ocean experiences the chart as a sequence of small differences between one position and the next.

That is a **field**, one value per place, and the word covers nothing else. The temperature in a room is a field. So is the wind over Halley's Atlantic, which hangs a speed and a heading at each place rather than a single number, and so is the deviation he charted, which needs only one.

On the twelve-port carrier the places are the ports. Give every port a number and each seam reads the gap between its two ends, thirty readings nobody had to supply separately. Chapter twelve had the same twelve numbers standing for something else, the private labeling each observer keeps at its own port. Twelve values and thirty differences: no more furniture than Halley put on his chart, and considerably less ocean.

## Twelve numbers that add to nothing

Chapter twelve left a law sitting at the ports: what a port pushes out along its five seams is the charge it holds, which is Gauss's law with the geometry taken out of it. By itself that law picks no field, since a current going round in a circle leaves every port's books balanced, so one list of charges comes back with a whole family. One member of the family has the least energy, the one with no circulation in it.

Energy on this carrier is a sum over the thirty seams. Each seam contributes the square of the number sitting on it. Squaring is what makes circulation expensive: a current running around a loop delivers nothing to any port and puts a value on every seam it crosses, and every one of those values gets squared and added in. So the arrangement of least energy is the one that has stopped moving anything it does not have to move, and on a wiring diagram with nineteen independent loops there is exactly one such arrangement for any list of charges.

Adding the same amount at all twelve ports changes no difference anywhere, so the twelve numbers are pinned only up to a common shift, exactly as chapter five's meridians were pinned only up to a choice of zero. The shift is fixed here by requiring the twelve values to add to nothing.

And the surface is closed. Twelve ports wired into a sphere have no edge, no outside and nowhere for a net outflow to go, so the twelve charges have to sum to zero or the law has no solution at all. A single unit of source at one port therefore arrives with its own balance: one twelfth taken off every port, the same subtraction everywhere, including at the port the source is sitting on. That subtraction is the only balance that treats the twelve ports alike.

Put one unit of source at one port, then, and ask what the other eleven read.

## Four fractions

Walk the wiring outward from the source and count hops. The eleven others sort into three groups: a ring of five, a second ring of five past it, and the port directly opposite, three hops out and the only one chapter fourteen leaves at that distance.

Four distances, therefore, and four numbers. Here they are, exact, in the order the walk produced them:

$$\frac{7}{36}, \qquad \frac{1}{90}, \qquad -\frac{7}{180}, \qquad -\frac{1}{18}$$

The first is the potential at the source itself, the second is what every port in the first ring reads, the third is what every port in the second ring reads, and the fourth is the reading at the far port.

Put them over a common denominator of a hundred and eighty and they read thirty-five, two, minus seven and minus ten. Add all twelve of them up, taking five at the second value and five at the third: thirty-five and five twos on one side, five sevens and a ten on the other, forty-five against forty-five. The twelve values add to nothing, which the closed surface required of them.

They fall. Each ring reads lower than the ring inside it, and the ordering was not put in anywhere. The hop counts came out of walking the wiring. The four values came out of the charge law and the least-energy condition. One of those was a walk and the other was a piece of algebra, and they agree about which port lies further out.

Subtract each value from the one inside it, in one-hundred-and-eightieths: thirty-five take two is thirty-three, two take minus seven is nine, and minus seven take minus ten is three. Over sixty those are eleven, three and one. The fall is steep near the source and shallow far from it, which is the shape of the thing Coulomb's thread was weighing.

The potential crosses zero somewhere between the first ring and the second. A potential around an isolated charge in open space never does that, and this one has no way to avoid it: twelve numbers that must sum to nothing, with the largest of them sitting at the source, have to go negative somewhere further out.

Compare what the same law costs in the ordinary setting. To write down the potential of a point charge in a laboratory you need the distance between two points, and a distance needs a coordinate system to measure in, and a coordinate system needs a space to be laid over. Three assumptions, stacked, before the first symbol goes on the paper, and not one of the three is a component anybody has found on the board.

The input to that calculation was a list of thirty pairs saying which port touches which. Nothing else went in. No coordinate, no distance, no direction, no angle, no space for anything to be at, and the hop counts that organize the answer were themselves read off the same list of thirty pairs. Coulomb's shape came out of a table of who touches whom, and it came out in exact fractions with denominators small enough to check by hand.

One force, worked to the last digit, out of a wiring diagram with thirty entries in it. The same thirty entries have to say how many forces there can be.

## Dials

Chapter five's group was eight moves of a paper square, and eight is a menu. You pick one.

Take instead a dial on the front of an amplifier. Turn it through any angle you like, including a millionth of a degree, and then turn it through another. The result is a setting of the dial, so the moves close. Every turn is undone by turning back, and leaving the dial alone is on the list. Closed, undoable, do-nothing included: a group, by chapter five's three conditions, with the difference that the moves come in a continuous supply rather than in a list of eight.

One thing about a dial that a menu of eight moves does not have: you can turn it by almost nothing. Set it a hair off the do-nothing position and you have a nudge, and nudges behave far better than the moves they generate. Two nudges in different directions can be added. A nudge can be doubled or halved. Every large turn is a small nudge repeated, so the whole continuous supply of moves is determined by the nudges available at the do-nothing setting, which is a much simpler object: a flat space of directions, with as many independent directions as the dials it came from.

The useful number about a collection of dials is how many independent dials it has. That number is its **dimension**. The amplifier has one. So does a sheet of paper spun about a pin stuck through it: every rotation of a plane is that one dial at some setting. A rigid object turned in three dimensions has three, which is why aircraft have exactly three words for it: pitch, roll and yaw.

Counting dimensions is a way of settling questions, and chapter six ran one without calling it that. The independent loops of a wiring diagram are a collection with a dimension of its own, and edges less dots plus one reads it off: thirty less twelve plus one, nineteen. Nineteen was then the number of conditions that had to hold at once for the twelve observers to agree, and two of them refused, both facts settled before anybody drew a field on the wiring.

Two more collections, the ones physicists reach for. Take an object with two slots in it, each holding a complex number, and ask for every reversible change that leaves the total size alone. Then ask the same of an object with three slots. Both counts are short enough to do on your fingers.

Each slot has a dial of its own, which turns that slot and leaves the others where they were. Each pair of slots has two more, which lean the two members of the pair against each other. Two slots: two single dials, one pair, two pair dials, four. Three slots: three single dials, three pairs, six pair dials, nine. Singles and pairs together come to the number of slots multiplied by itself, every time.

One dial in each collection turns every slot together by the same amount, which is a change of overall convention and says nothing about how the slots stand relative to each other. Divide it out. Four less one is three, nine less one is eight, and the count that survives needs nothing but how many slots there are:

$$n^2 - 1$$

Here n is the number of slots and n squared is that number multiplied by itself. Four slots give fifteen. Five give twenty-four. Between eight and fifteen this family gives nothing at all.

## A book on a table

Take a hardback off the shelf and lay it flat, cover upward, spine to your left.

Turn it a quarter turn clockwise on the tabletop. Then flip it away from you, end over end. Note where the cover faces and where the spine points.

Put it back the way it started and do the same two moves in the other order: flip it away from you end over end first, then turn it a quarter turn clockwise on the tabletop.

The cover is face down both times. The spine is nearest you in one case and furthest from you in the other, and the two results are a half turn apart.

A group in which every pair of moves gives the same result in either order is **abelian**, after Niels Henrik Abel, and one containing even a single pair that does not is non-abelian. The amplifier dial is abelian: thirty degrees then fifty is fifty degrees then thirty, and there is no third way for it to come out. Turns of a solid object in three dimensions are non-abelian, and the book on the table is the demonstration. The three-dial and eight-dial collections are non-abelian for the same reason, because three of the dials in each of them are turns of exactly that kind. The single dial is the only abelian collection here.

## Twelve dials on a closed box

Back to the carrier, which is a box with twelve ports and no lid.

Each port holds one reading. A response is whatever the box does when the readings are disturbed: twelve numbers go in, twelve come back, and the box is the thing in between. Reversible means the twelve that came back are enough to recover the twelve that went in, so nothing was thrown away in the middle. Ask for every response of that kind the ports admit. Those responses compose, they undo, and they come in a continuous supply, which makes them dials.

Count them. Twelve, one for each reading the ports hold, and no dial among the twelve leaves all twelve readings alone. Turn any one of them by any amount and some port reads differently afterward, which is the condition that welds the dials to the ports instead of to machinery standing behind them. Every direction of response does something, and there is nothing left over for a thirteenth to act on.

Twelve dials, and how many forces there can be is settled by how those twelve fall apart into independent pieces.

## The swap and the center

One move of the twelve ports came out of chapter fourteen indifferent to order: the swap that trades every port for its opposite. Do it before any of the sixty rotations or after, and the ports finish in the same arrangement either way. Twelve ports can be shuffled 479,001,600 ways, and of the 479,001,599 that move something, the swap is the only one that gives the same result in either order with all sixty rotations.

Collect, in any group, the moves that commute with every move in the group. That collection is the **center** of the group. Chapter nineteen used the same word for the questions inside an algebra that commute with every other question, which is where records live, and the property is identical, commuting with everything in sight; the difference is what is being commuted, an algebra of questions there and a group of moves here.

Groups differ enormously in how much of themselves lives in the center. Settings of the single amplifier dial all commute with each other, so that group is its own center, every last move of it. Turns of a solid object in three dimensions go the other way: the book on the table showed that a quarter turn and a flip disagree about order, and if you work through the possibilities the only turn that commutes with every other one is the turn that does nothing. Its center holds a single move, and the move is nothing happening, which is zero dials wide.

Go back to the clause chapter four spent a section on. Two parties who cannot look a rate up anywhere outside the room have to manufacture one out of what both of them hold. The same clause binds the twelve ports harder, since nothing anywhere holds a preferred labeling of them. Any change in how the readings are labeled has to be a change the box itself can perform.

Take that concretely. The twelve ports carry no names of their own. Numbering them is something an observer does on a piece of paper, a second observer numbers them differently, and the wiring is indifferent to both. What the box can do is turn its own dials. So the relabelings with any standing are the ones the dials produce, and applying one is a three-step operation built out of moves the box has: turn the dials to the new labeling, do whatever was going to be done, turn them back.

Feed a central move into that. It commutes with the relabeling, so the relabeling and its undoing slide together and cancel, and the move comes back out untouched. Every change of description leaves the center exactly where it was.

The sixty rotations of the wiring are changes of description, and each of them is an inside job by that same clause. So the center has to sit inside the part of the twelve readings that all sixty rotations leave fixed.

That part is easy to find, because chapter fourteen showed the sixty rotations carry every port to every other port. A set of twelve readings that survives all sixty unchanged has to read the same at every port. One number, repeated twelve times, and one number is one dial.

The center of the twelve is at most one dial wide.

## One, three, eight

A continuous group of this kind comes apart into its center and a remainder with no center in it at all. The center is one dial or none. The remainder is eleven dials or twelve.

A group that cannot be broken into independent pieces is called **simple**, and the sizes a simple group is allowed to have are a fixed list rather than a matter of taste: three, eight, ten, fourteen, fifteen, twenty-one, and on upward. Below twelve, the list holds exactly three entries. There is no simple group of size four, or five, or six, or seven, or nine, or eleven, or twelve.

So the count of forces is a sum to be hit exactly. The parts come off that list, the total to be reached is twelve less whatever sits in the center, and nothing may be left over at the end. Two targets, then, twelve and eleven. Take a pencil to both, using threes, eights and tens, as many of each as you like.

Suppose the center is nothing at all. Then twelve dials have to be assembled out of pieces of size three, eight and ten. Ten leaves two to find, and the list starts at three. Eight leaves four, and four is not on the list. Threes run three, six, nine, twelve, and the fourth one lands exactly. Three plus three plus three plus three is the only way to do it: four copies of the smallest simple group there is.

Check whether the wiring has room for four copies. Pair every port with the port opposite it, which gives the six axes of chapter fourteen. Each axis carries two readings, splitting into the part that is the same at both ends and the part that is equal and opposite at the two ends: six even readings and six odd ones. Among the six even readings sits the uniform one, the same number everywhere, and taking it out leaves five. The six odd readings fall into two separate blocks of three. So the twelve come apart as one, three, three and five, and no rotation of the wiring mixes any block with any other.

Four blocks of three are needed. Two exist. The rotations cannot get around that by shuffling four copies among themselves either, because a shuffle of four things would require the sixty rotations to fall into four families of fifteen, and no fifteen of the sixty close up into a group of their own.

So the center is exactly one dial, the remainder is eleven, and eleven has to be built out of three, eight and ten. Ten leaves one, and one is smaller than anything on offer. Two eights overshoot to sixteen. Threes run three, six, nine and then jump to twelve, stepping over eleven without touching it. Eight leaves three, and three is the first entry on the list. Three plus eight, and no other way at all.

$$1 + 3 + 8 = 12$$

One abelian dial, a three-dial piece and an eight-dial piece, with nothing left over and nothing borrowed. Two sums were tried, and the largest number that came up in trying them was sixteen.

Every simple group in the catalogue died at the same step, and all of them died before a single particle was looked at. A simple group's center is zero dials wide, by the definition of simple. This wiring hands over one central dial that no relabeling can touch, and it hands it over before anybody asks what the forces act on. Howard Georgi and Sheldon Glashow proposed in 1974 that the three forces sit inside one larger simple symmetry, and a carrier holding a central dial has nowhere to put them.

The allocation is fixed too, and not by size alone. The uniform reading, the one thing all sixty rotations leave alone, is the central dial. One of the two odd blocks of three carries the three-dial piece. The remaining block of three and the block of five carry the eight-dial piece between them, and they carry it because the eight dials of the three-slot group turn among themselves the same way those eight readings do under the rotations, which is a stronger match than eight equals three plus five.

Three pieces, and the world has names for them. The single abelian dial is the phase that chapter twelve's relabelings turned at each port, and its field is the one whose potential falls off as seven thirty-sixths, one ninetieth, minus seven one-hundred-and-eightieths and minus one eighteenth. That is electromagnetism. The three-dial piece is the weak force, which is why beta decay happens and why the Sun burns slowly enough to be worth living near. The eight-dial piece is the strong force, and its eight dials are the eight gluons that hold a proton together.

Look at the full list of what went into that. Twelve ports and thirty seams, which came out of chapter fourteen's counting. Sixty rotations, which came out of the same wiring. The requirement that every relabeling be an inside job, which is chapter four. The splitting of twelve readings into blocks of one, three, three and five, which is arithmetic on the same sixty rotations. And the list of sizes a simple group is allowed to have, which was settled by Wilhelm Killing in four papers between 1888 and 1890, and made rigorous by Élie Cartan in his thesis of 1894, for reasons that had nothing whatever to do with forces.

The count consulted no mass, no coupling strength, no decay rate and no accelerator run. Every input to it was a whole number obtained by counting something on a wiring diagram or looking something up in a classification finished before the electron was discovered. The forces come out one, three and eight.

Which opens a gap, because a force with nothing to push on is a dial wired to nothing. Each of the three pieces turns something, and the machine as counted holds nothing for them to turn. Whatever carries those charges has to fit inside the same twelve readings the dials were counted from, since there is nowhere else on the carrier to put it, and has to fit without disturbing the count that has just been made.

The catalogue of things that do exist makes that look difficult. The electron carries one unit of electric charge. The down quark carries one third of a unit and the up quark carries two thirds, and every quark carries one of three color charges besides, and the whole arrangement is repeated three times over at rising masses. Thirds, in a world assembled out of a table of whole numbers. And three of everything.