# 5. Why Do You and I See the Same Table?

Bob's kitchen table is oak, scratched, and a hundred and twelve centimeters along its long edge. At two in the morning he crosses the kitchen in the dark to get a glass of water, and he does it by memory and by hand: four steps from the doorway, the cold edge of the worktop under his left palm, then the corner of the table, which he takes wide because he has caught his hip on it before. He never sees it. Nothing reaches his eyes at all.

Alice, coming in behind him, turns on the light and sees the whole thing at once.

Two accounts of the table, and the physical signals behind them have nothing in common. Alice receives photons, scattered off the oak in a band of wavelengths a few hundred nanometers wide, focused onto a patch of light-sensitive cells at the back of each eye. Bob receives mechanical deformation of the skin of his palm, transduced by pressure receptors that would not respond to a laser shone directly at them. No signal either of them receives is available to the other. The two channels do not overlap anywhere in the physics.

They agree about the table instantly and without discussion, and neither of them experiences the agreement as an achievement.

Chapter four left them agreeing and did not say what they were agreeing about. The comparison threw away both descriptions. Something came through it anyway.

## What survives

Start with what does not survive.

The oak is brown. That belongs to Alice, and Bob has no version of it, because brown is what her visual system does with a distribution of wavelengths, and no amount of running your hand along a table delivers it. The oak is slightly waxy, and cooler near the window. That belongs to Bob. Neither of them is wrong, and neither of them can hand their version across.

What both of them have is that the edge runs straight for a hundred and twelve centimeters and then turns through a right angle. That the far corner stands a hundred and twelve by seventy from the near one. That the surface is rigid, so the distance between any two points on it stays what it was while you walk around it. That there is a hard flat obstruction between the doorway and the window, and that walking into the space it occupies will hurt.

Every item on that list is a relation between one part of the table and another. None of them is a property of the table on its own. The rigidity is a statement about pairs of points, the right angle is a statement about two edges, and the hundred and twelve centimeters is a comparison with a metal bar in a vault outside Paris, or with whatever has since taken over that bar's job. Two channels with no physics in common delivered the relations and dropped everything else.

A quantity that comes out the same under a change in how you look is called an **invariant** of that change.

The order matters here, and it is the order most treatments get backwards. The invariant comes first. You find out what survives, and only then do you have any business asking what it survived.

## A square and a pencil

Take a sheet of paper and draw a square on it, big enough to handle, five or six centimeters. Mark one corner with a pencil dot, and mark the same corner on the reverse side so you can find it after the square has been face down.

Cut it out. Then find every way to pick it up and put it back down so that the square occupies exactly the same patch of paper it did before. The outline has to land on the outline. Where the pencil dot ends up is your business.

Most people find four and stop. The four are the rotations: leave it alone, turn it a quarter turn, turn it a half turn, turn it three quarters. Then they remember they are allowed to turn it over, and there are four more, one flip about the vertical line through the middle, one about the horizontal line, and one about each of the two diagonals.

Eight. There are exactly eight, and it is worth spending the minute it takes to be sure of that, because a count you have made with your hands is a different kind of possession from a count you have been told.

Three things about those eight can be checked with the paper in your hand.

Do one move, then do another. Whatever you did, the result is one of the eight. A quarter turn counterclockwise followed by a flip about the diagonal that runs up to the right is a flip about the horizontal line, and you can watch the dot arrive to confirm it. There was nothing else it could have been, because the square is back on its outline, and there are only eight ways to be back on your outline. The set is closed. Nothing you can do with these moves takes you outside them.

Every move can be undone by another one on the list. The quarter turn is undone by the three-quarter turn. Each of the four flips is undone by itself. No move in the set strands you.

And leaving it alone is on the list. It looks like a formality, and it is load-bearing: without it, "undo the quarter turn" has no result to name.

Checking the closure by hand means checking sixty-four pairs of moves, and an eight-by-eight grid holds them all at once, the eight moves along the top and the same eight down the left edge. Each cell answers one question: do the move on the left, then the move along the top, and write down which of the eight you have ended up with. Sixty of the sixty-four are filled in below.

Eight short names first, because writing "the flip about the diagonal that runs up to the right" sixty-four times is worse than learning a letter. Each rotation gets the number of quarter turns it is, counted counterclockwise: 0 for leaving the square alone, 1 for the quarter turn, 2 for the half turn, 3 for three quarters. Each flip gets the line it is a flip about: V for the vertical line down the middle, H for the horizontal line across it, S for the diagonal that leans like a slash, running from the bottom left corner up to the top right, and B for the diagonal leaning the other way. The four lines belong to the table rather than to the square, so they stay put while the square moves. Rows are the first move, columns the second.

|then|0|1|2|3|V|H|S|B|
|-|-|-|-|-|-|-|-|-|
|0|0|1|2|3|V|H|S|B|
|1|1|2|3|0|S| |H|V|
|2|2|3|0|1|H|V|B|S|
|3|3|0|1| |B|S|V|H|
|V|V|B|H|S|0| |3|1|
|H|H| |V|B|2|0|1|3|
|S|S|V|B|H|1|3|0|2|
|B|B|H|S|V|3|1|2|0|

Four cells are empty. Pick the square up and fill them in.

Two of the four hold the same pair of moves in opposite orders. Row 1 under column H is a quarter turn and then a flip about the horizontal line. Row H under column 1 is the flip first and the quarter turn after. Both answers are diagonal flips, and they are flips about different diagonals, and the pencil dot settles which is which: one order leaves the dot in the corner it started from, and the other sends it to the corner diagonally opposite.

The other two blanks are quieter. Row 3 under column 3 is a three-quarter turn done twice, which is six quarter turns, which is a half turn. Row V under column H is a flip about the vertical line followed by a flip about the horizontal, and the square, turned over twice and lying face up again, comes back rotated through a half turn. Two mirrors make a rotation.

Count along any row of the finished grid. Each of the eight moves appears exactly once, and the same holds down every column. The row for the quarter turn shows why. Two different second moves cannot land the square in the same place there, because a three-quarter turn stuck in front of both sequences cancels the quarter turn and leaves the second move standing on its own, so those two second moves would have to have been one move. Eight second moves, eight different answers, and only eight moves available for the answers to be. All eight are in the row, once each, and nothing anywhere in the grid produces a ninth thing.

Closed, undoable, and containing the do-nothing move. A collection of changes with those properties is called a **group**, and mathematicians have been calling it that since the 1830s. The word arrives last on purpose. The object was in your hands before it had a name, and it would have been there if nobody had ever named it.

The three properties fit on one line of symbols. Write G for the collection of the eight moves, g and h for any two moves in it, and a dot between them for doing g and then doing h. Write e for leaving the square alone, and g with a raised minus one for the move that undoes g. The rounded E is read as "is one of". The conditions come in the order closure, do-nothing move, undo, because the undo is defined by what it gets you back to.

$$g \cdot h \in G \qquad e \cdot g = g \cdot e = g \qquad g \cdot g^{-1} = g^{-1} \cdot g = e$$

The first says the answer to any pair of moves is on the list. The second says that e changes nothing whichever side of another move it sits on, and the third that every move has an undo which works from either side. The do-nothing move is called the **identity**, the undo is the **inverse**, and one further condition is usually listed alongside them and costs nothing here: with three moves in a row, grouping the first two or the last two makes no difference, because the square gets picked up three times in the same order either way.

A group is a set with one operation obeying those conditions, and nothing in them mentions paper. The eight moves of the square are the **dihedral group of order eight**, order eight for the number of moves in it and dihedral from the Greek for two-faced, which here is a square with a front and a back. The rotations of an ammonia molecule satisfy the same three conditions, and so does every way of reordering a deck of cards, so anything proved from those conditions alone holds for the molecule, the deck and the paper square without being proved three times.

One more thing to see with the paper. The dot moves. Under the eight moves it visits all four corners, and every corner is reachable from every other. The set of places a thing can be sent by the moves in a group is its **orbit**, so the orbit of a corner is the four corners. The orbit of the midpoint of an edge is the four edge midpoints, which is a different set of four, and no move carries a corner to an edge midpoint.

The center of the square has an orbit of one. It sits where it sits under every one of the eight, and nothing you can do to the square touches it.

That is an invariant, in your hand, on a scrap of paper. Something the moves leave alone.

## Chapter four, from the other side

Look again at the three conditions.

Chapter four asked what has to hold between two record-keepers before they can be said to share a world, and it came back with conditions on the translations between their records. A translation had to be computable from inside, using only what the observer holds. It had to be undoable, so that going from Alice's account to Bob's and back returns Alice's account. Translations had to chain, so that a route through a third party gives the same answer as the direct one. And every closed loop of them had to come back where it started.

Those are the conditions on the paper square. Undoable is the inverse. Chaining is closure. Coming back where you started is the do-nothing move sitting at the end of every loop.

The overlap rule of chapter four and the symmetry of a square are one thing seen from two directions. The translations between observers form a group, in the way the moves of a square do, and the question "what do Alice and Bob agree about" is the question "what does that group leave alone".

Which gives the sentence this book runs on.

A symmetry is a change that leaves a named invariant alone. So a symmetry is the claim that part of your description was never physical.

Brown is not physical in this exact sense: it is what changes when you switch from Alice's channel to Bob's, and the world does not change when you make that switch. The distance from corner to corner does not change, and that one is the table.

## Nought degrees

The clearest case of a quantity everybody treated as a fact about the world, until it was subtracted, is the one printed on every map.

Latitude has an origin the planet supplies. The equator is where the spin axis says it is, and any surveyor anywhere can find it without asking permission. Longitude has no such thing. The Earth is very close to symmetric about its spin axis, which means that turning the whole planet through some angle about that axis leaves every measurable relation between places exactly as it was. Rotating everything is a symmetry, so the number you attach to a meridian is a labeling, and no experiment can read it off the ground.

For three centuries this was treated as an unsettled question rather than an empty one. The French measured from Paris, the Spanish from Cádiz or Tenerife, the Americans from Washington, and captains carried tables for converting between them. In October 1884 forty-one delegates from twenty-six nations sat down in Washington to fix it, and on the thirteenth they voted: twenty-two in favour of Greenwich, one against, which was San Domingo, and two abstentions, which were France and Brazil. France went on printing charts from Paris for another quarter of a century.

What the conference settled was a convention, and everybody in the room knew it. The tell is in the mechanism, because there is no other way to settle it. You cannot measure which meridian is the real zero, and you cannot be wrong about it either. What the vote bought was the ability to compare two charts without arithmetic, which is worth having, and which is a different kind of thing from a fact about the planet.

Distances between places survived the vote unchanged. Distance is the invariant, the labeling was never physical, and the world's response to twenty-two votes was to carry on exactly as before.

Every symmetry anybody has found in physics is that same operation run on a bigger description. Here is a change you can make to how the world is written down. Here is what survives it. Whatever did not survive was never a fact about the world, and the fact that people wrote it down for two hundred years is a fact about people.

## Vienna, orbit, and the Cretaceous

There is a symmetry so large that almost nobody notices they are assuming it.

A physicist in Vienna measures how fast a ball falls. A physicist on the International Space Station, four hundred kilometers up and moving at seven and a half kilometers a second, does the same experiment in a way that accounts for the local conditions. They get different readings, because the conditions differ, and the same law relating the readings to the conditions.

Do it again a hundred and fifty kilometers to the west. Same law. Do it in the Cretaceous, if you can arrange the transport: same law, and a version of this can actually be checked, because light arriving from a quasar eleven billion years old carries the spectral fingerprint of hydrogen behaving exactly as hydrogen behaves in a lamp in a basement in Vienna.

The laws do not depend on where you are. They do not depend on when you are. Those are two symmetries. Shifting everything sideways in space leaves the laws alone, and shifting everything along in time leaves the laws alone. Galileo's cabin is a third change of the same kind: put the whole laboratory into smooth motion and the laws do not notice that either, which is why the flies stayed where they were.

Strip them out and watch what happens to the idea of an experiment. If the law of falling depended on where you stood, a result obtained in Vienna would say nothing about a ball in Graz, and the Graz result would say nothing about the next street. If it depended on when, a measurement made on Tuesday would say nothing about Wednesday, which is to say it would say nothing at all, because a measurement is finished before anybody can use it. Repeating an experiment would be pointless. Publishing one would be a form of gossip. The word "law" would have no work to do.

Every experimental result anybody has ever quoted is leaning on those two symmetries. They are the reason a number obtained once is worth writing down.

And they pay. Emmy Noether proved in 1918 that a continuous symmetry of the laws always comes with a conserved quantity, exactly one, determined by the symmetry. Shifting in space gives you a quantity that never changes, and the quantity is momentum. Shifting in time gives you another, and it is energy.

Momentum is conserved because the laws do not care where you are. Energy is conserved because they do not care when. Those are two statements of one fact, about two directions. The procedure that turns a symmetry into its conserved quantity is exact and mechanical, and running it takes more mathematics than a paper square supplies.

## The letter to Chevalier

The paper square is the whole of what Évariste Galois was working on, run on something harder.

He was asking which polynomial equations can be solved by a formula built from the coefficients using addition, multiplication and the extraction of roots. The quadratic has such a formula, and every schoolchild is made to memorize it. The cubic and the quartic have them, found in Bologna in the 1500s, ugly and real. The fifth degree had resisted for two hundred and fifty years.

Galois stopped looking at the equation and looked instead at the symmetries of its solutions: the ways you can permute the roots among themselves without disturbing any relation between them that the coefficients can see. Those permutations form a group, in exactly the sense the square's eight moves do. He then showed that whether a formula exists is a property of that group, and that from the fifth degree onward the group is the wrong shape.

The submissions went badly. Cauchy took the first memoir in 1829 and it never appeared. Fourier took the second in 1830, for the Academy's Grand Prix, and died in May of that year with the manuscript among his papers, where it was not found. Poisson took the third in 1831 and reported in July that the reasoning was neither sufficiently clear nor sufficiently developed for the Academy to judge it.

On the night of 29 May 1832 Galois wrote to his friend Auguste Chevalier, setting out the results and asking him to make a public request: that Jacobi or Gauss be asked to give an opinion, not on whether the theorems were true, but on whether they were important. In the margin of one page, beside an argument he had not finished, he wrote that there was something to complete in the demonstration and that he did not have the time.

He was shot in a duel the following morning and died on 31 May, aged twenty. The letter was published in 1846. What was in it is taught to undergraduates in their first year, which is the ordinary fate of an idea that turns out to be the right one, and it is why the eight moves of a piece of paper have a name.

## The size of an atom

A notion that applies to everything explains nothing. Symmetry has a boundary, and the boundary can be found from where you are sitting.

Shift everything in space: symmetry. Shift everything in time: symmetry. Rotate everything: symmetry. Change the size of everything, so that every length in the universe is doubled and nothing else about the arrangement is altered: the world notices.

The evidence is in the room. Atoms have a definite size. A hydrogen atom is about one ångström across, a ten-billionth of a meter, and that number is the same in Vienna, in orbit and in the Cretaceous. It is a length built into the physics, and a world without scale symmetry is exactly a world in which such a length can exist. If doubling every length left the laws alone, no size could be preferred, and nothing would distinguish the atom you have from one twice as big. Atoms would come in every size, and chemistry, which depends on atoms of a given element being interchangeable, would not be available.

The consequences are structural rather than decorative. Galileo worked one of them out in 1638, in the *Two New Sciences*, and illustrated it with a drawing of a bone. Scale an animal up by a factor of ten in every direction and its weight goes up by a factor of a thousand, while the cross-section of the bones that carry the weight goes up by a hundred. The scaled animal is ten times worse off than the original, and its legs snap. Which is why you cannot design a bridge by building a small bridge and then being more ambitious. Every engineer who has tested a model in a wind tunnel has had to do arithmetic to correct for the fact that the model is not a small copy of the world, and that arithmetic is a bill for a symmetry the world does not have.

So symmetry is a claim with content, and the way you can tell is that it is false about scale.

## What the table was

Alice and Bob agree about the table because there is a change, from her channel to his, that leaves something alone, and what it leaves alone is the geometry: the edges, the right angle, the hundred and twelve centimeters, the rigidity. Everything else in either account belongs to the observer that produced it.

That is the same operation as the paper square, where the eight moves leave the center alone. It is the same operation as chapter four's overlaps, where the translations between records leave the shared content alone. It is the operation physics has been running since the 1830s, and each time it has been run, something everybody had taken for a property of the world has turned out to be a property of the description, and has been subtracted.

There is a question left over, and Alice and Bob cannot answer it between them, because two people who agree have very little to disagree about.

Put a third person in the room. Alice agrees with Bob. Bob agrees with Charlie. Charlie agrees with Alice. Every pair compares records, finds a translation, and comes away satisfied. Each pair shares one invariant, and each pair is content.

That has the sound of a settled matter, and three coins in three separate rooms are enough to take it apart.