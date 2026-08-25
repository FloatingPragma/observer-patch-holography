# 12. Why Is the World Made of Relations and Not of Things?

In the spring of 1918 Hermann Weyl sent Albert Einstein a paper arguing that the universe has no fixed unit of length. Einstein disagreed with it, presented it to the Prussian Academy of Sciences in Berlin anyway, and his objection was printed along with it, so that the proposal and the reason it could not be true went to press a few pages apart.

Let every point of space set its own unit of length, the way every workshop sets its own gauge blocks, and comparing two lengths at different places then needs a rule for carrying a unit from one to the other. Weyl's rule was the electromagnetic potential. Electricity and magnetism would be the bookkeeping of a freedom that had been sitting unexamined inside geometry, the freedom to recalibrate the ruler at every point, and the German word for setting the scale on an instrument is Eichung, calibration. The word that came back through the translation is gauge.

A workshop's gauge blocks are the object the word points at. In 1896 a Swedish machinist named Carl Edvard Johansson worked out that a set of 102 hardened steel blocks, finished flat enough that two of them cling together when wrung face to face, will combine to give any length between 0.5 and 300 millimeters in steps of a hundredth of a millimeter. He patented the set in 1901, and within twenty years it was what the world's engineering shops measured against. Every shop had its own box, every box was checked against a better box somewhere else, and Weyl's proposal was that the universe skips the last step: every point keeps its own box and no box anywhere is the master.

Einstein's objection was short. If carrying a unit of length depends on the route it takes, two identical objects sent by different routes arrive with different sizes. Atoms are objects of that kind and the light they emit is their scale marking, so hydrogen that spent a billion years in the arm of a distant galaxy would emit at frequencies displaced from hydrogen in a Berlin discharge tube by an amount depending on where it had been. Astronomy consists in large part of recognizing elements in starlight by matching their lines against laboratory samples. The lines match.

Weyl argued back for years and gave the scale geometry up in the early 1920s, and the word outlived it. In 1927 Fritz London noticed that the factor Weyl had written works exactly as intended if it is made imaginary, at which point it stops rescaling the size of anything and starts rotating the phase of a quantum wave function, which no instrument reads off. Weyl rewrote the theory in those terms in 1929, in a paper on the electron and gravitation. The phrase gauge invariance came through the rewrite with its modern meaning attached. What had to change was the quantity chosen freely at each point, which has to be one that no measurement reaches. The length of a ruler fails that condition in any laboratory with a spectrograph.

Twelve observers wired into thirty seams have a freedom of exactly that kind, for reasons with no geometry in them. Each keeps its state in its own labels, and nothing outside the twelve fixes the labeling, because there is nothing outside the twelve. Watch what the settling then produces. The repairs run, the disagreement drains out of the network, every one of the thirty seams reports itself content, and if you walk round the twelve afterwards and read off what each of them is holding, nothing in the run has pushed those numbers together. The seams converged and the labels did not. What is it that the twelve of them agree about?

## Three teams, one algorithm

Three teams are handed the same specification and told to implement it. One writes Python, one writes Java, one writes Rust. When the programs are done, the only comparison anybody can perform is to feed all three the same inputs and check the outputs, because that is the only place the three implementations touch. Everything else about them, variable names, memory layout, which loop runs backwards, is invisible to that test. Each team can choose it without consulting the others.

Chapter five made a symmetry the claim that part of the description was never physical. Here the choice gets made three times, independently, by parties who never see each other's source. A symmetry with one global choice is a claim about the whole world. A symmetry whose choice can be made separately in every place is a **gauge** symmetry, and twelve observers make twelve of them.

The picture has a client in it. The client is the part to delete. Somebody outside the three teams wrote the specification and knows what the program is supposed to do. Nothing plays that role here. The input-output behavior of the twelve observers is the only thing about them that exists. Each of them wrote itself. What is left when the client goes is three implementations whose agreement has nothing to appeal to but each other, and no envelope arriving in the morning with the requirements in it.

Local choice has a price and it is charged at the seams. Two patches sharing a seam compare port readings written in two private conventions. The comparison needs a translation. Chapter four supplied the shape of it, with two firms settling one account in two currencies and reading their exchange rate off the transaction they both recorded. Put a dictionary of that sort on every seam. The consistency condition stops being "my reading equals yours" and becomes "my reading equals this seam's dictionary applied to yours."

A rule that lives on an edge and translates one end into the other is called a **connection**. Chain the dictionaries around a closed loop and you come back holding a translation from your own conventions into your own conventions, which is chapter six's holonomy, and the one thing it tells you is whether the trip returned you unchanged.

Notice where the information went. Each observer's own labels are its private business. The thirty dictionaries are the only content the arrangement holds that anybody can check. In the large runs the freedom each patch holds is one of the six ways of permuting three labels, and every seam's dictionary is one of the same six. If a patch decides to swap the names of two of its labels, the dictionary reading into it picks up that swap, the dictionary reading out of it picks up the reverse, and every round trip through the patch comes back with what it came back with before. Redenominate one city's currency and each rate into and out of it changes, while the round-trip factor does not move.

## One number for the whole arrangement

Take the wiring itself, and give each observer the freedom London's rewrite left standing, an angle on a dial rather than one of six permutations. Twelve ports, thirty seams, five seams meeting at each port. A field on that graph is a number on each of the thirty seams. A relabeling is a number at each of the twelve ports, and what it does to a seam is add the difference between the numbers sitting at the seam's two ends.

Some fields can be manufactured that way out of nothing and some cannot. Every field splits into those two parts. The part no relabeling could have manufactured is the collection of loop readings. That part is the field strength, because it is the part a measurement can reach.

The split shows up in what two observers can argue about. A seam value is something each of them wrote in private. The gap between their two numbers closes if either one changes a convention. A loop reading is the result of a round trip. No convention either of them can change will move it. Relabel one port and all five seam values touching it move, while the round trips through it come back exactly as before.

Then a source. Charge moves and the source has to move with it, so the source is a current: one number on each seam, saying how much runs along it and in which direction.

Score the arrangement with a single number: the energy held in the field strength, less the thirty products of each seam's value with the current running along it. A number scored over a whole arrangement like that is called the **action**.

Apply a relabeling and work out what the score does. The first term does nothing, since the field strength was built to ignore relabelings. The second term moves, and the amount it moves by is one line of arithmetic: for each of the twelve ports, take the number that port chose and multiply it by the imbalance in the current there, meaning what flows out along its five seams less what flows in, then add the twelve products.

Both directions of the theorem are sitting in that sum. If the current balances at all twelve ports, every product is zero and no choice of labels moves the score, whatever numbers the observers pick. If the current fails to balance at even one port, the observer at that port can pick a number that makes its product nonzero, and the score moves. Send three units out of one port along its five seams and two units back in. The imbalance there is one. Let that observer add seven to its own labels, and the score of the whole arrangement moves by seven while the eleven other products sit at zero.

The second direction is the one doing work. A current that fails to balance somewhere gives the arrangement a score that depends on conventions the twelve observers picked privately, and any one of the twelve could then move that score over lunch by rewriting its own labels, having changed nothing whatever about the world. An arrangement like that has no score at all.

Charge conservation and gauge symmetry are one statement read in two directions. Both directions came out of the same twelve products.

Look at what went into that and what came out. In went a wiring diagram of twelve ports and thirty seams, a rule saying which changes count as relabelings, and a source. Out came a restriction on what a source is allowed to be, with a conservation law as its content. Nothing about electricity was assumed anywhere in the argument. The only geometry in it is a list of which port is wired to which.

## Nineteen loops

Count what the labels can reach. There are twelve dials to turn, one at each port, and turning all twelve by the same amount changes no seam, so eleven of the twelve directions do anything at all. Thirty seam values, eleven directions of pure convention, nineteen left over.

Nineteen turned up in chapter six as well, off this same wiring, out of a rule with three pieces of arithmetic in it: take the edges, subtract the dots, add one. Thirty seams, twelve ports, and thirty minus twelve plus one is nineteen.

Run the rule on two drawings small enough to check by eye. Sketch six observers in a row, which is five seams and six dots. Sketch two triangles sharing an edge, which is five lines and four dots. Do both subtractions in your head and keep the second answer.

Count the loops in the second drawing and the eye finds three: the left triangle, the right triangle, and the four-sided circuit around the outside. Walk the left triangle clockwise and then the right one clockwise. The shared edge gets crossed once in each direction and cancels. What survives of the two trips is the circuit around the outside. The third loop is the first two performed in sequence, two of the three are independent, and the drawing gives no hint as to which two. The eye overcounted by one on a figure with five lines in it, which is the argument for doing the arithmetic on a figure with thirty.

Three drawings, one rule, and the arithmetic for all of them fits on a line.

$$5 - 6 + 1 = 0 \qquad 5 - 4 + 1 = 2 \qquad 30 - 12 + 1 = 19$$

The row of six has nothing to walk around, the two triangles carry two independent loops where the eye found three, and twelve observers wired into thirty seams carry nineteen.

Nineteen does two jobs. It counts the sums that have to come out zero before the twelve can settle on anything, and it counts how much of the wiring is physics rather than paperwork.

The two nineteens are the same nineteen. Strip every label choice out of a field and what survives is nineteen numbers, one for each independent loop, and nothing else about the thirty seam values survives at all. Two fields with the same nineteen loop readings are one field written down twice.

Those nineteen numbers are also the only things any two of the twelve can both compute and both check. A round trip begins and ends at one port, so the observer sitting there can work it out from the dictionaries alone, without being told what anybody else calls their own state and without taking anybody's word for anything. The dictionaries sit on the seams where both ends can read them, and the shortest of the round trips chains three.

## The laziest field

Give the numbers on the seams a different job. They stop being anybody's private labeling and become the field itself. A field costs energy: square each seam's number and add the thirty of them. The law tying charge to that field reads at the ports: the net outflow along a port's five seams equals the charge sitting there. Give it a list of twelve charges and it answers with a family of fields rather than one, because adding a pure circulation around any loop leaves every port's outflow alone. The family is nineteen-dimensional, one dimension for each independent loop.

A law that answers with a nineteen-dimensional family has not finished the job. Two further conditions each pick out one member of it. They pick the same member. Exactly one field in the family has zero reading on every loop. Exactly one has the least energy. They are one field for a single reason: the energy of any solution equals the energy of that one plus the energy of the difference between the two, and the difference is pure circulation, which costs energy unless there is none of it.

The field around a charge is the solution that circulates nowhere. Every other member of the family is that field with a circulation laid on top of it, which costs energy and transports nothing. The operator that selects against it is chapter nine's repair rule, my value against the average of my neighbors, running on the twelve ports of this same graph, and the only reading that rule leaves alone is the one that is the same at every port.

On a twelve-port graph all of this is exact arithmetic rather than approximation. Put one unit of positive charge at a port and one unit of negative charge at a neighbor. Send the whole unit straight down the seam that joins them and the arrangement costs one, one squared being one. Let it settle into the field that circulates nowhere and the same two charges cost eleven thirtieths, and the nineteen thirtieths that went missing was circulation delivering nothing to either port. Eleven over thirty and nineteen over thirty, off the same wiring that gave eleven directions of convention and nineteen loops.

## One hundred and seven thousand nine hundred and fifty-nine seams

The same behavior is visible in the large runs, where nothing about it was arranged. In the sixty-five-thousand-patch run of chapter ten, 65,536 observers wired into 390,924 seams, the raw labels began by disagreeing across 206,910 seams and ended disagreeing across 314,869, while the transport-corrected disagreement began at 326,047 and ended at zero, with 27,250 seam dictionaries rewritten on the way.

The number of seams whose two ends hold different labels went up by 107,959. The disagreement anybody inside could detect went to nothing. Both counts came off the same run.

Count the other way at the end of that run. The seams whose two ends held identical labels numbered 76,055, and the seams whose dictionary was the permutation that does nothing numbered 76,055, and they were the same seams. Two ends of a seam carry the same label exactly when the translation between them happens to be the do-nothing one, which is a coincidence of conventions with no content in it. The other 314,869 seams were in perfect agreement and said so in different words.

What settles is the relation, never the values. It could not have been the values. The labels are each observer's private business, no procedure anywhere in the arrangement reads them, and a process that drove them together would be enforcing a convention that nothing can check.

## The residue in the loop

Chapter six ended with twelve observers, one flipped record, and a leftover of exactly one unit of disagreement that no repair removed and no ordering avoided. Repairing in a different order moved it somewhere else and left the total at one wherever it came to rest. It sat in the loops.

The nineteen loop readings are where it lives. A relabeling cannot reach it, because relabelings leave round trips alone, the way redenominating a currency leaves the arbitrage factor where it was. Repair cannot flush it, because repair works one seam at a time and no single seam holds it. A loop that comes back carrying a residue is a loop with something inside it. No rewriting of anybody's conventions disposes of it, because the conventions have been divided out of the reading before the residue is read.

That residue is charge. The quantity a current has to conserve for the arrangement to have a score at all, the quantity the law at the ports reads as an outflow, and the quantity trapped in a loop that fails to close are one quantity. A loop carrying a residue cannot be talked away as somebody's bookkeeping, because the bookkeeping was divided out of the reading before the reading was taken, and what survives that division has to be paid for in energy and balanced at every port. Weyl proposed in 1918 that electromagnetism is the bookkeeping of a freedom chosen separately at every point. He had the wrong freedom, and the bookkeeping was right.

Thirty seams. Twelve ports, five seams at each. Eleven directions of pure convention, nineteen loops, and one unit of charge sitting in the loops where no relabeling reaches it. Every number in that list was read off a wiring diagram. Nobody has opened one of the twelve observers to see what a port is.