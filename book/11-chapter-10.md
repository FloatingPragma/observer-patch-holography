# 10. Why Doesn't It Matter Who Goes First?

Among the crosswords in the May 1979 issue of *Dell Pencil Puzzles and Word Games* there was a nine by nine grid with some of its squares filled in and the rest blank. It was called Number Place. Every row, every column and each of the nine three by three boxes had to end up holding each digit from one to nine exactly once. The man who devised it was Howard Garns, a retired architect in Indianapolis, seventy-four years old, and Dell did not print his name on it. He died in October 1989. Five years before that, the Japanese publisher Nikoli had picked the puzzle up and given it the name it travels under, which is sudoku, and in 1997 a retired Hong Kong judge named Wayne Gould found one in a Tokyo bookshop, spent six years writing a program that could make them by the thousand, and talked The Times of London into printing one on 12 November 2004. Within weeks every British newspaper had a grid on its puzzle page.

Give the same grid to two people on the same train. One of them works the sevens first, because four of them are showing and they constrain each other. The other goes box by box from the bottom left. They do not speak, they do not look at each other's paper, and neither has the faintest idea what the other has written down. Forty minutes later both grids are full and all eighty-one squares match.

Nothing arranged that. Two different procedures ran in two different orders on different parts of the grid. The agreement at the end was neither negotiated nor checked. The two of them get off at different stops without ever finding out.

The reason is that the finished grid was fixed by the clues before either of them picked up a pen. There are 6,670,903,752,021,072,936,960 ways to fill a nine by nine grid legally, a figure Bertram Felgenhauer and Frazer Jarvis computed in 2005. A printed puzzle is a set of clues that cuts that number to one. In 2012 Gary McGuire, Bastian Tugemann and Gilles Civario finished an exhaustive search and established how few clues can manage it. Seventeen. Sixteen never suffices, on any grid, in any arrangement. Below seventeen the clues stop picking out a single completion. The two people on the train can hand in different grids and both be right.

The clues do half of the work and the pencils do the other half. Every move either of them makes is forced. A square gets a five because five is the only digit its row, its column and its box have left over for it, and once that is true it goes on being true whatever else gets filled in, because filling in another square removes possibilities and never restores one. A forced entry cannot be un-forced by a later forced entry. So the order the two of them work in decides which square gets filled at which minute. It has no way to reach the digit that square ends up with.

Repairs in the network of the last chapter happen when a mismatch is noticed. Which mismatch gets noticed first is whatever the world happens to hand over. There is no schedule, in the sense that there is no fact anywhere about the order. If the endpoint depends on that order, then the world an observer finds is a fact about scheduling. No observer inside can see a schedule.

## Two librarians

Take a shelf of books that ought to be in order and are not. The rule is one line: if you find two neighboring books out of order, swap them. Any pair, whenever you like.

Photograph the shelf. Let one librarian sort it. Then put the books back exactly as the photograph shows and hand the shelf to a second librarian who has been out of the room. Both finish with the shelf in order, which is no surprise, because there is one way for a shelf to be in order. Count their swaps. Each of them will have made exactly the same number. Neither of them chose it.

Count, at any moment, every pair of books in the wrong relative order, wherever on the shelf the two of them happen to sit. A book that belongs before another one and stands to its right is one such pair, whether the two are touching or at opposite ends of the run. Call each of those a wrong pair. Swapping two neighbors that are out of order fixes exactly one wrong pair, because it reverses the relative order of those two books and leaves the relative order of every other pair untouched. So the count drops by exactly one at every swap. It drops to zero. The number of swaps is the number of wrong pairs the shelf started with, which is a property of the shelf and not of the librarian. This is the descent argument of the last chapter, on furniture: a quantity that strictly falls and cannot fall below zero cannot go on falling.

The shelf does something else on the way down. Take a shelf where two swaps are available: books three and four are out of order, and so are books seven and eight. Do the first one. Do the second one instead. The two shelves you get are different from each other. From the first, swap seven and eight; from the second, swap three and four. Both land on the same shelf.

Draw those four shelves with arrows between them and the picture is a diamond: one at the top, two down the sides, one at the bottom where the paths meet. The **diamond property** is the name for that shape. The property it names is that any two moves out of a state can each be followed by one move that brings them back together. Allow the branches any number of moves to find each other again and the condition is called **local confluence**. Allow the branches themselves to be any length, so that two long histories out of one starting state have to be reconcilable, and the name for that is **confluence**.

Three books show what the looser condition buys. Put Chesterton, Cervantes and Bunyan in positions three, four and five. The pair three and four is out of order, the pair four and five is out of order, and the two swaps share a book, so they are not independent the way the first example's were. Swap three and four and the shelf reads Cervantes, Chesterton, Bunyan. Swap four and five instead and it reads Chesterton, Bunyan, Cervantes. Those two shelves are further apart than the first example's pair. No single swap carries either one to the other. From the first, swap four and five, then three and four. From the second, swap three and four, then four and five. Both arrive at Bunyan, Cervantes, Chesterton. Two moves on each branch rather than one, and three swaps down each route, which is the number of wrong pairs those three books started with.

Those last two conditions are separate statements. Local confluence says any two first steps can be reconciled. Confluence says any two histories can be. A history is a hundred moves of drift, and knowing that its opening move was harmless says nothing about where the other ninety-nine took it.

## A snowflake on Triple Divide Peak

Triple Divide Peak stands 8,020 feet up in Glacier National Park in Montana, on the one summit where two continental divides cross. Water leaving that summit on one face reaches the Pacific by way of the Columbia. On another face it reaches the Gulf of Mexico by way of the Missouri and the Mississippi. On the third it reaches Hudson Bay by way of the Saskatchewan. Two snowflakes landing a hand's width apart on the summit finish in different oceans.

Inside any one of those three basins the outcome stops depending on the route. A raindrop falling anywhere in the Columbia basin can take any of a billion routes down, through gullies and side creeks and a reservoir or two, and arrives at the same Pacific as every other drop in the basin. Different paths, one sea. That is confluence, drawn at the scale of half a continent.

Water is also guaranteed to stop. It goes downhill, its height is a quantity that strictly decreases, there is a floor, and every drop that lands on that mountain reaches sea level and stays. Termination promises that the process ends. It promises nothing about ending in one place. A mountain with three basins meeting on its summit is what the failure looks like.

The failure runs the other way as well. Loosen the librarians' rule so that any two neighboring books may be swapped whether or not they are out of order, and the branches go on being reconcilable, since enough swapping turns any shelf into any other shelf. What goes is the stopping. The two librarians can push Bunyan past Cervantes and back again forever, able to agree at every step and never arriving at a shelf to compare.

The mountain misleads in one place. Water runs downhill because there is a downhill, and it stops at the sea because the sea is sitting there waiting to be arrived at. The network has neither. Its settled configuration is the arrangement that satisfies every seam at once, picked out by the constraints in the way seventeen printed clues pick out one grid, before a single repair happens. Repairs delete the arrangements that fail. Nothing travels anywhere. Nobody who finishes a sudoku believes the grid was assembled by the pencil. The same restraint is owed to a network that ends up consistent.

Max Newman supplied the missing piece in 1942, in the *Annals of Mathematics*, in a paper on theories with a combinatorial definition of equivalence. He was a Cambridge topologist whose lectures had sent Alan Turing off to write "On Computable Numbers", and in the year that paper appeared he arrived at Bletchley Park, where he was given a section that everyone called the Newmanry and told to break the German high command's teleprinter cipher with machines. The first Colossus was delivered to him on 18 January 1944.

His result fits in a sentence. If a process always stops, and if every pair of single steps out of a state can be brought back together, then every pair of histories can be brought back together, and every starting state has exactly one endpoint. That is Newman's lemma.

The argument runs from the bottom of the descent upward. Take a state, and suppose that everything reachable below it has a single endpoint. That supposition is safe, because the process descends and cannot descend forever, so you can start at the bottom and work upward. Two moves lead out of your state. Local confluence brings those two moves back together somewhere below. The two states they led to both sit below the state you began at, so each of them has a single endpoint, and so does the state where they rejoined, which forces all three endpoints to be the same one. Your state inherits it. The induction runs over every level of the descent. The descent is finite because the disagreement it counts down is finite.

## Two things to check

Newman's first condition, that the process stops, was settled in the last chapter. Every accepted repair strictly lowers the total disagreement across the network, that total is a sum of terms none of which can go negative, and so it has a floor at zero. A quantity that drops by a positive amount every time and cannot pass zero runs out. The network stops.

The second condition is the one that does the work. Two repairs happen at two different seams. If those seams share no patch, each repair reads and writes numbers the other never touches, and the two orders give identical results, because the repairs are not interacting in the first place. That is the shelf's first example, where books three and four had nothing to do with books seven and eight.

The case that has to be checked is the shelf's second one. Draw three patches in a row. There is a seam between the first and the second, another between the second and the third, and the middle patch sits on both of them. A repair on the left seam moves the middle patch's reading. A repair on the right seam moves the same reading. The middle patch has only the one.

A repair does one thing. It takes what the two sides of a seam expose to each other and moves them together. What each side exposes is overlap data, information about a region both patches can see, and assembling overlapping regions is associative: glue the first two and then the third, or glue the last two and then the first, and the object you finish holding is the same object. Two repairs sharing a patch are two such gluings. The order in which overlapping regions were glued is not recorded anywhere in the result. Whichever repair went first, one further move on each branch brings both branches to the same configuration.

The network has both of Newman's conditions. Every configuration of it settles into exactly one state. That state carries no trace of the order in which the repairs happened.

Schedule independence is what objectivity means here. The twelve-observer universe of chapter six settled onto one arrangement from all 2,048 of its states, under sixteen orders of repair apiece. Newman's two conditions cover the configurations nobody enumerated, and past twelve observers that is all of them.

In most of engineering the order does change the answer. It changes the answer here too, on the smallest arrangement that can show it, which takes two observers and one seam.

Give each of them a single bit and one rule: when you disagree with your neighbor, adopt your neighbor's value. Nothing about that is unreasonable. It is the politest repair rule available, it strictly reduces disagreement, and it touches only the seam it is applied at.

Start them disagreeing. The left one holds false and the right one holds true.

Let the left one move. It takes the right one's value, both hold true, and since they agree neither will ever move again. Go back and let the right one move instead. It takes the left one's value, both hold false, and neither will ever move again.

Two endings, both final, out of one start. Whoever moves first wins, permanently. The two worlds are easy to tell apart, and that is the trouble rather than the consolation: both holding true means the left one moved first, both holding false means the right one did, so the contents of the world are a transcript of a schedule that no rule chose and nothing required.

Confluence belongs to the repair rule rather than to locality. It is not a courtesy extended to whoever asks for it. Leave those two patches wired exactly as they are and change only the rule, from adopt your neighbor's value to something that defers to a fixed criterion, say that whichever of them holds true gives way. From the same disagreeing start there is one move available instead of two, both orders arrive at false, and the fork is gone. Courtesy was the whole problem: a rule in which each party defers to the other leaves nothing to settle the case.

The twelve-observer network has confluence, and chapter six settled it by starting every configuration the arrangement admits, rather than by arguing that the alternative sounded unphysical.

The engineered protocols of chapter eight buy order-independence with sampling and probability. The guarantee they buy is a guarantee about how often. What the observers get instead is two conditions and an induction over a finite descent, which is why their version has no failure rate attached to it.

## What the network is made of

The arrangement in this chapter is one the OPH research team runs as a program, and the program is called OPH-FPE. It is not a numerical model of a continuous world. Every quantity inside it is a whole number or one of the six ways of permuting three labels, every move is exact, and no step anywhere rounds anything off. A run either meets a condition or fails to, and it returns the same verdict on every machine that executes it.

Each published run is archived entire: the configuration it was given, the seeds it drew from, the arrays it started with, the arrays it finished with, the full trace, and the state of all twelve ports on every patch. Beside them sits a verifier that borrows none of the simulator's own code. It builds the six permutations from scratch, recomputes the starting disagreements from the raw arrays, works out the finishing state, and recomputes every fingerprint. A reader with Python and a copy of the archive runs one command and either gets the same answer or finds out that they do not.

The run read here is the level-six rung of a tower built by subdividing an icosahedron over and over. It holds 81,920 patches and 122,880 seams between them, and every patch carries twelve local port slots, three of which route seams while nine stay exposed. At the start, 102,415 of those seams disagree, which is 83.3 percent of them.

## Two cycles

The trace is kept in cycles. A cycle is a round in which every patch with something to repair gets a turn. That count belongs to whoever reads the trace afterward. Nothing in the network counts cycles, nothing waits for the next one, and there is no tick anywhere in the machinery.

A run of a process that failed to be confluent would be an anecdote: what happened to one arrangement on one pass, with no warrant to say anything about the next pass, which would take a different order and could finish somewhere else entirely. Newman's two conditions are what make the last line of this trace worth more than the run that produced it. The lines in between are the accidents of one order. The line at cycle 1 would have read zero under any of them.

Cycle 0 accepts 81,920 repairs and leaves 20,495 seams in disagreement. Cycle 1 accepts the remaining 20,495 and leaves nothing. Of 122,880 seams, every single one agrees, and the trace runs on to cycle 15 without another mismatch posted against it.

Each repair removes exactly one disagreement, which is the librarians' shelf at the scale of a network, and the number of moves the run accepts is 102,415, the number of disagreements it began with. Neither figure has anything in it about the order the repairs came in.

## Sixteen shuffles

The last line of one trace is the last line of one trace. What settles the question is running the same arrangement sixteen times and shuffling the order before every pass.

All sixteen accept 102,415 moves. All sixteen finish with no seam in disagreement. All sixteen arrive at the same finishing state, and the check on that last point is a fingerprint: the whole terminal state is boiled down to one line of characters, and the sixteen lines are the same line.

Sixteen orders. One destination.

The same replay takes 512 places where two repairs share a patch and confirms that the two ways around the diamond meet, takes 512 places where two repairs have nothing in common and confirms that they commute, and relabels the local frames sixteen times over to confirm that nothing in the outcome was resting on the labels. None of those checks turns up a violation.

The records arrive after the arguing stops. Nothing is committed through cycle 6. Cycle 7 posts 47,412 of them, a little under three in five. Cycle 8 posts the rest, 81,920, one for every patch in the network, and from there to the end of the run nothing changes at all.

The spread of those records at the finish is 11.3134. The natural logarithm of 81,920 is 11.3135, and the whole of the gap between those two figures is the price of four coincidences: four pairs of patches hold matching readings, and each of the other 81,912 holds a reading that no other patch in the network holds. Spread counts how many different things a collection holds, on a scale where doubling the number of distinct items adds a fixed amount instead of doubling the figure.

## Agreement first

Records lag agreement by six cycles, and they lag it because of what a record is. A record is something you can go back to and find unchanged. Write one into a neighborhood whose seams disagree, and a repair comes through and revises what you wrote, which makes what you wrote a draft.

The run puts a number on that. A patch's reading becomes a record when it has held still for eight cycles together. Anything in the middle of being repaired fails that test, because being repaired is what failing to hold still consists of.

Read the rule closely and notice what it does not say. It does not tell a patch to look at its neighborhood, or to wait for anybody, or to ask whether the argument nearby has finished. It asks one question, and the question is about the patch itself: has this reading stopped changing. The neighborhood enters only through the back door, because a patch whose seams are in dispute is a patch whose reading keeps getting revised out from under it.

Agreement comes first and memory second. The order is forced by what a record is rather than by any instruction to observe it, and it runs the opposite way to the picture most people carry, in which observers write down what they see and then compare notes and argue about the differences. Through six cycles of this run, every patch in the network is in the middle of a disagreement or fresh out of one, and not one of them has written anything down, and there is nothing anywhere to compare.

## The last line

The last line of the trace reads like this. One hundred and twenty-two thousand eight hundred and eighty seams, every one of them in agreement. Eighty-one thousand nine hundred and twenty patches, each holding a record of its own. Whatever order the repairs ran in, that is where the network lands, and Newman's two conditions are the reason.

Everything settled. Nothing that settled is a number any two patches share. Each patch holds its own reading, taken through its own ports, of the piece of the world it can see. Walk up to two neighboring patches at the end of that run and ask each of them what the world is like, and you get two different answers from two observers in complete agreement.

The network has a settled state. It is written down in none of the eighty-one thousand nine hundred and twenty places where things get written down.