# 36. Why Does the Loop Have No Beginning?

The room was emptying. It was the afternoon of 7 September 1930, the last day of a three-day conference on the epistemology of the exact sciences that Kurt Reidemeister had organized at the University of Königsberg, and the closing item was a roundtable on the foundations of mathematics with Hans Hahn, Rudolf Carnap, Arend Heyting, John von Neumann, Arnold Scholz and Reidemeister himself at the table. Everyone had said their piece over the previous two days: Carnap for Russell's programme, Heyting for Brouwer's, von Neumann for Hilbert's. The proceedings went into the journal *Erkenntnis* the following year, volume 2, pages 87 to 190, with the roundtable transcript at pages 135 to 151, so it is possible to read what was said and in what order.

Near the end, a twenty-four-year-old Viennese logician remarked that one can give examples of statements of arithmetic that are true and that no proof inside the system reaches. Kurt Gödel had come to Königsberg to present a completeness theorem, the result he had built his reputation on, and he left having mentioned in passing the result that closed the question the three days had been about.

The discussion moved on. Almost nobody in the room reacted. One person did: von Neumann cornered him afterward, wanted the construction, went home and thought about it for two months, and wrote to Gödel on 20 November with a further consequence that Gödel had reached first.

The next morning the city of Königsberg made David Hilbert an honorary citizen. He was sixty-eight, he had been born there, he had retired from Göttingen that year, and he addressed the Society of German Natural Scientists and Physicians on knowledge of nature and logic. He closed with six words that were later cut into his gravestone: we must know, we will know. A recording of the ending survives, 373 words read for local radio and running four minutes. The answer had been supplied the previous afternoon, in a smaller room, at a table of seven, one of whom was listening.

What Gödel had was a construction. Hand it a system that can talk about its own formulas and it hands back a sentence about that system, assembled out of the system's own material and pointing at nothing outside it. A world made of observers comparing records leaves a bill that no amount of further comparing pays. It has to settle that bill the same way, out of its own material. Something has to be doing the comparing. Where did the network come from, what is it running on, and what started it?

## A program that reads its own text

Every programmer has wanted the same tool at some point. Something you can point at a piece of code that tells you whether the code finishes or runs forever, and answers by deciding rather than by waiting, correctly, on any program you hand it.

Suppose you have it. Call it the oracle. It takes a program's text along with an input, answers finishes or runs forever, and is never wrong.

Then write the following, which takes about four lines. It reads a piece of program text. It asks the oracle what that program does when it is fed its own text as its input. If the oracle says finishes, it goes into an infinite loop. If the oracle says runs forever, it stops immediately.

Hand that program its own text.

If it finishes, the oracle must have answered runs forever, since that is the only answer that makes it stop, and the oracle was wrong. If it runs forever, the oracle must have answered finishes, and the oracle was wrong the other way. There is no third case, so there is no oracle. Alan Turing published that argument in 1936.

Look at what did the work, because the contradiction is the disposable part. The four-line program is one whose behavior is the opposite of what is predicted about it. Hand somebody the transformation "do the opposite of the prediction" and the construction hands back a program satisfying it. Feeding a description into the thing it describes is a way of solving an equation. The equation is chapter two's: return what you were given.

The calculator in chapter two hunted for a number the operation left alone and found 1.6180339887. The same demand, made of descriptions instead of numbers, has the same kind of answer: a description the transformation returns unchanged. That is a fixed point, and self-reference is the construction that builds one. Chapter two's warning travels with it, because a fixed point can be a place nothing ever arrives at.

## Four rows and a diagonal

The move underneath the oracle is older than computers and fits on an envelope. Georg Cantor published it in 1891, four pages in the first yearbook of the German mathematical society, of which he was the first president. It runs on a numbered list.

Write four strings of four bits, one per line.

row 1:   1 1 0 0 row 2:   0 0 1 0 row 3:   1 1 1 0 row 4:   0 1 0 1

Walk the diagonal: the first bit of row one, the second bit of row two, the third of row three, the fourth of row four. That reads 1, 0, 1, 1. Flip every bit you landed on and you have built 0, 1, 0, 0.

Check it against the list. It differs from row one in the first place, which is where you flipped row one's own bit. It differs from row two in the second place, from row three in the third, from row four in the fourth. Four rows, four disagreements, one per row, each one forced by where the bit was taken from.

Nothing in the walk cared that the list was four long. Run it down a numbered list of strings that never end, one string per whole number, and what you build differs from the seventeenth string in the seventeenth place and from the millionth in the millionth, so it sits on no line of the list. A minute of bit flipping. No table of any size holds every string.

Turing's four lines are that walk with programs on the list and the oracle's prediction as the bit. Gödel's sentence is the identical construction pointed at proof instead of at halting. Number the formulas, so that a statement about numbers is also a statement about formulas, and build the sentence whose content is that it has no proof in the system. If the system proves it, the system proves a falsehood. If the system does not, the sentence says something true that no proof reaches. The lever was never a pun on the word this.

The numbering is the engineering and the self-reference is what falls out of it. Give every symbol a number. A formula is a list of symbols, so code the list by prime powers: symbols numbered three and one give two cubed times three, which is twenty-four, and twenty-four factors into primes one way only, so the code hands back three and one. Proofs are lists of formulas and get their numbers the same way, so the relation between a proof and what it proves is arithmetic, and "no number codes a proof of the formula coded by twenty-four" is a statement about numbers that a system able to do arithmetic can write down.

Getting a sentence to carry its own code with no regress is Quine's device, and works in English with no numbers in it at all. Take a phrase, write it once inside quotation marks, then write it again outside them.

"is twenty words long when appended to its own quotation" is twenty words long when appended to its own quotation.

Count. The quoted phrase is ten words and the phrase after it is the same ten, so the sentence is twenty words long, which is what it says about itself. It does that without the word this and without a second sentence propping it up.

Gödel's sentence is built that way. It gets called a paradox. It is not one. The liar sentence, the one saying that it is false, has no consistent reading: true makes it false, false makes it true, and the machinery seizes. Take Gödel's instead and suppose the system proves it. Then a proof of it exists, so what the sentence says is false, and the system can read that proof and derive the opposite of what it just proved. It contradicts itself, which is chapter two's disease, buying everything and ruling nothing out. So a system that contradicts nothing does not prove that sentence, which is what that sentence says about itself. The sentence is true.

That is what the table of seven was handed: a consistent system with a truth inside it that its own proofs do not reach. Hilbert wanted consistent and complete, and complete was the half that went.

## Two hands and a sheet of paper

In January 1948 M. C. Escher printed a lithograph, 28.2 by 33.2 centimeters, showing a sheet of paper lying on a surface. Two shirt cuffs are drawn on the sheet in flat gray pencil. Out of each cuff a hand rises off the paper in full modeling and shadow, and each hand holds a pencil, and each pencil is drawing the cuff of the other.

Chapter six borrowed the impossible staircase, which this artist put into a lithograph twelve years later. The staircase fails at the closing: every individual step up is unobjectionable and the loop cannot be assembled. *Drawing Hands* fails nowhere. Every line in it is a line an artist could draw. Both hands are drawn and both hands are drawing, and cover any part of it you like and what is left is consistent. What the picture lacks is a ground floor. Neither hand is the real one. There is no corner of the composition to stand in and ask which was drawn first.

Douglas Hofstadter named the arrangement in *Gödel, Escher, Bach* in 1979 and spent seven hundred pages on it. A **strange loop** is a hierarchy that climbs through its levels and arrives back at the level it started from, with no level underneath holding the rest of them up. The staircase is what a strange loop looks like when it fails. The hands are what one looks like when it holds. Everything remaining in this chapter is the second picture with something else in the cuffs.

## The description and the thing described

Chapter eleven put a name on what a network of observers holds in common. The record map takes the full state and returns the part that survives comparison, which is every seam entry and nothing that only one observer could have known. That list has a fiber, the set of consistent states carrying it, and when the fiber holds exactly one member, that member is the observable normal form. Objective reality is that state, up to changes nothing can see.

Every object in that sentence is a description or an operation on one. The record list is a description. The normal form is what the description specifies. And the operation running from the one to the other, take the description and return the world it names, is defined on any description at all, including a description of everything.

So point it at everything. Take the whole structure, records included, as one object, and apply the operation to it. What comes back is what that object's own description specifies. The demand is that the two agree.

**The universe is the observable normal form of its own description.**

One line covers the calculator, the four-line program, the sentence carrying its own code and the demand just made of the whole structure, and the same unknown stands on both sides of it.

$$x = f(x)$$

Here x is whatever is being solved for and f is the operation applied to it. On the calculator, x was a number and f was divide one by it and add one. Applied to everything, x is the whole structure with its records and f is take the description and return the world it specifies. Self-reference is that equation, and closure is that equation with the universe in it.

Chapter eleven's three verdicts apply to that description the way they apply to any other. The fiber could be empty. A description with an empty fiber names no consistent state, which is the reading a witness produces when no arrangement of the world satisfies the testimony. It could hold several, in which case more than one structure answers to the description and the description is of none of them in particular. Closure is the middle verdict, taken at the top: the fiber holds exactly one member, and the member is the structure the description came out of.

That single requirement is severe, in chapter two's sense, and does most of the work of existing. Let a structure's description specify something other than that structure, and the description is of the other one, so the structure it came out of is not described at all, which is a way of failing to be anything rather than a way of being wrong. A structure whose observers could never recover its architecture from their records fails the same demand from the other side, because then the description its observers hold is not a description of it, and the operation returns something else.

Nothing outside supplies the description. Asking what switched the universe on assumes an outside to switch it on from. The demand that picks out this object is the demand that there is no outside. Anybody asking has to supply the outside themselves. Every outside anybody has supplied has needed one of its own.

There is a dividend in that which nobody went looking for. A description obliged to specify the structure it came out of cannot leave a number loose, because a loose number is a second structure answering to the same description, and the fiber has to hold one member. Chapter twenty-eight is where that gets cashed. The grain, the fine-structure constant, the capacity and the cosmological constant are written in terms of each other, the loop admits one solution, and a quantity that every other physics has to measure and enter by hand is computed instead. Constants are what a world with no outside is not free to leave open.

## Three orders

The word loop invites a picture of something traveling backward into its own past. Three different orderings are in play. Every paradox that gets built out of the loop comes from confusing two of them.

The first is the fixed point. It contains no order whatsoever. It is a relation among the parts of one structure, holding or failing to hold. Nothing steps it, nothing iterates it, and there is no count of how many times it has been applied, because an equation is not a process. A loop with no beginning is what an equation looks like when nothing solves it first.

The second is an observer's record time: the chain of committed records that observer holds, in the order it holds them. The three things the word time has meant live inside this one. Order, from which repairs commit before which. Direction, from the irreversibility of a commit. Flow, the rate the observer's own state supplies. Memories, clock readings, plans, and the experience of working something out and then building it are events in this chain and nowhere else.

The third is the descent of repair: the order on authorized updates given by the disagreement count, which falls strictly at every accepted step and cannot fall forever. It is what proves the repairs terminate and what lets two schedules be compared. A repair step is not a second on anybody's clock. Two incomparable schedules do not make two parallel histories. Two repairs at opposite ends of the network, with no seam joining the patches that made them, stand in no order whatsoever: no fact settles which came first, and the question never enters anybody's records because there is no record it could enter.

Set the three side by side and the paradox has nowhere to live. First they work out the architecture, later they build it, is a statement about the second order, holding inside the records of particular observers. It says nothing about the first, because the first has no first. The recovered specification, the built machine and the inhabited structure are three readings of one object, and only an observer's own chain of records puts them in an order.

## A thousand machines and one machine

Run the network on one computer. One process, one memory, one loop that picks a seam, checks it, repairs it if the two sides disagree, and goes round again. It ends where chapter ten said it ends. The endpoint does not depend on which seam it picked first.

Then run the same network across a thousand computers, which is what anybody does in practice with eighty-two thousand patches. Cut the patches into a thousand pieces, one per machine. Each machine holds its own patches plus a copy of the seam data along its border, refreshed by messages that arrive late and out of order. Machines take checkpoints. One dies and restarts from its last checkpoint, replaying updates it had applied before it fell over. A supervisor notices another machine falling behind and re-cuts the entire partition to move work off it.

Anybody who has run a large computation knows what the second version looks like from the operator's chair. The logs differ. The timings differ. Two runs of the identical setup interleave their operations differently. Neither interleaving matches anything the single machine does.

Whether the two are the same universe takes one bookkeeping move to settle: say exactly what counts as the state of the world and what is scaffolding around it. Keep the authoritative record each patch holds. Quotient by the changes nothing can see, which chapter eleven's group names. Throw away the message queues, the retry counters, the labels saying which machine owns which patch, and any wall clock the machines happen to be reading. What is left is a state of the plain single-machine network. Every readout an observer inside can perform is computed from that.

Then classify what a distributed run is permitted to do, and it comes to three things.

A commit, which projects onto a legal repair path of the single-machine network. Something disagreed, a patch fixed it, and the fix corresponds to a repair the one machine could have made.

A stutter, which projects onto nothing at all. Delivering a message, refreshing a border copy, taking a checkpoint, preparing a transaction, aborting one, starting a machine, stopping a machine, replaying an update that has been applied once, and every byte of metadata describing how the work was cut up. The projected state before and after is the same state.

A rollback to an earlier committed state, carrying the name of the committed state it returns to.

That exhausts the list. The rest is an induction with three cases. Chapter ten's result is that the endpoint of the repair process is constant across everything reachable: if legal repairs get you from one state to another, both have the same endpoint. A commit moves along legal repairs, so the endpoint is unchanged. A stutter moves nowhere, so the endpoint is unchanged. A rollback lands on a state whose endpoint was the same one. Start the distributed run and the single-machine run from the same initial state, and after every event in the run, including the last, the endpoint is the same.

So the thousand-machine run settles into the state the one-machine run settles into. Every observer readout that goes through the projection agrees. Partition, worker count, schedule, restart history and repartition metadata appear in no observable.

That is a symmetry, in chapter five's sense and in no weaker one. Name the change: divide the work differently. Name the invariant it leaves alone: everything anybody can measure. Chapter five also supplied the deflation. It applies here at full strength. A symmetry is the claim that part of your description was never physical, and how the work was divided was never physical.

Nothing else in the architecture moves as much. Rotating the carrier moves twelve ports. A gauge change moves a convention chosen separately in every place. This one moves the entire execution, and goes uncounted as a symmetry because it looks like a fact about implementation rather than a fact about the world, which is what the labeling of the twelve ports and the choice of gauge in every place both looked like until somebody counted them.

A rendering loop has a frame number, and every object in the scene shares it, which is what makes it readable: ask any two objects what frame it is and they agree. Ask which step a given patch is on here, and two runs that agree on every measurement anybody can take come back with different numbers. A quantity that differs between runs nobody can tell apart belongs to the bookkeeping. No frame, no tick, no rendering loop, and no place to stand and watch the next frame get drawn. What there is instead is a state picked out by the constraints and arrived at by descent from wherever the system happened to be. The descent has no rate, because a rate would be a step count.

So a description of the universe, running inside the universe, is a distributed presentation of it. Its partition and its schedule reach no observable. It settles where the thing it is running inside settles.

## Five cards

Chapter six's observers have no clock, no leader and no queue. Can a thing built like that do anything besides agree?

Deal five cards face up in a row and number the places one to five. A shuffle is a rule saying where the card in each place goes. Call the plain one C: every card moves one place along, and the card in place five comes round to the front.

Two more shuffles, each of which is C wearing different clothes. Shuffle A sends place one to two, two to four, four to five, five to three, and three back to one. Shuffle B sends one to four, four to three, three to five, five to two, and two back to one. Trace either of them and you visit all five places once and return to where you began, which is exactly what C does. They differ in the route and in nothing else.

Chapter nineteen measured how far two operations are from being swappable by doing them in both orders and subtracting. In a group there is nothing to subtract, so you cancel instead: undo one, undo the other, then do each in turn. If the two get along, every doing cancels its undoing and the deck is as you found it. Try it with these two. Perform, in this order, A backwards, B backwards, A, and B.

Read the row off. Every card has moved one place along and the last has come round to the front. Four shuffles that were supposed to cancel have delivered C.

Five pieces of card on a table are enough to check that identity.

Write a program as a list of instructions, where an instruction reads one input bit and says: if that bit is one, do this shuffle, otherwise do that one. Run the instructions in order on the deck. Say the program computes a function of the input bits when the deck ends up cycled by C exactly on the inputs where the function is true, and untouched on the others.

A single instruction computes a single input bit. Building everything else out of that is the identity above, used once per gate. Given a program for one statement and a program for another, assemble four blocks in the order the identity wants them: the first in A's clothing run backwards, the second in B's clothing run backwards, the first in A's clothing, the second in B's clothing. If either statement is false, its two blocks do nothing and the surviving pair cancels, so the deck is untouched. If both are true, the four blocks are the four shuffles above and the deck comes out cycled, by a route that depends on the clothes and can be steered to any cycle you want by choosing them, which is what lets the step be applied again to its own output. That is the AND gate, and with the complement it is every gate there is. David Barrington published the construction in 1989 in the *Journal of Computer and System Sciences*.

Every function of bits is a tree of gates, so every function of bits is a fixed list of shuffles of five cards. The deck never holds more than five things in a row. The price sits in the length: four blocks per level of the tree, so a formula of depth ten compiles into a word of about a million shuffles, and depth is expensive.

Those five cards are hardware. Chapter fourteen counted the rotations of the twenty-faced carrier, got sixty, and found five families of half turns that those sixty rotations shuffle among themselves. Sixty permutations of five objects is the group that leaves the fifth-degree equation without a solution formula, and the carrier's own group of moves. The deck is the carrier. The shuffles are things the carrier does to itself.

The property making the construction work is the property Galois found. The leftover from those four steps, the thing that came out as C instead of as nothing, is itself a move of the group. Collect every leftover the group produces that way and run the same operation on the collection. In a group with the shape a formula in radicals requires, the leftovers shrink at every round and die out. In this group of sixty they never shrink at all: the leftovers are the whole group again, round after round. The five-card identity is that refusal to shrink, written out on one case. The absence of a formula for the quintic and the presence of universal computation are one property of one group, read from two ends.

Compile a function into a federation of observers instead of into a deck: one patch per gate, each patch holding a few registers, each patch checking the seams it touches and repairing its own registers and nothing else. The input arrives as patches that hold their bits and accept nothing else. Set every other register in the network to anything at all. Let the patches repair in any order, with no clock and no leader appearing anywhere in the argument. The network settles, and the value of the function is sitting in a designated register of the settled state, from every starting state, and any two settled states of that network agree there whatever order the repairs ran in.

The small cases can be checked on every input rather than on a sample. A full adder's sum bit, the three-input circuit that adds two bits and a carry and sits underneath every binary addition anybody has ever performed, comes out right on all eight inputs from every starting state. So do the two bit-mixing functions called choose and majority, which take three inputs each and sit inside the hash function that Bitcoin's ledger is built on. A two-input NAND, the gate that answers false only when both its inputs are true, runs end to end on both engines, the deck and the federation, and they agree on all four inputs.

Patch-local repair, with nothing global anywhere inside it, computes every function of bits.

## The hardware

Set the two results next to each other.

The first says the division of labor is invisible from inside. A description of the world, running inside the world, has a partition, a worker count, a schedule and a restart history, and not one of those reaches any observable. The second says the world's local repair is a general-purpose computer: the same rules that make two neighbors agree about a seam will, wired up, evaluate any function you can write down.

Together they close the loop. Observers inside the structure read their records, recover the architecture that produced them, and build it, in their own record time and at their own expense. What they build emits records of the kind they recovered it from. There is no step at which the built thing and the inhabited thing are two objects, because the only differences between them were partition, schedule and worker count, and those are the coordinates the first result deletes. The hardware the observers reverse engineered is the hardware they are running on.

Which leaves the observers themselves, who have been doing all the work and getting none of the explanation. A loop that computes has pieces. Some of those pieces spend everything they have on continuing to be pieces: holding a boundary, keeping a record, repairing themselves faster than they come apart. Every rule in the architecture is indifferent to which patterns survive. Some of the patterns are not.