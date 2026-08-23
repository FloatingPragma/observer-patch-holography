# 2. Why Is There Anything Rather Than Nothing?

Alfred North Whitehead and Bertrand Russell spent ten years trying to build arithmetic out of nothing, and they meant it literally: no numbers assumed, no quantities taken as given, no facts about counting smuggled in from the outside world. They began with a handful of logical marks and rules for shuffling them, and the plan was to show that everything true about numbers followed from the shuffling alone.

The result, *Principia Mathematica*, ran to three volumes. The proposition that one plus one equals two turns up in volume two, at proposition 110.643, after several hundred pages of preparation. Underneath it the authors added a line of commentary, which reads in full: "The above proposition is occasionally useful."

Russell wrote later that his intellect never quite recovered from the effort, and that he had been definitely less capable of difficult abstraction ever since. He was thirty-eight when they started. Cambridge University Press estimated the work would lose money and asked the authors to contribute fifty pounds each toward publishing it, which they did, and the book brought them no income.

It is easy to read that as an anecdote about eccentric mathematicians and miss what the project was. Whitehead and Russell were not being pedantic about two. They were testing whether a body of knowledge could stand on itself: whether you could put down marks and rules with no borrowed meaning, turn the handle, and have arithmetic come out. If it worked, then nothing about numbers had to be assumed, because everything about numbers would be a consequence. Ten years, three volumes, and fifty pounds each was what it cost to find out for arithmetic alone.

What they were building has a name. A **formal system** is a set of marks you are allowed to write, a starting position, and a set of moves that turn one arrangement of marks into another. Meaning is absent by construction. The marks stand for nothing, the moves have no reasons, and the system has no idea what it is about. Hand the rules to somebody who shares no language with you and they will apply them correctly.

Which makes it the only kind of object anybody can examine without assuming something first. Chapter one ended by throwing away everything: no space, no fields, no symmetry group, no numbers, no place for events to be. The question is whether anything at all survives that, and the question has a famous form. Leibniz put it in 1714, in a short piece written for a prince: why is there something rather than nothing? Nothing, he pointed out, is simpler and easier.

## Two answers that do not work

The traditional replies come in two shapes.

The first says: something outside made it. Something caused the world, and that something is the reason there is anything at all. This is the oldest answer on record and a serious one. Aristotle got to it by asking what moves a thing, and then what moves that, and noticing that the question does not stop on its own.

The trouble it runs into is arithmetic. Whatever you nominate as the cause, you can ask the same question about it, and you have two ways out, both bad. Either the chain runs backward forever, in which case you never arrive at a reason and have only postponed the problem an infinite number of times, or you stop at some point and declare that this one thing needs no cause. But if one thing can need no cause, the question is why that particular thing gets the exemption, and no answer to that is available except that we stopped there. The item nominated as needing no cause has an almost perfect historical record of being the item the arguer arrived with.

A first cause is also an assumed component, put on the board because the argument needs it rather than because anybody found it, which is the one rule a reverse engineer cannot break.

There is a second objection, and it is the one that does the damage here, because it does not depend on the regress at all.

Something outside. Something before.

Outside is a spatial word. It means separated by a distance, on the far side of a boundary, in a different region. Before is a temporal word. It means earlier along an ordering that both things sit in. Chapter one did not suggest that space and time might turn out to be built from something else. It took them off the starting inventory altogether, along with the arena and the parameter that things happen along. So there is no distance for a cause to sit at and no earlier for it to have acted in.

Which makes "something outside made it" a sentence assembled out of the two components a reverse engineer has not found on the board. The regress was never the real complaint. A regress is what turns up when a question has quietly helped itself to a place to stand and a moment to stand in, and then goes looking for somewhere to put the answer.

The second reply says it is just here. A brute fact, with nothing behind it, and asking why is a mistake. This has the virtue of stopping the regress immediately and it has one crippling defect. It explains nothing and it also forbids nothing. If the world is here for no reason, then a world with four spatial dimensions would also be here for no reason, and so would a world with none, and so would a world where the electron weighs twice what it does. You have not explained the world, you have declared that no explanation exists, and in the same motion you have thrown away any hope of understanding why it has the specific shape it has rather than one of the others. The twenty-six measured numbers of chapter one stop being a debt and become the natural state of affairs.

There is a way to feel the size of what that concedes. The mass of the electron is known to about eleven significant figures. Treat it as a brute fact and you have accepted that eleven digits arrived for no reason, and that the same is true of the other twenty-five numbers, and that the total quantity of unexplained precision in the world runs to a couple of hundred digits that simply are what they are. That is a great deal to put down to nothing in particular.

Both replies point away from the thing they are explaining, one to a cause and one to a shrug, and pointing away is the move that has just been withdrawn. With no outside and no before, there is nowhere for an account to come from except the thing itself. So what is left is not a third option chosen over two others. It is the only one still standing, and it stands by elimination: a structure accounts for itself, holding together in a way that requires nothing to be added to it, so that asking what put it there has the same status as asking what makes a triangle have three sides.

That also does something to the question in the title. "Why is there anything" wants a because, and a because wants somewhere to come from. Take away the outside and the before and the why has no room to operate in, while the question underneath it survives intact and turns out to be a question about structure: what does a thing have to be like, to be the kind of thing that can hold itself up? What follows is the beginning of that answer, and the beginning is small.

Stated like that it is a phrase rather than a proposal, and phrases of that shape have a long history of sounding like wisdom and meaning nothing. It gets edges from three things: what it means for a system to have consequences nobody chose, what it costs to contradict yourself, and what happens when you demand that something survive its own description.

## Two rules and a pencil

Here is a formal system small enough to run by hand. You will need about ninety seconds.

The marks are the letters a and b. The starting position is the string `ab`. There are two moves:

**Move one.** Take any string you have and wrap it: put an a on the front and a b on the back. So `ab` becomes `aabb`.

**Move two.** Take any string you have and write it twice. So `ab` becomes `abab`.

That is the whole system. Start from `ab` and apply the moves in any order you like, as many times as you like. From `ab` you can reach `aabb` and `abab`. From `aabb` you can reach `aaabbb` and `aabbaabb`. From `abab` you can reach `aababb` and `abababab`. Go on for a while and you will generate a large and reasonably complicated family of strings, and you can generate as many as you have patience for.

Here is a question about that family. Can you reach `aab`?

Try it for a minute and you will fail. By itself, failing to find something is not the same as showing it is not there. But look at the counts. The starting string has one a and one b. Move one adds one a and one b. Move two doubles both. So every string you can possibly reach has exactly as many a's as b's, in every case, forever, without anybody having written that down as a rule. `aab` has two a's and one b, so it is unreachable, and you know this without searching. You have proved something about an infinite collection of strings by looking at four lines of rules. The proof was sitting there before you began searching, which is the standard experience in this subject and does not become less irritating with practice.

The decisive quantity is the difference between the number of a's and the number of b's. It starts at zero, and neither move can change it: one move adds one to each side of the subtraction, the other doubles both. A quantity the rules cannot touch, however long you run them, is the most useful thing anybody can know about a system.

There is a second thing about that little system that's easy to overlook. The strings it makes are not about anything. `aabb` is not a description of some other, realer `aabb` sitting elsewhere in a warehouse of genuine strings. Ask whether the one in front of you is the authentic article or a copy and the question dies on the way out, because there is nothing available for it to be a copy of. You ran the rules, and what came out is what there is.

That is the first useful property of a system with no meaning in it. Its consequences are not limited to what its author had in mind. Whitehead and Russell wrote down their marks and their moves, and then spent a decade finding out what they had committed themselves to, and the answer took three volumes because the answer was not up to them. Nothing was consulted. Nothing was measured. The counts were fixed the moment the rules were written, and the ten years were spent reading them off.

## The price of contradictions

Suppose a bank's ledger says that an account holds nothing and also says that the same account holds a hundred pounds. Both entries, on the same page, at the same time.

A clerk looking at that does something sensible. They decide one entry is a mistake, or they pay out the hundred and open an investigation, or they go and find the deposit slip. What a clerk does not do is conclude that the account holds nine million pounds. People reason around contradictions all day without much difficulty: two entries that cannot both hold get quarantined, and the rest of the ledger goes on working.

A formal system has no clerk. It has marks and moves, it applies every move that applies, and it has no faculty for deciding which of its own lines is the suspicious one. Applied honestly, the moves take it somewhere no clerk would follow.

Watch. The ledger says the account holds a hundred pounds. From that alone I may write down a weaker sentence: either the account holds a hundred pounds, or the moon is a biscuit. That step costs nothing, because attaching an alternative to a sentence you already have cannot make it less true. Anybody who accepts the hundred pounds is obliged to accept it.

But the ledger also says the account holds nothing, so the first half of my weaker sentence is out. An either-or with its first half ruled out leaves the second half standing.

The moon is a biscuit.

No step in that is a trick, and a clerk would sign off on each one taken by itself. Run it again with a different second half and the account holds nine million pounds. Run it again and you get whatever you please, because the second half was never constrained by anything: it was chosen freely at the first step and it survived to the last. Two entries on one page, and the ledger certifies every sentence in the language.

An inconsistent description asserts every world at once, and a description that permits everything distinguishes nothing. It is not a false theory. It is not a theory.

Somebody found this out the expensive way, in June 1902, and the letter that did it came from Russell. Gottlob Frege had spent two decades building arithmetic out of logic, and the second volume of his *Grundgesetze der Arithmetik* was at the printer when Russell wrote to him on the sixteenth with one short observation about the class of all classes that are not members of themselves. Frege replied six days later. He had seen it at once, and what he had seen was not that some theorem of his was wrong. It was that his system now proved everything, which is the same as proving nothing. He added an appendix while the volume was in press, and it opens by saying that hardly anything more unwelcome can befall a scientific writer than to have one of the foundations of his edifice shaken after the work is finished.

Which is why *Principia Mathematica* runs to three volumes. The several hundred pages of preparation, and one plus one arriving at proposition 110.643, are not eccentricity. They are the price of a system that Russell's letter could not be written to.

Chapter one gave you the tool for seeing why that matters here. A theory with infinitely many free parameters fits anything and forbids nothing, which was the fatal problem with quantizing gravity by brute force. A contradiction is the same disease in its terminal form: not many knobs but all of them, not a poor prediction but every prediction at once. Consistency is what buys a description the right to rule anything out.

So the demand this book starts from can be stated. **A structure counts as existing exactly insofar as it accounts for itself with nothing left over.** No external cause is allowed, because a cause is one more thing standing in need of an account. No leftover assumptions are allowed, because an assumption is a thing you would have to explain. And it must not contradict itself, because a contradiction buys everything and therefore purchases nothing.

The obvious objection is that this demand is so severe that nothing could satisfy it, and the answer to that objection is a number.

## One demand, and a calculator

Take a calculator, or a phone, or a piece of paper if you are feeling nineteenth century. Pick any positive number at all. Any at all: your age, the year, seven.

Then do this. Divide one by your number, add one, and keep the answer. Then do it again to the answer. Then again.

Start with seven and you get 1.142857, then 1.875, then 1.533333, then 1.652174, then 1.605263, then 1.622951, then 1.616162, then 1.618750, then 1.617761. The numbers are hunting, overshooting one way and then the other, and closing in. Keep going and they settle on 1.6180339887.

Start with a thousand instead. Or with 0.001. Or with the number of stairs in your building. Same destination, every time, to as many decimal places as you have patience for.

There is a small bonus in the arithmetic if you look at the fractions rather than the
decimals. Starting from seven you get 8/7, then 15/8, then 23/15, then 38/23, then
61/38, then 99/61, then 160/99. Read the numbers down the page: 7, 8, 15, 23, 38, 61, 99,
160. Each one is the sum of the two before it. You did not ask for that and the rule you
were applying says nothing about addition.

What is 1.6180339887, and why does everything land on it?

It is the number that the operation leaves alone. If you feed in a number and get the same number back, then that number satisfies x equals one plus one over x, and there is only one positive number that does. Multiply through by x and it reads x squared equals x plus one, which is a quadratic, and its positive root is one plus the square root of five, all over two. Which is 1.6180339887, and which people have called the golden ratio since long before anybody could write it down that way.

The hunting has a mechanism. When your number is too big, one over it is too small, so the answer comes out below the target. When your number is too small, one over it is large, so the answer overshoots above. Every step lands on the far side of the destination from where you started, and lands closer. The error does not merely shrink, it changes sign each time and shrinks, which is why the decimals settle from both directions at once.

A number that a process leaves alone is called a **fixed point**, and this is the first of many in this book. The number itself turns up in enough seaside-gift-shop literature to have earned a bad reputation. Where it came from is the part that matters. Nothing was measured. No experiment was performed and no data existed at any point in the procedure. A single demand, that the operation return what it is given, reached into the continuum of real numbers and picked out exactly one.

That is the shape of the answer to the question in the title. A requirement severe enough to have exactly one solution.

That quadratic has two roots, and the other one is minus 0.618034. Check it: one divided by minus 0.618034 is minus 1.618034, and adding one gives you minus 0.618034 back. It is a fixed point, exactly as much as the other one is. But start your calculator anywhere near it and the numbers do not settle there, they flee. Being a point that a process leaves alone is one property. Being a point that a process moves toward is a different property, and a structure can perfectly well have fixed points that nothing ever arrives at.

## What one requirement is worth

Look at what the demand bought. Before it there was a continuum of positive numbers available, uncountably many, with no reason to prefer any one of them. After it, there was one. One number, to as many decimal places as anybody cares to compute.

That exchange rate between what goes in and what comes out is the engine of this book, and it fails in both directions.

Require nothing and you get nothing. This is the part that sounds wrong, because requiring nothing sounds like permitting everything, and permitting everything sounds generous. It is empty. A structure with no requirements on it does not have many possible forms; it has no form, because a form is exactly a set of things that are not allowed. If every arrangement is as good as every other, then nothing distinguishes anything, and there is no fact about what the structure is.

Which is the same failure the bank ledger had, arrived at through the opposite door. A contradiction permits every conclusion and therefore asserts nothing. A structure with no constraints permits every configuration and therefore is nothing. Total inconsistency and total freedom are emptiness, and a great deal of confused metaphysics comes from noticing only one of them.

Between those two failures is a window, and the window is far narrower than anybody expects. Put a small number of requirements on a structure and what survives is not a comfortable family of possible worlds to choose among. It is very nearly one world, with a definite number of directions in it, a definite list of things it is allowed to contain, and definite numbers attached to them, none of which anybody selected. The golden ratio is that phenomenon at the smallest scale it can occur at: one requirement, one survivor, out of a continuum.

This cuts both ways and the second way is the uncomfortable one. If adding a requirement collapses a continuum to a point, then a requirement that slips in without being noticed produces a world that arrives without being earned. Chapter one is three centuries of exactly that going wrong. So the constraints have to be counted, they have to be few, and every one of them has to be visible.

## Three requirements

Everything here comes out of three requirements, taken one at a time over the next three chapters.

The first says that whatever the world is made of, each piece of it sees only a part and keeps a record of what it saw. The second says that where two pieces overlap, their records have to agree, and that any relabeling of a piece has to be a move the system can carry out on itself, with no hand reaching in. The third says that beyond what those two force, nothing is assumed: whatever is not pinned down stays as unconstrained as it is possible to be.

There is very little in there. Nothing about space, because there is none yet. Nothing about time, particles, forces or numbers. No mention of what the pieces are made of, how many there are, or what the records say. The first two are constraints on bookkeeping and the third is a refusal to add anything.

Three, and not more, and if that sounds too few to build a universe from, the history of the last comparable attempt is instructive.

Euclid began his geometry with five postulates. Four of them are the sort of thing nobody argues with: you can draw a line between two points, you can extend a line, you can draw a circle with any center and radius, all right angles are equal. The fifth is about parallel lines and it is longer and clumsier than the other four put together. Euclid appears to have been embarrassed by it, because he avoids using it for as long as he possibly can. For two thousand years mathematicians tried to prove it from the other four and get rid of it, and failed, and assumed they had merely failed rather than attempted something impossible.

In 1733 an Italian Jesuit named Giovanni Saccheri published a book with the confident title *Euclid Freed of Every Flaw*, in which he tried the last remaining tactic. Assume the fifth postulate is false, grind out the consequences, and wait for a contradiction. He ground out a great many consequences. Triangles whose angles sum to less than two right angles. Lines that approach each other forever without meeting. An entire coherent body of theorems, dozens of them, all following perfectly from four postulates and the denial of the fifth, and not one contradiction anywhere in it.

What Saccheri had in his hands was a complete and consistent geometry of curved space, a century before anybody else found one. What he concluded was that the results were repugnant to the nature of a straight line, which he offered as the contradiction he had spent the book hunting for. He published, and died that year, having discovered a new kind of space and declined it on the grounds that he did not care for the look of it.

A hundred and eighty years later Einstein went looking for a mathematics of gravity and found he needed geometry without the fifth postulate: space whose triangles do not add up to two right angles, whose parallel lines do not stay parallel, and whose shape is settled by what is sitting in it. The postulate Saccheri was trying to rescue is false in the universe he wrote the book in.

Two things come out of that. Four postulates and a decision about a fifth had been quietly fixing the shape of space for two thousand years, and nobody found out by inspection, because you cannot see a postulate by staring at it. They found out by working: by taking the requirements seriously and grinding out what they force, which took twenty centuries and was being got wrong as late as 1733. And the consequences of a small set of requirements are not up for negotiation. Saccheri disliked his and they were true anyway.

That is the whole starting position, and there is a hole in it.

Look at what the three tools above have in common.

A formal system produces consequences, and those consequences are strings: arrangements of marks that stand one way rather than another and can be read back. Consistency is a property only certain objects can fail to have. A rock cannot contradict itself. It sits there being whatever shape it is, and no arrangement of atoms is a contradiction, because contradicting yourself means holding two claims about the same thing and a rock holds no claims. And a fixed point is whatever satisfies a condition, where a condition is a demand, and a demand has to be made of something and stated somewhere before it can do any demanding.

Every one of those is a property of a description rather than of a lump of world.

Which puts a clause into the founding demand that nobody wrote there. "Accounts for itself" has no application to a structure containing no accounts. Ask whether a universe with nothing in it that keeps a record is consistent, and the question slides off: there is nothing in such a place to be right or wrong about anything, so there is nothing for consistency to be a property of. The word does not become false. It stops reaching.

So the founding demand turns out to require something nobody put into it and nobody would have wanted. Somewhere in that structure there has to be a thing that holds a description, sets it against something, and is capable of getting it wrong. Otherwise the demand has nothing to bite on, and a structure the demand cannot bite on fails it, and by its own terms it does not exist.

A calculator, two letters of the alphabet, and a bank ledger with two entries on it. What comes out at the bottom of them is that the universe has to contain something capable of being wrong.