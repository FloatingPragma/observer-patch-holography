# 6. Why Isn't Agreeing With Your Neighbors Enough?

In 1785 the Marquis de Condorcet published an essay on applying analysis to decisions reached by majority vote, and buried in it is an observation that has been ruining committees ever since.

Take three people and three options. The first person ranks them A, B, C. The second ranks them B, C, A. The third ranks them C, A, B. Ask the committee to choose between A and B, and A wins by two votes to one. Ask it to choose between B and C, and B wins two to one. Ask it to choose between C and A, expecting a formality, and C wins two to one.

Every one of those votes is a real majority. Nobody was manipulated, nobody abstained, and each vote can be repeated as often as you like with the same result. There is no winner, and there is no winner in a way that no amount of further voting will fix, because the information that the three results are jointly impossible is not present in any of them.

Condorcet's committee is the first warning that pairwise agreement is a weaker thing than it sounds. Local data can be flawless, checkable and reproducible, and refuse to be assembled into anything. The voting version has no arithmetic attached to it that helps much. The version that does is simpler, and it can be checked on a table.

Give Alice, Bob and Charlie a coin each and put them in three separate rooms. Each of them looks at their own coin and sees heads or tails. Nobody can see anybody else's.

They are allowed to meet in the corridor two at a time. When a pair meets they compare coins and write down one word: **same**, if the two coins match, or **different**, if they do not. The pair does not write down which way either coin is facing, and this is the only restriction in the whole arrangement. What gets recorded is the relation between two coins, which is the only thing either of them could check together.

Three meetings happen. Alice and Bob meet, and record different. Bob and Charlie meet, and record different. Charlie and Alice meet, and record different.

Set the coins.

Suppose Alice's is heads. Bob's has to be tails, because they recorded different. Charlie's has to be heads, because Bob and Charlie recorded different. And Charlie and Alice recorded different, so Charlie's has to be tails. It was heads two sentences ago.

Try the other branch. Suppose Alice's is tails. Bob's is heads, Charlie's is tails, and Charlie and Alice are both tails, having recorded that they differ.

That is all the branches there are. There are eight ways to set three coins, and every one of them contradicts at least one of the three records.

Nobody lied. No coin was misread. Each of the three meetings produced a record that is perfectly satisfiable on its own, and any two of the three records can be satisfied together without the slightest difficulty. It is the third one that has nowhere to go, and which of the three you call the third one is your choice, because the three are symmetric and none of them is the culprit.

Notice what the three of them can do about it, which is nothing. They can meet in pairs for the rest of their lives. Every meeting will confirm what the previous meeting found. No pairwise comparison anywhere in the arrangement contains the information that the three records are jointly impossible, because that information is not held by any pair. The three of them can go on being individually correct indefinitely.

Chapter five put a third person in the room and left the three of them agreeing. This is the bill.

## Same and different, as arithmetic

The argument above is short enough to check by hand, and it will not survive being scaled up by hand. Twelve people generate sixty-six possible pairings, and staring at them is not a method. What is needed is a way of computing the answer, and it turns out to need the smallest arithmetic there is.

Write **0** for same and **1** for different.

Then chain two records together. If Alice and Bob are the same, and Bob and Charlie are the same, Alice and Charlie are the same: 0 and 0 give 0. If Alice and Bob are the same and Bob and Charlie differ, then Alice and Charlie differ: 0 and 1 give 1. And if Alice and Bob differ and Bob and Charlie differ, then Alice and Charlie are back to being the same, because two flips return you to where you started: 1 and 1 give 0.

That last line is the whole of it. Ordinary addition says one plus one is two. Here there is no two, because there are only two states in the world being described, and the count that matters is whether the number of flips is even or odd. So one plus one is zero.

This is called **arithmetic modulo 2**. It is the arithmetic of light switches: flick the switch twice and the room is as you found it. It is the arithmetic of turning a sock inside out. Every fact in it fits on two lines, and those two lines are the whole of what the coins require.

Take the three records again. Different, different, different: 1, 1, 1. Add them. One and one is zero, and zero and one is one.

The three records around the loop sum to 1.

For the coins to be settable at all, that sum has to be 0, because going all the way around the loop returns you to Alice's own coin, and Alice's coin does not differ from itself. An odd number of flips around a closed path is a demand that something be different from itself, and the coins are not able to arrange that.

Eight cases checked by hand, replaced by adding three ones.

## Loops

The picture generalizes as soon as it is drawn.

Put a dot for each observer and a line between two dots whenever those two have compared records. A collection of dots and lines like that is a **graph**, and it is one of the few pieces of mathematics that looks exactly like the thing it is about. The dots are called vertices, the lines edges, and two dots joined by a line are neighbors.

Walk along the edges from dot to dot and you have a **path**. Come back to the dot you started from without retracing your steps and you have a **cycle**. The triangle of Alice, Bob and Charlie is the smallest cycle there is.

Every edge carries a record, 0 or 1. Every cycle therefore carries a total, which is the sum of the records around it, computed in the arithmetic that has no two in it. The quantity you pick up by going all the way around a loop and arriving back where you started is called the **holonomy** of that loop.

The condition on the whole arrangement fits in one sentence.

**A consistent global assignment exists exactly when every cycle has holonomy zero.**

One direction of that is the coin argument again: if the coins can be set, then walking around any loop and coming home must find the coin you left, so the flips have to cancel.

The other direction is constructive, and it tells you what to do rather than what to hope.

Suppose every loop sums to zero. Pick any dot and set its coin however you like. Walk outward. Every time you cross an edge you know whether to keep the value or flip it, because the edge says same or different, so every dot you can reach gets a value. The only way this can go wrong is if two different routes to the same dot disagree, and two routes to the same dot form a loop, and every loop sums to zero. So they cannot disagree. The assignment is forced, everywhere, from one free choice at the start.

A route through the dots that reaches all of them without ever closing a loop is called a **spanning tree**, and it is the skeleton of that construction. For any connected graph, a spanning tree has exactly one edge fewer than it has dots. Every remaining edge, the ones the tree did not use, closes exactly one independent loop when you put it back.

Which gives you the count of conditions, from nothing but arithmetic. Take the number of edges, subtract the number of dots, add one. That is how many independent loops the graph has, and therefore how many separate sums have to come out zero. The triangle has three edges and three dots: three minus three plus one is one. One loop, one condition, and it failed.

## Ascending and descending

Roger Penrose went to the International Congress of Mathematicians in Amsterdam in 1954, where there was an exhibition of prints by M. C. Escher, and came away preoccupied. He and his father Lionel, a psychiatrist, spent some time afterward constructing figures that behave locally like ordinary drawings and cannot be assembled into an object, and published them in the *British Journal of Psychology* in 1958. One of them was a staircase.

Every flight goes up. Take any corner of the drawing, cover the rest, and there is nothing to object to: a step, then another step, each one higher than the last. Follow all four flights around and you arrive back at the landing you started from, having climbed the whole way.

The Penroses sent a copy to Escher, who put the staircase into a lithograph in March 1960 and populated it with a file of monks trudging around the roof of a monastery forever.

The staircase is the coin problem drawn. Each edge of the loop carries a record, in this case "up", every local record is fine, and the sum around the loop is not zero. What is impossible is not any of the steps. It is the closing.

## A magnet that cannot decide

Physics met this in the 1950s and gave it a name.

A magnet is many small magnetic moments, each of which can point one of two ways, and which push their neighbors around. In an antiferromagnet, every pair wants to be opposite: if this one points up, the one beside it should point down. On a square grid that is easy to arrange, and the result is the tidy alternating pattern of a chessboard.

Put the same rule on a triangular grid and it comes apart at the first triangle. Set the first moment up and the second down, which satisfies that pair. The third has to be opposite to both, and there is no third direction. Whatever it does, one of its two bonds is unhappy, and the choice of which bond to disappoint is arbitrary.

Gregory Wannier worked out the consequences for a whole triangular lattice in 1950 and found something that had no business being there: the system keeps a finite amount of disorder as the temperature goes to absolute zero. Ordinarily, cooling something to zero freezes it into its one best arrangement. Here there is no one best arrangement, only an enormous number of equally disappointing ones, and the system goes on having a choice after the energy has been taken away.

Gerard Toulouse gave the phenomenon its name in 1977, while working on spin glasses, and gave it the same test the coins gave it. Go around a closed loop of bonds and multiply out what each one demands. If the loop comes back consistent, the region can settle. If it does not, the loop is **frustrated**, and no arrangement of the moments inside it satisfies every bond, however long you wait.

Frustration is the coin problem in a material, and it is measurable. It changes the heat capacity, it shows up in neutron scattering, and it is the reason a class of magnets never orders properly no matter how carefully they are cooled. The obstruction lives in the loops, exactly where the arithmetic said it would.

## Twelve patches

The triangle is small enough to be dismissed as a puzzle. So take the smallest thing that deserves to be called a universe and check every state of it.

Twelve observers. Wire them as the corners of an icosahedron, the twenty-sided solid, so that each observer has exactly five neighbors and there are thirty overlaps in all. Twelve corners, thirty edges. The count above says the graph has thirty minus twelve plus one, which is nineteen independent loops, so nineteen conditions have to hold at once for the thing to be settable.

Each observer holds one bit. That is four thousand and ninety-six arrangements, and flipping all twelve at once changes no comparison anywhere, so there are two thousand and forty-eight distinct states.

Two thousand and forty-eight is a number you can exhaust. So exhaust it. Start the system in every one of those states and let the observers repair: find a neighbor you disagree with, and fix it. Do that in every order. The order matters enormously in principle, because these observers have no clock, no leader, and no queue, so the sixteen different schedules used here are sixteen different worlds as far as anyone inside can tell. Across the whole sweep the machine executes eleven thousand two hundred and sixty-four repairs.

Two thousand and forty-eight states, sixteen schedules apiece, and eleven thousand two hundred and sixty-four repairs is not a sample. It is the whole space. Nothing here rests on a well-chosen starting configuration or a lucky order, because there was no choosing and no luck: every configuration this universe admits was started, and every order these observers could repair in was run.

Two results come out, and the first one is the boring one that has to be true before the second one means anything.

Exactly one of the two thousand and forty-eight states is globally consistent, and every single run lands on it. Every starting point, every schedule, every order of repair. The answer does not depend on who went first. Nineteen loops, nineteen holonomies, all zero. The leftover at the end is zero, which is to say the observers end up in complete agreement and the agreement was never a matter of who spoke loudest.

## One edge

Flip a single edge. Change one record, out of thirty, from same to different. Change nothing else. Every observer has five neighbors, every overlap works, and the machine is untouched.

Two of the nineteen loops sum to 1.

Of the two thousand and forty-eight states, the number that are globally consistent is zero. There is no way to set twelve bits that satisfies all thirty records, and no amount of repair produces one, because the obstruction is a property of the loops rather than of the bits.

What happens instead is the interesting part, and it is not what a reader expects at this point. The machine does not thrash, and it does not fail to halt, and it does not give different answers on different days. It converges, exactly as before, from every one of the two thousand and forty-eight starting states and under every one of the sixteen schedules, onto one terminal arrangement.

That arrangement has a leftover of exactly one.

Every repair the observers can perform has been performed. Every disagreement that could be pushed somewhere else has been pushed somewhere else. And when the pushing stops, one unit of disagreement sits in the system. It cannot be removed by working harder, it cannot be removed by choosing a better order, and it cannot be pushed off the edge of the world, because a system of twelve observers wired into a closed surface has no edge to push it off.

You can move it. Repair in a different order and it sits somewhere else. What you cannot do is get rid of it, and the total is one wherever it happens to be sitting.

## The leftover

Look at what that leftover does.

It has a definite quantity, and the quantity is exactly one, not approximately one and not one on average. It has a location, in the sense that the two loops carrying holonomy tell you where it is. It survives every process the world can apply to it. It is not a property of any single observer, because every observer is behaving correctly and repairing whatever it finds. And it is not a property of any pair, because every pair is satisfiable and always was.

It is a property of the way the whole thing is wired, and it is as real as anything else in the arrangement, in the only sense of real that has been available since chapter four: it is what survives comparison between every point of view.

A conserved quantity, indivisible, located, indifferent to how you approach it, that no amount of local work will destroy.

The word for that costs twenty more chapters to earn, so it stays unsaid.

## The order of the meetings

Alice, Bob and Charlie met in a corridor, two at a time, and the account above knows which meetings happened before which. The twelve observers ran sixteen schedules, and calling them schedules is a claim that there is an order to run things in. Every one of those sentences assumed that comparisons are events, that events happen one after another, and that there is a fact about which.

Take that away.

Three observers with no clock, no shared present, and no way of telling a neighbor who has gone quiet from one that has stopped existing.

They have to agree anyway.