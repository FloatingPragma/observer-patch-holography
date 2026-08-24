# 11. What Does Everybody Actually End Up Holding?

In 1974 Elizabeth Loftus and John Palmer sat forty-five students down at the University of Washington and showed them seven short films of traffic accidents. Afterward each student was asked how fast the cars had been going. It was the same question every time except for one word. Some were asked how fast the cars were going when they smashed into each other, and the others got collided, bumped, hit, or contacted.

Everybody saw the same films. The smashed group put the cars at 40.8 miles an hour on average, the hit group at 34.0, and the contacted group at 31.8.

A second experiment used a hundred and fifty students and one film of a multiple-car crash. Fifty of them were asked about smashing, fifty about hitting, and fifty were asked nothing about speed at all. A week later, all of them were asked whether they had seen any broken glass. Sixteen of the fifty in the smashed group said they had, against seven in the hit group and six in the group that had never been asked about speed. There is no broken glass in the film.

Loftus and Palmer called the paper "Reconstruction of Automobile Destruction". Nobody in those rooms was lying, and nobody was playing back a stored recording. Each of them was assembling an account out of the material available, and the material available included the verb.

Courts worked this out considerably earlier and built an institution around it. A trial does not ask what is inside a witness's head, partly because there is no instrument for getting at it and partly because the contents turn out to depend on how you ask. What a trial produces instead is a record: the words actually spoken in the room, written down as they are spoken by somebody employed to write down words. Exhibits go into it. Private recollection does not, and neither does anything a juror happens to know from outside. Continental procedure has a maxim for the arrangement. Quod non est in actis, non est in mundo. What is not in the record is not in the world.

The record is a physical object made by a person under load. The Registered Professional Reporter certification, from the National Court Reporters Association, requires taking down two speakers at 225 words a minute with 95 percent accuracy, on a keyboard with 22 keys, which is fewer keys than there are letters and is the reason the output is phonetic rather than alphabetical. Everything a court later does with the case runs through what came off those 22 keys.

No authority outside the room supplies the verdict. Twelve people who sat through the trial compute it from the record. When they cannot agree about what a witness said they ask for the transcript and it is read back to them.

Chapter ten put sixty-five thousand patches through the repair process and the last of their 326,047 disagreements went at cycle 96, whatever order the repairs ran in. That settles the schedule. It says nothing about what any patch is holding at cycle 97. Each finished the run holding its own state, visible to nobody else, connected to the rest through what crosses its seams. So what do they end up holding in common, and how much of the world does it pin down?

## The list everybody can check

When Alice and Bob met in the corridor with their coins, they wrote down one word, same or different, and never which way either coin was facing. That was the only restriction in the whole arrangement. Two observers cannot compare what only one of them holds; they can compare the relation between them, and the relation is what gets written down.

Collect every one of those written comparisons into a single list. That list is what anybody can check, in the sense that any observer can walk up to the seam it names and find out whether the entry is right. Everything else in the network is somebody's private state.

The **record map** takes the full state of the network and returns the part of it that survives comparison: every seam entry, and nothing that only one observer could have known. Two observers who disagree about anything are disagreeing about an entry, because an entry is the whole of what either of them can put on the table.

What the record map deletes is real. Alice's coin is a feature of Alice, she can look at it whenever she likes, and it does not appear in the list because there is no procedure by which she and anybody else could jointly establish which way it faces. She has been deleting it at every meeting since the corridor, and so has everybody she has ever compared with.

A state of the network is **consistent** when every entry on that list is true of the patches it names. Twelve coins and thirty written records: consistent means all thirty records say what the pairs of coins they name actually do. The collection of consistent states is the **consistent set**.

Chapter six counted that set for the icosahedral net. Thirty records that pass the loop test, twelve coins, four thousand and ninety-six ways to set them, and exactly one setting satisfies every record, along with the setting you get by turning all twelve coins over at once.

## Three verdicts

Fix a list. Ask which consistent states carry it.

The answer is a set: all the ways the world could be that would produce exactly this record and no other. In court it is the set of stories the testimony fits. The prosecution has one of them in mind, the defense usually has another, and the size of the set is what the argument is about. The set of consistent states carrying a given reading is called the **fiber** of that reading.

Wire three patches to each other and to nobody else. Three coins, three seams carrying one entry each, and eight ways to set the coins.

Read all three entries as same. Heads all round satisfies them and tails all round satisfies them, and those two differ by turning over every coin in the world, which no entry sees, so they count once. One world.

Read all three entries as different instead. Chapter six wrote different as 1, added the three of them in the arithmetic where one and one make zero, got one rather than zero, and no setting of three coins survives a loop that fails to close.

Let the third patch go quiet. Its two seams go unread, the record shrinks to the single entry between the first two, and that entry reads same. Count this one yourself. Four of the eight settings have the first two coins matching: heads heads heads, heads heads tails, tails tails heads, tails tails tails. All four satisfy the record, because the record has nothing to say about the third coin. Pair each setting with its all-over flip, the way the first reading paired heads all round with tails all round, and the four collapse to two: one world in which the third coin agrees with the other two, one in which it differs. The record names both and picks neither.

Three readings, one question put to each, and the fibers came back holding one world, none, and two. One recipe produced all three counts. It takes three letters to write down. Write r for a reading, s for a state of the network, and R for the record map that strips a state down to its entries. The fiber of r is then every consistent s whose R(s) comes out r: keep the states that satisfy the entries, discard the ones whose entries read differently, count what is left over. Fix a map, name one of its outputs, collect the inputs that land on it, and that operation is called taking a preimage. A fiber is a preimage.

A fiber can have three sizes, counting a state and its all-over flip as one. The three sizes are three different situations rather than three points on a scale.

It can be empty. No consistent state produces this record, so the testimony describes nothing that could have happened, and no amount of further questioning will produce a story that fits. The reading is **unrealizable**.

It can hold exactly one. The record pins down the world completely. Any two observers who hold the record hold the same world whether or not they ever compare anything else. The reading is **reconstructing**.

It can hold several. Every state in it satisfies every record, and nothing anybody wrote distinguishes them. The reading is **ambiguous**. This is the fiber a jury is looking at when it acquits somebody it privately suspects, and the acquittal is the correct output: an ambiguous fiber has more than one member and the record does not name one.

Advocates argue about which of the three they are in. The defense wants the fiber shown to hold more than one member, and needs only one alternative story that fits every piece of testimony. The prosecution wants it shown to hold exactly one. A witness caught in an unrealizable reading has produced an account that no arrangement of the world satisfies, which is a result rather than a doubt.

Only the middle case convicts.

The twelve-patch net supplies all three. Its consistent record list is reconstructing: two coin settings satisfy it, they differ by flipping every coin in the world, and no meeting anyone could arrange would tell them apart. Flip one edge of the record, changing a single entry out of thirty from same to different, and the reading goes unrealizable: two of the nineteen independent loops carry a residue, and the number of coin settings satisfying all thirty records drops from two to zero.

For the ambiguous case, take away a patch's ability to talk. Its five seams go unread, so the record shrinks to twenty-five entries. The eleven patches that go on comparing are pinned exactly as before, and the silent patch's coin is not pinned at all: four settings of the twelve coins satisfy every record that exists, in two pairs that differ by the all-over flip. The reading leaves one bit undetermined, and the bit has an address.

Repair cannot fix that. Repair works by reducing disagreement, and there is no disagreement anywhere: every record is satisfied. A rule that picked one of the two would be picking on grounds no record contains, which makes it a choice rather than a repair.

## Four demands

Run the machine on the flipped-edge record anyway. Chapter six did, from all two thousand and forty-eight distinguishable starting states and under all sixteen schedules. It converged every time onto one terminal arrangement with a leftover of exactly one. The machine halts. It halts in the same place from everywhere. And the state it halts in violates a record, which is to say that halting and being consistent are two separate properties. The second one has to be checked by somebody.

The check is a procedure with a four-line specification. Take any state at all, read off its record, and return the consistent state that record names. It has to meet four demands. What comes back carries the record it was handed, since a procedure that revised the testimony on its way past would be answering a different question. A state that was consistent to begin with comes back untouched. Two states with the same record come back with the same answer. And when no consistent state carries the record at all, what comes back is a failure value, written bottom, rather than an answer. There is exactly one procedure meeting all four. It exists precisely when no reading is ambiguous.

The third demand is the one doing the work, and the one an ordinary repair procedure fails. Returning the same answer for two states with the same record forbids the procedure from consulting anything the record does not contain, including which state it was handed and which repairs it happened to run. Hand it an ambiguous reading and no procedure meets all four, because the third requires one answer and the second requires two.

Bottom is a proof that the fiber is empty, which is not what a search reports when it runs out of patience. On the twelve-patch net you can hold the proof in your hand. Add up the entries around each of the nineteen loops in the arithmetic where one and one make zero. Two of the sums come out one, and a loop summing to one demands that some coin differ from itself. Nineteen sums settle four thousand and ninety-six settings without a coin being set, and they name the two loops the obstruction sits in, which no amount of running through the settings would.

A normal form in chapter nine's sense was what a process stops at, which made it a statement about the process. This one is a statement about the reading: the settled state is determined by the record rather than by the path, so it is an **observation-determined normal form**. Two observers who agree on the record agree on the world. They can be handed the record by anybody, in any order, having watched none of the repairs.

Which makes the settling two functions run one after the other. Write N for the whole procedure, the one handed a state and returning the world it settles into. The third demand says N consults what it was handed only through the record, so it comes apart into two steps:

$$N(s) \;=\; \text{world}\big(R(s)\big)$$

Here world is the lookup that takes a reading to the consistent state that reading names, and returns bottom when the fiber is empty. The state s appears on the right inside the record map's brackets. The left-hand side reaches it through nothing else. Two networks holding different coins and writing down the same thirty entries settle into one world, and whatever separated them at the start is a difference neither can carry forward, because carrying it forward would take an entry, and the entries match.

Nineteen sums decided the twelve-patch case because its constraints are linear. For constraint systems at large, deciding whether any consistent state carries a given reading is NP-complete: hand somebody a proposed state and they can check it against the reading in seconds, and no known method of producing one beats searching. The ambiguity question comes out the same shape. Two consistent states sharing one record can be put in front of anybody in a moment, and ruling out every such pair is the expensive half. The record determines the world. Working out which world is a second problem, its cost set by the shape of the constraints rather than by the length of the record.

## Flip all twelve coins

Take the settled state of the twelve-patch net and turn every coin over. Every same stays same. Every different stays different. All thirty records are untouched, so the record map returns exactly what it returned before, and the two states are indistinguishable to every observer in the net, to every observer who joins later, and to any sequence of meetings anybody cares to arrange. That is why chapter six counted two thousand and forty-eight states rather than four thousand and ninety-six: the two members of each pair are one thing.

The silent patch's undetermined bit looks similar and is not. Flipping that coin also leaves every existing record alone, but the moment the patch's seams come back into use, somebody writes down a comparison that comes out differently. It is a difference this record fails to see. The all-over flip is a difference nothing could see, under any arrangement of comparisons, ever.

Changes nothing could see form a group, in chapter five's sense: doing two of them in a row is another one, doing nothing is one of them, and every one of them can be undone. The group has a job here that chapter five could only describe. **The symmetry group is exactly the set of changes the record map cannot see.** On the twelve-patch net that group has two members: leave the coins alone, or turn all twelve over.

That fixes what physics is allowed to talk about. Only quantities that survive a change of description may be called physical. Chapter five arrived at that sentence by watching what two accounts of one table had in common, and the settled answer enforces it: relabel every patch's internal state however you like, run the repair, and the result is the relabeled version of the result you would have got. The answer is a statement about whole families of descriptions rather than about any one of them, which is called **presentation invariance**. It means no internal labeling choice anybody makes can reach the physics.

So the word objective has an exact referent here, and it is a state. Objective reality is the observable normal form of the network: the consistent state that the shared record picks out, up to changes nothing can see. There is no further fact about the world hiding behind it, because a further fact would be a difference no record distinguishes, and those are exactly what the group deletes.

## Paths and fibers

Schedule independence and reconstruction sound alike. They are statements about different objects, and the first does not give the second.

Confluence is a property of paths. It takes one starting state, lets two repair sequences run out of it, and says the branches come back together. Reconstruction is a property of fibers. It takes two starting states that happen to carry the same record and says they settle to the same place. The first says nothing whatever about the second, because the first never mentions a second starting state.

The silent patch is the counterexample, small enough to check. Everything about it is confluent: every repair sequence terminates, every schedule from a given start lands in the same place, and Newman's lemma applies without amendment. Two observers who start from different coin settings and hold identical records land in different worlds anyway, because the silent coin is never touched by any repair and never recorded by anybody. Schedule independence in full, and a coin nobody can name.

Run the two together and you get a network that passes every confluence test there is and has no objective world in it. Any two observers can compare records forever without turning up a discrepancy. The thing they are converging on differs between them in a way no comparison reaches. The repairs are finished, the schedules agree, and there are two worlds.

The traffic runs the other way with a condition attached. A repair moves coins and leaves entries alone, so two branches out of one start carry that start's record all the way down, and if both branches stop, a record that pins down one consistent state hands them the same world. Schedule independence falls out, provided every state has somewhere to settle at all. Three properties, then: branches out of one state rejoin, endpoints from different states with one record agree, and every state has an endpoint to reach. The first and third leave the second free to fail, which is the silent patch. The second and third force the first.

## One over n

The trichotomy counts the members of a fiber and says nothing about how far apart they are. The distance is what decides whether the answer survives a closer look.

Instruments quote a number for exactly that, called sensitivity: how much of a difference in the world it takes to move the dial by a readable amount. A thermometer whose scale marks crowd together at the top of its range reads temperatures perfectly well and separates two nearby ones badly. A buyer who checks only that each mark has its own temperature has checked the wrong thing. A reconstruction is an instrument pointed backward, from the dial to the world, and needs its sensitivity quoted in that direction.

A reading can be refined. Resolve more seams and the record gets longer, and throwing away the extra entries hands back the shorter record you started with. Stack those refinements and you have a tower. The question is whether the answers at each level are converging on anything.

Build one. At level n there are two consistent states a full unit apart, and the record separates them: it reads 0 on the first and 1 over n on the second. Two distinct readings picking out two distinct states is reconstruction, so every level of the tower reconstructs. Then watch the levels go by. The readings converge on each other, the states sit a unit apart wherever you look, and the factor by which a difference in the record gets magnified into a difference in the state is n. At level ten thousand the record separates two worlds by one part in ten thousand. The limit has one reading and two worlds.

So being determined by the record is not the same as being determined stably by the record. What separates them is the largest factor by which a difference in the record can be magnified into a difference in the state, which is the conditioning of the reconstruction, and a tower of readings has a limit worth speaking of exactly when the conditioning is bounded all the way up rather than level by level. Chapter ten's curve went to zero across all 390,924 seams in the network, and 390,924 seams is one level. The bound is what says a finer one would have gone to the same place.

This is why a continuum is available at all. Space that looks smooth at the scale of a laboratory is a limit of finite readings, and a limit taken through a tower with unbounded conditioning is an artifact of whichever level you stopped at. With the bound, different solvers, different refinements and different orders land on one limit. The smooth description is a thing about the world rather than about the grid.

## Hold the boundary fixed

An accepted repair leaves the record it was handed exactly as it found it, so a settled state carries the record it started with, and what it settles into depends on that record from the first move to the last.

The settled state is unique given the boundary data, which is the part of the record fixed from outside the region doing the repairing. Every start carrying that boundary data lands in the same world. Change the boundary data and the region settles into a different world, just as consistent, just as settled, and just as indifferent to who went first.

On the twelve-patch net the boundary data is easy to point at. Repair a region of the net and the entries on the seams that leave the region are handed to it from outside: the patches inside cannot change them, and every setting the region can settle into has to agree with them. Those entries are the region's boundary, and the settled interior is a function of them. Two neighboring regions with different boundary entries settle into different interiors, and both are correct.

There is no single endpoint that the whole thing is heading toward. The terminal state is unique relative to its boundary record, and boundary records are not unique. Two different readings on the edge of a region produce two different settled interiors. Both of them are objective in the only sense the word has here, which is that every observer inside either one agrees with every other.

Which is the whole of what a jury does. It is handed a record it did not choose, and the verdict follows from the record. Hand it a different record and a different verdict follows, with no less force. Twelve people, one transcript, and a verdict that belongs to the transcript rather than to any of them.

So look at what the settled world is made of. Thirty entries, each one a comparison between two patches, none of them a property of either patch on its own. The coins were the things in this world, in the sense that they were what the patches actually held, and the coins are the part that came out in the wash: two settings of them produce every entry identically, and the answer is stated about the pair rather than about either member.

Thirty entries survived the procedure. Twelve coins did not appear anywhere in the answer, and the record map never once returned a coin.