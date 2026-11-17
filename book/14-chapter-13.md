# 13. What Is One Piece of Reality Actually Made Of?

The standard way to test a circuit board used to be a bed of nails: a fixture with a spring-loaded probe for every pad, which pressed the board down onto several hundred needles at once and read off which connections had taken solder. It worked while the connections were where a needle could touch them.

Then the boards got dense. Traces went down into inner layers, and chip packages arrived with their contacts in a grid underneath the body of the part instead of along its edges, so the pins of a soldered-down chip sat sealed in a sandwich a probe cannot enter. The nails had nowhere to land.

In 1985 a group of electronics firms formed the Joint Test Action Group to deal with it. Their answer gave up on reaching the pins and asked the chip to report itself. Put a small register cell beside every pin, chain those cells into one long shift register running around the boundary of the die, and add wires that let an outsider clock the chain out a bit at a time. In 1990 the arrangement became IEEE Standard 1149.1, titled Standard Test Access Port and Boundary-Scan Architecture. The port takes four wires: a test clock, a mode select, data in, data out, with an optional fifth for reset. Behind those four wires sit several million transistors, and the standard says nothing about them.

Hand somebody an unfamiliar board and the first thing they look for is the port: four or five pads in a row on a tenth-of-an-inch grid, often with no header soldered into them. Those pads say the board expects to be interrogated and on how many lines, before anybody knows what a single chip on it is.

Chapter twelve left twelve observers wired into thirty seams, nineteen independent loops, and one unit of charge sitting where no relabeling reached it. All of that came off a wiring diagram. A wiring diagram says which port is joined to which port and says nothing whatever about what sits behind either one.

Open one. The finite repair example has six parts and six moves. They describe an observer service with specified rules. The common matter and field model adds quantum evolution between observations: fermionic occupations and electric flux interact with local registers for energy response, volume, capacity and memory. Its Hamiltonian specifies that evolution before a run. A repair step belongs to the preparation or observer apparatus when that step is explicitly installed. It is not an extra operation silently inserted between moments of free evolution.

## The parts list

Start with the thing that holds a state. A **patch** is a machine small enough to be described completely: a bounded amount of internal state, a set of interfaces, and a finite list of moves it is allowed to make. Bounded has been in force since the definition of an observer. A patch that could hold unlimited state would have an unlimited amount to compare. The comparison would never finish.

A **port** is one place where a patch exposes part of itself to whatever is outside it. What appears at a port is a **packet**: a reading, drawn from a finite set of possible readings, that anything on the other side can pick up. The packet is not the patch's state. It is what the state looks like from the outside at that one interface, in the way that the cell beside a pin holds what that pin is doing rather than the internal state of the die.

Suppose a port could show any value on a continuous dial. Reading it to a hundred decimal places takes a hundred digits of storage in the reader, reading it exactly takes an infinite amount, and a machine with bounded state cannot hold either. A neighbor who writes down the first four digits and drops the rest has not read that port. It has read a different, coarser port with a finite list of settings, which is the port the architecture had all along.

Two ports routed to each other make a **seam**. Chapter six drew these as edges and counted thirty of them. A seam is where two patches read the same thing twice, once from each side, the only place in the architecture where a disagreement can be detected at all.

Everything any patch can ever learn about its neighbor arrives through the seams they share. Collect all of it, every reading available at every interface an observer has, and what you are holding is that observer's **screen**. A screen is where an observer's accessible information lives, in the same sense that the boundary-scan chain is where a chip's pin states live. It shows nothing to anybody, it faces no direction, and there is no seat in front of it.

An observer assembled out of several patches has a screen too, the readings on the ports facing outward from the assembly. The ports that patches inside the assembly share with each other are internal wiring to everything beyond it. So the size of a screen is a count of readings rather than a width in meters, and asking what an observer's screen looks like is asking how many distinct things its outward ports can be doing at once.

Last, a **record**. A patch that completes a move writes the outcome into the part of its state that later moves are required to leave alone. Chapter three priced this and chapter nine identified what gets written, which is the accepted repair. A record is a write that survives the next read, and survives it because the rules of the machine protect that region of the state rather than because writing is somehow permanent.

Six words, and every one of them names something you could point at on a board: the machine, the interface, the reading at the interface, the join between two interfaces, the set of all the readings, and the part of the state that stays put.

The last of the six carries the most. Everything an observer will ever call a fact is a record that survived every comparison anybody could make of it, and the world observers share is assembled out of those records afterwards, as chapter eleven showed. The record comes before the world it describes.

## Six moves and no others

A patch does one thing, over and over.

It reads the packets exposed on its ports. It compares them against its own, through a rule saying which of its readings answers which of the neighbor's, since two patches keeping their state in private conventions have to be told what corresponds to what before a comparison means anything. That rule is chapter twelve's dictionary on the seam. The technical word for whether the comparison is possible at all is commensurability. Then, if the comparison came back unequal, the patch chooses a move from a finite menu of allowed local updates, one entry off a printed card. It writes a record of the move it made. It exposes the resulting packet at its ports, so the neighbor's next read sees the new state rather than the old one. And any claim about what it did is settled by an exact check run on the exposed data, by anything that has access to that data, and the patch's own account of its work does not count.

Read, compare, choose, write, expose, verify. Those are the operations of this finite repair service. No move reaches a patch on the far side of its network, none consults a clock, and its menu has an end to it. A complete quantum measurement also retains the conditional quantum states and every complementary outcome. A classical repair trace alone does not specify that instrument or the intervening matter and field dynamics.

Chapter twelve's large arrangement was running that loop at 122,880 seams at once. A patch reads what its neighbor has put on their shared seam. It applies the dictionary sitting on that seam, which tells it which of the neighbor's three labels answers which of its own. Suppose they do not answer each other. The menu on that patch holds the six ways of permuting three labels, and the entry that gets taken is the one the dictionary names: the neighbor's label carried across the seam and written down at this end as its own. The dictionary is not on the menu, and no move in the run rewrites it. The patch takes one entry, writes down that it took it, and puts the resulting packet back on the port for the neighbor's next read. Nothing in that sequence required knowing anything at all about the far side of the network, which is the property that lets 81,920 of these run at once with no coordination between them.

The sixth move is the one engineers recognize on sight. A patch can be wrong. It can have run its comparison against a stale packet, or chosen a move that was not available, or written a record that does not follow from what it read. What settles the matter is a check performed on the readings themselves, which either passes or does not, and which anybody holding the same data performs identically. Anything that runs that check is a verifier. The data it runs on is chapter eleven's list everybody can check, the seam entries and nothing private to either end of a seam, so a verifier needs no access to the inside of the patch whose work it settles.

## Adding up the seams

Each seam can score itself. Take the two packets, apply the seam's dictionary to one of them, and produce a number that is zero when they answer each other and positive when they do not. That number is local: it involves the two ends of one seam and nothing else in the world.

Sum the seam scores and accept a move only when it lowers that sum. The sum takes finitely many values, and a quantity that strictly falls through finitely many values can fall only finitely many times, so the moves terminate.

Chapter ten watched that number come down. Of the 122,880 seams in that run, 102,415 were scoring above zero at the start, and at cycle 1 out of 16 none of them were. The settling curve counted the seams contributing anything at all. The quantity the machine was working against was their sum.

So every sequence of repairs stops. Wherever it stops, no patch has a move available: each one has read its ports, compared, and found nothing its menu can improve. That is a local normal form, exactly as much as the descent argument buys.

Everybody stopping and everybody stopping in the same place are two different claims. The score establishes only the first.

## What the second claim costs

Three more conditions, and each is a property the interface either has or lacks.

The first is that the repair menu is complete: for every mismatch a seam is able to report, the menu contains a move that clears it. Chapter nine showed what a missing entry looks like. Twelve observers passing units across seams, and a rule that acts whenever two neighbors differ by two or more. In 302 of the 1,352,078 arrangements, no seam has a gap of more than one, no rule can fire, and the arrangement sits lopsided. Every patch content, every local check passing, and the world not settled. The gap of one had no entry on the card.

The second is about the inside of a patch rather than the seam. When a repair changes what a patch exposes at a port, the patch has to have an interior state consistent with the new reading. If some admissible boundary value has no interior that produces it, then clearing the seam breaks something inside, repairing the inside changes what the port exposes, and the two repairs chase each other around the boundary forever. The condition that rules this out is a gluing condition: every value a port can legally be asked to display extends to a consistent state of the patch behind it.

A third requirement comes out of the relabeling freedom of chapter twelve. Two repairs leaving one state have to be reconcilable, which is chapter ten's diamond, and the diamond has to close in the description with private conventions divided out. Two patches that relabeled themselves differently have landed in the same place as far as anything checkable goes, so closing the diamond in the raw labels is the wrong test in both directions: it can fail where the world agrees, and pass where the world does not.

With the menu complete, the gluing condition holding, and the diamond closing where conventions have been divided out, the terminal state of the network is the same whatever order the repairs ran in. Chapter ten proved the shape of that argument on a shelf of books and a mountain with three basins. This is what it costs in parts: a scoring rule on every seam, a printed menu with no gaps in it, and interiors that can support every boundary reading their ports allow.

Two patches on one seam gave chapter ten two endings out of one start. The score in that arrangement runs down to zero, every move taken is an entry on the card applied to a mismatch that was there, and the sixth move's check passes on every reading either patch exposed. The descent supplies all of that and leaves both endings standing. A gap in the menu and an interior that cannot support its boundary are both printed into the card and the wiring before the first repair runs.

Those three requirements are the **synchronization contract**: obligations on the parts. Any assembly that meets them synchronizes. The descent gets the machine to stop. The contract is what makes the place it stops one place rather than a place per schedule, which is to say that whether there is a single world at the end of the repairs is a question about the connector.

## Four things that have to travel

A host needs maintenance, so the live machine running on it is stopped, written down and started again on another host, with the programs inside carrying on through the move. Server operators do this on a schedule. What makes it work is a decision about which parts of the machine have to travel and which parts can be left behind and rebuilt, and getting the list wrong shows up as a program that resumes into a world it does not recognize.

A machine that is interrupted and resumed is the same machine when enough of it has been put back. Enough, for an observer, is four things. Its records, meaning the writes that later moves protect. Its accessible state, meaning the state of what it can actually reach. Its interfaces, meaning which ports it has and what is routed to them. And its future law, meaning which moves and which readings its machinery allows. Those four together are a **checkpoint**.

The fourth item is the one that gets forgotten. Restore an observer's records and its accessible state into a machine whose menu of moves is different, and the restored thing reads the same history and takes different steps from it. The law that says which moves are available has to travel with the records, or the records land somewhere that treats them differently.

Two checkpoints that agree exactly on all four give the same probabilities to everything the observer subsequently sees, for every event it is able to register. If the two restored states differ slightly, the two futures differ by at most as much, which is what makes a rebuilt observer usable when the copy is imperfect.

A great deal of what a running system holds is outside that list. Worker identifiers. Queue positions. Retry counters, repair-cycle indices, timestamps, and the latency of individual packets. All of it is real, all of it is needed to replay the machine's history, and none of it is part of the observer or of the observer's time.

## A code with no distance in it

Look at what the seams do to the space of configurations. Each one demands that its two packets answer each other under its dictionary. A configuration that satisfies every one of those demands at once is legal, and one that fails any of them is not. Chapter twelve's arrangement had 81,920 patches holding one of six labelings apiece, which is finitely many possibilities, though writing out their number takes 63,747 digits, and that is before any seam dictionary is counted.

An object of that description has a standard name in engineering: it is a **constraint code**, the same kind of thing as the codes that protect a hard drive or a deep-space link, where a handful of parity checks pick the legal words out of all the possible words. This carrier is a finite constraint code.

What a code usually arrives with is a distance. Send every bit three times and the legal words are 000 and 111, which differ in three places. Corrupt one symbol and you hold 010, which is one step from 000 and two steps from 111, so the nearer legal word is the one that was sent and the error is corrected rather than merely noticed. Corrupt two symbols and 011 is nearer to 111, and the correction produces the wrong bit. The number three is the distance of that code, and the distance is what says how much the code can fix: a distance of three corrects one error, a distance of five corrects two, and a distance of two catches an error without being able to say which symbol moved.

The architecture supplies the checks and does not supply the distance. Nor does it supply a bottleneck across which information has to squeeze, a rate at which disagreement mixes through the network, or a bound on how many seconds settling takes. Each of those is a property of a particular wiring diagram, and has to be computed for that diagram, from the diagram. The six moves say which configurations are legal. How far apart two legal configurations sit depends on which port is joined to which port.

## What one boundary can hold

Take a patch and ask what anything outside it can determine about it.

Everything arrives through the ports. What arrives at a port is one packet from a finite set. So two interior states that expose the same packets on every port are the same state to every neighbor, to every neighbor's neighbor, and to the whole network. They can differ internally in any way you like. The difference has nowhere to be read.

Which turns the amount that can be known about a region into a count of arrangements of its boundary. The number of values one port can show, multiplied by itself once for every port the patch has. In the runs of chapter twelve the freedom each patch held was one of the six ways of permuting three labels, so a reading in that arrangement is one of six things, and the bound is six multiplied by itself once per port.

Assemble many patches into one observer and the same argument applies to the assembly. Adding a patch deep inside adds state and adds nothing to what an outsider can read, because every port it has is routed to another patch inside the assembly. The count moves only when a patch puts a port on the outside. Two assemblies with the same outward-facing ports and very different interiors are one object to everything beyond them, with the same list of possible readings.

The interior does not appear in that product anywhere. However much a patch holds inside, and however finely it holds it, what an outsider can distinguish is bounded by the readings on its boundary. That is half of the area law: the information associated with a region of space goes as the area of the surface around it rather than as the volume inside, and black holes are where that was first noticed. The other half is a number, how much information one unit of boundary area holds, and fixing that number takes a horizon and a thermodynamic argument. Getting this far took a machine with a bounded interior and finitely many ports, and nothing about gravity at all.

So what can an outsider hold about the interior? Odds spread over the interiors its boundary readings allow, and nothing besides, because anything more would be a claim no port reports. The spread that adds nothing to the readings is the even one. That is chapter two's third requirement applied to a single patch, and it is the state an unread region settles into. Its entropy is the logarithm of the number of interiors an outsider can tell apart, and the boundary count caps that number. The entropy a region can carry is a fact about its edge.

Six, multiplied by itself once per port. The six is an accident of one arrangement's three labels. The exponent is a count of holes in a connector. The six moves do not fix it.

The connector is a closed surface, its ports have to be arranged so that every one of them sees the same thing, and the seams have to pair up with nothing left over. Impose them all at once and the arithmetic permits one count and no other.

Count the ports.
