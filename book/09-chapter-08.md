# 8. Why Can Strangers Agree Without a Boss?

At 22:52 UTC on 21 October 2018, maintenance work on a failing piece of 100-gigabit optical equipment cut the link between GitHub's East Coast network hub and its primary East Coast data center. The link came back forty-three seconds later.

In those forty-three seconds the software that manages GitHub's database topology did what it was built to do. It runs a leader-election protocol called Raft. It could not reach the East Coast machines, and Raft carries no test that separates a machine which has stopped answering from one that has stopped. What it does instead is count. The nodes that could still hear one another formed a majority, agreed among themselves that the East Coast primary was gone, and promoted a cluster on the West Coast to accept writes in its place.

The East Coast primary had lost nothing except the ability to talk to the other coast. Applications in its own building went on sending it work, and it went on taking it. One of the busiest clusters accepted 954 writes during the window.

Then the link came back, and each coast was holding writes the other had never seen. Neither copy contained the other, so neither could be made the primary without throwing away data that had been acknowledged to somebody. GitHub ran for 24 hours and 11 minutes in a degraded state, showing information that was out of date and inconsistent, unable to publish pages or deliver webhooks. The divergent writes had to be reconciled by hand.

There is no bug in that story. Every component behaved as specified, on the information available to it, and the information available to it was insufficient in a way that no additional care removes. Forty-three seconds of missing fiber cost a day of the largest code host on the planet. The mechanism was that two groups of machines, neither able to see the other, both concluded they were in charge.

So how does a crowd of parties who cannot see each other's internal state, who share no clock, and who answer to nobody end up holding the same answer, and what does it cost? Engineers buy that agreement with hardware, messages and delay, and they know the price to several decimal places, because they pay it every day. Alice, Bob and Charlie get theirs for nothing.

## Two ways to be broken

In 1982 Leslie Lamport, Robert Shostak and Marshall Pease published a paper in *ACM Transactions on Programming Languages and Systems* called "The Byzantine Generals Problem", about a group of commanders surrounding a city who must agree on whether to attack, and some of whom are traitors.

Lamport wanted a nationality for the traitors that no reader would object to, and settled on Albania, which at the time was sealed tightly enough that he doubted any Albanian would ever see the paper. His colleague Jack Goldberg pointed out that there were Albanians outside Albania and that Albania might not stay sealed forever. Byzantium had the advantage of having been out of business since 1453.

A **crash fault** stops. The machine dies, the process is killed, the power goes, and what everybody else sees is silence.

A **Byzantine fault** keeps talking. It sends messages that are wrong, or that are inconsistent with each other, or that say one thing to one neighbor and something else to another, and it does all of this with no outward sign that anything has happened. Malice is one way to get there and by far the least common. A memory cell that flips, a disk returning data from before the last write, a program that applied half of an update and then fell over: each of these produces a participant confidently answering questions with garbage, and the confidence is the whole of the problem, because the other participants have no way to weigh it.

On 18 May 2003 a vote-counting machine in Schaerbeek, in Brussels, credited one candidate in the Belgian general election with 4,096 votes that nobody had cast. No fault was ever found in the software. The official finding was a spontaneous inversion of a single bit in memory at the thirteenth position, a position worth two to the twelfth, which is 4,096. The machine reported the result the way it reported every other result, with no error, no warning and no sign of distress. The reason anyone caught it was that the municipality had returned more votes than it had voters.

Crash faults are cheap to survive. You wait, you notice the silence, you carry on without the missing party. Byzantine faults are expensive. The exact price can be counted on a napkin.

## Four dots

Draw four dots. They are four participants. You would like the group to survive one of them being Byzantine, without knowing which one.

Decide that any three of them, agreeing, may commit a decision on behalf of the group. Call a set that large a **quorum**.

Take any two quorums. Each has three members drawn from four, so the two of them together name six memberships among four participants, which means they share at least two members: three plus three minus four is two. At most one participant of the four is faulty. Therefore at least one member of the shared pair is sound, and a sound participant does not vote for two different answers to the same question. So no two quorums can ratify conflicting decisions, whatever the faulty participant says to whom.

That is the entire safety argument, and the same subtraction runs at any size. Surviving three broken participants out of ten takes quorums of seven. Any two quorums of seven out of ten share four members, at most three of them faulty. Surviving four out of thirteen takes quorums of nine, which share five, at most four of them faulty. Either way at least one sound participant sits in both. The property is called **quorum intersection**, the reason the ratio is two thirds rather than one half. Ordinary majority voting survives crashes. Surviving liars costs an extra third of the hardware.

Read the ratio as a purchase order. To survive one broken participant you buy four machines and use them to do the work of one. Every decision costs a round of messages from three of them. To survive two you buy seven. The bill is linear in the number of liars you are willing to tolerate and it never gets cheaper, because the subtraction that makes the argument work is doing so at the tightest possible margin: two quorums of three out of four share two members, and if either quorum shrank to two, they could share one, and that one could be the liar.

Notice also what the four dots require you to have before any of this starts: a list of who the participants are, a way to reach any of them, and a fixed count of how many there are. Take away the roster and the arithmetic has nothing left to subtract.

## What Fischer, Lynch and Paterson proved

In April 1985 the *Journal of the ACM* published a paper by Michael Fischer, Nancy Lynch and Michael Paterson titled "Impossibility of Distributed Consensus with One Faulty Process".

Suppose messages always arrive eventually, but with no bound on how long they take. Suppose one single participant may crash, and only crash: no lying, no corruption, one machine that stops. Then there is no deterministic protocol that guarantees both that everybody who decides decides the same thing and that everybody decides at all.

The result does not say that no efficient protocol exists. It says that no protocol exists. The shape of the proof is that any such system passes through configurations from which both answers are reachable, and that a schedule delivering one message at a time can always steer it out of one such configuration and into another. Every message gets delivered. Nobody is starved. The system is kept poised forever by nothing more aggressive than a choice about the order of the mail.

Faster hardware does not touch this. Better engineering does not touch it. What the result attacks is the combination of three things: asynchrony with no bound on delay, a deterministic rule, and the demand that agreement and termination both hold always. Give up any one of them and consensus comes back.

Two of those surrenders are the ones the industry actually made.

The first gives up asynchrony. Assume that after some unknown point, message delays are bounded by some unknown constant, an assumption Cynthia Dwork, Nancy Lynch and Larry Stockmeyer formalized in the *Journal of the ACM* in 1988 and called partial synchrony. The protocol never violates agreement, whatever the network does; it promises to finish only during the stretches when the network behaves. Every leader-based system in production is standing on this hatch, which is why every one of them has a timeout in it, and why the timeout is the parameter that operators fight about.

The second gives up determinism. Let participants flip coins, which Michael Ben-Or showed in 1983 is enough: the protocol terminates with probability one, meaning that the runs in which it goes on forever have a total weight of zero, while carrying no promise about any particular deadline.

A leader with a timeout, or a coin. Every consensus protocol running in production is one of those two, or both.

The failover on 21 October 2018 was the first hatch working as designed. The nodes that could see each other formed a quorum, elected a new topology and never disagreed with one another about it, which is exactly what a protocol standing on partial synchrony promises: agreement always, progress when the network behaves. What diverged was the data underneath, because the databases accepting writes on the East Coast were not participants in the election and had made nobody any promises at all. Agreement held, and it held about which machine was in charge rather than about the 954 writes.

## A leader and a clock

Since September 2022 Ethereum has run its consensus on a schedule. Time is cut into slots of 12 seconds, one validator is drawn for each slot to propose a block, and a committee attests to what the proposer produced. Thirty-two slots make an epoch of 6.4 minutes. A checkpoint at an epoch boundary is finalized once the following epoch has been agreed, which puts final on the record 12.8 minutes after the fact and backs it with stakes that can be destroyed.

What that machinery buys is a single history that every participant agrees on, in one order, with a definite moment after which nothing in it can move. What it costs is on the face of the design. There is a leader for every slot. There is a wall clock that everybody reads, because a 12-second slot means nothing to a participant with no idea what time it is. Every full node holds the entire state, because checking the next block means having the previous one in full. And every transaction is placed in the one order relative to every other transaction, whether or not the two have anything whatsoever to do with each other. The system spends real work establishing which of them came first. Bitcoin buys the same single order with a lottery in place of a schedule, a race to find a number whose hash falls below a target that is retuned every 2,016 blocks to hold those 2,016 blocks to 1,209,600 seconds. Ten minutes a block is what that works out to, and 1,209,600 seconds is two weeks. Whoever wins a round is that round's leader, and the chain the winner extends is one line.

## Twenty strangers

In May 2018 a pseudonymous group calling itself Team Rocket posted a paper to a file sharing network, "Snowflake to Avalanche: A Novel Metastable Consensus Protocol Family for Cryptocurrencies", describing a way of reaching agreement in which nobody is in charge and nobody holds the tally.

Each participant keeps two things: its current preference, and a counter. That is the whole of its state. There is no vote count, no view of the network, no record of who said what.

Each round, a participant picks twenty others at random and asks each one what it prefers. If fifteen of the twenty give the same answer, the asker adopts that answer and increases its confidence in it. If twenty rounds in a row come back the same way, it decides, and the decision does not reverse. Twenty as the sample size, fifteen as the threshold, twenty as the number of consecutive rounds: those are the documented defaults on Avalanche's primary network, settings rather than laws. The protocol built out of the sample and the counter is called Snowball.

If a bag holds a hundred marbles and thirty of them are red, the probability of drawing a red one is thirty out of a hundred: the fraction of the cases that come out red, which is the fraction a long run of draws settles toward. That is all the sampling needs. When a sample of twenty comes back fifteen to five, that fraction is evidence about the fraction in the whole population, and it is evidence that can mislead, because a sample of twenty drawn from a population split sixty to forty will occasionally come back fifteen for the minority, in the way that a fair coin occasionally lands heads six times running. The counter is what handles this. Being misled twenty rounds consecutively is enormously less likely than being misled once. The decision threshold is the dial that pushes the chance of a wrong decision below one in a billion, a bound the published parameters hold with a fifth of the network lying.

Each participant sends twenty questions per round regardless of how large the network is. Twenty when there are a hundred participants, twenty when there are a hundred thousand. The number of rounds needed grows like the logarithm of the population, which is the slowest growth on offer, and there is no proposer, no miner, no round leader and no view change, so there is no single machine an attacker can take off the air to stop the thing.

And the dynamics are metastable. Start the network at an even split and it sits on a knife edge, where the sampling itself is the only thing capable of breaking the tie. As soon as one answer runs slightly ahead by accident, participants sampling the network get that answer slightly more often, adopt it, and become part of the majority that the next participant samples. The tilt feeds itself. Past a certain point the outcome is fixed, and the only thing that broke the tie was the order in which twenty names came out of a hat.

## Where the analogy inverts

A crowd of parties, none of them in charge, none of them holding the global picture, converging on one answer: that is what a network of overlapping observers does, and the engineered version of it settles real payments in production.

Watch an Avalanche node take its sample. It picks twenty participants uniformly at random from the whole validator set, which requires it to hold a roster of every participant in the network and to be able to open a connection to any of them. That roster is a global object, and holding one is the price of the sample being uniform. The locality in that protocol lives in the size of the sample, twenty out of however many, rather than in who the twenty are allowed to be.

An observer in this world has neighbors, fixed, and communicates through those and through nothing else. The twelve of chapter six each had five neighbors, and the five came with the wiring: the same five in every round of every schedule, with no round of repair ever producing a sixth. To reach a non-neighbor an observer has to go through a neighbor. There is no address book and no directory. A channel that bypassed the wiring would itself be an overlap, at which point it would be in the wiring.

So the engineered protocol is the loose version. What an observer can be affected by is settled by the diagram of who overlaps with whom. That diagram is finite and fixed.

Snowball's safety is probabilistic. It holds except with an error probability that the parameters drive down. One in a billion is a small number rather than zero. The order-independence of a repairing network of observers is of a different kind. Where it settles is fixed by the disagreements and the wiring, and running the repairs in a different order lands on the same place, in every run rather than in overwhelmingly many of them. The twelve-observer sweep of chapter six is what that looks like from below, where every state and every schedule was enumerated rather than sampled. The count of exceptions was zero because there was nowhere for an exception to hide.

An engineer would take that guarantee over Snowball's in an instant and cannot have it, because getting it requires the thing engineers are never given, which is control over how the participants are wired to each other.

## What a fixed neighbor list buys

Fix an observer and ask a question that the engineered protocols cannot even pose: what could the reading it holds possibly depend on?

Before any repair, its own state. After one repair, its own state and whatever its neighbors held, because a repair compares readings across a shared boundary and nothing else was consulted. After two repairs, the neighbors of those neighbors, and no further. The set of observers capable of having influenced what this one holds grows by exactly one step of the wiring per round of repair. At every stage it is a finite set you could write out by name.

The twelve-observer arrangement makes the counting easy, because the wiring is regular enough that everybody is in the same position. Take any observer of the twelve. Five of the others are its neighbors, five more are neighbors of those, and exactly one, the observer directly opposite it on the solid, is three steps away. So after one round of repair, six of the twelve can have touched what it holds. After two rounds, eleven. The twelfth arrives in the third round and not before, and there is nothing it can do to arrive sooner.

Call the set the **dependency cone**: the observers within reach, where reach is counted in steps along overlaps rather than in kilometers.

Everything outside the cone is absent from the reading in the strongest sense available, which is that there exists no route by which it could have entered. A signal too weak to detect is a different situation, and a better instrument fixes that one. If a distant observer repairs, or refuses to repair, or is deleted from the arrangement entirely, nothing about what this observer holds changes until the consequence has walked in along the overlaps, one step per round, the whole way.

Two distant readings in this world can be correlated, and chapter six built an arrangement where they are correlated by force, since the twelve observers had exactly one consistent state available to them and every one of them was pinned by it. Correlation is allowed and is ordinary. What no arrangement permits is a distant party choosing what you read. There is no dial they can turn, because the only channel between you is the wiring and the wiring has a pace, so whatever they do arrives as something that walked, or does not arrive.

Correlation without communication is the shape of the Delft result from the prologue, where two laboratories 1.3 kilometers apart got answers that matched too often and could not use the matching to send each other a single bit. In an arrangement whose only channels are its declared overlaps, that combination is arithmetic about which observers are within how many steps of which. The name for it is **no-signalling**, arriving here before there are any particles to attach it to.

## Only the conflicts

One feature of the engineered protocol carries over exactly. It looks at first like an efficiency trick.

Avalanche does not order transactions. It organizes them into conflict sets, where a conflict set is the group of transactions that try to spend the same coin, and runs its sampling only within those sets. Two transactions that touch nothing in common are never compared, never ordered, and never made to wait for each other. There is no fact in the system about which of them came first, and no work is done to manufacture one.

The leader family does the opposite, and pays for it in the queue. Every transaction goes into one line whether or not it interacts with anything else in the line, which is why a payment in Lisbon waits behind a payment in Osaka that has nothing to do with it.

A network of observers repairing their overlaps is in the first camp, there by construction rather than by choice. A repair happens where a comparison disagrees. Two disagreements in unrelated parts of the arrangement are not competing for anything, do not touch a common quantity, and have no order between them, which is chapter seven's partial order doing work in this world: the comparisons that conflict are ordered with respect to each other, and the rest of the arrangement is under no obligation to have an opinion about when they happened.

The protocol is mostly a list of absences: no leader, no clock, no roster beyond the immediate neighbors, no global tally, and no order imposed on things that do not interact. What survives is a set of observers, a set of overlaps, and a rule that fires wherever two readings of the same overlap fail to match.

Alice and Bob share an overlap. Alice reads twenty across it, Bob reads twenty-two, in whatever units that boundary is measured in. Both are careful, neither has malfunctioned, and both readings are of the same thing, which is why two numbers where there should be one is a problem rather than a curiosity. Something has to move. The protocol says nothing about which number moves, or by how much, or why that answer rather than one of the infinitely many others that would also leave the two of them matching.