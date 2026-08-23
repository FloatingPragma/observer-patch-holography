# 4. Why Do Two Descriptions Ever Agree?

Alice saw the food truck. Bob did not.

They were walking the same street in Barcelona, on opposite pavements, thirty meters apart, conducting the standard argument about where to eat. A delivery van had stopped on the corner of Carrer de Sepúlveda with its back doors open and a man unloading crates of bottled water, and from where Bob was walking the van filled the gap between two buildings completely. From where Alice was walking the van was a van, and behind it, in the gap, was a truck selling grilled sandwiches with two people queueing and a third buying two at once.

At the corner they compared accounts. Alice said there had been a food truck. Bob said there had been no food truck, an empty stretch of pavement, and a man with a great many crates of water.

Both accounts were accurate. That is the difficulty.

Alice and Bob have been doing this since February 1978, when Rivest, Shamir and Adleman needed two people to exchange a message in a paper on public-key cryptography and preferred names to letters. They have since been assigned to nearly every problem in science that turns on two parties who cannot see what the other one knows, a category that takes in more of physics than anybody expected in 1978, and in close to fifty years of continuous employment neither of them has ever been told what the message says.

The comfortable picture of what agreement means says that there is one street, that the street has a definite content, that each of them writes down what the street contains, and that if both of them do the job properly the two descriptions come out the same. Two people did the job properly and the descriptions came out different, and no amount of care by either of them would have fixed it, because the van is opaque and Bob was on the wrong side of it.

The question is harder than it looks. Take two record-keepers of the kind chapter three described, each seeing a bounded part, neither able to inspect the other's records, with no third party anywhere holding both. What would it even mean for their two accounts to agree? Something has to be demanded of them, or there is no shared world. Whatever is demanded has to be something they can check from inside their own windows, or it might as well be demanded of the weather.

## What can be demanded

Nothing can require Alice and Bob to match on the food truck. Only one of them was in a position to see it. A rule that penalised them for that would be a rule against having a point of view, and every observer in this book has one by construction.

What can be required is that they match where they both looked. Both of them saw the van. Both of them saw the corner, the traffic, the time of day, the number of lanes. That shared region is called the **overlap**, and it is the only place where a demand can land. Outside it, two accounts differing is two windows differing, which carries no information about anybody's competence.

Even inside the overlap the two accounts will fail to be identical, and this is worth being precise about, because it is where the naive picture does its second piece of damage. Alice wrote down that the van was on her left. Bob wrote down that the van was on his right. Those are different sentences. They are not competing claims about the van, and nothing is wrong with either of them.

What holds between them is a rule. To turn Alice's account into Bob's, swap left and right, because they were facing opposite ways. Apply the rule to Alice's whole record and you get Bob's record for everything they both saw. So the demand that survives contact with two points of view asks for a rule, and asks two things of it. That it carry one account onto the other everywhere both of them looked. And that each of them be able to work out what it is without being told.

Where the rule fails to line up, there is a leftover. If Alice has the van two car-lengths from the corner and Bob has it four, no relabeling of left and right repairs that. The leftover is a quantity. You can measure it, you can add up all the leftovers across everything two parties both looked at, and you get a single number that is zero exactly when their accounts are compatible and positive when they are not.

That number is the most important quantity in the book. Physics is what happens when something tries to make it smaller.

## Two detectors

Physics takes this problem seriously enough to spend money on it, and the money is the argument.

On 4 July 2012, two collaborations at CERN announced that they had found a new particle. ATLAS gave its mass as 126.0 billion electron-volts, with a statistical uncertainty of 0.4 and a systematic uncertainty of 0.4, at a significance of 5.9 standard deviations. CMS gave 125.3 billion electron-volts, give or take 0.6, at five standard deviations. Five is the bar particle physics sets before anybody uses the word discovery, and it means the odds of a fluctuation doing it are about one in three and a half million.

Those are different numbers. 126.0 is not 125.3, and a reader who was told that physics had converged on a value might reasonably ask which one it converged on. Neither, is the answer. What was announced was that two intervals overlapped, and overlapping is what agreement between measurements has always meant. A measurement is a claim about a range. Two measurements agree when the ranges have points in common, and the sharper the measurements the more that agreement is worth, because the ranges are smaller and there is less room to overlap by luck.

The reason there were two of them is the whole point. ATLAS and CMS sit at different places on the same ring and were built by different people to different designs, and the differences are not cosmetic. ATLAS bends particle tracks with a thin two-tesla solenoid and then bends them again in three enormous air-core toroids. CMS does it with one superconducting solenoid six meters across running at 3.8 tesla, which is a different way to spend the same money. Different magnets, different detector materials, different trigger electronics, different reconstruction software, different graduate students, different arguments in different corridors.

Each collaboration also blinded itself. The region of the data where the signal was expected to sit was masked while the analysis was being built, and it was unmasked only after the collaboration had committed to its procedures, so that nobody could adjust the method until the answer came out pleasing. Neither experiment tuned its analysis to the other's number.

Building the second detector cost roughly what building the first one cost, which is a remarkable sum of money to spend on a second opinion. It bought the one thing a single detector cannot produce at any price. One detector reporting 126 produces a number. Two detectors that share no components, no code and no staff, reporting 126.0 and 125.3 on the same afternoon, produce something that a number by itself is not: a fact about the world rather than a fact about an apparatus.

That is the street corner again, with a budget. Alice and Bob compared accounts on the overlap and found them compatible. The overlap of two collaborations is the particle, and the compatibility is the discovery.

## A clock face

It is quarter past three in the afternoon. The clock on the wall says three. A railway timetable says 15:15. Nobody is confused, and nobody thinks the clock and the timetable disagree about what time it is, because everybody knows the rule: fifteen and three are the same thing on a clock face.

Look at what that sentence does. Fifteen and three are different numbers. They stay different numbers. The clock declares that for its purposes, two numbers count as one thing when they differ by twelve, and having made that declaration it stops carrying the difference. A clock is a machine for forgetting the difference between three and fifteen, and forgetting it is the entire service the clock provides.

Do that systematically and you get the tool. Start with a collection of things. Choose a rule that says which of them count as the same. The rule has to be sane, in a way that takes ten seconds to check: everything counts as the same as itself; if this counts as that, then that counts as this; and if the first counts as the second and the second as the third, then the first counts as the third. Sweep up each family of things that count as one another into a single bundle. Each bundle is called an equivalence class, the collection of bundles is called a **quotient**, and the operation of passing from the original things to the bundles is the single most useful operation in the book.

Twelve bundles come out of the clock. One of them contains 3 and 15 and 27 and 39, and that bundle is what "three o'clock" names. There is no such thing as the true underlying number that three o'clock secretly is. There is a bundle, and 3 and 15 are two ways of writing it down.

Which brings the two friends back. Alice's record and Bob's record differ. Where they differ by something a permitted relabeling accounts for, they are two ways of writing down one thing, in the same sense that 3 and 15 are. Where they differ by more than that, the leftover is real and somebody's records have to change.

Everything hangs on the word permitted. The clock's rule is not a free choice made for convenience. If you declared that any two numbers count as the same you would get one bundle, which forgets everything, and a clock that says the same thing at every moment. The permitted relabelings are exactly the ones that no measurement made inside the system can detect, and no others, and working out which those are is a large part of the job of physics.

## The rate has to be computable in the room

One clause on all of this carries the rest of the book, and it is easiest to see with money.

Two firms keep books in different currencies and share an account. To settle it they need an exchange rate. In ordinary commercial life they look one up, which works because there is a market outside both firms that publishes rates all day.

Chapter three closed off that move. There is no outside. There is no venue where the correct description is kept, no referee holding both sets of books, and any procedure that requires one has smuggled in exactly the stage that chapter said cannot be assumed. The exchange rate has to be computable from inside the room, out of things both parties hold.

It can be. The two firms have both recorded the same transaction, one of them in each currency. The ratio of those two entries is the rate. It came out of the overlap, which is the one place both parties have data, and neither firm had to be told it by anybody.

Two properties come free with a rate obtained that way. It can be undone, because if the rate one way is known then the rate the other way is its reciprocal. And two of them can be chained, because a rate from the first currency to the second, followed by a rate from the second to the third, is a rate from the first to the third.

Chaining is where it stops being bookkeeping. Take three firms, each pair of which shares an account, so there are three rates. Convert from the first currency to the second, from the second to the third, and from the third back to the first. You have gone in a circle and arrived back where you started, and the amount you are holding had better be the amount you set out with.

If it is not, somebody can go round the loop repeatedly and end up with more money than they started with, having bought and sold nothing. Currency traders call this triangular arbitrage, they have machines watching for it, and when it opens it closes in milliseconds. The condition that closes it is the requirement that going around a loop of translations brings you back to yourself.

That condition is the whole content of what physicists call gauge structure. Accountants make it obvious and field theory makes it look hard, and it is doing the work behind electromagnetism twenty chapters from here. Descriptions related by a permitted relabeling say the same thing. The relabelings can be undone and chained. And every closed loop of them has to come back to the identity, or there is something to be extracted for free, and the world does not offer that.

## Ninety meters

In June 1792 two French astronomers left Paris in opposite directions to measure the world.

France wanted a unit of length that no king had defined, and the Academy of Sciences had settled on one ten-millionth of the distance from the north pole to the equator. To get it, somebody had to survey the meridian arc from Dunkirk to Barcelona by triangulation and scale up. Jean-Baptiste Delambre took the northern leg, Dunkirk down to Rodez. Pierre Méchain took the south, Rodez to Barcelona. The plan allowed a year. The survey took seven, most of them conducted through a revolution and a war, with both men repeatedly arrested by people who found a stranger on a hilltop with brass instruments and a telescope to be exactly what a spy would look like.

Méchain finished his southern latitude work at the fortress on Montjuïc, the hill above Barcelona harbour, and had earlier measured the latitude at his lodgings in the city. Both were careful measurements of the same quantity by the same excellent astronomer with the same instrument. In March 1794 he compared them and they disagreed by three seconds of arc.

Three seconds of arc is about a six-hundredth of the width of the full moon. On the ground, it is a little over ninety meters. It is a very small number, and it was fatal, because the two figures were both his and there was no third figure to adjudicate. War with Spain had closed the border behind him and he could not get back to the hill to remeasure.

He concealed it. When the international commission met in Paris in 1799 to fix the meter, the numbers Méchain supplied had been adjusted to remove the discrepancy. He spent the rest of his life trying to make the disagreement go away, went back into Spain in 1803 to redo the work, contracted yellow fever, and died at Castellón de la Plana on 20 September 1804 at the age of sixty. Delambre, preparing the official record of the survey afterwards, worked through the notebooks and found the alterations.

Set aside the concealment, which belongs to Méchain. The structure underneath it is the one the food truck had. Two records of the same quantity failed to agree on their overlap, the leftover would not go to zero, and the one move that could have driven it to zero, going back to the hill with the instrument, had been made unavailable. He was in the position of Alice and Bob with the van permanently in the way. A mismatch that cannot be repaired does not resolve itself and does not fade. It sits there.

## What the meter is

The meter they produced is wrong.

The quadrant they were measuring, from pole to equator, is 10,002,290 meters. They were trying to make it exactly ten million. Their meter is therefore short of its own definition by about a fifth of a millimeter, which is roughly the thickness of two sheets of paper, and the error has nothing to do with Méchain's three seconds of arc. It comes from the flattening of the Earth, which they had to estimate and which was not then known well enough.

Nothing whatever has gone wrong.

Every meter stick in the world was made to match the standard, and every meter stick in the world therefore matches every other one. The definition was replaced in 1889 with a platinum-iridium bar, replaced in 1960 with a count of 1,650,763.73 wavelengths of the orange-red light of krypton-86, and replaced again in 1983 with the distance light travels in one 299,792,458th of a second. Three replacements, and the meter has never changed length, which is a peculiar thing to be able to say about a length.

It has never changed length because the meter was never the bar in the vault. The meter is the bundle. It is the equivalence class of every stick, wavelength and timing that counts as agreeing with every other one, and a definition is a way of pointing at the bundle by holding up one member of it. Point at it with a bar, point at it with krypton, point at it with the speed of light and a clock, and you have pointed at the same bundle three times. A better representative buys precision. It has no power to change what is being represented, because it was never carrying the identity in the first place.

This is why the survey's error is not embarrassing and why French commerce did not collapse in 1799. The unit's job is to let two parties translate their measurements into each other's terms. It does that job through what everybody adopts, which fixes the translation. The pole is not consulted.

## One world

The naive picture of agreement wanted one street with a definite content and two matching descriptions of it, and a delivery van was enough to take it apart. There is no master copy that anybody holds a copy of. There are records held by parties who cannot see each other's, and there are overlaps, and on the overlaps there are translations, and the translations are constrained: they must be computable from inside, they must be undoable, they must chain, and every closed loop of them must come back to where it started.

What makes it one world is that the leftovers on the overlaps go to zero. That is the only sense in which Alice and Bob live in the same city, the only sense in which ATLAS and CMS found the same particle, and the sense that failed for Pierre Méchain on a hill above Barcelona.

Which raises the question the next chapter exists to answer. Alice and Bob compared notes and their accounts were compatible. Compatible accounts are two labelings of one thing, in the way that 3 and 15 are two labelings of one bundle, and the labeling is thrown away by the comparison.

So something survived the comparison, and it was not either of their descriptions.

They agree. On what?