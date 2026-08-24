# 19. Why Can't You Ask Everything At Once?

At Petco Park in San Diego on the night of 24 September 2010, Aroldis Chapman of the Cincinnati Reds released a fastball that the stadium's tracking system clocked at 105.1 miles an hour, which is 169.1 kilometers an hour.

Two instruments were pointed at that ball. A radar gun read its speed off the shift in a returned microwave pulse. A camera read its position off the light coming back from it. Both worked, at the same moment, on the same object, and neither one spoiled the other's answer. Run the radar first and the camera second, or the camera first and the radar second, and the pair of numbers comes out the same. Nobody at the ballpark has to decide which instrument goes first, because the ball has a position and a speed and both are sitting there waiting to be read off.

That is two questions about one object, answered together, in either order, with neither answer costing the other anything. Every measurement made anywhere before 1925 was assumed to work that way.

The world declines to extend that courtesy to every pair of questions. Which pairs it allows, and what asking anyway costs the answer already in hand, are settled before any instrument is switched on.

## Socks then shoes

Put on your socks, then your shoes. Then try it the other way.

The two orders are not the same operation. The difference is visible from across the room. Other pairs are indifferent: left sock then right sock, or right sock then left sock, and you finish in the same state either way. Dressing contains both kinds of pair.

So take the operations themselves as the objects worth studying, rather than the states they act on. You can do one operation and then another, which composes them. You can scale an operation, or add two of them, in any setting where addition means something. A collection of operations that admits all three of those moves is an **algebra**. The word carries no more than that. When every pair in an algebra composes the same way in both orders, the algebra is **commutative**. When some pair does not, the algebra is non-commutative. Adding up a shopping list is commutative. Getting dressed is not.

Ordinary numbers are the commutative case. Add them, scale them, multiply them in whichever order suits you, and the product comes out the same. Relax that one clause, leave every other rule of the arithmetic where it was, and you have crossed the whole distance between the physics of the nineteenth century and the physics of an atom.

The failure has a size, and the size is what you get by doing both orders and comparing. Write A and B for two operations, write AB for doing B and then A, and the difference is the **commutator**:

$$[A,B] = AB - BA$$

That expression is zero exactly when the two operations can be swapped freely. Its size when it is not zero says how much difference the swap makes. Socks and shoes have a large one. Left sock and right sock have zero.

Chapter nine had a version of this with money in it: two edits landing on one account, each of them correct arithmetic, and the balance depending on which one was written first. Order mattering is a property that operations can have in a world with no quantum mechanics anywhere near it. The questions you can put to a physical system are operations of exactly this kind, most pairs of them fail to commute, and the failure belongs to the pair of questions rather than to anybody's equipment.

The radar gun and the camera commute. **Two questions can be answered together exactly when their commutator is zero.** A pair whose commutator is something other than zero is a pair you have to choose between.

## Ten days on Helgoland

In June 1925 Werner Heisenberg was twenty-three and his hay fever had swollen his face badly enough that he left Göttingen for Helgoland, a red sandstone island in the North Sea with almost nothing growing on it and therefore almost no pollen. He spent about ten days there.

The problem he took with him was the hydrogen spectrum. Heat hydrogen and it emits light at a set of sharp wavelengths, the same set every time, catalogued since Johann Balmer fitted a formula to four of them in 1885. The standard picture had electrons on orbits, which reproduced those wavelengths and collapsed for every atom with two electrons in it.

What Heisenberg did was to throw out the orbits. Nobody had ever seen one. What a laboratory actually produces is a table: for each pair of atomic states, how bright the transition between them is and at what frequency. So he built his scheme out of the tables, with one entry for each pair of states, and no statement anywhere about where an electron is.

Then he had to multiply two tables together, because the energy of a system involves the square of a quantity, and squares need multiplication. The rule he was forced into combines a row of one table with a column of the other. It produced an answer. It also produced a different answer when he did the two tables in the other order.

He sent the work to Max Born, who read the multiplication rule and recognized it, recalling it from student days as the way you multiply matrices. Heisenberg's paper reached the *Zeitschrift für Physik* on 29 July 1925. Born and Pascual Jordan submitted the sequel that rewrote the whole thing in matrix language on 27 September. The story that gets told about a lone twenty-three-year-old on a rock at dawn leaves out that the theory arrived through three people over four months, one of whom had to remember an undergraduate course to see what the other one had written down.

Position and momentum came out of that scheme as tables that fail to commute. Measure a particle's position, then its momentum, and you have performed a different physical operation from measuring its momentum and then its position, and the difference is a fixed quantity with Planck's constant in it. It does not go to zero as the apparatus improves, because there is no apparatus in the statement.

The uncertainty relation is that fixed quantity, read as a floor. The spread in position multiplied by the spread in momentum cannot go below it, and where the floor sits is set by the size of the commutator of the two questions. Confine an electron to a region the size of an atom, about a tenth of a nanometer across, and the spread in its speed cannot be smaller than roughly six hundred kilometers a second, against the two thousand two hundred kilometers a second at which the electron in a ground-state hydrogen atom actually travels. The electron is not being jostled by the measurement into a speed it would otherwise not have had. There is no state of the electron with a sharp position and a sharp speed for the measurement to disturb.

## Three magnets in a row

On the night of 7 February 1922, in the building of the Physikalischer Verein in Frankfurt, Otto Stern and Walther Gerlach sent a beam of silver atoms from an oven through a magnetic field arranged to be much stronger at one edge of the gap than at the other, and let the beam land on a glass plate.

A magnet in an uneven field gets pushed toward the strong side or the weak side depending on which way it points. Silver atoms come out of an oven pointing every which way, so the deposit should have been a smear.

It was two spots.

The deposit was too thin to see at all until Stern leaned over the plate to look, and the sulfur in his cheap cigars turned the silver into black silver sulfide, which developed the trace like a photograph. Two spots, nothing in between, and the beam had been split into precisely two by a question that could have been asked along any direction they chose to point the magnet.

That is what **spin** is, operationally. Point a magnet along a direction and you get a two-outcome measurement: the atom goes one way or the other, and the only thing you get to choose is the direction. There is no third outcome and no continuous dial. Nothing rotates, either. The name arrived in 1925, when two graduate students in Leiden, George Uhlenbeck and Samuel Goudsmit, proposed that the electron was a small spinning ball of charge, and the word outlived the picture by a century.

Chain three of these and the arrangement stops behaving like a set of properties being read off.

Send the beam through a magnet pointing up, and keep only the atoms that went up. Send those through a second magnet pointing sideways. They split fifty-fifty, and every atom in that beam had been certified up a moment earlier. Keep the ones that went sideways-left. Send those through a third magnet pointing up again, and ask the question whose answer is on record.

Half of them go down.

The atoms were up. The sideways question was asked. Half the certified-up atoms come out down. Skip the middle magnet and every atom comes out up, every time, so the middle magnet is doing it.

The temptation is to say that the sideways magnet knocked the atoms about, and that a gentler magnet would knock them about less. Gentleness has been tried, at every scale anybody can build, and the answer never improves. What is going on is that the up-down question and the sideways question do not commute, and the second answer overwrote the first because the two of them were never jointly available to be written.

Where the beam sits between the second magnet and the third is a **superposition**: the sideways-left beam, described in up-down terms, carries two numbers, one for up and one for down. Those numbers are **amplitudes** rather than probabilities. Probabilities are their squares, and the difference is that amplitudes can cancel and probabilities cannot. The dark bands in the prologue's three-month photograph are places where two amplitudes of opposite sign met and left nothing. A beam that was secretly half up atoms and half down atoms, with nobody having looked, produces no cancellation anywhere, and a laboratory tells the two situations apart by bringing the split paths back together and seeing whether the original beam comes back.

## Four multiplications

Heisenberg's tables carry one entry for each pair of states, and a two-outcome question has two states, so its table is four numbers in a square, which is a size a person can multiply on the back of an envelope.

Write U for the up-down question and S for the sideways one, both laid out in up-down terms, with the physical scale stripped off so the entries come out as plain numbers.

$$U = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix} \qquad S = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$$

U is diagonal because up and down are the answers it sorts into, plus one for the atoms that go up and minus one for the ones that go down, with zeros in the corners that would mix them. S has its entries in exactly those corners, since it exchanges up with down.

The rule is the one Heisenberg was forced into. For the entry in a given row and column of the product, run along that row of the left table and down that column of the right one, multiply the pairs and add.

Do US, which is the sideways question asked first, and then do SU. Both products have zeros down the diagonal. US comes out with a one in the top right and a minus one in the bottom left, and SU arrives with those two entries carrying the other sign. Subtract, and nothing cancels: a two in the top right, a minus two in the bottom left. The two orders are exact opposites, entry for entry.

Those four numbers are the middle magnet. Half of the atoms certified up came out down because a subtraction anybody can do in a minute refuses to vanish, and there is nothing in those four entries to be gentler with.

Try a pair that behaves. Keep U, and write a second table with any two numbers you like down the diagonal, five in the top left and two in the bottom right, zeros in the other corners. Multiply it both ways round with U. Both products have five in the top left, minus two in the bottom right and zeros in the corners, so the difference is zero in all four places, whichever two numbers you picked. Tables that sort the atom into the same two answers commute, always, and the radar gun and the camera are a pair of that kind, which is why nobody at Petco Park had to decide which went first.

## A question that squares to itself

Take the smallest well-behaved kind of question: one that admits yes or no and nothing else. Ask it, get an answer, and ask it again straight away. The second answer agrees with the first. Whatever operation the asking performs, doing it twice does what doing it once does.

An operation with that property is a **projector**. It has a name from geometry, because the operation that flattens an arrow onto a plane has the same property: flatten a flat thing again and it does not get flatter. Every yes-or-no question about a physical system is one of these, and every projector is a question. Its opposite is what is left of the identity operation when you remove it, so a question and its opposite add up to certainty.

A complete list of mutually exclusive alternatives is a set of projectors that never overlap and that add to the operation that does nothing at all. The three magnets gave two such lists on the same atom, one for up-down and one for sideways, and neither list is a refinement of the other. The list a question is asked inside of is its **context**. In richer systems than a single atom the same question sits inside many different complete lists.

To turn questions into odds you need one more instrument, and it is arithmetic. Add up the diagonal entries of a square array of numbers. That sum is the **trace**, and has the one property that makes it usable: rewriting the array in a different set of directions leaves the trace alone. The trace of a projector counts how many independent directions the question says yes to.

A **state** is the thing that assigns odds. It is a rule that hands every question a number between zero and one, hands the question with the certain answer a one, and adds up correctly across a complete list of exclusive alternatives. Probability enters this world as an additive weight on yes-or-no questions rather than as a count of how often something happened in a long run of trials. Every such rule can be written as an array of its own, called a **density matrix**, paired with the question through the trace:

$$p = \mathrm{Tr}(\rho P)$$

Here P is the question written as a projector, the Greek letter rho is the density matrix holding the state of whatever you are asking about, and p is the probability of the answer yes.

Run it on the beam between the second magnet and the third. In up-down terms that beam is the array with a half in all four places. Asking whether the atom went down is the array with a one in the bottom right and zeros elsewhere. Multiply the two, add the diagonal entries of the product, and out comes one half, which is the number the third magnet reported. Then set the two corner halves to zero, which turns the beam into one that really is half up atoms and half down atoms, and run the same multiplication. One half again. Every up-down question hands the two arrays the same number, and the corners are the only place their difference lives.

The arrows those arrays act on, with lengths and angles supplied by the state, form a **Hilbert space**. Physics courses usually start there, with the space as the given object and the algebra of questions built afterward out of operations on it. The order runs the other way. What an observer has is a bounded list of questions with a composition rule, and the space is a way of writing that list down: pick a state, and the algebra plus the state manufacture a space of arrows on which the questions act as operations. Different observers with different question lists manufacture different spaces, and they agree wherever their lists overlap, which is the only sense in which the wave function is a public object.

The manufacturing takes a few steps, and none of them begins with a space. The algebra acts on itself, so that every question is both an operation and a thing operated on. The state supplies the geometry: it hands a number to every question, so it hands one to each pair of them, and a rule that gives a number to each pair of things is what lengths and angles are made of. Throw away whatever comes out with zero length, and what is left is a collection of arrows carrying the algebra as operations on it.

Two questions can be answered together exactly when a single complete list of alternatives is fine enough to answer both. That condition is chapter five's, wearing different clothes: compatible measurements are measurements that leave the same thing alone.

## Averaging over what you do not hold

An observer holds a patch. The questions it can put are the ones on its side of the seam, and the full description of the world involves everything past it.

So the state it holds is manufactured by averaging: take the description of the whole thing, and for every question you can actually ask, work out the odds and throw away everything that only a question from the far side could have detected. The operation is called the **partial trace**. Its output is a density matrix, which is why density matrices exist at all. A density matrix is what a description looks like after somebody has been cut out of the picture.

The averaging cannot be run backward. Two different worldwide descriptions can average to the same patch-level density matrix, and once you hold the density matrix there is nothing on your side of the seam that distinguishes them. It is a shadow, and no cleverness turns it back into the object, because the object had features the averaging was built to delete.

A courtroom artist settles the same question every working day, with a few minutes, a partial view, and a drawing to hand in. The drawing has to be consistent with what was seen, and it must not commit to a detail nobody saw, because a detail that happens to be wrong is worse than a vague one. Edwin Jaynes wrote the rule down for physics in the *Physical Review* in 1957: among all the descriptions consistent with what has actually been measured, take the one with the largest entropy, which is the one that adds nothing else. The averaged description a patch holds is exactly that: the maximum-entropy state consistent with the questions it can put and the answers it has recorded, and every feature it lacks is a feature nothing on its side could have established.

## Sorting the mail

Watch what writing a record does to an algebra.

Take a complete list of exclusive alternatives, and take any question at all. Cut the question into the pieces that live inside each alternative, keep those, and delete every piece that straddles two of them. That operation is called **pinching**, a mail room: letters keep their contents, they go into labeled bins, and what is destroyed is every comparison that only made sense between two bins.

Pinching is what writing a record does. The usual name in physics for the process is decoherence. The name describes the deletion rather than explaining it. What comes out of the mail room is an algebra with a property the original did not have: every surviving question commutes with every alternative on the list. The record cannot be disturbed by anything else the observer can ask, which is the whole of what makes it a record. You can read it without changing it, it reads the same the second time, and a neighbor who reads it gets your answer.

The set of things in an algebra that commute with everything in it is the **center**. Something is in there whatever the algebra is. Take the table with ones down the diagonal and zeros in the corners, the operation that does nothing, multiply either magnet by it in either order, and the magnet comes back with its four entries where they were. Doing nothing commutes with everything. Records live in the center.

Everything in the center commutes with everything else in it, so the center is a commutative algebra, and a commutative algebra of this kind has one shape available to it. Strip a table of its off-diagonal entries and what is left is a list of numbers, one for each alternative: a value attached to each way things could have turned out, which is a function on a set of outcomes. A classical description has that in it and nothing else, and the reason the public world wears that shape is a theorem about multiplication rather than a story about what a measurement does to a particle.

Look at what an observer finds when it looks at its own center. Questions with definite answers. Answers that survive being read. Answers that survive being copied to a neighbor, and that come back the same. Order that never matters, because everything in the center commutes with everything. Which is a complete description of the world of tables and chairs and pointer readings, arriving as the commuting part of a structure that is not commuting, rather than as a separate layer of reality sitting underneath.

The center is also the invariant part. Chapter five defined a symmetry as a change that leaves a named quantity alone, and said that a symmetry is the claim that part of your description was never physical. Chapter eleven built the record map, which takes the state of a network and returns the part that survives comparison. The center of an observer's algebra is what every internal change leaves alone. Symmetry, the record map, and the classicality of anything anybody can write down are one structure looked at three ways.

The relation between the two layers is sharper than a metaphor about shadows, and runs one way only. Upward, the record layer sits inside the operator layer perfectly: there is a map from records into questions that is faithful, that preserves sums and products and opposites, that lands region by region exactly inside each region's own algebra, and that gives every record the expectation value it had. Records are questions, of a particularly docile kind.

Downward there is no such map at all. A map that respected products would have to hand a single number to a whole block of questions that do not commute with each other, and no single number does that job, because the block contains pairs whose two orders of composition differ and a number multiplied by a number does not care about order. The best available downward map is an average. It preserves positivity and totals, it leaves records alone, and it fails to respect multiplication: there is an explicit projector in the block whose average of its square differs from the square of its average.

So the classical description is a subalgebra plus an average, the embedding is exact, and the average loses multiplication. There is no projection of the world onto the world of records, and the reason is arithmetic rather than epistemology.

Condition on a record after it is written, which is the operation that makes records usable at all. The state you hold afterward is the one that gives that question probability one and keeps every compatible question's odds unchanged, which is the update rule Gerhart Lüders wrote down in 1951. Apply it twice and the second application does nothing. The operation is a retraction onto certainty: one step takes you into the set of states that agree with the reading, and further reading leaves you where you are. That is why a committed record is a stable classical fact inside a structure whose questions mostly refuse to commute.

Which inverts the usual complaint about quantum mechanics. **Non-commutativity is what creates constraint.** If every question could be answered alongside every other, every observer's answers would be one row of an enormous spreadsheet, any two rows could be filled in independently, and nothing would connect what one observer holds to what another one holds. A world of that kind has no seams in it, nothing to reverse engineer, and no structure to find. The limits are the content.

## The observer, upgraded

Chapter three gave three items, arrived at by asking what the job requires and finding the cheapest object that does it: a bounded view, somewhere to write, and a reading that makes a difference. The third unpacks into three of its own, since a reading that makes a difference needs a way to update, a way to compare, and enough persistence to still be there at the next comparison. Those five become four objects. First, a local algebra of questions, which is the bounded view, bounded because the interface is. Second, a record algebra sitting in the center of that one, which is the somewhere to write, and centrality is what makes the writing stay put. Third, an interface that reads and updates the same way each time, so that comparison gives the same verdict on Tuesday as it gave on Monday. Fourth, enough checkpoint data to fix the odds on every question the observer can put next. The thermostat on the wall has all four: a question list one item long, a record of the setting, a comparison it performs identically every time, and a state that says what it will do next. The five-item version tells you what to look for on a bench. The four-object version can be handed to somebody else, because each of the four is a mathematical object with a name, and two observers can be asked whether they have the same one.

## One generator

Chapter eighteen's screwdriver turned an eighth at a time and drove a screw the whole way into the wood, one small move repeated being what a smooth change is made of. Do the same to an algebra. A change that is reversible and that preserves the structure, meaning it respects sums, products and opposites, is the algebraic version of a rigid motion, and a continuous family of them, starting from the change that does nothing, is a flow: for every real number t, one such change, with the change at t plus s equal to doing the one at t and then the one at s.

On a full block of non-commuting questions, every flow of that kind is generated by a single time-independent quantity: one operator, fixed once at the start, whose repeated infinitesimal application produces every member of the family. Marshall Stone proved the correspondence in the *Annals of Mathematics* in 1932. The generator is what physics calls the Hamiltonian. Chapter eighteen's energy is that generator. A world with continuous reversible change in it has a Hamiltonian whether anybody wrote one down or not.

Ask how the family changes over a single instant rather than over a stretch and what comes back is the equation of motion: applied to the state it is the **Schrödinger equation**, applied to the questions it is Heisenberg's, the same statement written from opposite ends. That is the second of chapter one's three axioms, arriving here as a consequence of continuity and reversibility.

The flow fixes its generator to within one real number. Add a fixed amount to the Hamiltonian everywhere at once, which is to say add a multiple of the operation that does nothing, and that piece commutes with every question in the block and moves none of them, so the family of changes comes back exactly as it was. Two generators producing the same flow differ by that one amount and by nothing else. The gaps between the energies, the ratios of the gaps, the whole spectrum up to a rigid shift, all of it is pinned by the motion. Which is why a laboratory reads differences in energy. Balmer's four wavelengths are differences, and the one number the flow declines to fix cancels out of every one of them.

The mirror statement is stranger and it is about the record layer. Every change of a record algebra that preserves the structure is a permutation of the finitely many labels the records can take. There are finitely many such permutations. A continuous family of them that starts at the identity has nowhere to go, since it would have to jump, so it sits at the identity forever. A purely classical public record cannot flow continuously at all.

Put those two together. The world of pointer readings has no continuous dynamics of its own, in the exact sense that its only structure-preserving motions are discrete relabelings. Whatever flows is flowing in the part that refuses to commute, and the record layer changes by relabeling or it changes not at all.

## A rule that fits every count

Max Born wrote the odds rule down in 1926, in a paper on collisions, putting it in the body of the text as the amplitude and correcting it to the square of the amplitude in a footnote. That footnote is the **Born rule**, one of the three things chapter one said quantum mechanics arrived with and nobody derived. A century of physicists treated the squaring as a postulate, because nobody could get it out of anything more basic.

Andrew Gleason derived it in 1957, in the *Journal of Mathematics and Mechanics*, from almost nothing: any assignment of numbers to yes-or-no questions that lands in the unit interval, gives the certain question one, and adds up across every complete list of exclusive alternatives is the trace against exactly one density matrix. The argument assumes no continuity anywhere and no shape whatever for the rule. There was one hole. Gleason's argument needs three dimensions or more, it fails at two, and two is where a spin lives. Paul Busch closed it in *Physical Review Letters* in 2003 by widening the questions to include blunt ones, which answer yes with a weight instead of sorting the world cleanly in two. The wider version covers every finite dimension including two.

A theorem that short, reaching that far, gets waved through and not believed, since it appears to conjure the whole of quantum probability out of bookkeeping.

Run a two-outcome instrument through its settings and count what comes out: eight settings, a pair of exclusive answers at each, and a frequency for each answer. Then go looking for any odds rule that fits those frequencies, with the trace formula given no privileges in the search.

One exists, and it is explicit. It reproduces every counted frequency exactly, as an identity between whole-number ratios rather than as a fit. Every question it scores lands between zero and one. The question with the certain answer gets one. A question and its opposite get numbers that add to one. And on every sum of exclusive questions that the run's own counts contain, it adds up correctly.

It is the odds rule of no state whatsoever. There is no density matrix at all whose trace against those questions returns those numbers.

Where it breaks is a pair of questions that nobody asked. Take the question this rule scores at 143/512, blunt enough that two copies of it side by side are again a legal question. Adding up demands that the pair score twice 143/512, so double it: 143/256. The rule scores the pair 35/64, which over the same denominator is 140/256. The gap is exactly three parts in two hundred and fifty-six, and that single pair is enough to prove the rule belongs to nobody's state.

That locates the boundary. Additivity restricted to the sums the counted data happens to contain amounts to nothing more than the numbers in each setting adding to one, and a demand that weak leaves the fake alive. The work is done by requiring the adding-up to hold for compatible questions the instrument was never pointed at. Impose that, and the same counts leave exactly one state standing, without a continuity assumption and without any appeal to long runs.

So the counted frequencies do not select the probability rule, and what does select it reaches past every count anybody took. The demand is legitimate because two observers whose patches overlap have to agree about a question they share, and the weight either of them gives it cannot depend on which other questions they chose to ask alongside it, because the other observer made a different choice and the seam has to close. Context independence is the overlap condition, applied to odds.

## Constraint

The commutator says how badly two questions refuse to be asked together, and the uncertainty relation prices that refusal in meters and kilograms and seconds. Questions that refuse nothing collect in the center, and what collects there is the world of tables and pointer readings. Once the refusals hold across every context at once, the odds have one form available to them, and it is the trace formula.

An observer holding an algebra of questions and a state has, in that pair, everything there is to have. Constraint is the part of the pair a commuting world could not supply, and constraint is what makes an undocumented machine readable from a seat inside it.

Constraint has a size. Put two observers far enough apart that no signal can cross between them while they work, hand each of them a pair of questions that do not commute, and compare the two sets of answers afterward. The agreement runs past anything a sealed envelope of instructions could have produced, and then it stops. Where it stops is one fixed number, 2.8284, and no apparatus anybody has built has ever gone past it.