# 22. Why Does Writing Something Down Cost Heat?

In 1824 a French military engineer of twenty-eight paid a Paris publisher named Bachelier to print six hundred copies of a book about steam engines, and then waited for a scientific community that did not read it. Sadi Carnot's *Reflections on the Motive Power of Fire* runs to a hundred and eighteen pages and asks the question every engine builder in Europe wanted answered, which is how much work you can get out of a fire.

His answer has no pistons in it. The most work any engine whatever can deliver, per unit of heat it swallows, is fixed by two numbers: the temperature of the thing it takes the heat from and the temperature of the thing it dumps the remainder into.

$$\eta = 1 - \frac{T_{\text{cold}}}{T_{\text{hot}}}$$

The Greek letter names the efficiency, the fraction of the heat that comes out as work. The two temperatures are the hot source and the cold sink, both counted from absolute zero. Brass or iron, steam or air, a better valve, a cleverer linkage, a bigger boiler: the ceiling sits where those two temperatures put it and no engineering moves it.

Carnot got that while believing heat was a weightless fluid called caloric, which fell from hot to cold the way water falls down a millrace and did as much work on the way. Caloric does not exist. Every heat engine built since has come in under the number he got out of it.

Cholera reached Paris in the spring of 1832 and killed something near twenty thousand people in the city that year. Carnot caught it on 24 August and was dead within the day, at thirty-six. His belongings were burned as a precaution against contagion, and nearly all of his papers went with them. What survived was a set of working notes his brother Hippolyte published in 1878, in which the fluid had been quietly dropped.

Look at what the surviving formula says. An engine's ceiling depends on the temperatures and on nothing else, which is a peculiar thing for a fact about brass and steam to be. Chapter twenty-one reached the same quantity from the other end, out of an algebra and a state, with no furnace anywhere in the argument. So what is a temperature a count of?

The everyday version of the same question sits on your desk. A drive with no moving parts gets warm while it is being written to. There is no friction in it to blame. No amount of engineering has ever brought the warmth to zero, in any machine anybody has built.

## The doorkeeper

On 11 December 1867 James Clerk Maxwell wrote to Peter Guthrie Tait describing a way to break the second law of thermodynamics using a very small employee.

Take a box of gas at one uniform temperature, divide it with a wall, and put a door in the wall. Temperature is an average: some molecules in there are moving fast and some slow. Station a being at the door with good eyes and a light touch. When a fast molecule approaches from the left, it opens the door and lets it through to the right. When a slow one approaches from the right, it opens the door and lets it through to the left. The door is weightless and frictionless, so opening and closing it costs nothing.

Wait. The right side heats up and the left side cools down, out of a box that started uniform. Then run Carnot's engine off the difference, and run it again, forever, in a box that nobody has done any work on.

Maxwell called his employee a finite being and later grumbled that it was really more of a valve. William Thomson named it a demon in *Nature* in 1874, meaning the word in its Greek sense, a spirit working quietly in the background.

The demon survived every attempt to kill it for sixty-two years, because every attempt went looking in the wrong place. People examined the door, the hinges, the light the demon needed to see by. Leo Szilard moved the search in 1929, in a paper in *Zeitschrift für Physik* whose title says where he was looking: on the decrease of entropy in a thermodynamic system by the intervention of intelligent beings. He stripped the gas down to a single molecule in a box, with a partition dropped in the middle. The demon looks, records which side the molecule is on, and uses that one recorded fact to let the molecule push the partition outward and lift a weight. The work per cycle comes out as the temperature times the logarithm of two, out of a box at one temperature, from one bit of recorded information.

The books balance only if the bit costs at least what the bit earns. Rolf Landauer, at IBM, worked out in 1961 which operation carries the charge, in a paper called "Irreversibility and heat generation in the computing process". Looking is not the expensive part. Forgetting is. An operation that takes two distinguishable states of a memory to one state has thrown away the distinction. The distinction has to leave the device, because the only exit a memory has is heat. The floor he computed for that heat depends on nothing but the temperature outside the device.

$$E \geq k_B T \ln 2$$

E is the energy released as heat when one bit is erased, T is the temperature of the surroundings the erasing device sits in, and the constant converts a temperature in degrees into an energy in joules. The logarithm of two is chapter sixteen's one question, the one you ask to tell two possibilities apart.

Put a room temperature into it, three hundred kelvin, and the bill for erasing one bit comes to 2.87 divided by ten to the twenty-first, in joules. A single photon of green light at a wavelength of 550 nanometers carries 3.61 divided by ten to the nineteenth, in joules, which is enough to pay off about a hundred and twenty-six bits.

Charles Bennett finished the demon in a review of the thermodynamics of computation in
1982. The demon's memory is finite. It can run a few cycles free, filling that memory with
records of which side each molecule was on, and then it is full, and to take another
measurement it has to clear space. Clearing the space costs the temperature times the
logarithm of two per bit, which is exactly what the molecule paid it. The demon was never
an engine. It was a machine borrowing against its own memory. The loan came due at the
same rate every time.

Antoine Bérut and colleagues put the number on a bench in 2012 and published it in *Nature*. A silica bead about two micrometers across, held in an optical trap shaped into two wells, one well for zero and one for one. Push the bead into a known well and you have erased a bit. Measure the heat that leaves. As the erasure cycles were slowed down, the mean heat per erasure came down onto that value and did not go below it.

Chapter three listed three requirements for anything that counts as an observer, said that writing a record costs energy, and left the price standing as a bill. That is the price.

## Four laws, numbered out of order

Thermodynamics has four laws. They were not discovered in the order they are taught.

Start with the one everybody uses without noticing. Nobody has ever established that a bath and a bowl of soup are at the same temperature by putting the soup in the bath. You put a thermometer in the soup, wait, read it, then put the same thermometer in the bath, wait, read it, and compare two numbers taken minutes apart with an instrument that has been in contact with each and never with both at once.

That procedure is a bet on transitivity. If the thermometer settles with the soup, and the same thermometer settles with the bath, then the soup and the bath will settle with each other. The bet is a law: two systems each in equilibrium with a third are in equilibrium with each other. Ralph Fowler and Edward Guggenheim wrote in 1939 that it "could with advantage be known as the zeroth law of thermodynamics", the other three having been numbered decades before anybody noticed this one was doing work underneath them.

Transitivity is what makes temperature a number. A relation that is transitive, and symmetric, and holds between a thing and itself, sorts everything in the world into classes with nothing left over, and classes in a row can be labeled by numbers on a scale. Without transitivity you could say that this is hotter than that, pair by pair, and there would be no scale to write either of them on, no degrees, no thermometer, and no sentence of the form "the water is at 47 degrees".

The first law gets taught as bookkeeping, in two terms. The energy of a thing changes when heat flows into it and when work is done on it, and heat and work are the two channels. The exact accounting has three terms rather than two. Work moves the energies of the available arrangements while the odds over them stay put, which is what a piston does. Heat changes the odds while the energies stay put, which is what contact with a flame does. Do both at once and there is a third term, the product of the two changes, which vanishes when you change one thing at a time and does not vanish otherwise. The textbook version is exact for a protocol that moves the piston and then opens the flame, and approximate for one that does both together.

The second law arrived in chapter sixteen, in a laboratory in Canberra with a latex bead in a dish of water. Take the odds a settled arrangement assigns, call them the reference, and ask how far your own description sits from it. That distance falls under any process and never climbs.

The third law is the one popular accounts skip. Walther Nernst put it to the Göttingen Academy in December 1905: as a system is cooled toward absolute zero, its entropy approaches a fixed floor. The floor is set by how many arrangements the system has at its lowest energy. One such arrangement and the floor is zero. Many, and it is the logarithm of how many, and sits there down to the last fraction of a degree.

Ordinary ice does this. Each oxygen atom in ice has four hydrogens near it, two close and two further off, and there are many ways to satisfy that rule across a whole crystal without ever violating it locally, which is chapter six's triangle of frustrated neighbors wearing a different mineral. Linus Pauling counted the arrangements in 1935 and got a floor of 0.805 calories per degree per mole. William Giauque and J. W. Stout measured the heat capacity of ice down to fifteen degrees above absolute zero and reported 0.82, give or take 0.05, in a paper of 1936. That is 3.4 joules per kelvin per mole of entropy that stays in a block of ice at the bottom of the temperature scale.

The other half of the third law is a prohibition. You cannot reach absolute zero in a finite number of steps, however good your refrigerator, because each step takes a fraction of what is left rather than a fixed amount. A sodium gas at the Massachusetts Institute of Technology was cooled in 2003 to 450 picokelvin, give or take 80, which is four hundred and fifty trillionths of a degree above a floor nobody reaches.

## One rule, read twice

Four laws, arrived at over eighty years by people working on engines, gases, chemical affinities and cold. They come out of the repair rule of chapter nine, read twice.

Read it first as a rule about descriptions. Jaynes's instruction from chapter nineteen applies unchanged: among all the descriptions consistent with what has actually been measured, take the one with the largest entropy, the one that adds nothing else. Apply that with one quantity held fixed, the average energy, and the answer is forced. The odds fall off exponentially with energy, cheap arrangements common and expensive ones rare, with a single multiplier in the exponent setting how fast the fall is. Chapter twenty-one ran into that exponential from the far side, in a two-outcome system whose flow came out as the energy divided by the temperature.

That multiplier is one over the temperature. The temperature arrives as a conversion rather than as a property of anything. Fix the average energy of a description and the least assumption you can make about the rest of it carries exactly one number. That number says how many joules the world charges for one yes-or-no question. A temperature is a count of joules per bit. At three hundred kelvin the count is 2.87 divided by ten to the twenty-first of them.

The zeroth law follows from the arithmetic of that exponential rather than from a habit of laboratories. Two systems that have settled with each other carry one multiplier between them, and equality of numbers is transitive, so a thermometer that matches the soup and matches the bath has established that the soup matches the bath. There is one condition on the thermometer: the instrument must have two arrangements of different energy. A thermometer whose two states cost the same reads the same number in a furnace and in liquid helium, because the temperature shows up only in the ratio of odds between two different energies. An instrument with one energy has no such ratio to offer.

Read the rule the second way, as a rule about transitions. A repair step is handed some things that neighbors have settled between them and some things nobody has settled. What should it do? The same instruction applies: change nothing you do not have to, and assume nothing you have not been given. Leave every settled quantity exactly as it stands, and redraw everything else from the reference, inside the set of arrangements that agree with what is settled.

That instruction picks out one map. Applying it twice does nothing that applying it once did not do, because the second application finds the same settled facts and the same reference. It leaves the reference where it is, since redrawing from the reference cannot move it. And it preserves the average of every quantity that can be read off the settled part, which is the first law in the form chapter nine derived it: a repair that changed a total would have to know the total. The heat and the work come off the same map. Shifting the energies of the arrangements without touching the settled odds is the work channel, redrawing the odds at fixed energies is the heat channel, and the cross term is what a step that does both at once picks up on the way through.

Three of the four laws drop out of that map without further argument. Distance to the reference falls at every step and never climbs, which is the second law, the data-processing fact of chapter sixteen with the process named: repair is error correction, and running your description through an error-correcting step cannot carry it further from the truth.

Every schedule of repairs on chapter nine's twelve observers lands on the same arrangement, and things that land on one fixed point together are in equilibrium with each other, which is transitivity again, reached from the transition side this time instead of from the exponential, and with the thermometer taken out of it. And one repair step leaves every arrangement with some weight on it: the step redraws from a reference that gives everything a positive share, so a possibility that had a share before the step has one after it. Finitely many steps extinguish nothing, which is why the floor of the third law is approached and never arrived at.

Landauer's bill is a corollary of the same map rather than a separate discovery. A step whose entropy falls by some amount expels at least that amount divided by the same multiplier. Take the amount to be one bit and out comes the temperature times the logarithm of two, in joules, per bit.

## A one-way loop that obeys the law

Chapter sixteen's fluctuation relations lean on one extra condition: measured against the reference, the weight a step carries from one state to a second matches the weight it carries back. The condition has a name, **detailed balance**.

The second law does not need it. The cheapest way to see that is a machine with three states and a permanent circulation in it.

Label three states one, two and three. From each state, the rule is: stay where you are with probability one half, or step clockwise to the next state with probability one half. There is no counterclockwise move at all. Nothing in this machine can go from state two to state one, ever, except the long way through state three, which takes four steps on average.

Check what it does to the flat description that puts a third of the weight on each state. State one receives half of its own third, from staying, and half of state three's third, from the clockwise step. That is a third. So is every other state, by the same count. The flat description does not move, which makes it the reference.

Detailed balance fails outright. The traffic from state one to state two is a third of a half, one sixth per step. The traffic from state two to state one is zero, because that move does not exist. There is a current going around the loop that never dies down.

Watch the second law hold anyway. Start with all the weight on state one, which is a description sitting the logarithm of three away from the flat one, 1.585 bits. One step spreads it into halves on states one and two, and the distance falls to 0.585 bits. A second step gives a quarter, a half and a quarter, at 0.085 bits. A third step brings it to 0.024. The descent is exact, it never reverses, and it happens inside a machine with a one-way street in it.

Stationarity is what the second law asks for: that the reference be left alone. Detailed balance is a stronger and separate demand. The extra content it carries buys the detailed fluctuation relations of chapter sixteen and Lars Onsager's paired transport coefficients. A world could have the second law with none of those and this loop is what it would look like.

## The price of settling

Two observers hold readings of one shared quantity and the readings disagree. Before any repair, the description of that pair has two live arrangements in it and no fact distinguishing them, which is one bit of entropy, the logarithm of two. Repair them. Both readings go to one value. The description has one arrangement in it and its entropy is zero.

Settling a disagreement lowers entropy. Every comparison in the world that comes out resolved, and every record written into the public account behind it, takes a description with two possibilities in it and hands back one.

Nothing is violated by that. The second law is about the distance between a description and the reference. That distance is the quantity that never climbs. The bare entropy of a description is a different number. A settling process drives it down by one bit each time it settles a bit. A rule that only ever lowers entropy carries none of the second law's content, which is why the three-state loop had to be scored on its distance from the flat description rather than on its entropy.

The bit had to go somewhere. It leaves the two observers as heat, at least the temperature times the logarithm of two of it, at whatever temperature the pair's own state sets, which is the multiplier in its own exponential. Chapter sixteen showed that a commit discards which of seven arrangements the seam used to hold, and priced it at 2.807 bits of the world's ability to say what it used to hold. The same commit has an energy attached: 2.807 bits of erasure at three hundred kelvin costs 8.06 divided by ten to the twenty-first joules, expelled as heat, per seam, per commit, forever, everywhere in a universe that is making a public record out of private disagreements.

## Five instruments, one quantity

The same number has been read off five different instruments by five sets of people who did not know they were measuring the same thing.

An engine reads it as the heat it took in divided by the temperature it took it at. That ratio is what Carnot's ceiling is a statement about. A gas reads it as Boltzmann's count of the arrangements a coarse description fails to distinguish. A memory reads it as Shannon's count of the questions it takes to identify what is stored, and Landauer's price per question. A horizon reads it as an area, in square meters, which is the strangest of the five. And a public record reads it as what a settled disagreement threw away, which is the one this world is built out of.

Every record ever written has been paid for in heat, at the temperature of whatever wrote it. The heat leaves through the edge of the region holding the record, which puts a bill on the edge rather than on the volume inside it. Ask how many records a region of space can hold, and the answer comes back in square meters.