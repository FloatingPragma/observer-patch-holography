# 30. What Is a Lagrangian, and Why Is It the Last Thing?

In 1788 the widow Desaint published a book on mechanics in Paris that contained no pictures, and its author put that on the first page as though it were the attraction. Joseph-Louis Lagrange had come down from Berlin the year before. The Académie sent the manuscript to a committee of Laplace, Cousin, Legendre and Condorcet, who approved it, and Legendre read the proofs.

No diagrams will be found in this work, the preface opens. The methods set out in it, Lagrange goes on, call for neither constructions nor geometrical or mechanical reasoning, but only algebraic operations subject to a regular and uniform procedure.

Mechanics until then was drawn. Newton's *Principia* argues in figures: an arc, a tangent, a shaded triangle, a limit taken as two points slide together, and the reader follows by looking. For a century after 1687 doing mechanics meant getting the picture right first, and the pictures were different for a pendulum, a spinning top and a planet, so each of them was a separate accomplishment. What Lagrange had found was that somebody who cannot draw at all obtains every one of those results by writing down a single function of the system and then turning a handle, and that the handle turns the same way in all three cases.

The line working physicists write down when somebody asks what the world is made of, the one that fits on a mug, is an object of the same kind: a few terms wide, with a subject folded into each one. Two questions come with an object like that. What sort of thing is it, that a whole subject folds into one line at all? And what does the folding cost, since a book of mechanics with the figures taken out has thrown something away, whatever else it has gained.

## Proof of a supreme being

On 15 April 1744, before the Académie Royale des Sciences in Paris, Pierre-Louis Moreau de Maupertuis read a paper called "Accord de différentes lois de la nature qui avaient jusqu'ici paru incompatibles", agreement between different laws of nature which had until then seemed incompatible. The laws were the ones for how light bends and how it bounces. He proposed that a single quantity, which he named the action, comes out as small as it can be whenever anything in nature changes.

Two years later in Berlin, running Frederick's academy, he published "Les Loix du mouvement et du repos déduites d'un principe métaphysique", and there the principle acquired a second job. Maupertuis presented it as a demonstration of a supreme being, and ranked it above the older proofs drawn from the beauty and the order of the world, on the ground that an equation is harder to argue with than a wonder.

The inference is the obvious one: a rule that picks the smallest total out of every total available looks like an economy, and an economy has a treasurer.

Maupertuis had earned the standing to say such things by going out and getting a number. In 1736 he led an expedition to the Torne valley in Swedish Lapland to measure the length of a degree of the meridian close to the Arctic Circle. The party came home in August 1737 with a degree longer than the one measured in France, which meant the Earth is flattened at the poles, which meant Newton was right and the Cassinis were wrong. Robert Levrac-Tournières painted him afterward in a fur hat and reindeer skins with one hand pressed down on a globe, squashing it out of round. Voltaire, an admirer at that stage, addressed him as the flattener of the world and of Cassini.

The rest went badly. In 1751 Samuel König claimed that Leibniz had had the least-action principle first and produced a letter to prove it, and the Berlin academy, with Maupertuis presiding, declared the letter a forgery. Voltaire took König's side and wrote the *Diatribe du docteur Akakia, médecin du Pape*, in which a papal physician works through the published proposals of a self-styled academy president, most of them quoted from what Maupertuis had actually printed. Frederick had the pamphlet burned by the public executioner on Christmas Eve 1752, with the king and Voltaire both standing at the fire. Maupertuis left Berlin and died in Basel in 1759.

His principle survived him exactly as he stated it, and so did his reading of it. A rule that scores a whole path and then selects one looks as though the path was compared against paths that never happened.

## A ball, described twice

Throw a ball across a room. It leaves your hand at one place and arrives at another a second later.

Here is the ordinary description. At each moment the ball has a position and a velocity. Gravity pulls down, which changes the velocity a little, which changes the position a little, and repeating that as finely as you like traces out the arc. The description looks forward one step at a time and needs nothing except what is true where the ball is.

Here is the other one. Draw every arc you like between those two points taking the same second: the true one, one that goes twice as high and comes down in a rush, one that dawdles halfway across and then sprints, one shaped like a staircase. Score each of them with a single number. The arc the ball flies is picked out by what its score does when you nudge it.

The scoring rule has two parts, and the future appears in neither.

The first part is one function. At any instant, two things about the ball matter: where it is, and how fast it is going. Feed both into a single function and it hands back a number. For a thrown ball, that function is the part that grows with speed less the part that grows with height.

Both parts are countable for a ball you could actually throw. The height part is the **potential energy**, what the ball's position owes: the mass times the strength of gravity times the height. Gravity at the Earth's surface pulls at 9.81 meters per second per second, and rounding that to ten costs two percent and buys arithmetic you can do standing up, so a one-kilogram ball held a meter off the floor owes ten joules. The speed part is the **kinetic energy**, which is what the motion holds: half the mass times the square of the speed. That same ball moving at two meters per second holds two joules, and at four meters per second it holds eight, because the square is where the growth is. The function is the second of those less the first, joule for joule.

The second part is a running total. The number the function hands back is per second, so run along the arc from the first moment to the last and keep a running total of it, the way a taxi meter accumulates while the wheels turn. One arc in, one number out. Every arc on the list gets a total this way, including the ridiculous ones.

Selection uses slopes. Stand on a hillside and step one pace in any direction: your height changes by an amount proportional to the size of the step. Stand on the floor of the valley and step one pace: your height barely changes, because the ground under you is level, and the change you do get shrinks with the square of the step rather than with the step. The bottom of a valley is where the slope is zero.

The arc a thrown ball flies is where the slope is zero. Push the middle of it up by a hair, holding both ends where they are, and the running total does not change. Push it down by a hair and the running total does not change. Every neighboring arc scores the same to first order, and the score is flat there in the way the valley floor is flat.

Take the crudest version of that and check it by hand. Three moments and two hops: a one-kilogram ball starts on the floor, ends on the floor one second later, and the only freedom is how high it is at the halfway mark. Raising the midpoint costs in the speed part, because the ball has to climb further and fall further in the same second. It pays in the height part, because the ball spends the middle of the flight higher up. Each hop lasts half a second and there are two of them, so anything holding steady through both contributes its own value times one second to the running total, and the halves cancel out of the arithmetic before it starts.

A peak of one meter. The ball climbs a meter in half a second, so it travels at two meters per second and holds two joules of motion, and its height averages half a meter across each hop, which owes five joules. Two less five is minus three.

A peak of one and a half meters. Three meters per second, four and a half joules of motion, average height three quarters of a meter, seven and a half joules owed. Minus three again.

A peak of one and a quarter. Two and a half meters per second, three and an eighth joules of motion, average height five eighths of a meter, six and a quarter joules owed. Minus three and an eighth.

The arc that overshoots and the arc that undershoots score the same number, and the one between them scores below both, by an eighth of a joule-second in either direction. A single line says why, with h standing for the height of the halfway mark in meters and the total coming out in joule-seconds:

$$\text{score} \;=\; 2h^{2} \;-\; 5h$$

The cost grows with the square of the midpoint height and the payment grows in proportion to it, so the cost overtakes. At one midpoint height and no other, raising the peak by a hair costs exactly what it pays, and there the total has slope zero. That is the height the rule selects, one and a quarter meters, and it is where a ball airborne for a second peaks: it leaves the floor at five meters per second, sheds ten of those every second, and stops climbing after half of one. Three moments is the coarsest account of a flight with any freedom left in it, and it lands on the height exactly. Cut the second into four hops instead of two and the corners of the selected shape land on a parabola, and into forty and they land on the same one.

Look at what had to be in your hand to do any of that. All three sums used the landing. Fixing the speed on each hop took the ball being back on the floor at one second, and without the finish there was nothing to add up. The ball has the start. It leaves your hand with a place and a speed and no information whatever about where it comes down, and the rule says it flies the arc whose total sits at the bottom of a parabola whose shape depends on where it comes down. Something in that description knows the ending before it happens, and it is not the ball. Maupertuis put a treasurer in the gap, which at least has the merit of being a mechanism.

The one function of where you are and how fast you are going is the **Lagrangian**. Its running total along a path is the **action**, which is the word chapter twelve used for a single number scored over a whole arrangement, and this is the same object scored along a history instead of across a wiring. The rule that selects the path where wiggling changes nothing is **stationary action**.

Everything about a system is in its Lagrangian. Hand somebody that function and every quantity the system has gets its equation of motion, by the same handle-turning in each case: take the function's sensitivity to how fast the quantity is changing, ask how that sensitivity changes as the run proceeds, and set it equal to the function's sensitivity to the quantity itself. One line goes in and the entire behavior comes out. **A Lagrangian is a compression**: the thing you keep when you throw away the histories, chosen so that the histories can be regenerated from it.

Watch it come apart once. The Lagrangian of a thrown ball has two pieces, the speed piece and the height piece. Its sensitivity to how fast the ball is going comes from the speed piece and is the mass times the velocity, which is the ball's momentum. Its sensitivity to where the ball is comes from the height piece and is the mass times the strength of gravity, pointing down. Set the rate of change of the first equal to the second and out drops the statement that the ball's momentum gains mass times gravity downward every second, which was Newton's second law and which nobody wrote into the function. Two pieces of arithmetic and the handle, and the same handle applied to a spinning top gives its equations too.

## Göttingen, 1918

The second thing a Lagrangian compresses is the reason anything is conserved.

Emmy Noether arrived at Göttingen in the spring of 1915 at the invitation of Hilbert and Klein, and worked there for years with no salary and no post. The faculty blocked her habilitation. Her lecture courses appeared in the catalogue under Hilbert's name, with herself listed as the person who would actually be giving them. When she had the result that carries her name she could not read it to the Göttingen Scientific Society, not being a member of it; Klein read it out for her on 26 July 1918, and it went into the society's proceedings under the title "Invariante Variationsprobleme", pages 235 to 257.

In April 1933 the Law for the Restoration of the Professional Civil Service removed her from the university. She went to Bryn Mawr College in Pennsylvania, and once a week to Princeton to lecture. In the spring of 1935 she had surgery for a tumor, the operation went well, and four days afterward she developed a fever and lost consciousness. She died on 14 April 1935, at fifty-three.

Chapter five gave the correspondence and chapter eighteen used it: every continuous symmetry of the laws comes with exactly one conserved quantity, determined by the symmetry, and momentum and energy are what you get from the laws not caring where and not caring when. The machine that turns the one into the other is a recipe with three steps.

Take the Lagrangian. Find a change you can make to it continuously, by an amount you can dial down to nothing, which leaves its value alone: slide the whole system a hair sideways, turn it by a small angle, advance every clock by the same interval, add the same small amount to a phase everywhere. Each of those is a direction of change with a size attached to it.

Then, for every quantity the system has, take the Lagrangian's sensitivity to how fast that quantity is changing. That sensitivity is the momentum belonging to the quantity, and it falls out of the function rather than being defined off to one side. Multiply each sensitivity by how far the symmetry moves its own quantity, and add the products over the whole list. Advancing the clock takes one more move, since shifting the clock shifts the running total too: subtract the function's own value times the interval.

The number you get does not change along any history the Lagrangian selects. One symmetry in, one number out, by a recipe that does not care which symmetry it was fed.

Slide sideways and the number is momentum. Turn and it is angular momentum. Advance the clocks and it is energy. Add the same phase everywhere and it is electric charge, arriving by a second road. Chapter twelve reached it from the other end, by showing that a force whose charge leaks stops being a force at all, and the recipe here starts from the same phase and hands back the quantity you can total up.

The same statement holds on a discrete chain with no limit taken anywhere. A run is a list of states with one step between each neighboring pair, and a symmetry is a change that carries each state to another while leaving each step's contribution to the running total alone. Take the contribution's sensitivity to the change across one step and it comes out the same across every step of the run. Segment by segment, the number holds, and there is no continuum anywhere in the argument.

Bookkeeping in this world was never a separate law added on top of the dynamics. It is a property of the line, extracted by a fixed recipe, and a world has as many conservation laws as its line has continuous symmetries.

## A product of small numbers

None of that touches Maupertuis's reading. The recipe is exact, and a rule that scores whole paths and selects the flat one continues to look like a comparison, which needs the losers to be somewhere in order to lose.

Start again from the smallest assumption a world can be given.

Something is in some configuration. There is a rule saying what it does next, and the rule has three properties. It is local, so the next step depends only on the configuration it is standing in. It is memoryless, so the rule does not consult how that configuration was arrived at. And it is strictly positive, so nothing is flatly forbidden: every configuration it could move to has some chance, however small.

That is the whole of the input. There is no action in it, no Lagrangian, no scoring, and nothing that looks at a history.

Ask about histories anyway, since a history is a thing that has happened and can be asked about. Take a list of configurations, one after another, from a start to a finish. What are the odds of that particular list? Each step carries odds of its own and the steps do not consult each other, so the odds of the list are the odds of the first step times the odds of the second times the odds of the third, all the way to the end.

A product of a thousand numbers each below one is a number nobody can read. Chapter sixteen introduced the operation that fixes it: a logarithm converts multiplying into adding. Take the logarithm of each step's odds, flip the sign so the contributions come out positive, and add them along the history instead of multiplying. Written out, with the odds of the first step called p one, the odds of the second p two, and the history running to n steps:

$$-\log\left(p_{1}\,p_{2}\cdots p_{n}\right) \;=\; \left(-\log p_{1}\right) \;+\; \left(-\log p_{2}\right) \;+\; \cdots \;+\; \left(-\log p_{n}\right)$$

On the left is one number belonging to a whole history. Check it on the smallest case: two steps, each a coin flip, so the history has odds of one in four, which is two of chapter sixteen's yes-or-no questions, one from each step. On the right are n numbers, each belonging to a single step and worked out from the two configurations that step ran between, none of them consulting any other step, and the equals sign in the middle is the line of algebra that takes the foresight out.

Look at what that produces. It is a running total along a path, one contribution per step, and each contribution depends on where its step started and where it ended and on nothing else anywhere in the history. That is the shape of the thing the thrown ball was scored with, arrived at from a rule that has never heard of it.

A history with a large total took improbable steps. A history with the smallest total is the most probable history there is. And a history whose total does not move when you jiggle one of its interior configurations, holding the ends, is sitting on the valley floor in that direction, which is what it means to be locally most likely. **Least action and maximum probability are two readouts of one thing.**

Run the total backwards. Give each history a weight that shrinks as its total grows, one over the exponential of the total, and out comes the odds of that history again, because exponentiating undoes the logarithm that built the total. Nothing has to be divided by anything to bring the weights to one across all the histories, since the odds added to one before any logarithm was taken.

Then suppose somebody hands you a different running total, and a different overall scale to multiply it by, and their combination also reproduces the same odds. Their scale times their total equals the total above plus a constant, at every history, and there is no other freedom anywhere in the problem. You may slide the zero, and you may change the units. Everything else is forced.

A world whose step rule is local, memoryless and strictly positive has an action, then, whether or not anybody writes one down.

Nothing in that argument looks ahead, because the step rule has nothing to look ahead with. The endpoint turns up in the statement of the principle because the statement is about a completed history, and a completed history has a finish in it the same way a receipt has a total printed at the bottom. A machine repairing locally, with no view past its own neighbors, produces a selection that looks global, because that is what taking the logarithm of a product of local factors does.

Maupertuis had found something real. The quantity a world comes out smallest on is improbability, and a chain of local steps runs that tally by itself, without a manager and without a comparison ever being performed.

## Plateau's wire frames

Joseph Plateau lost the sight of both eyes by the end of 1843, at forty-two. Ghent promoted him to full professor a few months later, and for the next four decades his wife Fanny, and afterwards his son and son-in-law, did the seeing for him. They dipped wire frames of every shape into soapy water mixed with glycerine, drew them out, and described aloud exactly what surface had formed. In 1873 he published nearly a thousand pages on the statics of liquids held together by molecular forces alone, every word of it dictated. Jean Taylor proved his rules about how the films meet each other in 1976, a hundred and three years afterward, as consequences of minimizing area.

Dip a bent wire in soap solution and the film stretched across it is the surface of least area that wire will hold. It looks like the answer to a question about all possible surfaces. The film has no way of finding out what any other surface would have cost. Each patch of it pulls on the patches around it, the pulling is settled locally, patch against neighbor, and when no patch is being pulled harder one way than another the film stops moving. Where it stops has the least area available. Nothing in the film has ever evaluated an area.

The settling is also invisible, which is most of why the film looks clever. By the time the frame is steady in your hand the film has finished, so nobody has ever watched a soap film choose. What Plateau's family described to him across forty years was the finished surface, every time.

## The line

Physics as it is normally taught reaches a Lagrangian in the first week and spends the rest of the time taking it apart. The spacetime is there before the pen touches the page, the group is written down because the detectors reported it, the fields are listed because they were found, and the numbers are entered because laboratories measured them.

The line has four terms in it, written schematically, each standing for a block that opens out further.

$$\mathcal{L} \;=\; \mathcal{L}_{\text{carry}} \;+\; \mathcal{L}_{\text{matter}} \;+\; \mathcal{L}_{\text{couple}} \;+\; \mathcal{L}_{\text{select}}$$

Four things added together, and the chapters behind this one built every symbol in them.

The script L on the left is a Lagrangian in the sense taught above, with one change of scale. Where the thrown ball had a position and a speed, this one takes a field's value at a place and how fast that value is changing, and its running total runs over all of space and all of time rather than along a single arc.

The first term is the carriers. Chapter twenty-five split the twelve dials of the wiring into a central one, a block of three and a block of eight, with nothing else arithmetically available, and this term holds the energy stored in what those dials do around a loop. Chapter twelve's loop readings are the field strength, and squaring them and adding over the arrangement is the term.

The second term is matter, carrying chapter twenty-six's fifteen states, three families over. Inside it sits a comparison between a reading here and a reading next door, taken after the reading next door has been run through the dictionary sitting on the seam between them, which is chapter twelve's connection written with derivatives instead of exchange rates. That comparison is the entire content of the phrase "the force acts on matter".

The third term holds the couplings: the channels through which the vector-selecting field touches matter, taken from the same scan over subsets that produced the fifteen states, with strengths that chapter twenty-nine's balance condition fixes for the charged leptons, the third of them inside a window seventy-two electron-volts wide.

The fourth term belongs to the vector-selecting field itself, chapter twenty-seven's single spin-zero type, which carries no polarization index because a magnitude has no direction for a rotation to mix, and whose vacuum vector has a size of 246.22 billion electron-volts.

Alongside the four sits gravity's term, chapter twenty-four's curvature, whose coupling is Newton's constant, which is the area of one cell of a screen, 2.61 times ten to the minus seventy square meters. And beside that the constant term, the one every light ray reads as exactly zero, which is why the cosmological constant sits in the field equation as a quantity no local calculation was ever going to reach.

The prologue promised that by the last third of the book you would be able to read that line and say where every symbol in it came from, and it made no other promise. Chapter one gave the standing list of what a framework begins by assuming: a spacetime, a collection of fields living on it, a symmetry group chosen to match what was observed, and a page of numbers measured in laboratories. The four terms above are that list with none of the four assumed. The spacetime is chapter fifteen's completion of a tally of comparisons. The fields are chapter twenty-seven's five places a pattern can sit. The group is one, three and eight because twelve dials come apart into a center of one and a remainder of eleven, and eleven splits as three plus eight and in no other way. Where the wiring settles a number the chapters behind this one give it, and where the wiring stops, chapter twenty-nine printed the exact angle at which it stops.

Those are the same marks somebody else writes on line one, and the difference shows up in a single place. Written first, the line has three families in it because three families were found. Written here, it has three because an oriented face has three corners.

## What a compression costs

A compression is allowed to be lossy, and the useful ones generally are.

This one presumes a continuum. The running total runs over space and time, the fields are smooth, and the handle-turning that produces the equations of motion needs slopes to exist at every point of every path. Chapter fifteen was careful about that limit for a reason that bites here: an observer holds six whole numbers and sits at some finite number of passes, which is a scatter of isolated points with gaps between them that never close at any stage. Smoothness stands to that scatter as pi stands to 3.1416. The line is written where the tally is heading, and no observer's tally has ever been there.

The second loss is sharper, and it can be written down exactly.

Take the line and add a term chosen so that it comes out zero at every configuration the world realizes. On a chain, the realized configurations sit at the corners, and a term that vanishes at every corner adds nothing to the score of any realized history. So every realized history keeps its old score. Swap an interior configuration for another the world can realize, holding the ends, and the answer comes back what it was: the same run is the most probable run, and it sits at the bottom against every single-step swap.

The momentum changes. Momentum is the line's sensitivity to how fast a quantity is changing, and a term pinned to zero at every corner tilts as it passes through, in a direction no realized step ever moves. So the two lines assign different momenta to the same motion, and therefore different energies, since energy is what you get by trading the speeds for the momenta.

The freedom is one number. The added term bends with speed by any positive amount you like, and its bend with position is then forced, the two summing to one particular number the step odds themselves set. Under that constraint the most probable run stays at the bottom when an interior configuration slides continuously rather than jumping corner to corner, and the momentum stays a strictly increasing, invertible function of speed. Each setting is a different line carrying a different energy. Nothing in any history this world produces distinguishes them from each other. The histories are the same histories, and they were the only evidence there was.

So the line is the best compression anybody has written down, and a compression is not the thing compressed. Lagrange knew what he had taken out, which is why he mentioned it in the first sentence.

Hand the line to somebody with a computer and enough patience and out comes the arc of a thrown ball to as many decimal places as the ball has, the spectrum of hydrogen, the lifetime of a muon and the rate at which the sun burns. Then put the ball in their hand and push it, and watch it go. The line scores the path it takes. It says nothing whatever about why a pushed thing moves.