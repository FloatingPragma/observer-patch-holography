# The candidate laws

OPH builds its Universe as a program: a large network of very simple
observers, each holding a few records and comparing them with its
neighbours. The theory around that program is finished. What is missing is
the native repair law, the one local rule every observer follows. A search
over finite families of laws, each tested on what its simulated observers
record, has narrowed that law to three candidates. This page states them,
what has been proven about them and what is open. The
[simulator](../simulator/README.md) runs them, and its replay app shows every
experiment.

## What the three laws share

| Part | What it does | Where it is specified |
| --- | --- | --- |
| Reflected routing | One constituent lives on eight slots per address, paired along three spatial axes and one private cyclic register. Each update reflects the equal-slot part of the eight amplitudes and sends every slot one address along its axis into the partner slot. | [simulator spec, sections 1 to 3](../simulator/SPEC.md#3-update) |
| Contact | Where two constituents share an address in matched slots, a contact coin with phase φ acts on the eight matched states, and both then step together: a bound pair. One constituent carries two rails, with contact phases π/16 and π/16 + π/2048; their relative phase is the pair's clock hand. | [simulator spec, contact](../simulator/SPEC.md#3-update) |
| Capture | Where the two constituents meet head-on, a reversible vertex turns the collision, at the angle π/4, into a matched pair plus one marker that leaves on the junction's own wire. | [simulator spec, capture](../simulator/SPEC.md#3-update) |
| Maturation | A newborn pair carries an age counter. While young it takes a partial hop $e^{-i\theta}[\cos\theta\,I + i\sin\theta\,F]$ between the two ends of each routed edge, with $\theta = \tfrac{\pi}{2} f\big((j+1)/M\big)$ at age $j$; mature, it takes the full hop. The counter cycles after 8192 updates. | [simulator law registry](../simulator/law.py) |

## The three laws

They differ only in the opening schedule $f$ and its duration $M$.

| Law | Opening | Correct direction among late crossings | Separated after the marker's second return |
| --- | --- | ---: | ---: |
| Cosine 32 | $f(u) = (1 - \cos\pi u)/2$ over 32 updates | 95.24% | 1.06% |
| Cosine 128 | the same profile over 128 updates | 96.56% | 0.53% |
| Linear 128 | $f(u) = u$ over 128 updates | 95.97% | 6.67% |

Cosine 128 has the highest directional fraction and the smallest breakup on
the tested fixtures. No measured record selects one of the three for
nature: they share the mass shell of the mature pair and pass the same clock
tests.

## Proven

These are exact statements about the specified construction. The simulator
checks each one against the literal update.

- **A relativistic field from one local rule.** For one constituent released
  in the equal-slot state, $\langle e|W^n e\rangle = \cos(nE)$ with
  $\cos E$ the mean of the four cosines of the wave vector. Under refinement
  the field obeys the free Klein–Gordon equation
  $\varphi_{tt} = \tfrac14\nabla^2\varphi - \tfrac{Q^2}{4}\varphi$ in three
  space dimensions, with signal speed $c = 1/2$ and error of order $h^2$ at
  fixed chart time. [Continuum run](../simulator/README.md#what-the-app-shows)
- **Pairs stay together.** The matched coincident states of the contact law
  form an invariant subspace, and so does its complement: a bound pair never
  separates under the contact law, and without the capture vertex separated
  constituents never bind.
- **Moving clocks run slow.** The mature pair's band is
  $\cos E = a(k)\cos(\varphi/2)$, so at long wavelengths
  $E^2 \simeq (|k|^2 + \varphi^2)/4$ and the splitting between the rails falls
  with motion by the relativistic factor. At every finite momentum
  $1 - |v|^2/c^2 - r^2 \ge 0$: the law's clocks run slightly slower than the
  smooth formula.
- **The matter-wave relation.** On a plane wave the mature step obeys
  $\cos\omega = \cos m\,(2w_{\rm private} + 2w_{\rm spatial}\sum_i\cos k_i)$
  with $m = \varphi/2$, so
  $\omega^2 = m^2 + 0.24919629\,|k|^2 + O(|k|^4)$ in radar and home-clock
  units, against the scalar standard $c^2 = 0.24935023$.
- **A strict front.** Every update moves amplitude one address at most, and
  each forward hop beyond the smooth cone costs a factor 1/16 in probability.
- **Reversibility.** The capture vertex, the partial hops and the age counter
  have exact inverses.

## Measured on simulated observers

Each result comes from records that the simulated observers hold, with the
forecast sealed before the records were read. The
[journal](../docs/research/laws/JOURNAL.md) gives the protocols and numbers.

| Observable | Result | Journal |
| --- | --- | --- |
| Time dilation of a moving clock | Free clocks matched the flat-metric forecast within 0.035 on every long-wave route; two differently built clocks agreed on five journeys; the bound pair's clock kept proper time, and the finite law predicted its residual within 1.32 counting standard errors on six sealed flights. | [steps 1, 3, 4](../docs/research/laws/JOURNAL.md#step-1-clocks-and-spacetime-from-observer-records) |
| Three space dimensions | Nine stations predicted 120 held-out field values with RMS error 0.0011; dropping any one spatial derivative raised it to 0.31–0.32. | [step 1](../docs/research/laws/JOURNAL.md#step-1-clocks-and-spacetime-from-observer-records) |
| A message on the free wave | A compact free-wave launch carried no fresh bit ahead of the smooth cone; the detector probability before the early deadline stayed below $10^{-34}$. | [step 2](../docs/research/laws/JOURNAL.md#step-2-a-race-between-the-addressed-relay-and-the-free-wave) |
| Formation | Two separately launched constituents formed a persistent pair in 8,145 of 16,384 aimed launches; zero with the vertex disabled. | [step 5](../docs/research/laws/JOURNAL.md#step-5-a-pair-formed-in-an-encounter) |
| Direction | Ordinary incoming trains sent 95–97% of the crossing pairs the right way. | [step 6](../docs/research/laws/JOURNAL.md#step-6-three-laws-that-form-steer-and-keep-a-clock) |
| A steady clock | The whole receiving array read the formed pair's clock at one rate to 0.21% over six durations. | [step 6](../docs/research/laws/JOURNAL.md#step-6-three-laws-that-form-steer-and-keep-a-clock) |
| The matter-wave relation | All 14 held-out directions at periods 64, 128 and 256 obeyed $\omega^2 = m^2 + c^2\vert k\vert^2$ within 3.5% of the kinetic part; a law with half the spatial weight failed every one. | [step 7](../docs/research/laws/JOURNAL.md#step-7-the-matter-wave-relation-of-formed-matter) |

## Open

- **An individual clock.** The steady rate belongs to the whole array; widening
  the incoming train changed it by 31%. One stationary component has to keep
  one rate at different widths before it serves as a portable clock.
- **Frames.** A travelling laboratory that ranges by echoes and measures a
  persistent rod in its own frame.
- **Access.** An inhabitant who can trigger an addressed relay sends messages
  ahead of the smooth cone, and an equipped receiver decoded constituent
  signals early at short range. Which operations inhabitants can build is
  part of the law.
- **Formation in an ordinary history.** The tested births use prepared
  incoming beams; one qualified marker costs about $2^{53}$ incident pairs on
  the largest cube.
- **Everything else.** Charge, the gauge forces, gravity, cosmic expansion
  and the identification of particles are work in progress on these laws.

Supplied by the makers and listed with every result: the carrier and its
sizes, the contact phases, the capture angle and the marker wires, the
opening schedule and the age register, the preparations, the detectors and
the archived scalar standard.

## Run the laws

```sh
python -m model                      # the candidate laws' exactness checks
python -m simulator laws             # every law of the family and its standing
python -m simulator run all          # every experiment, and the offline replay app
```

## Theory packages

Independent finite theorems and their tests sit beside the candidate laws,
each with its own README: [support_refinement](support_refinement/README.md),
[matched_einstein_limit](matched_einstein_limit/README.md),
[defect_particle_classes](defect_particle_classes/README.md),
[carrier_selection_reduction](carrier_selection_reduction/README.md),
[finite_instrument_calibration](finite_instrument_calibration/README.md),
[maxwell_local_selection](maxwell_local_selection/README.md) and
[maxwell_weak_field](maxwell_weak_field/README.md). They keep their own
hypotheses and carriers.
