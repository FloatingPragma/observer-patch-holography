# The OPH law

OPH builds its Universe as a program: a large network of very simple
observers, each holding a few records and comparing them with its
neighbours. Every one of them follows the same local rule, the OPH law. This
page states the law, what is proven about it and what its simulated
observers measure. The [simulator](../simulator/README.md) runs it, and its
replay app shows every experiment.

## The law

| Part | What it does |
| --- | --- |
| Reflected routing | One constituent lives on eight slots per address, paired along three spatial axes and one private cyclic register. Each update reflects the equal-slot part of the eight amplitudes and sends every slot one address along its axis into the partner slot. |
| Contact | Where two constituents share an address in matched slots, a contact coin acts on the eight matched states, and both then step together: a bound pair. One constituent carries two rails, with contact phases π/16 and π/16 + π/2048; their relative phase is the pair's clock. |
| Capture | Where two constituents meet head-on, a reversible vertex turns the collision, at the angle π/4, into a matched pair plus one marker that leaves on the junction's own wire. |
| Opening | A newborn pair carries an age counter of 8192 states. At age $j < 128$ it takes the partial hop $e^{-i\theta}[\cos\theta\,I + i\sin\theta\,F]$ between the two ends of each routed edge, with $\theta = \tfrac{\pi}{2}\,\tfrac{1 - \cos(\pi (j+1)/128)}{2}$; from age 128 on it takes the full hop $F$. |

The [simulator specification](../simulator/SPEC.md) states every part
exactly, and [law.py](../simulator/law.py) executes it.

## Proven

These are exact statements about the law. The simulator checks each one
against the literal update.

- **A relativistic field from one local rule.** For one constituent released
  in the equal-slot state, $\langle e|W^n e\rangle = \cos(nE)$ with
  $\cos E$ the mean of the four cosines of the wave vector. Under refinement
  the field obeys the free Klein–Gordon equation
  $\varphi_{tt} = \tfrac14\nabla^2\varphi - \tfrac{Q^2}{4}\varphi$ in three
  space dimensions, with signal speed $c = 1/2$ and error of order $h^2$ at
  fixed chart time.
- **Pairs stay together.** The matched coincident states of the contact law
  form an invariant subspace, and so does its complement: a bound pair never
  separates under the contact law, and two constituents bind only through
  the capture vertex.
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
forecast sealed before the records were read.

| Observable | Result |
| --- | --- |
| Time dilation of a moving clock | Free clocks match the flat-metric forecast within 0.035 on every long-wave route, two differently built clocks agree on five journeys, and the bound pair's clock keeps proper time; the finite law predicts its residual within 1.32 counting standard errors on six sealed flights. |
| Three space dimensions | Nine stations predict 120 held-out field values with RMS error 0.0011; dropping any one spatial derivative raises it to 0.31–0.32. |
| A message on the free wave | A compact free-wave launch carries no fresh bit ahead of the smooth cone; the detector probability before the early deadline stays below $10^{-34}$. |
| Formation | Two separately launched constituents form a persistent pair through the capture vertex; with the vertex disabled none forms. |
| Persistence | A bound pair's circulation is unchanged by every update at every age, and a returning marker cannot release it: zero separations in five cohorts of 4,194,304 trials, where independent transport of the two constituents separates 75%. |
| Direction | Ordinary incoming trains send 96.56% of the crossing pairs the right way, and 0.53% of the pairs separate after the marker's second return on a closed wire. |
| A steady clock | The whole receiving array reads the formed pair's clock at 4.395–4.400 × 10⁻⁴ radians per home tick, steady to 0.21% over six durations. |
| The matter-wave relation | All 14 held-out directions at periods 64, 128 and 256 obey $\omega^2 = m^2 + c^2\vert k\vert^2$ within 3.5% of the kinetic part; a law with half the spatial weight fails every one. |

The [replay app](../simulator/README.md#what-the-app-shows) shows these
experiments.

## Work in progress

- **A portable clock.** The steady rate belongs to the whole receiving array;
  widening the incoming train changes it by 31%. One stationary component
  with one rate at every width gives an individual clock.
- **Frames.** A travelling laboratory that ranges by echoes and measures a
  persistent rod in its own frame.
- **Access.** Which operations inhabitants can build: an inhabitant who can
  trigger an addressed relay sends messages ahead of the smooth cone.
- **Formation in an ordinary history.** The tested births use prepared
  incoming beams.
- **Charge, forces and gravity.** Charge, the gauge forces, gravity, cosmic
  expansion and the identification of particles on this process.

Supplied with every result: the carrier and its sizes, the preparations,
the detectors and the scalar reference standard.

## Run the law

```sh
python -m model                      # the law's exactness checks
python -m simulator laws             # the law and the controls of its experiments
python -m simulator run all          # every experiment, and the offline replay app
```

## Theory packages

Independent finite theorems and their tests sit beside the law, each with
its own README: [support_refinement](support_refinement/README.md),
[matched_einstein_limit](matched_einstein_limit/README.md),
[defect_particle_classes](defect_particle_classes/README.md),
[carrier_selection_reduction](carrier_selection_reduction/README.md),
[finite_instrument_calibration](finite_instrument_calibration/README.md),
[maxwell_local_selection](maxwell_local_selection/README.md) and
[maxwell_weak_field](maxwell_weak_field/README.md). They keep their own
hypotheses and carriers.
