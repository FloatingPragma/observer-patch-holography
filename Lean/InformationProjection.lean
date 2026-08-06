import InformationProjection.PathGibbs
import InformationProjection.GlobalObjective
import InformationProjection.HistoryLaw

/-!
# Information-projection umbrella root

Conditional finite-history information-projection results: a supplied
exponential path tilt, exact least action for a supplied modal history, and
quantitative plus asymptotic zero-noise concentration.  Constructing the tilt
from OPH source data remains an attachment obligation.
`InformationProjection.GlobalObjective` represents the
four-law package's state and transition receipts as the two marginal
minimizations of one divergence to one product reference, with joint
minimum uniquely at the pinned pair and a control showing the objective
pins its reference.
`InformationProjection.HistoryLaw` derives the exponential
history law: the tilted path law is the unique KL minimizer at fixed
mean action, the committed Gibbs packet's law is the tilt with its
normalization forced, and the named corollary states the history law as
the information projection of the reference onto the mean-action
constraint.
-/
