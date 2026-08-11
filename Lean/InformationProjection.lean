import InformationProjection.PathGibbs
import InformationProjection.GlobalObjective
import InformationProjection.HistoryLaw
import InformationProjection.SourceHistoryPacket
import InformationProjection.LogTransitionAction
import InformationProjection.ReferenceNormalForm

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
`InformationProjection.SourceHistoryPacket` mirrors the B7
source-attached history-law packet as exact literals extracted from one
retained locally hash-pinned run: both path laws are strictly positive and
normalized, the reference law is the Markov product of the committed
mixing chain with the cleared-denominator identity kernel-decided, the
tilt at multiplier zero is the reference exactly, the `HistoryLaw`
interface is inhabited at the packet literals at every multiplier, and
an intermediate-value receipt locates a multiplier matching the
declared empirical mean action inside the payload's rational bracket.
The multiplier selection principle from source data remains the open
core of B7.
`InformationProjection.LogTransitionAction` proves that every
strictly positive Markov path law is the Gibbs tilt of the step-uniform
reference by the log-transition action, that any action reproducing the
law lies in the additive-constant and multiplier-rescale gauge orbit of
that action, and instantiates both statements with kernel-decided
receipts on the committed source chain, including the negative control
that the repair-count action reproduces the chain law at no multiplier.
-/
