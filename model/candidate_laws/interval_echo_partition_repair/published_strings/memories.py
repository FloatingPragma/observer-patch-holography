"""Finite material entropy and vacuum modular readings on one set of endpoints."""
from dataclasses import dataclass
import math


THETA = math.pi/6
P = 3/8
LOG_ODDS = math.log((1-P)/P)
ENDPOINT_GAP = LOG_ODDS/(2*math.pi)
SENSITIVITY = math.sqrt(3)/4


def _relative_term(x):
    """(1+x) log(1+x) - x near the reference, retaining the quadratic term."""
    if abs(x) < 0.001:
        return x*x*(1/2 + x*(-1/6 + x*(1/12 + x*(-1/20 + x/30))))
    return (1+x)*math.log1p(x)-x


def _sin_minus_angle(x):
    if abs(x) < 0.05:
        square = x*x
        return x*square*(-1/6 + square*(1/120 + square*(-1/5040 + square/362880)))
    return math.sin(x)-x


@dataclass(frozen=True)
class Pair:
    """Rotate the 00,11 plane of a pair with spectrum (5,1,1,1)/8."""
    angle_increment: float

    @classmethod
    def from_weighted_work(cls, work):
        # Calibrations are properties of the endpoint vacuum and its local gate.
        return cls(work/(ENDPOINT_GAP*SENSITIVITY))

    @property
    def probability_change(self):
        d = self.angle_increment
        return 0.5*math.sin(2*THETA+d)*math.sin(d)

    @property
    def cut_relative_entropy(self):
        d = self.probability_change
        return P*_relative_term(d/P) + (1-P)*_relative_term(-d/(1-P))

    @property
    def cut_modular_change(self):
        return LOG_ODDS*self.probability_change

    @property
    def cut_entropy_change(self):
        # The exact difference of the two binary entropies, evaluated stably.
        return self.cut_modular_change-self.cut_relative_entropy

    @property
    def entropy_beyond_linear(self):
        """Cut entropy minus its calibrated first variation, retaining O(d²)."""
        d = self.angle_increment
        # p(d)-p(0) = (sqrt(3)/8) sin(2d) + sin(d)^2/4.
        nonlinear_p = SENSITIVITY*_sin_minus_angle(2*d)/2 + math.sin(d)**2/4
        return LOG_ODDS*nonlinear_p-self.cut_relative_entropy

    @property
    def whole_modular_change(self):
        # Whole-pair entropy is invariant under the rotation. This is also its KL.
        return 0.5*math.log(5)*math.sin(self.angle_increment)**2
