"""Bounded observer histories for the declared coherent source program."""
from dataclasses import dataclass
from .program import Observer, State


@dataclass(frozen=True)
class Session:
    """Finite record horizon with complete unnormalized measurement branches.

    wave_steps counts autonomous W applications. Preparation and apparatus
    have separate labels; every requested operation consumes one history slot.
    The caller supplies max_events as an explicit capacity parameter.
    """
    observer: Observer
    state: State
    max_events: int
    history: tuple = ()
    wave_steps: int = 0

    def __post_init__(self):
        self.observer._validate(self.state)
        if type(self.max_events) is not int or self.max_events < 0:
            raise ValueError('max_events must be a nonnegative integer')
        if len(self.history) > self.max_events:
            raise ValueError('history exceeds the declared capacity')

    def _space(self):
        if len(self.history) >= self.max_events:
            raise OverflowError('coherent observer record capacity exhausted')

    def _append(self, state, record, wave_increment=0):
        self._space()
        entry = tuple(sorted(record.items()))
        return Session(self.observer, state, self.max_events,
                       self.history + (entry,), self.wave_steps + wave_increment)

    def advance(self):
        self._space()
        return self._append(self.observer.step(self.state), {'kind': 'wave_step'}, 1)

    def phase(self, address):
        self._space()
        state, record = self.observer.quarter_turn(self.state, address)
        return self._append(state, record)

    def detect(self, address):
        self._space()
        return tuple(self._append(state, record)
                     for record, state in self.observer.register_measurement(self.state, address))

    def detect_return(self, port, internal=0):
        self._space()
        return tuple(self._append(state, record)
                     for record, state in self.observer.embedded_return_measurement(self.state, port, internal))

    def prepare(self):
        self._space()
        return tuple(self._append(state, record)
                     for record, state in self.observer.preparation_filter(self.state))
