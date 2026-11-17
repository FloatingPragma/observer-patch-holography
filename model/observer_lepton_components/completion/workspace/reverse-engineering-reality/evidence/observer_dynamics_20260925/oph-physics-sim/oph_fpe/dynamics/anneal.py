from __future__ import annotations


def beta_at(schedule: dict, cycle: int, total_cycles: int) -> float:
    kind = schedule.get("kind", "geometric")
    beta_start = float(schedule.get("beta_start", 0.1))
    beta_end = float(schedule.get("beta_end", 5.0))
    if total_cycles <= 1:
        return beta_end
    frac = min(1.0, max(0.0, cycle / (total_cycles - 1)))
    if kind == "linear":
        return beta_start + frac * (beta_end - beta_start)
    if kind == "geometric":
        if beta_start <= 0 or beta_end <= 0:
            raise ValueError("geometric beta schedule requires positive endpoints")
        return beta_start * ((beta_end / beta_start) ** frac)
    if kind == "constant":
        return beta_start
    raise ValueError(f"unknown beta schedule: {kind}")
