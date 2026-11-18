"""Sequential native/observer oracle, independent of the compiled gates.

The inherited native event is the source component. The observer schedule is
evaluated directly from its declared blank preparation, without calling edge,
prefixes, preparation, or public_record from the compiled program.
"""

from math import sqrt

from model.native_memory_clock.native import native_event

from .program import Data


def sequential(config, source, address):
    """Complete two-event semantics for one prepared first-address branch."""
    first = native_event(source, address, config.background, config.law)
    observed = (first.target >> address) & 1
    pointer_a = (2 | observed) if config.write_a else 0
    courier = pointer_a if config.communication else 0
    pointer_b = courier
    second_address = 6 if config.feedback and pointer_b == 3 else 0
    second = native_event(first.target, second_address, config.background, config.law)
    return {
        "initial_source": source,
        "first_address": address,
        "first_source": first.target,
        "first_old_bit": first.old_bit,
        "observed_bit": observed,
        "pointer_a": pointer_a,
        "courier": courier,
        "pointer_b": pointer_b,
        "second_address": second_address,
        "second_old_bit": second.old_bit,
        "final_source": second.target,
    }


def expected_prefixes(config, source, address):
    """Every actual register after each prefix, including unread buffers."""
    row = sequential(config, source, address)
    first, pa, courier = row["first_source"], row["pointer_a"], row["courier"]
    old = row["first_old_bit"]
    return (
        Data((source, 0, 0, 0, 0), (address, 0, 0, 0), (0, 0, 0),
             0, 0, (0, 0), 0),
        Data((0, first, 0, 0, 0), (0, address, 0, 0), (0, 0, 0),
             0, 0, (old, 0), 0),
        Data((0, 0, first, 0, 0), (0, 0, address, 0), (0, courier, 0),
             pa, 0, (old, 0), 0),
        Data((0, 0, 0, first, 0), (0, 0, 0, address), (0, 0, courier),
             pa, 0, (old, 0), 0),
        Data((0, 0, 0, 0, row["final_source"]), (0, 0, 0, address),
             (0, 0, courier), pa, row["pointer_b"],
             (old, row["second_old_bit"]), row["second_address"]),
    )


def record(data):
    """Complete basis-record meaning, without a clock-position inference."""
    return (data.pointer_a >> 1, data.pointer_a & 1,
            data.pointer_b >> 1, data.pointer_b & 1, data.second_address,
            data.buffers, data.addresses, data.couriers, data.old_bits)


def expected_records(config, source, address):
    """All five complete register records for independent receipt replay."""
    return tuple(record(data) for data in expected_prefixes(config, source, address))


def archive_label(data):
    """Read all nonsource registers while retaining the quantum source.

    This complete diagnostic instrument also reads unused address registers.
    It is a separate terminal setting from reading every source buffer.
    """
    return (data.pointer_a, data.pointer_b, data.second_address,
            data.addresses, data.couriers, data.old_bits)


def branch_columns(config, sources, stage=4):
    """Sparse Kraus columns for the complete source-retaining instrument.

    The input columns are the supplied distinct source labels. Each outcome
    maps (output buffer tuple, input source) to its amplitude. Every address
    branch is retained, with no successful-arrival normalization.
    """
    if type(stage) is not int or not 0 <= stage <= 4:
        raise ValueError("stage must be an integer in [0, 4]")
    if not sources or len(set(sources)) != len(sources):
        raise ValueError("distinct nonempty source labels required")
    columns = {}
    for source in sources:
        for address in range(12):
            data = expected_prefixes(config, source, address)[stage]
            key = (data.buffers, source)
            block = columns.setdefault(archive_label(data), {})
            block[key] = block.get(key, 0.0) + 1 / sqrt(12)
    return columns
