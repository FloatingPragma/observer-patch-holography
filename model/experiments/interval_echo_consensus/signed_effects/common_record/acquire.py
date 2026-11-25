"""An addressed Interval Echo preparation with a persistent repair record."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from ..phase_table import acquire as instruments
from ..phase_table import measure as native
from ..sequential import require
from ..tagged_preparation import source

HERE = Path(__file__).resolve().parent
COIN = (144, 145)
TAG = (156, 157)
TRACE = HERE / "records.bin.gz"
RECEIPT = HERE / "results.json"
COLUMNS = ("source_ticket", "repair_before", "context_ticket", "instrument_raw",
           "x_readback_raw", "repair_after")
ABSENT = 255
SEED = 17911168
LAUNCHES = 65536


def sha(data):
    return hashlib.sha256(data).hexdigest()


def calibrations():
    census = json.loads(source.ORIGINAL.read_bytes())["context_web"]["diagonal_context"]
    counts = census["realized_outcome_counts"]
    require(counts == [111, 68] and sum(counts) == census["block_restricted_mass"],
            "historical calibration differs from the theorem operand")
    return counts, [slot["context"] for slot in native.context_table()]


def material():
    """Compute the coin law and the source/tag columns from the full native gate."""
    board = native.federation(0)
    for a, b in (COIN, TAG, (0, 1)):
        require(b in board.neighbours[a], "missing addressed native contact")
    pair = [6, 2]  # flattened ternary words (2,0) and (0,2)
    half = source.native_gate("contact", 0.5)
    coin = abs(half[np.ix_(pair, pair)])**2
    require(np.max(abs(coin - 0.5)) < 3e-13, "coin depends on its prior word")
    full = source.native_gate("contact", 1)
    idle = np.eye(9)
    rails, tags = [], []
    for gate in (idle, full):
        rail = gate[[7, 5], 7]  # excess on rail zero or rail one
        tag = abs(gate[pair, pair[0]])**2
        require(abs(np.linalg.norm(rail) - 1) < 3e-13 and abs(tag.sum() - 1) < 3e-13,
                "source or tag escaped its material code")
        rails.append(rail)
        tags.append(tag)
    require(np.max(abs(np.abs(rails)**2 - np.eye(2))) < 3e-13 and
            np.max(abs(np.array(tags) - np.eye(2))) < 3e-13, "paired write failed")
    return coin, rails, np.array(tags)


def reader_kernels(contexts, rails):
    """Joint distribution of the sharp pointer and a subsequent real X read."""
    # Verify the physical copy and address-word action used to factor the branches.
    native.address.local_check()
    instruments.cleared_record_check("sharp")
    x_word = instruments.x_analyzer_word()
    analyzer = native.source_unitary(tuple(x_word), 1)
    copy = native.record_word("sharp")
    address_copy = native.address.pulse_word(native.address.schedule())
    kernels = np.zeros((len(contexts), 2, 17, 17))
    protected = set(COIN + TAG + (native.address.REFERENCE,))
    operation_hashes = {}
    for c, context in enumerate(contexts):
        front = native.context_word(context) + native.dilation_word("sharp") + copy
        matrices = native.raw_branch_maps(context, "sharp", sign=1)
        effect = matrices[0].conj().T @ matrices[0]
        require(np.linalg.norm(effect - native.projector(context)) < 4e-12,
                "native sharp head disagrees with archived context")
        operation = []
        for raw, branch in matrices.items():
            words = front + native.post_word(context, "sharp", raw)
            words += instruments.reset_word(raw) + x_word + address_copy
            require(all(a not in protected and b not in protected for _, a, b, _ in words),
                    "the apparatus touched the coin, repair record or balancing charge")
            operation.append(instruments.serialized(words))
            for rail, state in enumerate(rails):
                kernels[c, rail, raw, :12] = abs(analyzer @ branch @ state)**2
        operation_hashes[context] = sha(json.dumps(operation, separators=(",", ":")).encode())
    require(np.max(abs(kernels.sum(axis=(2, 3)) - 1)) < 8e-12, "lost public branch")
    phase = native.raw_branch_maps("phase", "sharp", 1)[0]
    orientation = [float(np.linalg.norm(phase @ instruments.input_vector(name, 1))**2)
                   for name in ("plus_i", "minus_i")]
    require(np.max(abs(np.array(orientation) - [1, 0])) < 4e-12,
            "phase head has wrong orientation")
    return kernels, [round(p, 12) for p in orientation], operation_hashes


def fair_word(rng, coin, initial, width):
    previous = initial
    word = np.zeros(len(initial), dtype=np.uint8)
    for _ in range(width):
        previous = (rng.random(len(initial)) < coin[1, previous]).astype(np.uint8)
        word = 2 * word + previous
    return word


def acquire(shots, seed):
    counts, contexts = calibrations()
    coin, rails, tags = material()
    kernels, orientation, operation_hashes = reader_kernels(contexts, rails)
    rng = np.random.Generator(np.random.PCG64(seed))
    records = np.full((shots, len(COLUMNS)), ABSENT, dtype=np.uint8)
    tickets = fair_word(rng, coin, np.zeros(shots, dtype=np.uint8), 8)
    records[:, 0] = tickets
    admitted = np.flatnonzero(tickets < sum(counts))
    rail = (tickets[admitted] >= counts[0]).astype(np.uint8)
    # The controller writes both materials before taking any selector readings.
    tag = (rng.random(len(rail)) < tags[rail, 1]).astype(np.uint8)
    records[admitted, 1] = tag
    selected = fair_word(rng, coin, tickets[admitted] & 1, 3)
    records[admitted, 2] = selected
    for c in range(len(contexts)):
        for b in range(2):
            rows = admitted[(selected == c) & (rail == b)]
            law = kernels[c, b].ravel()
            # Both stages are retained. The joint Kraus law is equivalent to
            # drawing the first record followed by its conditional readback.
            encoded = rng.choice(len(law), size=len(rows), p=law / law.sum())
            records[rows, 3], records[rows, 4] = np.divmod(encoded, 17)
    # Every intervening operation acts as identity on the tag tensor factor.
    # Its second sharp occupation read acts on the collapsed orientation.
    records[admitted, 5] = (rng.random(len(tag)) < np.eye(2)[tag, 1]).astype(np.uint8)
    return records, kernels, orientation, operation_hashes


def summarize(records, kernels, orientation, operation_hashes, seed):
    """Read the public stream and check joint custody against the physical branches."""
    counts, contexts = calibrations()
    ticket, tag, selected, raw, readback, final_tag = records.T
    admitted = ticket < sum(counts)
    require(np.all(records[~admitted, 1:] == ABSENT), "refusal carries instrument data")
    expected_tag = (ticket[admitted] >= counts[0]).astype(np.uint8)
    require(np.array_equal(tag[admitted], expected_tag) and
            np.array_equal(final_tag[admitted], expected_tag), "ticket and tag disagree")
    require(np.all(selected[admitted] < len(contexts)) and
            np.all(raw[admitted] < 17) and np.all(readback[admitted] < 17),
            "malformed public record")
    broken = int(np.count_nonzero(np.roll(tag[admitted], 1) != expected_tag))
    banks = []
    for c, name in enumerate(contexts):
        for b in range(2):
            inside = admitted & (selected == c) & (tag == b)
            n = int(inside.sum())
            require(n > 0, "missing selected preparation")
            event = 17 * raw[inside].astype(int) + readback[inside]
            observed = np.bincount(event, minlength=17**2).reshape(17, 17)
            require(not np.any(observed[kernels[c, b] < 1e-18]), "forbidden public event")
            radius = np.sqrt(np.log(2 * 16 * 17**2 / 0.01) / (2 * n))
            require(np.max(abs(observed / n - kernels[c, b])) <= radius,
                    "public frequencies outside simultaneous 99 percent bound")
            bins = [[int(i), int(j), int(observed[i, j])] for i, j in zip(*np.nonzero(observed))]
            banks.append({"context": name, "repair": b, "launches": n, "nonzero_joint_counts": bins})
    inputs = list(instruments.source_hashes())
    inputs += [str(p.relative_to(native.ROOT)) for p in
               (source.ORIGINAL, Path(source.__file__), Path(__file__),
                Path(native.fixed_pointer.__file__), Path(native.sequential.__file__))]
    return {"word_sha256_by_context_including_all_raw_branches": operation_hashes,
            "candidate": "Interval Echo", "seed_pcg64": seed, "launches": len(records),
            "columns_uint8": list(COLUMNS), "refusal_sentinel": ABSENT,
            "source_ticket_bits": 8, "selector_bits_after_admission": 3,
            "field_bit_order": "most significant bit first", "coin_ports": COIN, "repair_ports": TAG,
            "contexts_in_bit_order": contexts, "admitted_ticket_fibres": counts,
            "refusal_tickets": 256 - sum(counts), "refused_launches": int((~admitted).sum()),
            "admitted_tag_counts": np.bincount(tag[admitted], minlength=2).tolist(),
            "shifted_tag_disagreements_preserving_marginal": broken,
            "native_phase_yes_for_plus_i_and_minus_i": orientation,
            "omitted_joint_bins_are_zero": True, "banks": banks,
            "uncompressed_trace_sha256": sha(records.tobytes()),
            "source_sha256": {k: sha((native.ROOT / k).read_bytes()) for k in sorted(set(inputs))}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    fresh, kernels, orientation, operation_hashes = acquire(LAUNCHES, SEED)
    if args.write:
        TRACE.write_bytes(gzip.compress(fresh.tobytes(), mtime=0))
    stored = np.frombuffer(gzip.decompress(TRACE.read_bytes()), dtype=np.uint8).reshape(-1, len(COLUMNS))
    require(np.array_equal(stored, fresh), "public trace differs from source replay")
    receipt = summarize(stored, kernels, orientation, operation_hashes, SEED)
    receipt["run_sha256"] = sha(json.dumps(receipt, separators=(",", ":"), sort_keys=True).encode())
    encoded = json.dumps(receipt, indent=2) + "\n"
    if args.write:
        RECEIPT.write_text(encoded)
    require(RECEIPT.read_text() == encoded, "receipt differs from independently reread records")
    print(f"{len(stored)} launches, {receipt['refused_launches']} refusals, "
          f"{len(receipt['banks'])} context/repair banks, verified")


if __name__ == "__main__":
    main()
