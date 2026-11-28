"""Certify the rank of the joint spectator-vacuum transfer on written fan paths."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, factorial, isqrt, prod
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
FIELD = HERE.parent
ORDER_LIMIT = 72
START_PRIME = 1_000_003


def integer_generator(path):
    """Integer matrix A = 2 H in the fixture's complete physical basis."""
    matrix = (2 * path.h).tocsr(copy=True)
    matrix.sum_duplicates()
    matrix.sort_indices()
    matrix.eliminate_zeros()
    integer = matrix.astype(np.int64)
    np.testing.assert_array_equal(matrix.data, integer.data)
    assert (integer - integer.T).nnz == 0
    row_bound = int(abs(integer).sum(axis=1).max())
    assert row_bound > 0
    return integer, row_bound


def integer_port_columns(path, c):
    """Use the same local code and every spectator vacuum at port zero of fan c."""
    block = path.blocks[c, 0]
    code = np.rint(2 * path.code).astype(np.int64)
    np.testing.assert_array_equal(2 * path.code, code)
    assert set(code.ravel()) == {-1, 1}
    columns = np.zeros((len(path.basis), 2), dtype=np.int64)
    for spectator, words in enumerate(block.words):
        if spectator == c:
            assert len(words) == len(code)
        else:
            neutral = path.field.ring_laplacian(words)
            np.testing.assert_array_equal(neutral @ np.ones(len(words)), 0)
    for i in range(2):
        row = np.ones(1, dtype=np.int64)
        for spectator, words in enumerate(block.words):
            local = (
                code[:, i] if spectator == c else np.ones(len(words), dtype=np.int64)
            )
            row = np.kron(row, local)
        columns[block.rows, i] = row
    normalization = prod(block.shape)
    np.testing.assert_array_equal(
        columns.T @ columns, normalization * np.eye(2, dtype=int)
    )
    np.testing.assert_allclose(
        columns / np.sqrt(normalization), path.source(c), rtol=0, atol=3e-16
    )
    # Coherently projecting each spectator onto its own uniform vacuum is the
    # uniform linear combination of the word-resolved receiving amplitudes.
    listening = path.decode(columns, c)
    projected = listening.sum(axis=0) / np.sqrt(listening.shape[0])
    np.testing.assert_allclose(
        projected, columns.T @ columns / np.sqrt(normalization), rtol=0, atol=1e-14
    )
    return columns, normalization


def determinant_residues(matrix, sender, receiver, row_bound, p, order):
    """Ordinary-time Taylor coefficients with the common n! removed."""
    assert p > order
    assert max(row_bound, np.count_nonzero(receiver[:, 0])) * p < 2**63
    receive_rows = np.flatnonzero(receiver[:, 0])
    receive = receiver[receive_rows].T
    propagated = sender % p
    moments = []
    for k in range(order + 1):
        moments.append(((receive @ propagated[receive_rows]) % p).tolist())
        if k != order:
            propagated = (matrix @ propagated) % p
    return [
        sum(
            comb(n, j) * (heads[0][0] * tails[1][1] - heads[0][1] * tails[1][0])
            for j, (heads, tails) in enumerate(
                zip(moments[: n + 1], reversed(moments[: n + 1]), strict=True)
            )
        )
        % p
        for n in range(order + 1)
    ]


def primes():
    number = START_PRIME
    while True:
        if all(number % divisor for divisor in range(3, isqrt(number) + 1, 2)):
            yield number
        number += 2


def model_fingerprint(matrix, sender, receiver):
    digest = sha256()
    for array in [
        matrix.shape,
        matrix.indptr,
        matrix.indices,
        matrix.data,
        sender,
        receiver,
    ]:
        digest.update(np.asarray(array, dtype="<i8").tobytes())
    return digest.hexdigest()


def exact_series(matrix, sender, receiver, norm_sender, norm_receiver, row_bound):
    """Certify the signed integer coefficients by the Chinese remainder theorem."""
    sequence = primes()
    prime = next(sequence)
    residues = determinant_residues(
        matrix, sender, receiver, row_bound, prime, ORDER_LIMIT
    )
    first_residue = next((n for n, value in enumerate(residues) if value), None)
    cut = ORDER_LIMIT if first_residue is None else first_residue
    bound = 2 * norm_sender * norm_receiver * (2 * row_bound) ** cut
    product = prime
    reconstructed = residues[: cut + 1]
    used = [prime]
    while product <= 2 * bound:
        p = next(sequence)
        residues_here = determinant_residues(
            matrix, sender, receiver, row_bound, p, cut
        )
        modular_inverse = pow(product, -1, p)
        reconstructed = [
            x + product * (((y - x) * modular_inverse) % p)
            for x, y in zip(reconstructed, residues_here, strict=True)
        ]
        product *= p
        used.append(p)
    reconstructed = [x if 2 * x < product else x - product for x in reconstructed]
    assert max(map(abs, reconstructed)) <= bound
    first_exact = next((n for n, value in enumerate(reconstructed) if value), None)
    if first_residue is not None:
        assert first_exact is not None and first_exact <= first_residue
    result = {
        "determinant_residues_at_first_prime": residues,
        "first_nonzero_at_first_prime": first_residue,
        "coefficient_absolute_bound_at_certified_order": str(bound),
        "crt_primes": used,
        "crt_modulus_product": str(product),
        "exact_integer_coefficients": [str(value) for value in reconstructed],
        "first_exact_nonzero_order": first_exact,
    }
    if first_exact is not None:
        n = first_exact
        Taylor = Fraction(
            reconstructed[n], (2**n) * factorial(n) * norm_sender * norm_receiver
        )
        # The higher Taylor terms obey 2 B^m / m! in absolute value.
        # On this punctured interval the total remainder is at most half the
        # magnitude of the exact leading term.
        radius = min(
            Fraction(1, row_bound),
            abs(Taylor) * Fraction(factorial(n + 1), 12 * row_bound ** (n + 1)),
        )
        result.update(
            {
                "leading_taylor_factor_before_minus_i_power": str(Taylor),
                "leading_taylor_minus_i_power": n % 4,
                "certified_open_punctured_interval_radius": str(radius),
            }
        )
    return result


def reflection_certificate(path, matrix, sender, receiver, norm_receiver):
    """Check a terminal-only involution for the same detector and sender."""
    rim = (0, 5, 4, 3, 2, 1)
    ports = path.cells[-1]
    assert [path.board.degree[p] for p in ports[1:]] == [5] * 5
    port_image = {p: ports[rim[j]] for j, p in enumerate(ports)}
    word_image = {}
    for word in {words[-1] for _, words in path.basis}:
        image = 0
        for (a, b), bit in path.field.bits.items():
            head = a if word & bit else b
            reflected = tuple(sorted((rim[a], rim[b])))
            if rim[head] == reflected[0]:
                image += path.field.bits[reflected]
        word_image[word] = image
    permutation = np.array(
        [
            path.ids[port_image.get(port, port), words[:-1] + (word_image[words[-1]],)]
            for port, words in path.basis
        ]
    )
    np.testing.assert_array_equal(permutation[permutation], np.arange(len(path.basis)))
    assert (matrix[permutation][:, permutation] - matrix).nnz == 0
    np.testing.assert_array_equal(sender[permutation], sender)
    projector_image = receiver.T @ receiver[permutation]
    assert not (projector_image % norm_receiver).any()
    parity = projector_image // norm_receiver
    np.testing.assert_array_equal(receiver[permutation], receiver @ parity)
    assert np.array_equal(parity, [[0, 1], [1, 0]]) or np.array_equal(
        parity, [[0, -1], [-1, 0]]
    )
    # The symmetry acts only on the terminal word, so this equation also
    # holds in each spectator word read separately.
    rows = path.blocks[path.size - 1, 0].rows
    terminal_words = len(path.blocks[path.size - 1, 0].words[-1])
    first_row = np.arange(rows.start, rows.stop)
    np.testing.assert_array_equal(
        (permutation[first_row] - rows.start) // terminal_words,
        (first_row - rows.start) // terminal_words,
    )
    return {
        "terminal_local_port_permutation": list(rim),
        "terminal_code_action": parity.tolist(),
        "invariant_sender": True,
        "fixed_spectator_words_at_readout": True,
        "each_accepted_receiver_image_dimension_at_most": 1,
        "determinant_identically_zero": True,
    }


def read_path(path, obstruction=False):
    matrix, row_bound = integer_generator(path)
    sender, norm_sender = integer_port_columns(path, 0)
    receiver, norm_receiver = integer_port_columns(path, path.size - 1)
    result = {
        "native_board_cells": path.board.cells,
        "path_fans": path.size,
        "bridges": path.bridges,
        "ordered_fan_ports": path.cells,
        "fan_degrees": [[int(path.board.degree[p]) for p in fan] for fan in path.cells],
        "prepared_charge_sign": path.sign,
        "configuration_count": matrix.shape[0],
        "integer_matrix_maximum_abs_row_sum": row_bound,
        "joint_source_and_receiver_squared_normalizers": [norm_sender, norm_receiver],
        "configuration_matrix_and_read_columns_sha256": model_fingerprint(
            matrix, sender, receiver
        ),
    }
    if obstruction:
        result["exact_reflection"] = reflection_certificate(
            path, matrix, sender, receiver, norm_receiver
        )
    else:
        result["determinant_series"] = exact_series(
            matrix, sender, receiver, norm_sender, norm_receiver, row_bound
        )
    return result


def run():
    sys.path.insert(0, str(FIELD))
    from fan_path.fixture import FanPath, archived_fixture

    output = []
    with archived_fixture() as archived:
        for sign in (-1, 1):
            output.append(read_path(FanPath(archived, 2, sign)))
            output.extend(
                read_path(FanPath(archived, size, sign, asymmetric=True))
                for size in (3, 4)
            )
            output.extend(
                read_path(FanPath(archived, size, sign), obstruction=True)
                for size in (3, 4)
            )
    return {
        "max_determinant_order_for_initial_prime": ORDER_LIMIT,
        "first_prime": START_PRIME,
        "archive_manifest_sha256": sha256(
            (FIELD / "records.sha256").read_bytes()
        ).hexdigest(),
        "configurations": output,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    content = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x", encoding="utf-8") as file:
            file.write(content)
    elif json.loads(content) != json.loads((HERE / "transfer_series.json").read_text()):
        raise ValueError(
            "Transfer Taylor certificates differ from the archived receipt."
        )
    else:
        print("Reproduced exact transfer rank and terminal reflection certificates.")


if __name__ == "__main__":
    main()
