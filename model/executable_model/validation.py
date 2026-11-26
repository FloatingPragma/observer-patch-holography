"""Live checks on the labelled federations supplied to candidate laws."""

from . import carrier


def require_three_gateways(board):
    """Count distinct external ports on every owner, including owners with zero external ports."""
    counts = [0] * board.cells
    for port in {port for seam in board.external for port in seam}:
        counts[port // 12] += 1
    for owner, count in enumerate(counts):
        if count != 3:
            raise ValueError(f'{board.name}: owner {owner} has {count} external ports, expected 3')


def run(levels=(0, 1)):
    """Return the carrier receipt after checking every owner of each requested federation."""
    result = {'cell': carrier.cell_checks()}
    for level in levels:
        board = carrier.federation(level)
        require_three_gateways(board)
        entry = carrier.federation_checks(board)
        entry['seam_list_sha256'] = carrier.seam_digest(board)
        ends = carrier.EMITTER_SEAM[level]
        entry['emitter_seam'] = list(ends)
        entry['emitter_seam_end_has_external_seam'] = any(p in seam for seam in board.external for p in ends)
        result[f'level{level}'] = entry
    return result
