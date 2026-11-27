"""Live checks on the labelled federations supplied to candidate laws."""

from . import carrier


# The canonical sorted labelled seams, shared by the registered candidate receipts.
FEDERATION_SEAM_SHA256 = (
    '4dadf98983c609e5d305099dc18c19c4f668f9425f12ba75511dbd148c2d60b8',
    '52825a1e46714d3487d4411c20d5ac879750db0ccea164e71f0ad7c7b7576f71',
    'ee40a97c9f9405c00de156d98edabe73c4e04f59b787d80abc2308769c5032f8',
)


def require_three_gateways(board):
    """Count distinct external ports on every owner, including owners with zero external ports."""
    counts = [0] * board.cells
    for port in {port for seam in board.external for port in seam}:
        counts[port // 12] += 1
    for owner, count in enumerate(counts):
        if count != 3:
            raise ValueError(f'{board.name}: owner {owner} has {count} external ports, expected 3')


def registered_seam_digest(board, level):
    require_three_gateways(board)
    digest = carrier.seam_digest(board)
    if digest != FEDERATION_SEAM_SHA256[level]:
        raise ValueError(f'{board.name}: labelled seams differ from registered level {level}')
    return digest


def run(levels=(0, 1)):
    """Return the carrier receipt after checking every owner of each requested federation."""
    result = {'cell': carrier.cell_checks()}
    for level in levels:
        board = carrier.federation(level)
        digest = registered_seam_digest(board, level)
        entry = carrier.federation_checks(board)
        entry['seam_list_sha256'] = digest
        ends = carrier.EMITTER_SEAM[level]
        entry['emitter_seam'] = list(ends)
        entry['emitter_seam_end_has_external_seam'] = any(p in seam for seam in board.external for p in ends)
        result[f'level{level}'] = entry
    return result
