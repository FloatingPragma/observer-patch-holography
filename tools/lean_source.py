"""Position-preserving views of Lean comments and source for repository checks."""
from __future__ import annotations

import re
from collections.abc import Iterator


# Match a quote character literal before treating its payload as a string.
# Lean identifiers can end in a prime, so a lone prime does not start a literal.
_OPENER = re.compile(r'''/-|--|'"'|\br(\#*)"|"''')
_BLOCK = re.compile(r"/-|-/")
_NON_NEWLINE = re.compile(r"[^\n]+")


def _hidden_spans(text: str) -> Iterator[tuple[int, int, bool]]:
    """Yield spans for comments (True) and literals (False) in source order."""
    i = 0
    while match := _OPENER.search(text, i):
        token = match[0]
        start, i = match.span()
        comment = token in ("--", "/-")
        if token == "--":
            end = text.find("\n", i)
            i = len(text) if end < 0 else end
        elif token == "/-":
            depth = 1
            while depth and (boundary := _BLOCK.search(text, i)):
                depth += 1 if boundary[0] == "/-" else -1
                i = boundary.end()
            if depth:
                i = len(text)
        elif match[1] is not None:
            # Raw strings close at a quote followed by the opening hash count.
            delimiter = '"' + match[1]
            end = text.find(delimiter, i)
            i = len(text) if end < 0 else end + len(delimiter)
        elif token == '"':
            while i < len(text):
                char = text[i]
                i += 1
                if char == '"':
                    break
                if char == "\\" and i < len(text):
                    i += 1
        yield start, i, comment


def _blank(text: str) -> str:
    return _NON_NEWLINE.sub(lambda match: " " * len(match[0]), text)


def code_view(text: str) -> str:
    """Mask comments and quoted literals, preserving token boundaries and positions."""
    pieces = []
    pos = 0
    for start, end, _ in _hidden_spans(text):
        pieces.extend((text[pos:start], _blank(text[start:end])))
        pos = end
    pieces.append(text[pos:])
    return "".join(pieces)


def comment_view(text: str) -> str:
    """Keep comments only, with their original offsets and line breaks."""
    pieces = []
    pos = 0
    for start, end, comment in _hidden_spans(text):
        if comment:
            pieces.extend((_blank(text[pos:start]), text[start:end]))
            pos = end
    pieces.append(_blank(text[pos:]))
    return "".join(pieces)
