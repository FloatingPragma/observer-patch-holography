# Author Contribution Tracking

This directory keeps a reproducible git-based contribution snapshot for paper authorship decisions.

Use it when you need to answer questions like:

- who materially touched which paper file
- how broad each co-author's footprint is across papers, book, README, and supporting files
- how the current paper-order recommendation compares across contributors

## Scope

The tracker is intentionally conservative:

- it uses non-merge git commits only
- it tracks the canonical paper repo only
- it reports raw git metrics first and lets authorship decisions stay human

This means it will not capture off-Git work such as private notes, oral theorem input, review chat, or draft comments that never landed as commits.

## Inputs

- Author/alias mapping: [authors.json](./authors.json)
- Report generator: [generate_report.py](./generate_report.py)

## Regenerate

From the repo root:

```bash
python3 tracking/author_contributions/generate_report.py
```

This updates:

- [latest.json](./latest.json)
- [latest.md](./latest.md)

## Recommended Use

For paper ordering, read the report in this order:

1. paper rankings for the exact target paper files
2. recent paper commits for theorem-bearing or proof-bearing changes
3. broader category footprint as a tie-breaker

Do not use raw commit count alone as the sole authorship rule.

## License

This tracking surface is part of the OPH public repository. See the main
[LICENSE](../../LICENSE).
