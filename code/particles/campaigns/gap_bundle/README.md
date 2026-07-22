# Particle Gap Bundle Campaign

This is the code-repo-owned launch surface for the remaining particle blockers.
The work is bundled, not one issue at a time.

Bundle packets:

- `particles.gapbundle.electroweak-root`
- `particles.gapbundle.spectrum-source`
- `particles.gapbundle.qcd-thomson-backend`

The integration gate runs only after those packets return.

Dry run from `oph-meta/`:

```bash
./reverse-engineering-reality/code/particles/campaigns/gap_bundle/run_bundle.sh --max-workers 3 --dry-run
```

Live run:

```bash
./reverse-engineering-reality/code/particles/campaigns/gap_bundle/run_bundle.sh --max-workers 3
```

The script intentionally reuses the established Oracle runner in
`automation/oph_completion/` and the established worker homes under
`~/.oracle-oph/`.

## License

This particle campaign surface is part of the OPH public repository. See the
main [LICENSE](../../../../LICENSE).
