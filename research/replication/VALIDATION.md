# Artifact validation — 2026-10-23

- Audited all 48 historical real-project executions against their raw phase reports and recomputed summary statistics.
- Audited the historical 20-step synthetic trace; retained its records with the documented empty-snapshot limitation.
- Prepared new Python 3.11.13 environments for all four pinned projects in a clean source copy. Dependency compatibility checks passed for all four. Checkouts came from local copies of the exact upstream commits; most dependencies came from the local cache, with an online fetch for Rich's remaining dependencies. This validates installation without claiming a new GitHub download test.
- Executed the Click smoke command: collection plus four complete executions, each collecting 611 tests, with 588 passed, 22 skipped, and one xfailed. The independent audit passed.
- Executed and audited the corrected synthetic workload: warmup plus 20 baseline/FAST pairs, finishing with 140 passing tests. Each FAST report records the current snapshot file list; the audit checks it against the executed tests. Full reports and the environment freeze are in `data/corrected-synthetic/`.
- Executed the snapshot workload through all 40 mutations after the seed. Every detected added/modified/deleted count matched the generated mutation. Saved 41 exact-byte records under `data/snapshot-validation/`: working files grew from 100 MiB to 114.5 MiB; snapshot repository size grew from 185,281 to 268,295 bytes.
- Checked that the auditor rejects incorrect summary statistics, an omitted test report, and an empty synthetic snapshot.
- Regenerated both measured figure PDFs and PNG previews; rendered the PDFs for visual inspection. Adjusted the case-study value labels to avoid overlap.
- Checked Python syntax, local documentation links, and whitespace errors.

The new runs validate the reproduction path. They do not replace the archived real-project timing study or establish performance across machines. No tool runtime code or manuscript text was changed as part of this artifact preparation.
