# Experiment replication

Run commands below from the repository root. This directory contains the experiment code and recorded evidence independently of the manuscript sources.

## Audit the recorded results (no installation or network needed)

Requires Python 3.8+:

```bash
python3 research/replication/audit_results.py
```

The audit checks SHA-256 digests, reads the raw pytest phase reports directly from `data/full-run-reports.zip`, compares collected and executed test identities and outcomes across conditions, and recomputes the real-project means and sample standard deviations. It also checks the historical and corrected synthetic workloads' test counts, outcomes, and measured FAST components, including the corrected run's snapshot file lists. Missing reports, failed runs, mismatched outcomes, or inconsistent statistics cause a nonzero exit.

The archive contains stdout, stderr, commands, process times, test identities, setup/call/teardown outcomes, environment freezes, and hardware information. Original absolute paths in these historical records are provenance only; the scripts do not use them. Synthetic process wall times are recorded in `synthetic-measurements.json`; the independent audit verifies synthetic component timings against the observer reports, rather than reconstructing process wall times from test durations.

## Rerun the experiments

Prerequisites: Git, [uv](https://docs.astral.sh/uv/), network access to GitHub and Python package indexes, and a Unix-like system. Preparation installs Python 3.11.13 and creates separate environments. Each dependency freeze comes from the recorded run; `appnope` is restricted to macOS. Other platforms may collect or skip different tests. The original measurements used macOS 26.5.2, Apple M4 (10 cores), 16 GiB RAM; see `data/hardware.json` and `data/verified-summary.json` for provenance and source hashes.

```bash
python3 research/replication/prepare_projects.py
python3 research/replication/verify_projects.py
python3 research/replication/measure_synthetic.py
output/research/replication/envs/rich/bin/python research/replication/measure_snapshot.py
python3 research/replication/audit_results.py --results-dir output/research/replication
```

Generated checkouts, environments, caches, logs, and JSON results go to `output/research/replication/`, which is ignored by Git. The scripts preserve the bundled paper results. To repeat an experiment, use a new `--output-dir` consistently for preparation and measurement; measurement scripts refuse to overwrite existing project results. Allow roughly 30 minutes for all real-project runs, plus dependency installation. Tests execute third-party project code; use a disposable environment.

Generated synthetic and snapshot projects have their own empty Git repositories so the enclosing repository's `output/` ignore rule does not hide their files from snapshot detection. Snapshot working-file sizes exclude `.git/` and `.fast/`; the script checks each detected mutation against the generated change counts.

For a shorter check:

```bash
python3 research/replication/prepare_projects.py --projects click --output-dir output/research/smoke
python3 research/replication/verify_projects.py --projects click --repetitions 1 --output-dir output/research/smoke
python3 research/replication/audit_results.py --results-dir output/research/smoke
```

`measure_synthetic.py --python /absolute/path/to/python` can use any prepared environment containing this package and pytest. Otherwise it uses the Rich environment in the chosen output directory. Run the snapshot script with a prepared environment's Python, as shown above. It needs about 120 MiB for the working files and performs 40 mutations after the seed snapshot.

## Experiment definitions

| Experiment | Code | Recorded evidence |
| --- | --- | --- |
| Real-project case study | `prepare_projects.py`, `verify_projects.py`, `fast_verification_audit.py` | `data/verified-summary.json` and `runs/` inside the ZIP |
| Executed incremental workload | `measure_synthetic.py`, `fast_verification_audit.py` | Corrected run in `data/corrected-synthetic/`; historical run retained in `data/synthetic-measurements.json` and the original ZIP |
| Snapshot storage growth | `measure_snapshot.py` (standalone version of the notebook) | Original rounded output in [snapshot_cache_growth.ipynb](../snapshot_cache_growth.ipynb); new runs save exact byte counts in `snapshot/measurements.json` |

### Real-project case study

[projects.json](projects.json) pins Requests 2.31.0, Rich 13.7.0, Click 8.1.7, and attrs 23.2.0 to exact Git commits. Preparation verifies each commit and installs the dependency versions in `requirements/<project>.txt`, this tool, and the target project. The measured FAST source revision and recursive Python source hashes are in the recorded summary. Reruns use the tool source in your checkout and record its revision and hashes; later versions need not produce identical timings.

For each project, the harness collects the whole `tests/` suite, executes one baseline screening run and one cold FAST run, then five baseline/FAST pairs with alternating order. That is 12 complete executions per project, 48 total. Collection-only output is excluded from execution timings. Both conditions load the same observer and FAST plugin; only the FAST condition enables prioritization. No early-stop flag or test-file exclusion is used. Project `addopts` and unrelated pytest plugin autoloading are disabled; required plugins are listed explicitly in `projects.json`. Each process has a 600-second timeout.

The cold run follows snapshot initialization with no persisted signatures. Warm repetitions reuse the snapshot and signatures without source changes. Means and sample standard deviations use only the five warm runs per condition. Expected skips and xfails are retained and compared by test identity, including xfails reported during setup.

| Project | Collected | Passed | Skipped | Xfailed |
| --- | ---: | ---: | ---: | ---: |
| Requests | 596 | 582 | 13 | 1 |
| Rich | 826 | 803 | 23 | 0 |
| Click | 611 | 588 | 22 | 1 |
| attrs | 1410 | 1405 | 4 | 1 |

Process wall time includes the full pytest process. FAST driver time wraps `_run_fast_blackbox`; partition, preparation, and ranking times come from its returned measurements. The wall-time difference between conditions is not a direct measurement of prioritization overhead. These observations cover one machine and four pinned suites, with natural test-order effects and system noise.

### Executed incremental workload

Start with ten files containing five passing tests each. At every one of 20 steps, add one five-test file and append a comment to one original file; delete an earlier added file at steps 10 and 20. Execute baseline and FAST after each mutation and require all collected tests to pass. A FAST warmup precedes the measured steps. The final suite has 140 tests. This is a generated workload executed by pytest; its timings are measured. It is not a fault-detection or industrial-scale evaluation.

**Correction found while preparing this artifact:** the historical synthetic project lived under the parent repository's ignored `output/` directory. Its snapshot tree was empty. Those timings came from real passing test executions, but cannot substantiate incremental snapshot change detection. The corrected harness isolates the project, checks that every current test file is in the snapshot, and records that list in the raw observer reports. The corrected measurements and full reports are retained separately in `data/corrected-synthetic/`. The trace renderer uses the corrected run by default; a manuscript using the earlier trace should update that figure and any derived numbers. The four real-project suites have their own Git checkouts and are unaffected by this enclosing-ignore issue.

### Snapshot storage growth

Seed 50 files of 2 MiB each, then run 40 mutations with random seed 42: rewrite 18–32 base files per step, add 2–6 temporary files of 0.5 or 1 MiB, and remove temporary files after six steps. Every payload repeats a label, making it highly compressible. Record change counts, working-file bytes, snapshot repository bytes, and loose-object counts after snapshot processing. This measures this compressible workload, rather than a general storage bound for arbitrary test data. Git version and repository metadata can affect byte counts.

The original notebook did not export its exact byte-count JSON. Its saved, rounded output remains available; we do not treat reconstructed rounded values as raw measurements. The standalone script retains the notebook's workload and saves every new record automatically. A new 41-record validation run is included in `data/snapshot-validation/`, separately from the notebook's historical output.

## Regenerate the measured figures

Create a plotting environment:

```bash
uv venv output/research/plot-env --python 3.11.13
uv pip install --python output/research/plot-env/bin/python -r research/replication/requirements/figures.txt
output/research/plot-env/bin/python research/replication/render_case_study.py
output/research/plot-env/bin/python research/replication/render_measured_trace.py
```

This writes cropped vector PDFs and PNG previews (`case_study_bar` and `benchmark_grid`) under `output/research/replication/figures/`, using the bundled real-project and corrected synthetic measurements. Ready-to-use copies are included in [figures/](figures/). Each renderer also accepts `--input PATH` and `--output-dir PATH` to plot a new run. The case-study renderer expects all four projects. To reconstruct the historical trace explicitly, pass `--input research/replication/data/synthetic-measurements.json`. The legacy `research/figures/generate_figures.py` contains preliminary and constructed values and is not part of this reproduction path.

For the snapshot chart, the original notebook includes the plotting code. New exact-byte measurements are saved separately so they can be compared with its rounded historical output.

## Files and preservation

`data/` contains the original verified summary, measured trace, full reports archive, hardware record, and `SHA256SUMS.json`. The checksum file detects changes relative to the bundled records; it is not an independent authenticity guarantee. `audit_results.py` does not extract the archive or modify any evidence.

The scripts and recorded results can be distributed with the same GitHub release as the tool. A permanent archive/DOI can be added when that release is deposited; no DOI has been assigned here.
