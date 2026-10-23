# fast-tcp

`fast-tcp` packages FAST-style similarity-based test case prioritization for real projects. It includes a standalone CLI, a pytest plugin, and helpers for Go `testing` and Ant/JUnit workflows.

## Highlights

- FAST-family prioritization (`FAST-pw`, `FAST-one`, `FAST-log`, `FAST-sqrt`, `FAST-all`)
- Incremental runs backed by a Git-based snapshot cache under `.fast/`
- Native pytest integration plus helper setup for Go and Ant
- Python API for driving black-box prioritization directly

## Installation

### TestPyPI

Version `1.3.1` is available on [TestPyPI](https://test.pypi.org/project/fast-tcp/1.3.1/). Install from TestPyPI with PyPI as an additional index for dependencies:

```bash
python -m pip install \
  --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple/ \
  fast-tcp==1.3.1
```

The package has not yet been published to production PyPI.

### From this repository

```bash
git clone https://github.com/titi-lima/FAST.git
cd FAST
python -m pip install -e .
```

Python `3.8+` is the supported baseline for the packaged tool.

## Pytest quick start

```bash
python -m pip install pytest
fast-tcp init pytest
pytest --fast-tcp --fast-tcp-debug
```

Useful pytest flags:

- `--fast-tcp`
- `--fast-tcp-algo FAST-pw`
- `--fast-tcp-r 1`
- `--fast-tcp-b 10`
- `--fast-tcp-k 5`
- `--fast-tcp-budget 0`
- `--fast-tcp-debug`

## Other integrations

### Vitest (experimental)

The package includes a Vitest scaffold; successful execution of this integration has not been verified.

```bash
fast-tcp init vitest
```

### Go `testing`

```bash
fast-tcp init gotest
bash .fast/tools/gotest/run-fast.sh
```

### Ant / JUnit

```bash
fast-tcp init ant
ant fast-tcp
```

## Standalone CLI

The direct CLI works on plain-text black-box suites.

```bash
fast-tcp \
  --test-dir .fast/in \
  --algo FAST-pw \
  --entity bbox \
  --file-naming entity-suffix \
  --output-dir .fast/out
```

Key arguments:

- `--test-dir`: directory containing the input `.txt` suite
- `--algo`: one of `FAST-pw`, `FAST-one`, `FAST-log`, `FAST-sqrt`, `FAST-all`
- `--entity`: currently `bbox`
- `--repetitions`: repetition count passed to the legacy prioritization driver
- `--pattern`: file glob used when discovering suites
- `--debug`: prints timing and partition information

## Python API

```python
from fast_tcp import run_blackbox_file

partition_time, prep_time, prio_time, prioritized_order = run_blackbox_file(
    algo="FAST-pw",
    input_file=".fast/in/suite-bbox.txt",
    signature_dir=".fast/signatures",
    k=5,
    r=1,
    b=10,
    budget=0,
    project_root=".",
    debug=True,
)
```

`prioritized_order` contains the prioritized 1-based test indices from the input suite.

## Project layout

Runtime data is written under `.fast/`:

- `.fast/signatures/`: persisted MinHash signatures
- `.fast/snapshot.git/`: Git-based snapshot cache
- `.fast/in/`: generated integration inputs
- `.fast/out/`: prioritization outputs

## Research background

This implementation is derived from:

> Breno Miranda, Emilio Cruciani, Roberto Verdecchia, and Antonia Bertolino. 2018. FAST Approaches to Scalable Similarity-based Test Case Prioritization. ICSE 2018. DOI: [10.1145/3180155.3180210](https://doi.org/10.1145/3180155.3180210)

Repository-level documentation is in the [GitHub README](https://github.com/titi-lima/FAST/blob/master/README.md). Experiment scripts, pinned environments, recorded full-suite reports, and reproduction instructions are in [research/replication/](https://github.com/titi-lima/FAST/tree/master/research/replication).
