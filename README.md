# fast-tcp

`fast-tcp` is a practical implementation of FAST-style similarity-based test case prioritization. This repository contains the Python package, ready-to-run integrations for common test frameworks, example projects, and the accompanying research/demo artifacts.

## What is in this repository

- `fast_tcp/`: the installable Python package and framework integrations
- `examples/`: runnable example projects for pytest, Vitest, Go `testing`, and Ant/JUnit
- [research/](research/README.md): experiment scripts, pinned environments, recorded results, paper draft, and demo materials
- [PUBLISHING.md](PUBLISHING.md): PyPI and TestPyPI release instructions

## Installation

### From source

```bash
git clone https://github.com/titi-lima/FAST.git
cd FAST
python -m pip install -e . datasketch==1.6.5
```

### From TestPyPI

Version `1.3.0` is available on [TestPyPI](https://test.pypi.org/project/fast-tcp/1.3.0/). Install the dependencies from PyPI, then install the package from TestPyPI:

```bash
python -m pip install datasketch==1.6.5 xxhash
python -m pip install --index-url https://test.pypi.org/simple/ --no-deps fast-tcp==1.3.0
```

The package has not yet been published to production PyPI. Maintainer release instructions are in [PUBLISHING.md](PUBLISHING.md).

These commands pin `datasketch==1.6.5` because `fast-tcp` 1.3.0 is incompatible with `datasketch` 2.0.

## Quick start with pytest

```bash
python -m pip install pytest
fast-tcp init pytest
pytest --fast-tcp --fast-tcp-debug
```

`fast-tcp init pytest` creates a `.fast/` workspace and updates `pytest.ini` so the package can persist signatures and filesystem snapshots between runs.

## Other supported integrations

- `Vitest`: `fast-tcp init vitest` then `npm run test:fast`
- `Go testing`: `fast-tcp init gotest` then `bash .fast/tools/gotest/run-fast.sh`
- `Ant/JUnit`: `fast-tcp init ant` then `ant fast-tcp`
- Direct CLI: `fast-tcp --test-dir <dir> --algo FAST-pw --entity bbox ...`

## Direct CLI example

The standalone CLI prioritizes plain-text black-box suites where each line represents one test command or identifier.

```bash
fast-tcp \
  --test-dir .fast/in \
  --algo FAST-pw \
  --entity bbox \
  --file-naming entity-suffix \
  --output-dir .fast/out
```

## Python API example

```python
from fast_tcp import run_blackbox_file

partition_time, prep_time, prio_time, order = run_blackbox_file(
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

## Research context

The implementation is based on:

> Breno Miranda, Emilio Cruciani, Roberto Verdecchia, and Antonia Bertolino. 2018. FAST Approaches to Scalable Similarity-based Test Case Prioritization. ICSE 2018. DOI: [10.1145/3180155.3180210](https://doi.org/10.1145/3180155.3180210)

Experiment scripts, pinned project revisions and dependencies, complete test reports, and instructions to regenerate the measured figures are in [research/replication/](research/replication/README.md). They run from this repository; the Overleaf project is not required. Start with the offline audit of the recorded results, or follow the instructions to rerun the experiments.

The paper draft, exploratory notebooks, and demo materials are indexed in [research/](research/README.md).

## Repository notes

- Supported Python baseline in this repo is `3.8+`
- The package currently exposes black-box prioritization through the public CLI
- Snapshot data and generated integration assets live under `.fast/`
