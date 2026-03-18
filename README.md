# fast-tcp

`fast-tcp` is a practical implementation of FAST-style similarity-based test case prioritization. This repository contains the Python package, ready-to-run integrations for common test frameworks, example projects, and the accompanying research/demo artifacts.

## What is in this repository

- `fast_tcp/`: the installable Python package and framework integrations
- `examples/`: runnable example projects for pytest, Vitest, Go `testing`, and Ant/JUnit
- `research/`: the paper draft, figures, notebooks, and demo materials
- [PUBLISHING.md](PUBLISHING.md): PyPI and TestPyPI release instructions

## Installation

### From source

```bash
git clone https://github.com/icse18-FAST/FAST.git
cd FAST
python -m pip install -e .
```

### PyPI release target

The package is prepared for PyPI/TestPyPI publication. Until the first release is published, use the source install above and the release workflow documented in [PUBLISHING.md](PUBLISHING.md).

## Quick start with pytest

```bash
python -m pip install -e . pytest
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

The current tooling, paper draft, and artifact notes for this repository live under [research/](research/).

## Repository notes

- Supported Python baseline in this repo is `3.8+`
- The package currently exposes black-box prioritization through the public CLI
- Snapshot data and generated integration assets live under `.fast/`
