# Pytest integration example (FAST TCP)

This example shows how to integrate fast-tcp with pytest to prioritize the order in which tests run.

It implements a small pytest plugin in `conftest.py` that:
- Collects the test nodeids
- Runs FAST in black-box mode (string shingling over nodeids)
- Reorders the collected items before execution

## Prerequisites

- Python 3.8+
- Install packages:
  
```bash
# initialize virtual environment
python3 -m venv virtual
source virtual/bin/activate
# install deps (if in pre-release, make sure to add --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple so you'll get the testpypi version)
pip install -r requirements.txt 
```

## Run
You can run this from any directory, but we'll assume you're doing it from the root of the repository.
```bash
pytest -q examples/python/pytest-fast
```

Or to see plugin options:

```bash
pytest -q examples/python/pytest-fast --help | sed -n '/fast-tcp:/,/^$/p'
```

## Plugin options

- `--fast-tcp` (bool): enable FAST-based ordering (disabled by default)
- `--fast-tcp-algo {FAST-pw,FAST-one,FAST-log,FAST-sqrt,FAST-all}`: choose variant (default: `FAST-pw`)
- `--fast-tcp-r INT`: number of rows (default: 1)
- `--fast-tcp-b INT`: number of bands (default: 10)
- `--fast-tcp-k INT`: shingle size for black-box (default: 5)
- `--fast-tcp-budget INT`: optional budget B (0 means all tests)

Examples:

```bash
# Enable FAST with defaults
pytest -q examples/python/pytest-fast --fast-tcp

# Use FAST-log variant
pytest -q examples/python/pytest-fast --fast-tcp --fast-tcp-algo FAST-log

# Increase shingle size and bands
pytest -q examples/python/pytest-fast --fast-tcp --fast-tcp-k 7 --fast-tcp-b 20
```

## Notes

- This example uses black-box prioritization over test nodeids. It does not require coverage.
