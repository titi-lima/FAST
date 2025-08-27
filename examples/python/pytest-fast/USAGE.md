## FAST TCP for Pytest — Quick Start

Make Pytest run tests in an order prioritized by FAST. The plugin is auto-discovered when `fast-tcp` is installed.

### 1) Initialize your project

From your project root:

```bash
fast-tcp init pytest --project-dir .
```

This will create/update `pytest.ini` to enable FAST by default with sensible parameters. To keep FAST disabled by default, run:

```bash
fast-tcp init pytest --project-dir . --no-enable
```

### 2) Run tests

If FAST was enabled in `pytest.ini`:

```bash
pytest
```

If not enabled by default, turn it on per run:

```bash
pytest --fast-tcp
```

### Options

You can override parameters at the CLI:

```bash
pytest --fast-tcp \
       --fast-tcp-algo FAST-pw \
       --fast-tcp-r 1 \
       --fast-tcp-b 10 \
       --fast-tcp-k 5 \
       --fast-tcp-budget 0
```

- `--fast-tcp` enables prioritization
- `--fast-tcp-algo` one of: `FAST-pw`, `FAST-one`, `FAST-log`, `FAST-sqrt`, `FAST-all`
- `--fast-tcp-r` rows (default: 1)
- `--fast-tcp-b` bands (default: 10)
- `--fast-tcp-k` shingle size for black-box over nodeids (default: 5)
- `--fast-tcp-budget` number of tests to prioritize strictly (0=all)

### Notes

- The plugin works in black-box mode using test nodeids; no coverage required.
- It preserves tests omitted by a budget by appending them in original order.

