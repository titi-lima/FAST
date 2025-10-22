### FAST TCP + Go (testing) Integration

This scaffolding lets you run Go tests in a FAST-prioritized order.

### Quick start

1) Initialize at your project root:

```bash
fast-tcp init gotest
```

This creates:
- `.fast/in` and `.fast/out` directories
- `.fast/tools/gotest/` helper scripts
- `.fast/USAGE.md` with these instructions

2) Run prioritized tests:

```bash
bash .fast/tools/gotest/run-fast.sh
```

What happens:
- Discover test function names (`TestXxx`) across `*_test.go`
- Generate bbox tokens per test case
- Run FAST prioritization via the Python CLI
- Map prioritized IDs back to Go test names
- Build package test binaries and execute tests in that order via `-test.run`

### Customize

Re-run init with flags to set algorithm and repetitions:

```bash
fast-tcp init gotest --algo FAST-pw --repetitions 3
```

You can also edit `.fast/tools/gotest/run-fast.sh` to pass extra flags or adjust discovery.

### Requirements

- Python 3.8+
- Go 1.20+
- `fast-tcp` installed in your Python environment



