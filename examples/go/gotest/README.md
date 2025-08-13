## Go (standard `testing`) integration with FAST TCP

This example shows how to integrate a Go project that uses the standard `testing` package with FAST TCP for headless test case prioritization.

Scope: black-box (bbox) prioritization only. We extract tokens from `_test.go` functions and let FAST TCP compute a prioritized order of tests. Then we run the compiled test binary in that order.

### Prerequisites

- Go 1.20+ on PATH
- Python 3.8+
- FAST TCP installed or available to run from this repository root, e.g.:
  - Development install from repo root: `pip install -e .`
  - Or run via module path: `python -m fast_tcp.cli --help`

### Layout

```
examples/go/gotest/
├── go.mod
├── hello/
│   ├── hello.go
│   └── hello_test.go
├── scripts/
│   ├── generate-bbox.py         # extract tokens from Test* functions
│   ├── tests-map.py             # list discovered Test* names in discovery order
│   ├── build-prioritized-tests.py
│   └── run-prioritized-tests.sh # end-to-end: generate → prioritize → run
├── .fast/
│   ├── in/   (created by scripts)
│   └── out/  (created by scripts)
└── README.md
```

### 1) (Optional) Run tests normally

```bash
# From repo root
cd examples/go/gotest
go test ./...
```

### 2) Generate black-box input from Go tests

The script extracts bodies of `func TestXxx(t *testing.T)` in `*_test.go` files and writes one line per test case.

```bash
# From repo root
cd examples/go/gotest
python3 scripts/generate-bbox.py \
  --tests-dir . \
  --out .fast/in/gotest-bbox.txt
```

### 3) Compute prioritized order and run tests accordingly

This one-shot script generates selectors, runs FAST TCP, maps IDs → test names, compiles the test binary, and executes tests one-by-one in the prioritized order using `-test.run`.

```bash
# From repo root
cd examples/go/gotest
./scripts/run-prioritized-tests.sh
```

Notes:
- The Go standard `testing` package does not provide a built-in way to impose an execution order within a single run. We emulate prioritization by compiling the test binary (`go test -c`) and invoking it multiple times, once per prioritized test name, using the `-test.run ^<TestName>$` filter.
- Outputs are stored under `.fast/out/<dataset>/prioritized/` and the prioritized test list under `.fast/in/prioritized-tests.txt`.



