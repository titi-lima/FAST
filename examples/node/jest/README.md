### Node + Jest example (bbox-only)

This example shows how to integrate FAST TCP as a headless test case prioritizer with a Node project using Jest.

It follows the bbox-only example conventions: we extract tokens from individual test bodies, run the FAST prioritizer to compute an order, map prioritized IDs back to test names and file paths, and finally run Jest in that prioritized order.

### Layout

```
examples/node/jest/
├── __tests__/                # sample tests
├── src/                      # sample code
├── scripts/
│   ├── generate-bbox.py      # tokens-from-tests extractor
│   ├── tests-map.py          # discovery order of tests (file + test name)
│   ├── build-prioritized-tests.py  # map CLI pickle to file+name list
│   └── run-prioritized-tests.sh    # end-to-end runner
├── .fast/
│   ├── in/                   # generated bbox inputs
│   └── out/                  # CLI outputs
├── package.json
└── README.md
```

### Prereqs

- Python 3.9+
- Node.js 18+ and npm

### Run end-to-end

From repo root (or this directory), run:

```bash
bash examples/node/jest/scripts/run-prioritized-tests.sh
```

What it does:

- Generates `.fast/in/test-names.tsv` (file + test names)
- Generates `.fast/in/jest-bbox.txt` (tokens, one line per test case)
- Runs prioritization via the FAST TCP CLI
- Maps prioritized IDs to `(file, testName)` pairs
- Installs dev deps and runs Jest once per test in priority order

Notes:

- This runs Jest one test at a time using `-t` against a specific file for clarity. For larger suites, you can batch or aggregate by file using a custom Jest Sequencer.
- bbox is the only supported mode in examples per the current scope.



