## Java (Ant) integration with FAST TCP

This example shows how to integrate a Java 21 project built with Ant to produce inputs for FAST TCP and run prioritization in a headless way.

The example focuses on black-box prioritization (bbox), where each test case is represented by text tokens (e.g., test method bodies). White-box per-test coverage requires per-test coverage tooling, which is not included here.

### Prerequisites

- Java 21 (JDK 21) available on PATH
- Apache Ant 1.10+
- Python 3.8+ and this repository installed or available to run:
  - From repository root:
    - Development install: `pip install -e .`
    - Or run via module path: `python -m fast_tcp.cli --help`

### Layout

```
examples/java/ant/
├── build.xml
├── lib/
├── scripts/
│   ├── get-junit.sh
│   ├── generate-bbox.py
│   └── prioritize-bbox.sh
├── data/                  # Generated FAST TCP inputs (created by script)
├── out/                   # FAST TCP outputs (created by script)
└── src/
    ├── main/java/com/example/App.java
    └── test/java/com/example/AppTest.java
```

### 1) (Optional) Fetch JUnit and run tests

This project includes a minimal JUnit 5 test. If you want to compile and run tests via Ant:

```bash
# From repo root
cd examples/java/ant
./scripts/get-junit.sh
ant clean compile test
```

Note: Running tests is not required to generate black-box inputs; the bbox input is extracted from test sources.

### 2) Generate black-box input from tests

The script extracts tokens from `@Test` methods in `src/test/java` and writes one line per test case:

```bash
# From repo root
cd examples/java/ant
python3 scripts/generate-bbox.py --tests-dir src/test/java --out .fast/in/ant-sample-bbox.txt
```

This produces `data/ant-sample-bbox.txt` where each line is the tokenized content representing a test case.

### 3) Run FAST TCP prioritization (black-box)

Run FAST TCP using the generated input file. The helper script wraps the CLI for convenience:

```bash
# From repo root
cd examples/java/ant
./scripts/prioritize-bbox.sh
```

This will:
- Point FAST TCP to `examples/java/ant/data` as the test directory
- Use the `FAST-pw` algorithm on `bbox` entity
- Save results under `examples/java/ant/out/<dataset>/prioritized`

You can customize the algorithm and repetitions. For example:

```bash
python -m fast_tcp.cli \
  --test-dir $(pwd)/examples/java/ant/data \
  --algo FAST-sqrt \
  --entity bbox \
  --repetitions 5 \
  --output-dir $(pwd)/examples/java/ant/out
```

