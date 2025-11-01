# FAST TCP: Scalable Similarity-based Test Case Prioritization

A Python library for test case prioritization using similarity-based algorithms including LSH-based approaches and competitor algorithms.

This library is based on the research paper:

> Breno Miranda, Emilio Cruciani, Roberto Verdecchia, and Antonia Bertolino. 2018. FAST Approaches to Scalable Similarity-based Test Case Prioritization. In *Proceedings of ICSE'18: 40th International Conference on Software Engineering, Gothenburg, Sweden, May 27-June 3, 2018 (ICSE'18)*, 11 pages. DOI: [10.1145/3180155.3180210](http://dx.doi.org/10.1145/3180155.3180210)

## Installation

### From PyPI (when published)
```bash
pip install fast-tcp
```

### From Source
```bash
git clone https://github.com/icse18-FAST/FAST.git
cd FAST
pip install .
```

### For Development
```bash
git clone https://github.com/icse18-FAST/FAST.git
cd FAST
pip install -e .
```

## Quick Start

### Command Line Interface

After installation, you can use the `fast-tcp` command:

```bash
# Run test case prioritization
fast-tcp prioritize flex_v3 bbox FAST-pw 50

# Run scalability experiments
fast-tcp scalability 1000 small FAST-pw

# Generate synthetic input for scalability testing
fast-tcp generate-input 1000 small

# Plot scalability results
fast-tcp plot-results small prioritization FAST-pw FAST-one

# Clean preprocessed files
fast-tcp clean
```

### Python API

You can also use the library directly in Python:

```python
from fast_tcp import run_blackbox_file

input_file = "path/to/test_suite.txt"
signature_dir = ".fast/signatures"
r, b, k = 1, 10, 5

# FAST-pw prioritization
prep_time, prio_time, prioritized_tests = run_blackbox_file(
    algo="FAST-pw",
    input_file=input_file,
    signature_dir=signature_dir,
    k=k,
    r=r,
    b=b,
    budget=0,
)

# Alternate FAST variants (FAST-log, FAST-sqrt, FAST-one, FAST-all)
_, _, prioritized_log = run_blackbox_file(
    algo="FAST-log",
    input_file=input_file,
    signature_dir=signature_dir,
    k=k,
    r=r,
    b=b,
    budget=0,
)
```

## Available Algorithms

### FAST Algorithms (LSH-based)
- **FAST-pw**: Pairwise comparison with LSH candidate filtering
- **FAST-one**: Fixed sample size of 1
- **FAST-log**: Logarithmic sample size
- **FAST-sqrt**: Square root sample size  
- **FAST-all**: Full comparison (no sampling)

### Competitor Algorithms
- **STR**: String-based similarity (black-box only)
- **I-TSD**: Information Theory-based Test Suite Dissimilarity (black-box only)
- **GT**: Greedy Total (white-box only)
- **GA**: Genetic Algorithm (white-box only)
- **GA-S**: Genetic Algorithm with Seeding (white-box only)
- **ART-D**: Adaptive Random Testing - Distance (white-box only)
- **ART-F**: Adaptive Random Testing - Failure (white-box only)

## Command Reference

### Prioritize Command

```bash
fast-tcp prioritize <dataset> <entity> <algorithm> <repetitions>
```

**Parameters:**
- `dataset`: Test suite to prioritize
  - Options: `flex_v3`, `grep_v3`, `gzip_v1`, `make_v1`, `sed_v6`, `closure_v0`, `lang_v0`, `math_v0`, `chart_v0`, `time_v0`
- `entity`: Type of prioritization
  - Options: `bbox` (black-box), `function`, `branch`, `line` (white-box)
- `algorithm`: Algorithm to use
  - Options: `FAST-pw`, `FAST-one`, `FAST-log`, `FAST-sqrt`, `FAST-all`, `STR`, `I-TSD`, `ART-D`, `ART-F`, `GT`, `GA`, `GA-S`
- `repetitions`: Number of prioritization runs (positive integer)

**Note:** Some algorithms are specific to black-box or white-box prioritization:
- Black-box only: `STR`, `I-TSD`
- White-box only: `ART-D`, `ART-F`, `GT`, `GA`, `GA-S`

### Scalability Command

```bash
fast-tcp scalability <tssize> <tcsize> <algorithm>
```

**Parameters:**
- `tssize`: Number of test cases in the test suite (positive integer)
- `tcsize`: Size category of test cases
  - Options: `small` (avg 100 elements), `medium` (avg 1K elements), `large` (avg 10K elements)
- `algorithm`: Algorithm to benchmark

### Generate Input Command

```bash
fast-tcp generate-input <test_suite_size> <test_case_size>
```

Generates synthetic test suites for scalability experiments.

### Plot Results Command

```bash
fast-tcp plot-results <test_case_size> <time_type> <algorithm1> [<algorithm2> ...]
```

**Parameters:**
- `test_case_size`: Size category (`small`, `medium`, `large`)
- `time_type`: Type of time measurement (`prioritization`, `total`)
- `algorithms`: One or more algorithms to include in the plot

### Clean Command

```bash
fast-tcp clean
```

Removes preprocessed input files for a clean environment.

## Directory Structure

When running experiments, the following directories are used:

```
./
├── input/          # Input datasets and coverage information
├── output/         # Prioritization results and performance data
└── scalability/    # Scalability experiment data
    ├── input/      # Generated synthetic test suites
    ├── output/     # Scalability results
    └── plots/      # Generated plots
```

## Python API Reference

Use the `run_blackbox_file` helper to execute any FAST variant in Python:

```python
from fast_tcp import run_blackbox_file

prep_time, prio_time, prioritized = run_blackbox_file(
    algo="FAST-pw",
    input_file="input/flex_v3/flex-bbox.txt",
    signature_dir=".fast/signatures",
    k=5,
    r=1,
    b=10,
    budget=0,
)
```

The helper caches MinHash signatures in the provided directory, making
subsequent runs faster if the test suite is unchanged.

## Examples

See `test.py` in the repository root for a minimal, self-contained example that
demonstrates preparation and prioritization flows using the new FAST core.

## Contributing

This library is based on research code. For issues or improvements, please refer to the original repository at https://github.com/icse18-FAST/FAST.

## License

This software is distributed under the GNU General Public License v3 (GPLv3). See the LICENSE file for details.

## Citation

If you use this library in your research, please cite:

```bibtex
@inproceedings{miranda2018fast,
  title={FAST approaches to scalable similarity-based test case prioritization},
  author={Miranda, Breno and Cruciani, Emilio and Verdecchia, Roberto and Bertolino, Antonia},
  booktitle={Proceedings of the 40th international conference on software engineering},
  pages={222--232},
  year={2018}
}
``` 