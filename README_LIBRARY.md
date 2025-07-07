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
import fast_tcp

# Load a test suite and run FAST-pw prioritization
from fast_tcp import fast_pw, fast_

# Run FAST-pw algorithm
input_file = "path/to/test_suite.txt"
r, b = 1, 10  # LSH parameters
prioritized_tests = fast_pw(input_file, r, b, bbox=True)

# Run other FAST variants
def sqrt_selsize(x):
    return int(x**0.5) + 1

prioritized_tests_sqrt = fast_(input_file, sqrt_selsize, r, b, bbox=True)
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

### Core Functions

```python
from fast_tcp import fast_pw, fast_

# FAST-pw: Pairwise comparison with LSH
stime, ptime, prioritized = fast_pw(
    input_file,     # Path to test suite file
    r,              # LSH rows parameter
    b,              # LSH bands parameter
    bbox=False,     # True for black-box, False for white-box
    k=5,            # k-shingle size (for black-box)
    memory=True,    # Keep signatures in memory
    B=0             # Budget (0 = unlimited)
)

# FAST with custom sample size function
stime, ptime, prioritized = fast_(
    input_file,     # Path to test suite file  
    selsize,        # Sample size function
    r,              # LSH rows parameter
    b,              # LSH bands parameter
    bbox=False,     # True for black-box, False for white-box
    k=5,            # k-shingle size (for black-box)
    memory=True,    # Keep signatures in memory
    B=0             # Budget (0 = unlimited)
)
```

### Competitor Algorithms

```python
from fast_tcp import str_, i_tsd, gt, ga, ga_s, artd, artf

# String-based similarity (black-box)
stime, ptime, prioritized = str_(input_file)

# Information Theory-based (black-box)  
stime, ptime, prioritized = i_tsd(input_file)

# Greedy Total (white-box)
stime, ptime, prioritized = gt(input_file)

# And so on for other algorithms...
```

### Utility Functions

```python
from fast_tcp import apfd

# Calculate APFD metric
apfd_score = apfd(prioritized_tests, fault_matrix, java_flag=False)
```

## Examples

### Basic Usage

```python
import fast_tcp

# Run FAST-pw on a black-box test suite
input_file = "input/flex_v3/flex-bbox.txt"
stime, ptime, prioritized = fast_tcp.fast_pw(
    input_file, r=1, b=10, bbox=True, k=5
)

print(f"Signature time: {stime:.2f}s")
print(f"Prioritization time: {ptime:.2f}s") 
print(f"Prioritized order: {prioritized[:10]}")  # First 10 tests
```

### Custom Sample Size

```python
import math
from fast_tcp import fast_

def logarithmic_sample_size(x):
    return int(math.log(x, 2)) + 1

# Run FAST with logarithmic sample size
stime, ptime, prioritized = fast_(
    "input/test_suite.txt",
    logarithmic_sample_size, 
    r=1, b=10, bbox=False
)
```

### Batch Processing

```python
import fast_tcp

algorithms = ['FAST-pw', 'FAST-sqrt', 'GT']
datasets = ['flex_v3', 'grep_v3']

for dataset in datasets:
    for algorithm in algorithms:
        print(f"Running {algorithm} on {dataset}")
        # Run experiments programmatically
        # ... your experiment code here
```

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