# CLI Testing Guide

This directory contains comprehensive test scripts to ensure the FAST TCP CLI remains functional during development.

## Test Scripts

### `quick_test.py` ⚡ (Recommended for Development)
A fast, lightweight test script that validates core CLI functionality:

- **Help and argument parsing**: Tests `--help` command and required argument validation  
- **Error handling**: Tests missing required arguments
- **Basic execution**: Tests with real sample data if available
- **Runtime**: ~5-10 seconds

### `test_cli.py`
A comprehensive Python test suite that validates:

- **Help and argument parsing**: Tests `--help` command and required argument validation
- **Error handling**: Tests invalid algorithms, entities, and directories
- **File discovery**: Tests automatic file detection and custom patterns
- **Algorithm execution**: Tests core prioritization algorithms
- **Output generation**: Validates that results are properly generated
- **Runtime**: ~30-60 seconds (quick mode) or 5-10 minutes (all algorithms)

### `test_cli.sh`
A shell script wrapper that provides convenient testing commands:

```bash
# Quick tests (recommended for regular development)
./test_cli.sh

# Comprehensive tests (tests all algorithms - slower)
./test_cli.sh --all

# Manual testing mode (opens interactive shell)
./test_cli.sh --manual
```

## Usage

### Quick Testing (Recommended)
For regular development work, run the lightweight quick test:

```bash
python quick_test.py
```

This tests core CLI functionality with minimal overhead (~5-10 seconds).

For more comprehensive testing without testing all algorithms:

```bash
python test_cli.py
```

This tests core functionality plus file discovery and some algorithms (~30-60 seconds).

### Comprehensive Testing
Before major releases or when making algorithm changes:

```bash
python test_cli.py --all-algorithms
```

This tests all available algorithms: FAST-pw, FAST-one, FAST-log, FAST-sqrt, FAST-all, STR, I-TSD, ART-D, ART-F, GT, GA, GA-S.

### Integration with Development
Add this to your development workflow:

1. **During development**: Run `python quick_test.py` for fast feedback
2. **Before making changes**: Run `python test_cli.py` to establish baseline  
3. **After making changes**: Run `python quick_test.py` then `python test_cli.py` to ensure no regression
4. **Before committing**: Run `./test_cli.sh --all` for comprehensive validation

## Test Coverage

The test suite covers:

### Core CLI Functionality
- [x] Help command (`--help`)
- [x] Required argument validation
- [x] Invalid algorithm/entity handling
- [x] Non-existent directory handling

### File Discovery
- [x] Auto file detection (based on filename patterns)
- [x] Entity-suffix naming convention
- [x] Custom file patterns
- [x] Multiple file types (bbox, function, branch, line)

### Algorithm Testing
- [x] FAST variants (pw, one, log, sqrt, all)
- [x] Competitor algorithms (STR, I-TSD, ART-D, ART-F, GT, GA, GA-S)
- [x] Multiple repetitions
- [x] Different entity types

### Output Generation
- [x] Default output directory creation
- [x] Custom output directory
- [x] Result file generation
- [x] Proper directory structure

## Test Data

The test suite automatically generates synthetic test data:

- **Black-box data**: Command-line style strings with error/feature/category patterns
- **White-box data**: Space-separated numeric entity coverage patterns
- **Fault matrices**: Pickle files with simulated fault detection data

## Troubleshooting

### Common Issues

1. **Module not found error**:
   ```bash
   pip install -e .
   ```

2. **Permission denied**:
   ```bash
   chmod +x test_cli.sh
   ```

3. **Tests failing after changes**:
   - Check the detailed output for specific failure reasons
   - Ensure your changes maintain backward compatibility
   - Verify that required dependencies are still available

### Test Output

Successful test output should look like:
```
============================================================
FAST TCP CLI Test Suite
============================================================
Testing --help command...
✓ PASSED

Testing missing required arguments...
✓ PASSED

...

============================================================
Test Results: 10/10 tests passed
============================================================
```

### Adding New Tests

To add new test cases, extend the `CLITestSuite` class in `test_cli.py`:

```python
def test_your_new_feature(self) -> bool:
    """Test description."""
    print("Testing your new feature...")
    
    # Test implementation
    success = your_test_logic()
    
    self.test_results.append({
        "test": "your_new_feature",
        "success": success,
        "returncode": 0,
        "details": "Test details"
    })
    
    return success
```

Then add it to the `test_methods` list in `run_all_tests()`.

## CI/CD Integration

For continuous integration, use:

```bash
# In your CI script
python test_cli.py
exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo "CLI tests passed"
else
    echo "CLI tests failed"
    exit 1
fi
```

The test script returns exit code 0 for success and 1 for failure, making it suitable for automated testing pipelines.