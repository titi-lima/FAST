#!/usr/bin/env python3
"""
Quick CLI validation script - tests basic functionality without complex algorithms.
This is perfect for rapid development iteration.
"""

import subprocess
import sys
import tempfile
from pathlib import Path


def run_command(cmd, timeout=30):
    """Run a command and return success status and output."""
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return False, "", "Command timed out"
    except Exception as e:
        return False, "", str(e)


def test_help():
    """Test that --help works."""
    print("Testing --help command...")
    success, stdout, stderr = run_command(
        [sys.executable, "-m", "fast_tcp.cli", "--help"]
    )

    if success and "Test Case Prioritization Tool" in stdout:
        print("✓ --help works")
        return True
    else:
        print(f"✗ --help failed: {stderr}")
        return False


def test_missing_args():
    """Test that missing required args are caught."""
    print("Testing missing required arguments...")
    success, stdout, stderr = run_command([sys.executable, "-m", "fast_tcp.cli"])

    if not success and "required" in stderr.lower():
        print("✓ Missing args properly detected")
        return True
    else:
        print(f"✗ Missing args not detected properly")
        return False


def test_basic_execution():
    """Test basic execution with real data."""
    print("Testing basic execution with sample data...")

    # Use existing sample data if available
    sample_dirs = ["input/chart_v0", "input/math_v0", "input/lang_v0"]

    test_dir = None
    for sample_dir in sample_dirs:
        if Path(sample_dir).exists():
            test_dir = sample_dir
            break

    if not test_dir:
        print("ℹ No sample data found, skipping basic execution test")
        return True

    with tempfile.TemporaryDirectory() as temp_output:
        success, stdout, stderr = run_command(
            [
                sys.executable,
                "-m",
                "fast_tcp.cli",
                "--test-dir",
                test_dir,
                "--algo",
                "FAST-pw",
                "--entity",
                "function",
                "--output-dir",
                temp_output,
                "--repetitions",
                "1",
            ],
            timeout=60,
        )

        if success:
            print("✓ Basic execution works")
            return True
        else:
            print(f"✗ Basic execution failed")
            print(f"  stdout: {stdout[:200]}")
            print(f"  stderr: {stderr[:200]}")
            return False


def main():
    """Run quick tests."""
    print("FAST TCP CLI Quick Test")
    print("======================")

    tests = [
        test_help,
        test_missing_args,
        test_basic_execution,
    ]

    passed = 0
    total = len(tests)

    for test in tests:
        if test():
            passed += 1
        print()

    print("=" * 30)
    print(f"Results: {passed}/{total} tests passed")

    if passed == total:
        print("🎉 All quick tests passed!")
        return 0
    else:
        print("❌ Some tests failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
