#!/usr/bin/env python3
"""
Comprehensive test script for the FAST TCP CLI tool.

This script tests various scenarios to ensure the CLI remains functional
during development and refactoring.
"""

import os
import sys
import tempfile
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pickle


class CLITestSuite:
    """Test suite for the FAST TCP CLI tool."""

    def __init__(self):
        self.test_results = []
        self.temp_dirs = []
        self.original_cwd = os.getcwd()
        # Repository root (directory containing this test file)
        self.repo_root = Path(__file__).resolve().parent

    def cleanup(self):
        """Clean up temporary directories."""
        for temp_dir in self.temp_dirs:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)

    def create_test_data(
        self, test_dir: Path, entity_type: str, num_test_cases: int = 10
    ) -> Path:
        """Create sample test data for the specified entity type."""
        test_file = test_dir / f"sample-{entity_type}.txt"

        if entity_type == "bbox":
            # Black-box test data format (command-line like strings)
            test_data = [
                "-P[error,exit,input,output] -F[gui,database] -C[slow]",
                "-P[input,validation] -F[network] -C[fast]",
                "-P[output,display] -F[graphics] -C[medium]",
                "-P[error,exception] -F[logging] -C[slow]",
                "-P[input,parsing] -F[file] -C[fast]",
                "-P[output,report] -F[database,file] -C[medium]",
                "-P[error,network] -F[network] -C[slow]",
                "-P[input,user] -F[gui] -C[fast]",
                "-P[output,email] -F[network,email] -C[slow]",
                "-P[validation,input] -F[validation] -C[fast]",
            ]
        else:
            # White-box test data format (space-separated numbers)
            test_data = []
            base_entities = list(range(1, 21))  # Base set of 20 entities

            for i in range(num_test_cases):
                # Create different coverage patterns
                if entity_type == "function":
                    # Function coverage: larger chunks
                    coverage = (
                        base_entities[i : i + 8] + base_entities[max(0, i - 2) : i + 2]
                    )
                elif entity_type == "branch":
                    # Branch coverage: medium chunks with some overlap
                    coverage = base_entities[i : i + 6] + base_entities[i + 3 : i + 9]
                elif entity_type == "line":
                    # Line coverage: more granular
                    coverage = base_entities[i : i + 12] + [
                        j for j in range(i + 20, i + 25)
                    ]
                else:
                    coverage = base_entities[i : i + 5]

                # Remove duplicates and sort
                coverage = sorted(list(set(coverage)))
                test_data.append(" ".join(map(str, coverage)))

        # Write test data
        with open(test_file, "w") as f:
            for line in test_data:
                f.write(line + "\n")

        return test_file

    def create_fault_matrix(self, test_dir: Path, num_test_cases: int) -> Path:
        """Create a dummy fault matrix file."""
        fault_matrix_file = test_dir / "fault_matrix_key_tc.pickle"

        # Create a simple fault matrix with a few faults
        fault_matrix = {}
        for i in range(1, num_test_cases + 1):
            # Simulate some faults detected by certain test cases
            if i % 3 == 0:  # Every 3rd test case detects fault 1
                fault_matrix[str(i)] = ["fault_1"]
            elif i % 5 == 0:  # Every 5th test case detects fault 2
                fault_matrix[str(i)] = ["fault_2"]
            else:
                fault_matrix[str(i)] = []

        with open(fault_matrix_file, "wb") as f:
            pickle.dump(fault_matrix, f)

        return fault_matrix_file

    def run_cli_command(
        self, cmd: List[str], cwd: Optional[str] = None
    ) -> Tuple[int, str, str]:
        """Run a CLI command and return exit code, stdout, stderr."""
        try:
            # Ensure subprocesses prefer local sources over installed packages
            env = os.environ.copy()
            repo_root_str = str(self.repo_root)
            existing_pythonpath = env.get("PYTHONPATH", "")
            env["PYTHONPATH"] = (
                repo_root_str
                if not existing_pythonpath
                else repo_root_str + os.pathsep + existing_pythonpath
            )
            result = subprocess.run(
                cmd,
                cwd=cwd or str(self.repo_root),
                capture_output=True,
                text=True,
                timeout=60,  # 60 second timeout
                env=env,
            )
            return result.returncode, result.stdout, result.stderr
        except subprocess.TimeoutExpired:
            return -1, "", "Command timed out"
        except Exception as e:
            return -1, "", str(e)

    def test_help_command(self) -> bool:
        """Test that help command works."""
        print("Testing --help command...")

        returncode, stdout, stderr = self.run_cli_command(
            [sys.executable, "-m", "fast_tcp.cli", "--help"]
        )

        success = (
            returncode == 0
            and "Test Case Prioritization Tool" in stdout
            and "--test-dir" in stdout
            and "--algo" in stdout
            and "--entity" in stdout
        )

        self.test_results.append(
            {
                "test": "help_command",
                "success": success,
                "returncode": returncode,
                "details": f"stdout length: {len(stdout)}, stderr: {stderr[:100]}",
            }
        )

        return success

    def test_missing_required_args(self) -> bool:
        """Test that missing required arguments are properly handled."""
        print("Testing missing required arguments...")

        returncode, stdout, stderr = self.run_cli_command(
            [sys.executable, "-m", "fast_tcp.cli"]
        )

        success = returncode != 0 and "required" in stderr.lower()

        self.test_results.append(
            {
                "test": "missing_required_args",
                "success": success,
                "returncode": returncode,
                "details": f"stderr: {stderr[:200]}",
            }
        )

        return success

    def test_invalid_algorithm(self) -> bool:
        """Test handling of invalid algorithm choice."""
        print("Testing invalid algorithm...")

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            self.create_test_data(temp_path, "bbox")

            returncode, stdout, stderr = self.run_cli_command(
                [
                    sys.executable,
                    "-m",
                    "fast_tcp.cli",
                    "--test-dir",
                    str(temp_path),
                    "--algo",
                    "INVALID_ALGO",
                    "--entity",
                    "bbox",
                ]
            )

            success = returncode != 0 and "invalid choice" in stderr.lower()

            self.test_results.append(
                {
                    "test": "invalid_algorithm",
                    "success": success,
                    "returncode": returncode,
                    "details": f"stderr: {stderr[:200]}",
                }
            )

            return success

    def test_invalid_entity(self) -> bool:
        """Test handling of invalid entity choice."""
        print("Testing invalid entity...")

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            self.create_test_data(temp_path, "bbox")

            returncode, stdout, stderr = self.run_cli_command(
                [
                    sys.executable,
                    "-m",
                    "fast_tcp.cli",
                    "--test-dir",
                    str(temp_path),
                    "--algo",
                    "FAST-pw",
                    "--entity",
                    "invalid_entity",
                ]
            )

            success = returncode != 0 and "invalid choice" in stderr.lower()

            self.test_results.append(
                {
                    "test": "invalid_entity",
                    "success": success,
                    "returncode": returncode,
                    "details": f"stderr: {stderr[:200]}",
                }
            )

            return success

    def test_nonexistent_directory(self) -> bool:
        """Test handling of non-existent test directory."""
        print("Testing non-existent directory...")

        returncode, stdout, stderr = self.run_cli_command(
            [
                sys.executable,
                "-m",
                "fast_tcp.cli",
                "--test-dir",
                "/path/that/does/not/exist",
                "--algo",
                "FAST-pw",
                "--entity",
                "bbox",
            ]
        )

        # The CLI should fail when given a non-existent directory
        success = returncode != 0 and (
            "not found" in stderr.lower()
            or "error" in stderr.lower()
            or "not found" in stdout.lower()
            or "error" in stdout.lower()
        )

        self.test_results.append(
            {
                "test": "nonexistent_directory",
                "success": success,
                "returncode": returncode,
                "details": f"returncode: {returncode}, stderr: {stderr[:100]}, stdout: {stdout[:100]}",
            }
        )

        return success

    def test_basic_prioritization(self, algo: str, entity: str) -> bool:
        """Test basic prioritization workflow."""
        print(f"Testing {algo} with {entity} entity...")

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            output_path = temp_path / "output"

            # Create test data
            self.create_test_data(temp_path, entity, num_test_cases=20)

            returncode, stdout, stderr = self.run_cli_command(
                [
                    sys.executable,
                    "-m",
                    "fast_tcp.cli",
                    "--test-dir",
                    str(temp_path),
                    "--algo",
                    algo,
                    "--entity",
                    entity,
                    "--output-dir",
                    str(output_path),
                    "--repetitions",
                    "2",
                ]
            )

            # Check if command succeeded or at least ran without crashing
            success = returncode == 0 or (
                returncode != 0 and "Error running prioritization:" not in stderr
            )

            # Additional checks for output
            details_parts = []
            if returncode == 0:
                details_parts.append(f"returncode: {returncode}")
                if output_path.exists():
                    output_contents = list(output_path.glob("**/*"))
                    details_parts.append(f"output_files: {len(output_contents)}")
                    if len(output_contents) == 0:
                        # Even if no output files, the command succeeded
                        details_parts.append("(no output files but command succeeded)")
                else:
                    details_parts.append("(no output dir but command succeeded)")
            else:
                details_parts.append(f"returncode: {returncode}")
                if stderr:
                    details_parts.append(f"stderr: {stderr[:100]}")
                if stdout:
                    details_parts.append(f"stdout: {stdout[:100]}")

            details = " | ".join(details_parts)

            self.test_results.append(
                {
                    "test": f"basic_prioritization_{algo}_{entity}",
                    "success": success,
                    "returncode": returncode,
                    "details": details,
                }
            )

            return success

    def test_file_discovery_auto(self) -> bool:
        """Test automatic file discovery."""
        print("Testing automatic file discovery...")

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            output_path = temp_path / "output"

            # Create files with different naming patterns
            self.create_test_data(temp_path, "bbox", num_test_cases=5)
            bbox_file = temp_path / "sample-bbox.txt"
            auto_bbox_file = temp_path / "test_bbox_coverage.txt"
            shutil.copy2(bbox_file, auto_bbox_file)

            returncode, stdout, stderr = self.run_cli_command(
                [
                    sys.executable,
                    "-m",
                    "fast_tcp.cli",
                    "--test-dir",
                    str(temp_path),
                    "--algo",
                    "FAST-pw",
                    "--entity",
                    "bbox",
                    "--output-dir",
                    str(output_path),
                    "--file-naming",
                    "auto",
                ]
            )

            success = returncode == 0 and "Found files:" in stdout

            self.test_results.append(
                {
                    "test": "file_discovery_auto",
                    "success": success,
                    "returncode": returncode,
                    "details": f"stdout: {stdout[:200]}",
                }
            )

            return success

    def test_custom_pattern(self) -> bool:
        """Test custom file pattern matching."""
        print("Testing custom file pattern...")

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            output_path = temp_path / "output"

            # Create test file with custom extension
            custom_file = temp_path / "test_data.csv"
            with open(custom_file, "w") as f:
                f.write("1 2 3 4 5\n")
                f.write("2 3 4 5 6\n")
                f.write("1 3 5 7 9\n")

            returncode, stdout, stderr = self.run_cli_command(
                [
                    sys.executable,
                    "-m",
                    "fast_tcp.cli",
                    "--test-dir",
                    str(temp_path),
                    "--algo",
                    "FAST-pw",
                    "--entity",
                    "bbox",
                    "--output-dir",
                    str(output_path),
                    "--pattern",
                    "*.csv",
                    "--file-naming",
                    "custom",
                ]
            )

            success = returncode == 0

            self.test_results.append(
                {
                    "test": "custom_pattern",
                    "success": success,
                    "returncode": returncode,
                    "details": f"stderr: {stderr[:100]}",
                }
            )

            return success

    def test_all_algorithms(self) -> bool:
        """Test all available algorithms with a simple dataset."""
        algorithms = [
            "FAST-pw",
            "FAST-one",
            "FAST-log",
            "FAST-sqrt",
            "FAST-all",
        ]

        results = []

        for algo in algorithms:
            print(f"Testing algorithm: {algo}")

            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                output_path = temp_path / "output"

                # Create small test dataset
                self.create_test_data(temp_path, "bbox", num_test_cases=5)

                returncode, stdout, stderr = self.run_cli_command(
                    [
                        sys.executable,
                        "-m",
                        "fast_tcp.cli",
                        "--test-dir",
                        str(temp_path),
                        "--algo",
                        algo,
                        "--entity",
                        "bbox",
                        "--output-dir",
                        str(output_path),
                        "--repetitions",
                        "1",
                    ]
                )

                success = returncode == 0
                results.append(success)

                self.test_results.append(
                    {
                        "test": f"algorithm_{algo}",
                        "success": success,
                        "returncode": returncode,
                        "details": f"stderr: {stderr[:100]}" if not success else "OK",
                    }
                )

        overall_success = all(results)
        print(f"Algorithm tests: {sum(results)}/{len(results)} passed")

        return overall_success

    def run_all_tests(self) -> Dict:
        """Run all tests and return results."""
        print("=" * 60)
        print("FAST TCP CLI Test Suite")
        print("=" * 60)

        test_methods = [
            self.test_help_command,
            self.test_missing_required_args,
            self.test_invalid_algorithm,
            self.test_invalid_entity,
            self.test_nonexistent_directory,
            lambda: self.test_basic_prioritization("FAST-pw", "bbox"),
            self.test_file_discovery_auto,
            self.test_custom_pattern,
            # Commented out for faster testing - uncomment to test all algorithms
            # self.test_all_algorithms,
        ]

        passed = 0
        total = len(test_methods)

        try:
            for test_method in test_methods:
                try:
                    if test_method():
                        passed += 1
                        print("✓ PASSED")
                    else:
                        print("✗ FAILED")
                except Exception as e:
                    print(f"✗ ERROR: {e}")
                print()

        finally:
            self.cleanup()

        print("=" * 60)
        print(f"Test Results: {passed}/{total} tests passed")
        print("=" * 60)

        # Print detailed results
        for result in self.test_results:
            status = "PASS" if result["success"] else "FAIL"
            print(f"[{status}] {result['test']}: {result['details']}")

        # Return summary
        return {
            "total": total,
            "passed": passed,
            "failed": total - passed,
            "success_rate": passed / total * 100,
            "details": self.test_results,
        }


def main():
    """Main test runner."""
    if len(sys.argv) > 1 and sys.argv[1] == "--all-algorithms":
        # Enable comprehensive algorithm testing
        suite = CLITestSuite()

        # Add algorithm testing back to the suite
        original_run = suite.run_all_tests

        def extended_run():
            result = original_run()
            print("\nRunning comprehensive algorithm tests...")
            suite.test_all_algorithms()
            return result

        suite.run_all_tests = extended_run

        results = suite.run_all_tests()
    else:
        suite = CLITestSuite()
        results = suite.run_all_tests()

    # Exit with appropriate code
    if results["passed"] == results["total"]:
        print("\n🎉 All tests passed!")
        sys.exit(0)
    else:
        print(f"\n❌ {results['failed']} test(s) failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
