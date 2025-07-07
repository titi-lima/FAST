"""
Command-line interface for FAST: Scalable Similarity-based Test Case Prioritization
"""

import argparse
import sys
import os
from pathlib import Path

from . import prioritize, scalability


def create_parser():
    """Create the main argument parser."""
    parser = argparse.ArgumentParser(
        description="FAST: Scalable Similarity-based Test Case Prioritization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  fast-tcp prioritize flex_v3 bbox FAST-pw 50
  fast-tcp scalability 1000 small FAST-pw
  fast-tcp generate-input 1000 small
  fast-tcp plot-results small prioritization FAST-pw FAST-one
        """,
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Prioritize command
    prioritize_parser = subparsers.add_parser(
        "prioritize",
        help="Run test case prioritization",
        description="Run test case prioritization on specified datasets",
    )
    prioritize_parser.add_argument(
        "dataset",
        choices=[
            "flex_v3",
            "grep_v3",
            "gzip_v1",
            "make_v1",
            "sed_v6",
            "closure_v0",
            "lang_v0",
            "math_v0",
            "chart_v0",
            "time_v0",
        ],
        help="Test suite to prioritize",
    )
    prioritize_parser.add_argument(
        "entity",
        choices=["bbox", "function", "branch", "line"],
        help="BB or WB (function, branch, line) prioritization",
    )
    prioritize_parser.add_argument(
        "algorithm",
        choices=[
            "FAST-pw",
            "FAST-one",
            "FAST-log",
            "FAST-sqrt",
            "FAST-all",
            "STR",
            "I-TSD",
            "ART-D",
            "ART-F",
            "GT",
            "GA",
            "GA-S",
        ],
        help="Algorithm used for prioritization",
    )
    prioritize_parser.add_argument(
        "repetitions", type=int, help="Number of prioritizations to compute"
    )

    # Scalability command
    scalability_parser = subparsers.add_parser(
        "scalability",
        help="Run scalability experiments",
        description="Measure scalability of TCP approaches",
    )
    scalability_parser.add_argument(
        "tssize", type=int, help="Number of test cases in the test suite"
    )
    scalability_parser.add_argument(
        "tcsize",
        choices=["small", "medium", "large"],
        help="Size of the test cases (small, medium, large)",
    )
    scalability_parser.add_argument(
        "algorithm",
        choices=[
            "FAST-pw",
            "FAST-one",
            "FAST-log",
            "FAST-sqrt",
            "FAST-all",
            "STR",
            "I-TSD",
            "ART-D",
            "ART-F",
            "GT",
            "GA",
            "GA-S",
        ],
        help="Algorithm used for prioritization",
    )

    # Generate input command
    generate_parser = subparsers.add_parser(
        "generate-input",
        help="Generate synthetic input for scalability experiments",
        description="Generate synthetic test suites for scalability testing",
    )
    generate_parser.add_argument(
        "test_suite_size", type=int, help="Number of test cases to generate"
    )
    generate_parser.add_argument(
        "test_case_size",
        choices=["small", "medium", "large"],
        help="Size category for test cases",
    )

    # Plot results command
    plot_parser = subparsers.add_parser(
        "plot-results",
        help="Plot scalability results",
        description="Generate plots from scalability experiment results",
    )
    plot_parser.add_argument(
        "test_case_size",
        choices=["small", "medium", "large"],
        help="Test case size category",
    )
    plot_parser.add_argument(
        "time_type",
        choices=["prioritization", "total"],
        help="Type of time measurement to plot",
    )
    plot_parser.add_argument(
        "algorithms",
        nargs="+",
        choices=[
            "FAST-pw",
            "FAST-one",
            "FAST-log",
            "FAST-sqrt",
            "FAST-all",
            "STR",
            "I-TSD",
            "ART-D",
            "ART-F",
            "GT",
            "GA",
            "GA-S",
        ],
        help="Algorithms to include in the plot",
    )

    # Clean command
    clean_parser = subparsers.add_parser(
        "clean",
        help="Clean preprocessed input files",
        description="Remove preprocessed input files for a clean environment",
    )

    return parser


def run_prioritize(args):
    """Run the prioritize command."""
    # Change to the project root directory if needed
    original_cwd = os.getcwd()

    try:
        # Mock sys.argv for the prioritize module
        sys.argv = [
            "prioritize.py",
            args.dataset,
            args.entity,
            args.algorithm,
            str(args.repetitions),
        ]

        # Import and run the prioritize functionality
        from . import prioritize as prio_module

        # Create necessary directories
        os.makedirs("output", exist_ok=True)
        os.makedirs("input", exist_ok=True)

        # Call the main function from prioritize module
        prio_module.main()

    except SystemExit:
        # Handle sys.exit() calls in the original module
        pass
    finally:
        os.chdir(original_cwd)


def run_scalability(args):
    """Run the scalability command."""
    original_cwd = os.getcwd()

    try:
        # Mock sys.argv for the scalability module
        sys.argv = ["scalability.py", str(args.tssize), args.tcsize, args.algorithm]

        # Import and run the scalability functionality
        from . import scalability as scal_module

        # Create necessary directories
        os.makedirs("scalability/input", exist_ok=True)
        os.makedirs("scalability/output", exist_ok=True)

        # Call the main function from scalability module
        scal_module.main()

    except SystemExit:
        # Handle sys.exit() calls in the original module
        pass
    finally:
        os.chdir(original_cwd)


def run_generate_input(args):
    """Run the generate input command."""
    try:
        import subprocess
        import sys

        # Try to run the tools script
        script_path = (
            Path(__file__).parent.parent
            / "tools"
            / "generate-scalability-synthetic-input.py"
        )

        if script_path.exists():
            subprocess.run(
                [
                    sys.executable,
                    str(script_path),
                    str(args.test_suite_size),
                    args.test_case_size,
                ],
                check=True,
            )
        else:
            print(f"Warning: Could not find {script_path}")
            print("Please ensure you're running from the FAST project directory")

    except subprocess.CalledProcessError as e:
        print(f"Error running generate input script: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


def run_plot_results(args):
    """Run the plot results command."""
    try:
        import subprocess
        import sys

        # Try to run the tools script
        script_path = (
            Path(__file__).parent.parent / "tools" / "plot-scalability-results.py"
        )

        if script_path.exists():
            cmd = [
                sys.executable,
                str(script_path),
                args.test_case_size,
                args.time_type,
            ] + args.algorithms
            subprocess.run(cmd, check=True)
        else:
            print(f"Warning: Could not find {script_path}")
            print("Please ensure you're running from the FAST project directory")

    except subprocess.CalledProcessError as e:
        print(f"Error running plot results script: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


def run_clean(args):
    """Run the clean command."""
    try:
        import subprocess
        import sys

        # Try to run the tools script
        script_path = (
            Path(__file__).parent.parent / "tools" / "clean-preprocessed-input.py"
        )

        if script_path.exists():
            subprocess.run([sys.executable, str(script_path)], check=True)
        else:
            print(f"Warning: Could not find {script_path}")
            print("Please ensure you're running from the FAST project directory")

    except subprocess.CalledProcessError as e:
        print(f"Error running clean script: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


def main():
    """Main CLI entry point."""
    parser = create_parser()
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return

    if args.command == "prioritize":
        run_prioritize(args)
    elif args.command == "scalability":
        run_scalability(args)
    elif args.command == "generate-input":
        run_generate_input(args)
    elif args.command == "plot-results":
        run_plot_results(args)
    elif args.command == "clean":
        run_clean(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
