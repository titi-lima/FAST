#!/usr/bin/env python3
"""
Generate a large synthetic pytest suite to stress-test FAST TCP ordering.

This creates many `test_*.py` files under a target directory, each with a
configurable number of tests. Each test has a small sleep to simulate runtime
and some deterministic variations to avoid identical bodies.

Usage:
  python scripts/generate-large-suite.py --out-dir generated --num-files 200 --tests-per-file 10

After generation, run from repo root:
  pytest -q examples/python/pytest-fast/generated
  pytest -q examples/python/pytest-fast/generated --fast-tcp  # if not enabled via pytest.ini
"""

from __future__ import annotations

import argparse
import os
import random
from pathlib import Path


TEMPLATE = """
import time
import math


{tests}
"""


def _make_test_body(test_index: int, sleep_ms: int) -> str:
    sleep_seconds = max(0, sleep_ms) / 1000.0
    # Use math to slightly vary work to avoid identical bytecode across tests
    return (
        f"def test_generated_{test_index:04d}():\n"
        f"    time.sleep({sleep_seconds})\n"
        f"    values = [i * {test_index + 3} % 7 for i in range(10)]\n"
        f"    assert sorted(values)[0] in (0, 1, 2, 3, 4, 5, 6)\n"
    )


def generate_suite(out_dir: Path, num_files: int, tests_per_file: int, min_ms: int, max_ms: int, seed: int) -> None:
    rng = random.Random(seed)
    out_dir.mkdir(parents=True, exist_ok=True)

    test_counter = 0
    for file_index in range(num_files):
        tests_src_lines = []
        for _ in range(tests_per_file):
            duration_ms = rng.randint(min_ms, max_ms)
            tests_src_lines.append(_make_test_body(test_counter, duration_ms))
            test_counter += 1

        file_src = TEMPLATE.format(tests="\n\n".join(tests_src_lines))
        file_path = out_dir / f"test_mod_{file_index:04d}.py"
        file_path.write_text(file_src, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a large synthetic pytest suite")
    parser.add_argument("--out-dir", default="generated", help="Directory to write tests into (relative or absolute)")
    parser.add_argument("--num-files", type=int, default=100, help="Number of test files to create")
    parser.add_argument("--tests-per-file", type=int, default=10, help="Number of tests per file")
    parser.add_argument("--min-ms", type=int, default=5, help="Minimum sleep per test (ms)")
    parser.add_argument("--max-ms", type=int, default=20, help="Maximum sleep per test (ms)")
    parser.add_argument("--seed", type=int, default=1337, help="Random seed for reproducibility")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        # Resolve relative to this example directory to avoid polluting repo root
        script_dir = Path(__file__).resolve().parent
        example_root = script_dir.parent
        out_dir = (example_root / out_dir).resolve()

    generate_suite(
        out_dir=out_dir,
        num_files=args.num_files,
        tests_per_file=args.tests_per_file,
        min_ms=args.min_ms,
        max_ms=args.max_ms,
        seed=args.seed,
    )
    print(f"Generated tests: {args.num_files * args.tests_per_file} in {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


