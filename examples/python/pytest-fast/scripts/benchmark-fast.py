#!/usr/bin/env python3
"""
Benchmark pytest execution with and without FAST TCP ordering over a target test tree.

It runs pytest twice (baseline and FAST-enabled), captures durations, and prints a summary.
The target is a path (e.g., the generated suite). It supports overriding FAST options.

Usage (from repo root):
  python examples/python/pytest-fast/scripts/benchmark-fast.py \
      --target examples/python/pytest-fast/generated \
      --repeat 3 --fast-algo FAST-pw --fast-r 1 --fast-b 10 --fast-k 5 --fast-budget 0
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path


def run_pytest(cmd: list[str]) -> tuple[int, float]:
    start = time.perf_counter()
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    duration = time.perf_counter() - start
    return proc.returncode, duration


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Benchmark pytest with/without FAST")
    p.add_argument("--target", required=True, help="Path to tests directory or file")
    p.add_argument("--repeat", type=int, default=1, help="Number of repetitions per mode")
    p.add_argument("--pytest-bin", default="pytest", help="Pytest executable")
    p.add_argument("--quiet", action="store_true", help="Run pytest with -q to reduce noise")

    # FAST options
    p.add_argument("--fast-algo", default="FAST-pw")
    p.add_argument("--fast-r", type=int, default=1)
    p.add_argument("--fast-b", type=int, default=10)
    p.add_argument("--fast-k", type=int, default=5)
    p.add_argument("--fast-budget", type=int, default=0)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    target = Path(args.target)
    if not target.exists():
        print(f"Target not found: {target}")
        return 2

    base_cmd = [args.pytest_bin]
    if args.quiet:
        base_cmd.append("-q")
    base_cmd.append(str(target))

    fast_cmd = base_cmd + [
        "--fast-tcp",
        "--fast-tcp-algo",
        args.fast_algo,
        "--fast-tcp-r",
        str(args.fast_r),
        "--fast-tcp-b",
        str(args.fast_b),
        "--fast-tcp-k",
        str(args.fast_k),
        "--fast-tcp-budget",
        str(args.fast_budget),
    ]

    # Baseline
    base_durations: list[float] = []
    for i in range(args.repeat):
        code, dur = run_pytest(base_cmd)
        base_durations.append(dur)
        print(f"baseline run {i+1}/{args.repeat}: rc={code} time={dur:.3f}s")
        if code != 0:
            print("Warning: baseline pytest returned non-zero code")

    # FAST
    fast_durations: list[float] = []
    for i in range(args.repeat):
        code, dur = run_pytest(fast_cmd)
        fast_durations.append(dur)
        print(f"fast     run {i+1}/{args.repeat}: rc={code} time={dur:.3f}s")
        if code != 0:
            print("Warning: FAST pytest returned non-zero code")

    def stats(values: list[float]) -> tuple[float, float, float]:
        if not values:
            return (0.0, 0.0, 0.0)
        return (min(values), sum(values) / len(values), max(values))

    b_min, b_avg, b_max = stats(base_durations)
    f_min, f_avg, f_max = stats(fast_durations)

    print("\nSummary (seconds):")
    print(f"  baseline: min={b_min:.3f} avg={b_avg:.3f} max={b_max:.3f}")
    print(f"  fast    : min={f_min:.3f} avg={f_avg:.3f} max={f_max:.3f}")
    if b_avg > 0:
        delta = b_avg - f_avg
        pct = (delta / b_avg) * 100.0
        print(f"  delta   : {delta:+.3f}s ({pct:+.1f}%) vs baseline (avg)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


