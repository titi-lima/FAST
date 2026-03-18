#!/usr/bin/env python3
"""Profiling benchmark for fast-tcp to identify performance bottlenecks.

This script creates a synthetic test suite and profiles the fast-tcp
prioritization pipeline to identify the biggest time consumers.

Usage:
    python profile_benchmark.py [--tests N] [--cprofile] [--trace]
"""

from __future__ import annotations

import argparse
import cProfile
import io
import os
import pstats
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fast_tcp import fast
from fast_tcp import prioritize
from fast_tcp import snapshot_cache


@dataclass
class TimingResult:
    """Store timing breakdown for a single run."""
    total_tests: int
    snapshot_detect_time: float
    partition_time: float
    prep_time: float
    prio_time: float
    total_time: float
    
    def __str__(self) -> str:
        return (
            f"Tests: {self.total_tests}\n"
            f"  Snapshot detection: {self.snapshot_detect_time:.4f}s\n"
            f"  Partition:          {self.partition_time:.4f}s\n"
            f"  Preparation:        {self.prep_time:.4f}s\n"
            f"  Prioritization:     {self.prio_time:.4f}s\n"
            f"  Total:              {self.total_time:.4f}s"
        )


def generate_test_content(file_idx: int, test_idx: int, version: int = 0) -> str:
    """Generate unique test content for k-shingle diversity."""
    unique_hash = hash((file_idx, test_idx, version))
    return f"""def test_func_{file_idx:04d}_{test_idx:03d}():
    # Unique content: file={file_idx} test={test_idx} version={version} hash={unique_hash}
    x = {file_idx} + {test_idx} + {version}
    assert x >= 0
"""


def create_test_suite(playground_dir: Path, num_files: int, tests_per_file: int) -> List[Path]:
    """Create synthetic test files."""
    test_dir = playground_dir / "tests"
    test_dir.mkdir(parents=True, exist_ok=True)
    
    test_files: List[Path] = []
    for file_idx in range(num_files):
        test_file = test_dir / f"test_mod_{file_idx:04d}.py"
        content_lines = [f'"""Test module {file_idx}."""\n']
        for test_idx in range(tests_per_file):
            content_lines.append(generate_test_content(file_idx, test_idx))
        test_file.write_text("\n".join(content_lines))
        test_files.append(test_file)
    
    return test_files


def init_fast_tcp(playground_dir: Path) -> None:
    """Initialize fast-tcp for the project (manual initialization)."""
    # Create the .fast directory structure manually
    fast_dir = playground_dir / ".fast"
    fast_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize snapshot cache
    snapshot_cache.initialize_snapshot_cache(playground_dir)


def collect_tests(playground_dir: Path) -> str:
    """Collect test items and create input file."""
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "tests", "--collect-only", "-q"],
        cwd=playground_dir,
        capture_output=True,
        text=True,
    )
    
    # Parse collected tests
    test_items = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if "::" in line and "test_" in line:
            test_items.append(line)
    
    # Write to input file
    input_file = playground_dir / ".fast" / "test_suite.txt"
    input_file.write_text("\n".join(test_items))
    
    return str(input_file)


def profile_run_with_breakdown(
    playground_dir: Path,
    input_file: str,
    signature_dir: str,
) -> TimingResult:
    """Run prioritization with detailed timing breakdown."""
    
    # Configure FAST module
    prioritize._configure_fast(
        algo="FAST-pw",
        k=5,
        r=1,
        b=10,
        budget=0,  # All tests
        signature_dir=signature_dir,
    )
    
    # Load test suite
    new_suite = prioritize._load_test_suite(input_file)
    total_tests = len(new_suite)
    
    start_total = time.perf_counter()
    
    # Step 1: Snapshot detection (if applicable)
    start_snapshot = time.perf_counter()
    snapshot_root = prioritize._resolve_project_root(input_file, str(playground_dir))
    snapshot_diff = None
    if snapshot_root:
        try:
            snapshot_diff = snapshot_cache.detect_changes_preparation(snapshot_root)
        except Exception:
            pass
    snapshot_detect_time = time.perf_counter() - start_snapshot
    
    # Get old suite (empty for first run)
    old_suite: List[Tuple[str, str]] = []
    
    # Clean signature directory
    if os.path.exists(signature_dir):
        shutil.rmtree(signature_dir)
    os.makedirs(signature_dir, exist_ok=True)
    
    # Step 2: Partition
    start_partition = time.perf_counter()
    new_tests, old_tests, del_tests = prioritize.partition_test_suite(new_suite, old_suite)
    partition_time = time.perf_counter() - start_partition
    
    # Step 3: Preparation (signature generation)
    start_prep = time.perf_counter()
    fast.preparation(new_suite, del_tests)
    prep_time = time.perf_counter() - start_prep
    
    # Step 4: Prioritization
    start_prio = time.perf_counter()
    prioritized = fast.prioritization(new_suite, new_tests=new_tests, old_tests=old_tests)
    prio_time = time.perf_counter() - start_prio
    
    total_time = time.perf_counter() - start_total
    
    # Persist snapshot if needed
    if snapshot_root:
        try:
            snapshot_cache.snapshot_prioritization(snapshot_root)
        except Exception:
            pass
    
    return TimingResult(
        total_tests=total_tests,
        snapshot_detect_time=snapshot_detect_time,
        partition_time=partition_time,
        prep_time=prep_time,
        prio_time=prio_time,
        total_time=total_time,
    )


def run_with_cprofile(
    playground_dir: Path,
    input_file: str,
    signature_dir: str,
) -> Tuple[TimingResult, pstats.Stats]:
    """Run with cProfile and return stats."""
    
    profiler = cProfile.Profile()
    profiler.enable()
    
    result = profile_run_with_breakdown(playground_dir, input_file, signature_dir)
    
    profiler.disable()
    
    # Create stats (print to stdout by default)
    stats = pstats.Stats(profiler)
    stats.sort_stats("cumulative")
    
    return result, stats


def print_profile_summary(stats: pstats.Stats, top_n: int = 30) -> None:
    """Print a summary of the profiling results."""
    print("\n" + "=" * 80)
    print("PROFILING RESULTS - Top Time Consumers")
    print("=" * 80)
    
    # Print stats directly to stdout
    stats.print_stats(top_n)
    
    # Also show callers for the top functions
    print("\n" + "=" * 80)
    print("CALLER BREAKDOWN - Who calls the slow functions?")
    print("=" * 80)
    stats.print_callers(15)


def run_subprocess_trace(playground_dir: Path, input_file: str) -> None:
    """Run with subprocess tracing to see git command overhead."""
    print("\n" + "=" * 80)
    print("SUBPROCESS TRACE - Tracking external command calls")
    print("=" * 80)
    
    # Patch subprocess.run to trace calls
    original_run = subprocess.run
    call_times: List[Tuple[str, float]] = []
    
    def traced_run(*args, **kwargs):
        cmd = args[0] if args else kwargs.get("args", [])
        cmd_str = " ".join(str(c) for c in cmd) if isinstance(cmd, (list, tuple)) else str(cmd)
        start = time.perf_counter()
        result = original_run(*args, **kwargs)
        elapsed = time.perf_counter() - start
        call_times.append((cmd_str, elapsed))
        return result
    
    subprocess.run = traced_run
    
    try:
        # Run the actual benchmark
        signature_dir = str(playground_dir / ".fast" / "signatures_trace")
        if os.path.exists(signature_dir):
            shutil.rmtree(signature_dir)
        
        result = profile_run_with_breakdown(playground_dir, input_file, signature_dir)
        print(f"\nTiming breakdown:\n{result}")
        
        # Analyze subprocess calls
        total_subprocess_time = sum(t for _, t in call_times)
        git_calls = [(c, t) for c, t in call_times if "git" in c.lower()]
        total_git_time = sum(t for _, t in git_calls)
        
        print(f"\n\nSubprocess call summary:")
        print(f"  Total subprocess calls: {len(call_times)}")
        print(f"  Total subprocess time:  {total_subprocess_time:.4f}s")
        print(f"  Git calls:              {len(git_calls)}")
        print(f"  Git time:               {total_git_time:.4f}s")
        print(f"  Git % of total:         {100 * total_git_time / result.total_time:.1f}%")
        
        # Group by command type
        from collections import defaultdict
        cmd_groups = defaultdict(lambda: {"count": 0, "time": 0.0})
        for cmd, t in call_times:
            parts = cmd.split()
            if len(parts) >= 2 and parts[0] == "git":
                key = f"git {parts[1]}"
            elif len(parts) >= 1:
                key = parts[0]
            else:
                key = cmd
            cmd_groups[key]["count"] += 1
            cmd_groups[key]["time"] += t
        
        print("\n  By command type:")
        sorted_groups = sorted(cmd_groups.items(), key=lambda x: x[1]["time"], reverse=True)
        for cmd, data in sorted_groups[:15]:
            print(f"    {cmd:30s} count={data['count']:5d}  time={data['time']:.4f}s")
        
    finally:
        subprocess.run = original_run


def main():
    parser = argparse.ArgumentParser(description="Profile fast-tcp performance")
    parser.add_argument("--tests", type=int, default=200,
                        help="Number of test files to generate (default: 200)")
    parser.add_argument("--tests-per-file", type=int, default=3,
                        help="Tests per file (default: 3)")
    parser.add_argument("--cprofile", action="store_true",
                        help="Run with cProfile for detailed function-level profiling")
    parser.add_argument("--trace", action="store_true",
                        help="Trace subprocess calls to measure git overhead")
    parser.add_argument("--keep", action="store_true",
                        help="Keep the temporary playground directory")
    args = parser.parse_args()
    
    total_tests = args.tests * args.tests_per_file
    print(f"Fast-TCP Performance Profiler")
    print(f"=" * 60)
    print(f"Generating {args.tests} test files with {args.tests_per_file} tests each")
    print(f"Total tests: {total_tests}")
    
    # Create temporary playground
    playground_dir = Path(tempfile.mkdtemp(prefix="fast_tcp_profile_"))
    print(f"Playground: {playground_dir}")
    
    try:
        # Setup
        print("\n[1/4] Creating test suite...")
        create_test_suite(playground_dir, args.tests, args.tests_per_file)
        
        print("[2/4] Initializing fast-tcp...")
        init_fast_tcp(playground_dir)
        
        print("[3/4] Collecting tests...")
        input_file = collect_tests(playground_dir)
        print(f"       Input file: {input_file}")
        
        signature_dir = str(playground_dir / ".fast" / "signatures")
        
        # Run profiling
        print("[4/4] Running profiled benchmark...")
        
        if args.trace:
            run_subprocess_trace(playground_dir, input_file)
        
        if args.cprofile or not args.trace:
            result, stats = run_with_cprofile(playground_dir, input_file, signature_dir)
            
            print("\n" + "=" * 80)
            print("TIMING BREAKDOWN")
            print("=" * 80)
            print(result)
            
            if args.cprofile:
                print_profile_summary(stats)
            else:
                # Basic run - give summary
                overhead_per_test = result.total_time / total_tests * 1000
                print("\n" + "=" * 80)
                print("ANALYSIS SUMMARY")
                print("=" * 80)
                print(f"\nOverhead per test: {overhead_per_test:.2f}ms")
                print("\nRun with --cprofile for detailed function breakdown")
                print("Run with --trace for subprocess call analysis")
        
    finally:
        if args.keep:
            print(f"\nPlayground preserved at: {playground_dir}")
        else:
            shutil.rmtree(playground_dir)
            print(f"\nPlayground cleaned up")


if __name__ == "__main__":
    main()
