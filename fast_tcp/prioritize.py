"""Prioritization driver that bridges the legacy CLI with the new FAST module."""

from __future__ import annotations

import os
import pickle
import shutil
import sys
import time
from typing import Iterable, List, Sequence, Set, Tuple

import xxhash

from . import fast


usage = """USAGE: python py/prioritize.py <dataset> <entity> <algorithm> <repetitions>
OPTIONS:
  <dataset>: test suite to prioritize.
    options: flex_v3, grep_v3, gzip_v1, make_v1, sed_v6, closure_v0, lang_v0, math_v0, chart_v0, time_v0
  <entity>: BB prioritization.
    options: bbox
  <algorithm>: algorithm used for prioritization.
    options: FAST-pw, FAST-one, FAST-log, FAST-sqrt, FAST-all
  <repetitions>: number of prioritization to compute.
    options: positive integer value, e.g. 50
"""


ALLOWED_ALGOS = {
    "FAST-pw",
    "FAST-one",
    "FAST-log",
    "FAST-sqrt",
    "FAST-all",
}


def _configure_fast(
    *,
    algo: str,
    k: int,
    r: int,
    b: int,
    budget: int,
    signature_dir: str,
) -> None:
    """Mutate the global FAST module parameters for the current run."""

    fast.DEFAULT_ALG = algo
    fast.DEFAULT_K = k
    fast.DEFAULT_R = r
    fast.DEFAULT_B = b
    fast.DEFAULT_H = fast.DEFAULT_R * fast.DEFAULT_B
    fast.DEFAULT_BUDGET = budget
    fast.SIGNATURE_FOLDER = signature_dir


def _load_test_suite(path: str) -> List[Tuple[str, str]]:
    """Load a plain-text test suite where each line is a test command."""

    tests: List[Tuple[str, str]] = []
    with open(path, "r", encoding="utf-8") as fh:
        for idx, line in enumerate(fh, start=1):
            content = line.rstrip("\r\n")
            if not content:
                continue
            tests.append((str(idx), content))
    if not tests:
        raise ValueError(f"No tests found in {path}")
    return tests


def partition_test_suite(
    new_test_suite: List[Tuple[str, str]], old_test_suite: List[Tuple[str, str]]
) -> Tuple[Set[str], Set[str], Set[str]]:
    """
    Partition tests into unchanged, deleted, and new sets by comparing
    (test_id, hash(test)) pairs from old and new test suites.

    Parameters:
    - new_test_suite: iterable of (t_id, test_content) representing the current suite.
    - old_test_suite: iterable of (t_id, test_content) representing the previous suite.

    Returns:
    - new_tests: set of t_id present in new suite but not in old suite (or changed).
    - old_tests: set of t_id present in both suites with identical content.
    - del_tests: set of t_id present in old suite but not in new suite (or changed).
    """
    new_hash = {
        (t_id, xxhash.xxh64_intdigest(t.encode("utf8"))) for t_id, t in new_test_suite
    }
    old_hash = {
        (t_id, xxhash.xxh64_intdigest(t.encode("utf8"))) for t_id, t in old_test_suite
    }

    new_tests = {t_id for t_id, _ in new_hash - old_hash}  # new or changed
    old_tests = {t_id for t_id, _ in old_hash & new_hash}  # unchanged
    del_tests = {t_id for t_id, _ in old_hash - new_hash}  # deleted or changed

    return new_tests, old_tests, del_tests


def run_blackbox_file(
    *,
    algo: str,
    input_file: str,
    signature_dir: str,
    k: int,
    r: int,
    b: int,
    budget: int,
    old_suite: Iterable[Tuple[str, str]] | None = None,
    debug: bool = False,
) -> Tuple[float, float, float, List[int]]:
    """Run one prioritization repetition and collect timing metrics."""

    new_suite = _load_test_suite(input_file)
    old_suite = list(old_suite or [])

    if os.path.exists(signature_dir):
        shutil.rmtree(signature_dir)
    os.makedirs(signature_dir, exist_ok=True)

    _configure_fast(
        algo=algo,
        k=k,
        r=r,
        b=b,
        budget=budget,
        signature_dir=signature_dir,
    )

    if debug:
        print(f"[DEBUG] Starting partition phase")

    start_partition = time.perf_counter()
    new_tests, old_tests, del_tests = partition_test_suite(new_suite, old_suite)
    partition_time = time.perf_counter() - start_partition
    if debug:
        print(f"[DEBUG]   Partition time: {partition_time:.4f}s")

    if debug:
        print(f"[DEBUG] Starting preparation phase (k={k}, r={r}, b={b})")

    start_prep = time.perf_counter()
    fast.preparation(new_suite, del_tests)
    prep_time = time.perf_counter() - start_prep

    if debug:
        print(f"[DEBUG]   Preparation time: {prep_time:.4f}s")
        print(f"[DEBUG] Starting prioritization phase")

    start_prio = time.perf_counter()
    prioritized = fast.prioritization(
        new_suite, new_tests=new_tests, old_tests=old_tests
    )
    prio_time = time.perf_counter() - start_prio

    if debug:
        print(f"[DEBUG]   Prioritization time: {prio_time:.4f}s")

    # The FAST module returns a list of IDs in priority order.
    prioritized_ids = [int(t_id) for t_id in prioritized]

    return partition_time, prep_time, prio_time, prioritized_ids


def bbox_prioritization(
    name: str,
    prog: str,
    version: str,
    entity: str,
    *,
    k: int,
    r: int,
    b: int,
    repeats: int,
    budget: int,
    debug: bool = False,
) -> None:
    """Prioritize the specified dataset.

    This function expects the dataset to be in the format <prog>_<version>/<prog>-<entity>.txt.
    The dataset is expected to be a plain-text file where each line is a test command.
    It will prioritize the dataset and save the results to the output directory.

    The output directory is in the format <prog>_<version>/prioritized.

    The prioritized results are saved in the format <name>-<entity>-<run>.pickle.

    The output is saved in the format <name>-<entity>.tsv.

    Args:
        name: The name of the algorithm.
        prog: The name of the program.
        version: The version of the program.
        entity: The entity to prioritize.
        k: The k parameter.
        r: The r parameter.
        b: The b parameter.
        budget: The budget.
        debug: Whether to print debug information.
        repeats: The number of repetitions.

    Returns:
        None

    Raises:
        FileNotFoundError: If the input file is not found.
    """

    input_file = f"input/{prog}_{version}/{prog}-{entity}.txt"
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    output_root = f"output/{prog}_{version}"
    prioritized_dir = os.path.join(output_root, "prioritized")
    os.makedirs(prioritized_dir, exist_ok=True)

    signature_base = os.path.join(output_root, "signatures", name)

    prep_times: List[float] = []
    prio_times: List[float] = []
    partition_times: List[float] = []

    print(f"Running {name} prioritization on {entity} ({repeats} repetition(s))")
    for run in range(repeats):
        sig_dir = os.path.join(signature_base, f"run_{run + 1}")
        print(f"  Repetition {run + 1}/{repeats}")
        partition_time, prep_time, prio_time, prioritized = run_blackbox_file(
            algo=name,
            input_file=input_file,
            signature_dir=sig_dir,
            k=k,
            r=r,
            b=b,
            budget=budget,
            debug=debug,
        )
        partition_times.append(partition_time)
        prep_times.append(prep_time)
        prio_times.append(prio_time)

        write_prioritization(prioritized_dir, name, entity, run, prioritized)
        print(
            f"    Partition time: {partition_time:.4f}s | Preparation time: {prep_time:.4f}s | Prioritization time: {prio_time:.4f}s"
        )

    write_output(output_root, entity, name, partition_times, prep_times, prio_times)


def write_prioritization(
    path: str, name: str, entity: str, run: int, prioritization: Iterable[int]
) -> None:
    os.makedirs(path, exist_ok=True)
    fout = os.path.join(path, f"{name}-{entity}-{run + 1}.pickle")
    with open(fout, "wb") as fh:
        pickle.dump(list(prioritization), fh)


def write_output(
    outpath: str,
    entity: str,
    algo_name: str,
    partition_times: Sequence[float],
    prep_times: Sequence[float],
    prio_times: Sequence[float],
) -> None:
    os.makedirs(outpath, exist_ok=True)
    fileout = os.path.join(outpath, f"{algo_name}-{entity}.tsv")
    with open(fileout, "w", encoding="utf-8") as fout:
        fout.write("SignatureTime\tPrioritizationTime\n")
        for st, pt, pt2 in zip(partition_times, prep_times, prio_times):
            fout.write(f"{st}\t{pt}\t{pt2}\n")


def main() -> None:
    """Main function that can be called from CLI or directly."""

    if len(sys.argv) != 5:
        print("Wrong input.")
        print(usage)
        raise SystemExit(1)

    prog_v, entity, algname, repeats_str = sys.argv[1:]

    if "_" not in prog_v:
        print(
            "<dataset> input should be in format 'name_version' (e.g., 'flex_v3', 'custom_v1')."
        )
        print(usage)
        raise SystemExit(1)

    if entity != "bbox":
        print("<entity> input incorrect.")
        print(usage)
        raise SystemExit(1)

    if algname not in ALLOWED_ALGOS:
        print("<algorithm> input incorrect.")
        print(usage)
        raise SystemExit(1)

    try:
        repeats = int(repeats_str)
    except ValueError:
        print("<repetitions> must be a positive integer.")
        print(usage)
        raise SystemExit(1)

    if repeats <= 0:
        print("<repetitions> must be positive.")
        print(usage)
        raise SystemExit(1)

    prog, version = prog_v.rsplit("_", 1)

    # FAST parameters mirror the defaults from the configuration file.
    k = fast.DEFAULT_K
    r = fast.DEFAULT_R
    b = fast.DEFAULT_B
    budget = fast.DEFAULT_BUDGET

    # Check for debug mode via environment variable
    debug = os.environ.get("FAST_TCP_DEBUG", "").lower() in ("1", "true", "yes")

    bbox_prioritization(
        algname,
        prog,
        version,
        entity,
        k=k,
        r=r,
        b=b,
        repeats=repeats,
        budget=budget,
        debug=debug,
    )


if __name__ == "__main__":
    main()
