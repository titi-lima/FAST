"""Prioritization driver that bridges the legacy CLI with the new FAST module."""

from __future__ import annotations

import io
import os
import pickle
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Set, Tuple

import xxhash

from . import fast
from . import snapshot_cache as _snapshot_cache
from .snapshot_cache import (
    detect_changes_preparation,
    load_file_from_snapshot,
    snapshot_prioritization,
)


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


def _debug_print(enabled: bool, message: str) -> None:
    if enabled:
        print(f"[DEBUG] {message}")


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


def _parse_test_suite_lines(lines: Iterable[str]) -> List[Tuple[str, str]]:
    tests: List[Tuple[str, str]] = []
    for idx, line in enumerate(lines, start=1):
        content = line.rstrip("\r\n")
        if not content:
            continue
        tests.append((str(idx), content))
    return tests


def _suite_cache_base(root: Path) -> Path:
    return root / ".fast" / "cache" / "suites"


def _suite_cache_path(root: Path, relative_path: Path) -> Path:
    return _suite_cache_base(root) / relative_path


def _load_cached_suite(root: Path, relative_path: Path) -> List[Tuple[str, str]] | None:
    cache_file = _suite_cache_path(root, relative_path)
    if not cache_file.exists():
        return None
    try:
        with open(cache_file, "r", encoding="utf-8") as fh:
            return _parse_test_suite_lines(fh)
    except OSError:
        return None


def _store_cached_suite(
    root: Path, relative_path: Path, suite: Sequence[Tuple[str, str]]
) -> None:
    cache_file = _suite_cache_path(root, relative_path)
    try:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_file, "w", encoding="utf-8") as fh:
            for _, test in suite:
                fh.write(test + "\n")
    except OSError:
        pass


_TEST_FILE_EXTS = (
    ".py",
    ".pyw",
    ".js",
    ".cjs",
    ".mjs",
    ".ts",
    ".tsx",
    ".jsx",
    ".java",
    ".kt",
    ".go",
    ".rb",
    ".php",
    ".rs",
)


def _extract_dependency_paths(command: str, project_root: Path) -> Set[Path]:
    candidates: Set[Path] = set()

    def _consider(token: str) -> None:
        token = token.strip()
        if not token:
            return
        if "[" in token and token.endswith("]"):
            token = token[: token.index("[")]
        token = token.replace("\\", "/")
        if not any(token.endswith(ext) for ext in _TEST_FILE_EXTS):
            if "/" not in token:
                return
        path = Path(token)
        candidate = path if path.is_absolute() else (project_root / path)
        try:
            rel = candidate.resolve(strict=False).relative_to(project_root)
        except ValueError:
            return
        candidates.add(rel)

    segments = command.split()
    if segments:
        primary = segments[0]
        if "::" in primary:
            _consider(primary.split("::", 1)[0])
        else:
            _consider(primary)

    for token in re.split(r"[,\s:]+", command.replace("::", " ")):
        _consider(token)

    return candidates


def _filter_suite_by_changed_files(
    suite: List[Tuple[str, str]],
    changed_files: Set[str],
    project_root: Path,
    debug: bool,
) -> Tuple[List[Tuple[str, str]], Set[str]]:
    if not suite or not changed_files:
        return suite, set()

    changed_paths = {Path(path) for path in changed_files}
    kept: List[Tuple[str, str]] = []
    dropped_ids: Set[str] = set()

    for test_id, command in suite:
        deps = _extract_dependency_paths(command, project_root)
        if any(dep in changed_paths for dep in deps):
            dropped_ids.add(test_id)
        else:
            kept.append((test_id, command))

    if debug and dropped_ids:
        ordered = sorted(
            dropped_ids,
            key=int if all(t.isdigit() for t in dropped_ids) else str,
        )
        _debug_print(debug, "Invalidated cached tests: " + ", ".join(ordered))

    return kept, dropped_ids


def _load_test_suite(path: str) -> List[Tuple[str, str]]:
    """Load a plain-text test suite where each line is a test command."""

    with open(path, "r", encoding="utf-8") as fh:
        tests = _parse_test_suite_lines(fh)
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


def _resolve_project_root(
    input_file: str, explicit_root: Optional[str]
) -> Optional[Path]:
    if explicit_root:
        candidate = Path(explicit_root).expanduser().resolve()
        if (candidate / ".fast").is_dir():
            return candidate

    env_root = os.environ.get("FAST_TCP_PROJECT_ROOT")
    if env_root:
        candidate = Path(env_root).expanduser().resolve()
        if (candidate / ".fast").is_dir():
            return candidate

    path = Path(input_file).expanduser().resolve()
    for candidate in (path.parent, *path.parents):
        if (candidate / ".fast").is_dir():
            return candidate
    return None


def run_blackbox_file(
    *,
    algo: str,
    input_file: str,
    signature_dir: str,
    k: int,
    r: int,
    b: int,
    budget: int,
    project_root: str | None = None,
    old_suite: Iterable[Tuple[str, str]] | None = None,
    debug: bool = False,
) -> Tuple[float, float, float, List[int]]:
    """Run one prioritization repetition and collect timing metrics."""

    new_suite = _load_test_suite(input_file)

    snapshot_root: Optional[Path] = None
    snapshot_pending = False
    snapshot_should_persist = False
    previous_suite: List[Tuple[str, str]] = []
    relative_input_path: Optional[Path] = None

    if old_suite is not None:
        previous_suite = list(old_suite)
        _debug_print(
            debug,
            f"Using caller-provided previous suite ({len(previous_suite)} tests)",
        )
    else:
        snapshot_root = _resolve_project_root(input_file, project_root)
        if snapshot_root:
            _debug_print(debug, f"Detected snapshot root at {snapshot_root}")
            try:
                snapshot_diff = detect_changes_preparation(snapshot_root)
            except Exception as exc:
                _debug_print(debug, f"Snapshot diff failed: {exc}")
                _snapshot_cache.discard_pending_snapshot(snapshot_root)
                snapshot_root = None
            else:
                _debug_print(
                    debug,
                    (
                        "Snapshot diff resolved: "
                        f"old_tree={snapshot_diff.old_tree or '<none>'} "
                        f"new_tree={snapshot_diff.new_tree}"
                    ),
                )
                _debug_print(
                    debug,
                    (
                        "Filesystem changes: "
                        f"+{len(snapshot_diff.added)} "
                        f"~{len(snapshot_diff.modified)} "
                        f"-{len(snapshot_diff.deleted)}"
                    ),
                )
                if debug:
                    for path in sorted(snapshot_diff.added):
                        _debug_print(debug, f"  [ADDED] {path}")
                    for path in sorted(snapshot_diff.modified):
                        _debug_print(debug, f"  [MODIFIED] {path}")
                    for path in sorted(snapshot_diff.deleted):
                        _debug_print(debug, f"  [DELETED] {path}")
                    if (
                        not snapshot_diff.added
                        and not snapshot_diff.modified
                        and not snapshot_diff.deleted
                    ):
                        _debug_print(debug, "  (no tracked file changes)")

                rel_path: Optional[Path] = None
                try:
                    rel_path = (
                        Path(input_file)
                        .expanduser()
                        .resolve()
                        .relative_to(snapshot_root)
                    )
                    relative_input_path = rel_path
                except ValueError:
                    rel_path = None

                if rel_path is not None and snapshot_diff.old_tree:
                    blob = load_file_from_snapshot(
                        snapshot_root, rel_path.as_posix(), tree=snapshot_diff.old_tree
                    )
                    if blob is not None:
                        try:
                            previous_suite = _parse_test_suite_lines(
                                io.StringIO(blob.decode("utf-8"))
                            )
                            _debug_print(
                                debug,
                                (
                                    "Recovered previous test suite from snapshot "
                                    f"({len(previous_suite)} tests)"
                                ),
                            )
                        except UnicodeDecodeError as exc:
                            _debug_print(
                                debug,
                                f"Failed to decode snapshot suite for {rel_path}: {exc}",
                            )
                    else:
                        _debug_print(
                            debug,
                            "No previous snapshot available for input file; treating as new suite",
                        )
                if not previous_suite and rel_path is not None:
                    cached_suite = _load_cached_suite(snapshot_root, rel_path)
                    if cached_suite:
                        previous_suite = cached_suite
                        _debug_print(
                            debug,
                            (
                                "Recovered previous test suite from cache "
                                f"({len(previous_suite)} tests)"
                            ),
                        )
                if previous_suite:
                    changed_files = (
                        snapshot_diff.added
                        | snapshot_diff.modified
                        | snapshot_diff.deleted
                    )
                    previous_suite, dropped_ids = _filter_suite_by_changed_files(
                        previous_suite, changed_files, snapshot_root, debug
                    )
                    if debug and dropped_ids:
                        for test_id in sorted(
                            dropped_ids,
                            key=int if all(t.isdigit() for t in dropped_ids) else str,
                        ):
                            _debug_print(
                                debug,
                                f"  -> Invalidated cached test {test_id} (dependency changed)",
                            )
                snapshot_pending = True
        else:
            _debug_print(
                debug,
                "Snapshot root could not be determined; running without cached suite",
            )

    old_suite = previous_suite

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

    _debug_print(debug, "Starting partition phase")

    start_partition = time.perf_counter()
    new_tests, old_tests, del_tests = partition_test_suite(new_suite, old_suite)
    partition_time = time.perf_counter() - start_partition
    if debug:
        _debug_print(
            debug,
            (
                f"Partition summary: {len(new_tests)} new / "
                f"{len(old_tests)} unchanged / {len(del_tests)} deleted"
            ),
        )
        if new_tests:
            _debug_print(
                debug,
                "  New tests: "
                + ", ".join(
                    sorted(
                        new_tests,
                        key=int if all(t.isdigit() for t in new_tests) else str,
                    )
                ),
            )
        if del_tests:
            _debug_print(
                debug,
                "  Deleted tests: "
                + ", ".join(
                    sorted(
                        del_tests,
                        key=int if all(t.isdigit() for t in del_tests) else str,
                    )
                ),
            )
        if not new_tests and not del_tests:
            _debug_print(debug, "  No test-level changes detected")
        _debug_print(debug, f"  Partition time: {partition_time:.4f}s")

    _debug_print(debug, f"Starting preparation phase (k={k}, r={r}, b={b})")

    start_prep = time.perf_counter()
    fast.preparation(new_suite, del_tests)
    prep_time = time.perf_counter() - start_prep

    if debug:
        _debug_print(debug, f"  Preparation time: {prep_time:.4f}s")
        _debug_print(debug, "Starting prioritization phase")

    start_prio = time.perf_counter()
    prioritized = fast.prioritization(
        new_suite, new_tests=new_tests, old_tests=old_tests
    )
    prio_time = time.perf_counter() - start_prio

    if debug:
        _debug_print(debug, f"  Prioritization time: {prio_time:.4f}s")
        if prioritized:
            preview = ", ".join(prioritized[:10])
            if len(prioritized) > 10:
                preview += ", ..."
            _debug_print(
                debug,
                f"  Prioritized order preview ({len(prioritized)} total): {preview}",
            )
        else:
            _debug_print(debug, "  Prioritized order is empty")

    # The FAST module returns a list of IDs in priority order.
    prioritized_ids = [int(t_id) for t_id in prioritized]

    try:
        snapshot_should_persist = snapshot_pending
        return partition_time, prep_time, prio_time, prioritized_ids
    finally:
        if snapshot_root:
            if snapshot_should_persist:
                try:
                    snapshot_prioritization(snapshot_root)
                except Exception as exc:
                    _debug_print(debug, f"Snapshot persistence failed: {exc}")
            elif snapshot_pending:
                _snapshot_cache.discard_pending_snapshot(snapshot_root)
            if relative_input_path is not None:
                _store_cached_suite(snapshot_root, relative_input_path, new_suite)
                _debug_print(
                    debug,
                    f"Stored current test suite in cache for {relative_input_path.as_posix()}",
                )


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
    project_root: str,
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
        print(f"  Project root: {project_root}")
        partition_time, prep_time, prio_time, prioritized = run_blackbox_file(
            algo=name,
            input_file=input_file,
            signature_dir=sig_dir,
            k=k,
            r=r,
            b=b,
            budget=budget,
            debug=debug,
            project_root=project_root,
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
        project_root=Path(sys.argv[1]).parent.as_posix(),
        k=k,
        r=r,
        b=b,
        repeats=repeats,
        budget=budget,
        debug=debug,
    )


if __name__ == "__main__":
    main()
