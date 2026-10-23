"""Run the notebook's 100 MiB, 40-mutation snapshot experiment; save exact bytes."""
from __future__ import annotations

import json
import math
import os
import random
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import argparse
import platform
from common import DEFAULT_OUTPUT, provenance

from fast_tcp import snapshot_cache

PLAYGROUND_ROOT = DEFAULT_OUTPUT / "snapshot/project"
SUITE_DIR = PLAYGROUND_ROOT / "tests"
DELTA_DIR = SUITE_DIR / "delta"

TARGET_SUITE_MB = 100
FILE_CHUNK_MB = 2
FILE_CHUNK_BYTES = FILE_CHUNK_MB * 1024 * 1024
BASE_FILE_COUNT = TARGET_SUITE_MB // FILE_CHUNK_MB

ITERATIONS = 40
MODIFY_RANGE: Tuple[int, int] = (18, 32)
DELTA_ADDITION_RANGE: Tuple[int, int] = (2, 6)
DELTA_FILE_TTL = 6
DELTA_FILE_SIZES = [FILE_CHUNK_BYTES // 4, FILE_CHUNK_BYTES // 2]

DELTA_REGISTRY: Dict[Path, int] = {}

@dataclass
class SnapshotRecord:
    iteration: int
    tag: str
    added: int
    modified: int
    deleted: int
    repo_size_bytes: int
    loose_objects: int
    working_set_bytes: int


def format_bytes(num: int) -> str:
    if num <= 0:
        return "0 B"
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    current = float(num)
    for unit in units:
        if current < 1024 or unit == units[-1]:
            return f"{current:0.1f} {unit}"
        current /= 1024
    return f"{current:0.1f} PiB"


def ensure_playground_clean() -> None:
    if PLAYGROUND_ROOT.exists():
        raise SystemExit("Snapshot project already exists; use a new --output-dir to preserve the earlier run")
    SUITE_DIR.mkdir(parents=True, exist_ok=True)
    # Stop Git discovering the parent FAST repository and its ignored output/.
    subprocess.run(["git", "init", "--quiet", str(PLAYGROUND_ROOT)], check=True)
    DELTA_REGISTRY.clear()


def write_payload(path: Path, label: str, size_bytes: int = FILE_CHUNK_BYTES) -> None:
    block = (label.encode("utf-8") + b"|")
    repeats = (size_bytes // len(block)) + 1
    payload = (block * repeats)[:size_bytes]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def seed_suite() -> None:
    for idx in range(BASE_FILE_COUNT):
        target = SUITE_DIR / f"test_chunk_{idx:04d}.bin"
        write_payload(target, f"seed-{idx}")


def working_set_bytes() -> int:
    total = 0
    for path in PLAYGROUND_ROOT.rglob("*"):
        if path.is_file() and ".fast" not in path.parts and ".git" not in path.parts:
            total += path.stat().st_size
    return total


def snapshot_repo_path() -> Path:
    return PLAYGROUND_ROOT / ".fast" / "snapshot.git"


def snapshot_repo_size_bytes() -> int:
    repo = snapshot_repo_path()
    if not repo.exists():
        return 0
    return sum(p.stat().st_size for p in repo.rglob("*") if p.is_file())


def git_stats() -> Dict[str, int]:
    repo = snapshot_repo_path()
    if not repo.exists():
        return {"count": 0, "size_bytes": 0}
    env = os.environ.copy()
    env["GIT_DIR"] = str(repo)
    proc = subprocess.run(
        ["git", "count-objects", "-v"],
        cwd=PLAYGROUND_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    stats: Dict[str, int] = {}
    for line in proc.stdout.splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        try:
            stats[key.strip()] = int(value.strip().split()[0])
        except ValueError:
            continue
    stats["size_bytes"] = snapshot_repo_size_bytes()
    return stats


def age_delta_files() -> int:
    """Decrease TTL for adhoc files and delete those that expired."""
    expired = 0
    for path in list(DELTA_REGISTRY.keys()):
        DELTA_REGISTRY[path] -= 1
        if DELTA_REGISTRY[path] <= 0:
            if path.exists():
                path.unlink()
            expired += 1
            DELTA_REGISTRY.pop(path, None)
    return expired


def mutate_suite(iteration: int) -> Dict[str, int]:
    min_mod, max_mod = MODIFY_RANGE
    max_mod = min(max_mod, BASE_FILE_COUNT)
    min_mod = min(min_mod, max_mod)
    modifications = random.randint(min_mod, max_mod)
    indices = random.sample(range(BASE_FILE_COUNT), modifications)
    for idx in indices:
        path = SUITE_DIR / f"test_chunk_{idx:04d}.bin"
        write_payload(path, f"iter-{iteration}-file-{idx}")

    expired = age_delta_files()

    additions = random.randint(*DELTA_ADDITION_RANGE)
    DELTA_DIR.mkdir(parents=True, exist_ok=True)
    for add_idx in range(additions):
        size_bytes = random.choice(DELTA_FILE_SIZES)
        salt = random.randint(0, 9999)
        extra_path = DELTA_DIR / f"iter_{iteration:04d}_{add_idx:02d}_{salt:04d}.bin"
        write_payload(extra_path, f"delta-{iteration}-{add_idx}", size_bytes=size_bytes)
        DELTA_REGISTRY[extra_path] = DELTA_FILE_TTL

    return {
        "modified": modifications,
        "delta_added": additions,
        "delta_expired": expired,
    }


def run_snapshot_cycle(iteration: int, tag: str) -> SnapshotRecord:
    diff = snapshot_cache.detect_changes_preparation(PLAYGROUND_ROOT)
    snapshot_cache.snapshot_prioritization(PLAYGROUND_ROOT)
    stats = git_stats()
    record = SnapshotRecord(
        iteration=iteration,
        tag=tag,
        added=len(diff.added),
        modified=len(diff.modified),
        deleted=len(diff.deleted),
        repo_size_bytes=stats.get("size_bytes", 0),
        loose_objects=stats.get("count", 0),
        working_set_bytes=working_set_bytes(),
    )
    return record


def main():
    global PLAYGROUND_ROOT, SUITE_DIR, DELTA_DIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output_dir.resolve() / "snapshot"
    PLAYGROUND_ROOT = output / "project"
    SUITE_DIR = PLAYGROUND_ROOT / "tests"
    DELTA_DIR = SUITE_DIR / "delta"
    random.seed(42)
    ensure_playground_clean()
    seed_suite()
    records = []
    for iteration in range(ITERATIONS + 1):
        if iteration:
            mutation = mutate_suite(iteration)
        record = run_snapshot_cycle(iteration, "seed" if not iteration else f"iter-{iteration}")
        expected = (mutation["delta_added"], mutation["modified"], mutation["delta_expired"]) if iteration else (BASE_FILE_COUNT, 0, 0)
        if (record.added, record.modified, record.deleted) != expected:
            raise RuntimeError(f"Snapshot did not detect the generated changes at iteration {iteration}: {record}")
        records.append(record.__dict__)
        (output / "measurements.json").write_text(json.dumps(records, indent=2))
        print(json.dumps(record.__dict__), flush=True)
    metadata = {**provenance(), "python": platform.python_version(),
                "git": subprocess.check_output(["git", "--version"], text=True).strip(),
                "seed": 42, "iterations": ITERATIONS, "payload": "repeated labels, highly compressible"}
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
