"""Git tree snapshot cache for incremental FAST prioritization.

This module provides utilities to lazily snapshot a project's working tree into a
standalone bare Git repository living under ``.fast/snapshot.git``. The snapshot
is used to detect filesystem changes between prioritization runs without relying
on the user's own VCS history.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import platform
import stat
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

try:
    from . import __version__ as _FAST_TCP_VERSION
except Exception:  # pragma: no cover - fallback when package metadata unavailable
    _FAST_TCP_VERSION = "unknown"


__all__ = [
    "SnapshotDiff",
    "initialize_snapshot_cache",
    "detect_changes_preparation",
    "snapshot_prioritization",
    "discard_pending_snapshot",
    "load_file_from_snapshot",
]

_DEFAULT_IGNORE_DIRS: Tuple[str, ...] = (
    ".git",
    ".fast",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    ".nox",
    ".tox",
    "build",
    "dist",
    "node_modules",
    ".venv",
    "venv",
    ".idea",
    ".vscode",
)
_DEFAULT_IGNORE_FILE_SUFFIXES: Tuple[str, ...] = (
    ".pyc",
    ".pyo",
    ".pyd",
    ".so",
    ".dll",
    ".dylib",
    ".class",
    ".o",
    ".obj",
)
_DEFAULT_IGNORE_FILE_NAMES: Tuple[str, ...] = (".DS_Store",)


@dataclass(frozen=True)
class SnapshotDiff:
    """Represents a diff between the previous and current filesystem snapshot."""

    added: Set[str]
    modified: Set[str]
    deleted: Set[str]
    new_tree: str
    old_tree: Optional[str]


@dataclass
class _PendingSnapshot:
    tree: str
    ignore_hash: str
    mtime_hash: str


_PENDING_TREES: Dict[Path, _PendingSnapshot] = {}


# Cache for mtime hash to avoid recomputation within same run
_MTIME_HASH_CACHE: Dict[Path, Tuple[str, float]] = (
    {}
)  # path -> (hash, computation_time)


class _GitIgnoreHelper:
    """Utility to leverage the project's `.gitignore` patterns when available."""

    def __init__(self, root: Path) -> None:
        self._root = root
        self._cache: Dict[str, bool] = {}
        self._env = dict(os.environ)
        self._env.pop("GIT_DIR", None)
        self._env.pop("GIT_INDEX_FILE", None)
        self._enabled = self._detect_repo()

    def _detect_repo(self) -> bool:
        try:
            proc = subprocess.run(
                ["git", "-C", str(self._root), "rev-parse", "--is-inside-work-tree"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
                env=self._env,
            )
        except FileNotFoundError:  # pragma: no cover - depends on environment
            return False
        if proc.returncode != 0:
            return False
        return proc.stdout.strip().lower() in {b"true", b"1"}

    def is_ignored(self, path: Path, *, is_dir: bool) -> bool:
        if not self._enabled:
            return False
        try:
            rel = path.relative_to(self._root)
        except ValueError:
            return False
        key = rel.as_posix()
        if is_dir and not key.endswith("/"):
            key += "/"
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        try:
            proc = subprocess.run(
                [
                    "git",
                    "-C",
                    str(self._root),
                    "check-ignore",
                    "-q",
                    "--no-index",
                    key,
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
                env=self._env,
            )
        except FileNotFoundError:  # pragma: no cover - depends on environment
            self._enabled = False
            return False
        if proc.returncode not in (0, 1):
            # Disable further gitignore checks if Git errors (e.g., not a repository).
            self._enabled = False
            return False
        ignored = proc.returncode == 0
        self._cache[key] = ignored
        return ignored


def initialize_snapshot_cache(project_root: Path | str) -> None:
    """Ensure the `.fast/snapshot.git` bare repository exists for a project."""

    root = Path(project_root).resolve()
    fast_dir = root / ".fast"
    snapshot_dir = fast_dir / "snapshot.git"
    if snapshot_dir.exists():
        return

    fast_dir.mkdir(parents=True, exist_ok=True)
    manifest_dir = fast_dir / "manifests"
    manifest_dir.mkdir(parents=True, exist_ok=True)

    try:
        subprocess.run(
            ["git", "init", "--bare", str(snapshot_dir)],
            cwd=str(root),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError as exc:  # pragma: no cover - depends on environment
        raise RuntimeError(
            "Git executable not found; install Git to enable snapshot caching."
        ) from exc

    manifest_path = manifest_dir / "latest.json"
    if not manifest_path.exists():
        manifest_path.write_text(
            json.dumps(
                {
                    "created_at": _dt.datetime.now(tz=_dt.timezone.utc).isoformat(),
                    "tree": None,
                    "kind": "none",
                    "tool_version": _FAST_TCP_VERSION,
                    "python": platform.python_version(),
                    "ignore_hash": _compute_ignore_hash(),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )


def detect_changes_preparation(project_root: Path | str) -> SnapshotDiff:
    """Build a temporary tree for `project_root` and diff it against the snapshot.

    The resulting diff lists are returned to the caller. The snapshot ref is not
    mutated here; callers should invoke :func:`snapshot_prioritization_step2`
    after prioritization completes to persist the new tree.
    """

    root = Path(project_root).resolve()
    initialize_snapshot_cache(root)

    # Compute current mtime hash for short-circuit check
    current_mtime_hash = _compute_mtime_hash(root)
    current_ignore_hash = _compute_ignore_hash()

    # Check if we can short-circuit (no changes since last snapshot)
    manifest = _read_latest_manifest(root)
    env = _snapshot_env(root, None)
    old_tree = _resolve_latest_tree(root, env)

    if (
        manifest
        and old_tree is not None
        and manifest.get("mtime_hash") == current_mtime_hash
        and manifest.get("ignore_hash") == current_ignore_hash
    ):
        # No files have changed - reuse previous snapshot
        pending = _PendingSnapshot(
            tree=old_tree,
            ignore_hash=current_ignore_hash,
            mtime_hash=current_mtime_hash,
        )
        _PENDING_TREES[root] = pending
        return SnapshotDiff(
            added=set(),
            modified=set(),
            deleted=set(),
            new_tree=old_tree,
            old_tree=old_tree,
        )

    # Files have changed - need to build new tree
    index_path = _select_index_path(root)
    env = _snapshot_env(root, index_path)

    tracked_files: List[str] = []
    try:
        new_tree = _build_tree(root, env, tracked_files)
    finally:
        _cleanup_index(index_path)

    pending = _PendingSnapshot(
        tree=new_tree,
        ignore_hash=current_ignore_hash,
        mtime_hash=current_mtime_hash,
    )
    _PENDING_TREES[root] = pending

    # Invalidate ignore hash if it changed
    if manifest and manifest.get("ignore_hash") != current_ignore_hash:
        old_tree = None

    if old_tree is None:
        added = set(tracked_files)
        modified: Set[str] = set()
        deleted: Set[str] = set()
    else:
        env_for_diff = _snapshot_env(root, None)
        diff_proc = _run_git(
            ["diff-tree", "--no-commit-id", "--name-status", "-r", old_tree, new_tree],
            cwd=root,
            env=env_for_diff,
            check=True,
        )
        added, modified, deleted = _parse_diff_output(diff_proc.stdout.decode("utf-8"))

    return SnapshotDiff(
        added=added,
        modified=modified,
        deleted=deleted,
        new_tree=new_tree,
        old_tree=old_tree,
    )


def snapshot_prioritization(project_root: Path | str) -> str:
    """Persist the tree produced during Step 1 into the snapshot repository."""

    root = Path(project_root).resolve()
    initialize_snapshot_cache(root)

    pending = _PENDING_TREES.pop(root, None)
    write_env = _snapshot_env(root, None)
    if pending is None:
        index_path = _select_index_path(root)
        env = _snapshot_env(root, index_path)
        try:
            tree = _build_tree(root, env, tracked_files=None)
            ignore_hash = _compute_ignore_hash()
            mtime_hash = _compute_mtime_hash(root)
        finally:
            _cleanup_index(index_path)
    else:
        tree = pending.tree
        ignore_hash = pending.ignore_hash
        mtime_hash = pending.mtime_hash
        env = write_env

    previous_tree = _resolve_latest_tree(root, write_env)
    tree_changed = previous_tree != tree

    if tree_changed:
        _update_snapshot_ref(root, write_env, tree)

    _write_manifest(root, tree, ignore_hash, mtime_hash)

    # Clear mtime cache after snapshot is persisted
    _clear_mtime_cache(root)

    if tree_changed:
        _prune_snapshot_objects(root, write_env)
    return tree


def discard_pending_snapshot(project_root: Path | str) -> None:
    """Drop any staged tree for ``project_root`` without updating refs."""

    root = Path(project_root).resolve()
    _PENDING_TREES.pop(root, None)


def load_file_from_snapshot(
    project_root: Path | str, relative_path: str, tree: Optional[str] = None
) -> Optional[bytes]:
    """Load a file from the cached snapshot tree."""

    root = Path(project_root).resolve()
    snapshot_dir = root / ".fast" / "snapshot.git"
    if not snapshot_dir.exists():
        return None

    env = _snapshot_env(root, None)
    if tree is None:
        tree = _resolve_latest_tree(root, env)
        if tree is None:
            return None

    proc = _run_git(
        ["show", f"{tree}:{relative_path}"],
        cwd=root,
        env=env,
        check=False,
    )
    if proc.returncode != 0:
        return None
    return proc.stdout


def _select_index_path(root: Path) -> Path:
    pid = os.getpid()
    return root / ".fast" / f"tmp_index.{pid}"


def _cleanup_index(index_path: Optional[Path]) -> None:
    if index_path and index_path.exists():
        try:
            index_path.unlink()
        except FileNotFoundError:
            pass


def _snapshot_env(root: Path, index_path: Optional[Path]) -> Dict[str, str]:
    env = dict(os.environ)
    env["GIT_DIR"] = str(root / ".fast" / "snapshot.git")
    if index_path is not None:
        env["GIT_INDEX_FILE"] = str(index_path)
    elif "GIT_INDEX_FILE" in env:
        env.pop("GIT_INDEX_FILE")
    return env


def _compute_ignore_hash() -> str:
    payload = "\n".join(
        list(_DEFAULT_IGNORE_DIRS)
        + list(_DEFAULT_IGNORE_FILE_SUFFIXES)
        + list(_DEFAULT_IGNORE_FILE_NAMES)
    )
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    return digest


def _compute_mtime_hash(root: Path) -> str:
    """Compute a hash of all tracked file paths and their mtimes.

    This is used as a fast short-circuit check: if the mtime hash hasn't changed
    since the last snapshot, we can skip all git operations and reuse the
    previous snapshot tree.

    The hash includes:
    - File paths (sorted for determinism)
    - File mtimes (nanosecond precision where available)
    - File sizes (to catch edge cases where mtime doesn't change)
    """
    # Check cache first (valid for current process only)
    cached = _MTIME_HASH_CACHE.get(root)
    if cached is not None:
        return cached[0]

    hasher = hashlib.sha256()
    file_data: List[Tuple[str, int, int]] = []  # (path, mtime_ns, size)

    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        current = Path(dirpath)

        # Filter directories (same logic as _collect_files)
        filtered_dirnames: List[str] = []
        for d in dirnames:
            dir_path = current / d
            if not _should_ignore_path(dir_path):
                filtered_dirnames.append(d)
        dirnames[:] = filtered_dirnames

        for name in filenames:
            path = current / name
            if _should_ignore_path(path):
                continue

            rel_path = path.relative_to(root).as_posix()
            try:
                st = path.lstat()  # lstat to handle symlinks correctly
                file_data.append((rel_path, st.st_mtime_ns, st.st_size))
            except OSError:
                continue

    # Sort for deterministic hash
    file_data.sort(key=lambda x: x[0])

    for rel_path, mtime_ns, size in file_data:
        hasher.update(f"{rel_path}:{mtime_ns}:{size}\n".encode("utf-8"))

    result = hasher.hexdigest()
    _MTIME_HASH_CACHE[root] = (result, 0.0)
    return result


def _clear_mtime_cache(root: Path) -> None:
    """Clear the mtime hash cache for a project root."""
    _MTIME_HASH_CACHE.pop(root, None)


def _resolve_latest_tree(root: Path, env: Dict[str, str]) -> Optional[str]:
    proc = _run_git(
        ["rev-parse", "--verify", "refs/snapshots/latest"],
        cwd=root,
        env=env,
        check=False,
    )
    if proc.returncode != 0:
        return None

    object_id = proc.stdout.decode("utf-8").strip()
    type_proc = _run_git(["cat-file", "-t", object_id], cwd=root, env=env, check=True)
    object_type = type_proc.stdout.decode("utf-8").strip()
    if object_type == "commit":
        tree_proc = _run_git(
            ["rev-parse", f"{object_id}^{{tree}}"], cwd=root, env=env, check=True
        )
        return tree_proc.stdout.decode("utf-8").strip()
    return object_id


def _build_tree(
    root: Path,
    env: Dict[str, str],
    tracked_files: Optional[List[str]],
) -> str:
    """Build a git tree from the current working directory.

    This optimized version batches git operations to minimize subprocess calls:
    - Uses `git hash-object -w --stdin-paths` to hash all regular files in one call
    - Uses `git update-index --index-info` to update the index in one call
    """
    # Collect all files first
    file_entries = _collect_files(root, env)

    if not file_entries:
        # Empty tree
        proc = _run_git(["write-tree"], cwd=root, env=env, check=True)
        return proc.stdout.decode("utf-8").strip()

    # Separate regular files from symlinks (symlinks need special handling)
    regular_files: List[Tuple[str, str]] = []  # (rel_path, mode)
    symlink_entries: List[Tuple[str, str, str]] = []  # (rel_path, mode, sha)

    for rel_path, mode, is_symlink, symlink_target in file_entries:
        if tracked_files is not None:
            tracked_files.append(rel_path)
        if is_symlink:
            # Hash symlink targets individually (they're usually few)
            proc = _run_git(
                ["hash-object", "-w", "--stdin"],
                cwd=root,
                env=env,
                input=symlink_target.encode("utf-8") if symlink_target else b"",
                check=True,
            )
            sha = proc.stdout.decode("utf-8").strip()
            symlink_entries.append((rel_path, mode, sha))
        else:
            regular_files.append((rel_path, mode))

    # Batch hash all regular files in one call
    file_shas: Dict[str, str] = {}
    if regular_files:
        file_paths = [rel_path for rel_path, _ in regular_files]
        file_shas = _batch_hash_files(root, env, file_paths)

    # Build index entries for batch update
    index_entries: List[str] = []

    # Add regular files
    for rel_path, mode in regular_files:
        sha = file_shas.get(rel_path)
        if sha:
            index_entries.append(f"{mode} {sha}\t{rel_path}")

    # Add symlinks
    for rel_path, mode, sha in symlink_entries:
        index_entries.append(f"{mode} {sha}\t{rel_path}")

    # Batch update index in one call
    if index_entries:
        _batch_update_index(root, env, index_entries)

    proc = _run_git(["write-tree"], cwd=root, env=env, check=True)
    return proc.stdout.decode("utf-8").strip()


def _collect_files(
    root: Path, env: Dict[str, str]
) -> List[Tuple[str, str, bool, Optional[str]]]:
    """Collect all files to be tracked.

    Returns list of (rel_path, mode, is_symlink, symlink_target).
    """
    gitignore = _GitIgnoreHelper(root)
    results: List[Tuple[str, str, bool, Optional[str]]] = []

    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        current = Path(dirpath)
        filtered_dirnames: List[str] = []
        for d in dirnames:
            dir_path = current / d
            if _should_ignore_path(dir_path):
                continue
            if gitignore.is_ignored(dir_path, is_dir=True):
                continue
            filtered_dirnames.append(d)
        dirnames[:] = filtered_dirnames

        for name in filenames:
            path = current / name
            if _should_ignore_path(path):
                continue
            if gitignore.is_ignored(path, is_dir=False):
                continue

            rel_path = path.relative_to(root).as_posix()
            if path.is_symlink():
                target = os.readlink(path)
                results.append((rel_path, "120000", True, target))
            else:
                try:
                    st_mode = path.stat().st_mode
                except OSError:
                    continue
                if st_mode & stat.S_IXUSR:
                    mode = "100755"
                else:
                    mode = "100644"
                results.append((rel_path, mode, False, None))

    return results


def _batch_hash_files(
    root: Path, env: Dict[str, str], file_paths: List[str]
) -> Dict[str, str]:
    """Hash multiple files in a single git hash-object call.

    Uses `git hash-object -w --stdin-paths` which reads file paths from stdin
    and outputs one SHA per line in the same order.
    """
    if not file_paths:
        return {}

    # git hash-object --stdin-paths expects one path per line
    input_data = "\n".join(file_paths).encode("utf-8")

    proc = _run_git(
        ["hash-object", "-w", "--stdin-paths"],
        cwd=root,
        env=env,
        input=input_data,
        check=True,
    )

    shas = proc.stdout.decode("utf-8").strip().split("\n")

    # Map paths to their SHAs
    result: Dict[str, str] = {}
    for path, sha in zip(file_paths, shas):
        if sha:
            result[path] = sha

    return result


def _batch_update_index(root: Path, env: Dict[str, str], entries: List[str]) -> None:
    """Update git index with multiple entries in a single call.

    Uses `git update-index --index-info` which reads entries from stdin
    in the format: mode SP sha TAB path LF
    """
    if not entries:
        return

    input_data = "\n".join(entries).encode("utf-8")

    _run_git(
        ["update-index", "--index-info"],
        cwd=root,
        env=env,
        input=input_data,
        check=True,
    )


def _should_ignore_path(path: Path) -> bool:
    if path.name in _DEFAULT_IGNORE_FILE_NAMES:
        return True
    if "snapshot.git" in path.parts:
        return True
    if any(part in _DEFAULT_IGNORE_DIRS for part in path.parts):
        return True
    if path.is_file() and path.suffix in _DEFAULT_IGNORE_FILE_SUFFIXES:
        return True
    return False


def _parse_diff_output(
    diff_output: str,
) -> Tuple[Set[str], Set[str], Set[str]]:
    added: Set[str] = set()
    modified: Set[str] = set()
    deleted: Set[str] = set()

    for line in diff_output.splitlines():
        if not line:
            continue
        status, _, path = line.partition("\t")
        if not path:
            continue
        if status == "A":
            added.add(path)
        elif status == "M":
            modified.add(path)
        elif status == "D":
            deleted.add(path)
        else:
            modified.add(path)
    return added, modified, deleted


def _update_snapshot_ref(root: Path, env: Dict[str, str], tree: str) -> None:
    _run_git(
        ["update-ref", "refs/snapshots/latest", tree],
        cwd=root,
        env=env,
        check=True,
    )


def _write_manifest(root: Path, tree: str, ignore_hash: str, mtime_hash: str) -> None:
    manifest_dir = root / ".fast" / "manifests"
    manifest_data = {
        "created_at": _dt.datetime.now(tz=_dt.timezone.utc).isoformat(),
        "tree": tree,
        "kind": "tree",
        "tool_version": _FAST_TCP_VERSION,
        "python": platform.python_version(),
        "ignore_hash": ignore_hash,
        "mtime_hash": mtime_hash,
    }
    if _manifest_matches(root, manifest_data):
        return
    manifest_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = manifest_dir / "latest.json"
    manifest_path.write_text(
        json.dumps(manifest_data, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _manifest_matches(root: Path, manifest_data: dict[str, Any]) -> bool:
    current = _read_latest_manifest(root)
    if not current:
        return False
    comparable_fields = (
        "tree",
        "kind",
        "tool_version",
        "python",
        "ignore_hash",
        "mtime_hash",
    )
    return all(
        current.get(field) == manifest_data.get(field) for field in comparable_fields
    )


def _read_latest_manifest(root: Path) -> Optional[Dict[str, object]]:
    manifest_path = root / ".fast" / "manifests" / "latest.json"
    if not manifest_path.exists():
        return None
    try:
        return json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def _prune_snapshot_objects(root: Path, env: Dict[str, str]) -> None:
    _run_git(
        ["gc", "--prune=now"],
        cwd=root,
        env=env,
        check=True,
    )


def _run_git(
    args: Sequence[str],
    *,
    cwd: Path,
    env: Dict[str, str],
    check: bool,
    input: Optional[bytes] = None,
) -> subprocess.CompletedProcess[bytes]:
    try:
        proc = subprocess.run(
            ["git", *args],
            cwd=str(cwd),
            env=env,
            input=input,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    except FileNotFoundError as exc:  # pragma: no cover - depends on environment
        raise RuntimeError(
            "Git executable not found; install Git to use snapshot caching."
        ) from exc

    if check and proc.returncode != 0:
        raise subprocess.CalledProcessError(
            returncode=proc.returncode,
            cmd=["git", *args],
            output=proc.stdout,
            stderr=proc.stderr,
        )
    return proc
