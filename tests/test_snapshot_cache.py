"""Tests for the snapshot cache incremental Git repository."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

from fast_tcp import snapshot_cache


GIT_AVAILABLE = shutil.which("git") is not None

pytestmark = pytest.mark.skipif(
    not GIT_AVAILABLE, reason="Git is required for snapshot cache tests."
)


def _write_file(root: Path, rel_path: str, content: str) -> Path:
    path = root / rel_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def test_snapshot_cache_classifies_added_modified_deleted(tmp_path: Path) -> None:
    project_root = tmp_path / "proj"
    project_root.mkdir()

    # Initial set of files.
    _write_file(project_root, "tests/test_alpha.py", "print('alpha v1')\n")
    _write_file(project_root, "tests/test_bravo.py", "print('bravo v1')\n")
    _write_file(project_root, "tests/subpkg/test_charlie.py", "print('charlie v1')\n")

    first_diff = snapshot_cache.detect_changes_preparation(project_root)
    assert first_diff.added == {
        "tests/test_alpha.py",
        "tests/test_bravo.py",
        "tests/subpkg/test_charlie.py",
    }
    assert not first_diff.modified
    assert not first_diff.deleted
    snapshot_cache.snapshot_prioritization(project_root)

    # Mutate the suite: modify one file, delete one, create another.
    _write_file(project_root, "tests/test_alpha.py", "print('alpha v2')\n")
    (project_root / "tests/test_bravo.py").unlink()
    _write_file(project_root, "tests/subpkg/test_new.py", "print('new test')\n")

    second_diff = snapshot_cache.detect_changes_preparation(project_root)
    assert second_diff.modified == {"tests/test_alpha.py"}
    assert second_diff.added == {"tests/subpkg/test_new.py"}
    assert second_diff.deleted == {"tests/test_bravo.py"}
    # Ensure untouched file is not flagged as changed.
    assert "tests/subpkg/test_charlie.py" not in (
        second_diff.added | second_diff.modified | second_diff.deleted
    )
    snapshot_cache.snapshot_prioritization(project_root)
