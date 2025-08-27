"""Pytest plugin for FAST TCP test prioritization.

This plugin reorders collected tests using FAST in black-box mode over pytest nodeids.
It registers a `--fast-tcp` flag and related configuration options. When the flag is
enabled, the `pytest_collection_modifyitems` hook computes a prioritized execution
order and rewrites the collected items accordingly.

The plugin is auto-discovered via the `pytest11` entry point when the package is
installed. Users can enable it via the CLI flag or by setting `addopts` in pytest.ini.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import List

import pytest

from ... import fast as _fast


def _write_blackbox_input(nodeids: List[str], tmp_dir: Path) -> Path:
    input_file = tmp_dir / "suite-bbox.txt"
    with open(input_file, "w", encoding="utf-8") as f:
        for nodeid in nodeids:
            f.write(nodeid + "\n")
    return input_file


def _run_fast_blackbox(
    *, input_file: Path, algo: str, r: int, b: int, k: int, budget: int
) -> List[int]:
    if algo == "FAST-pw":
        _, _, order = _fast.fast_pw(
            str(input_file), r=r, b=b, bbox=True, k=k, memory=True, B=budget
        )
        return order

    def _sel_all(x: int) -> int:
        return x

    def _sel_sqrt(x: int) -> int:
        import math

        return int(math.sqrt(x)) + 1

    def _sel_log(x: int) -> int:
        import math

        return int(math.log(x, 2)) + 1

    def _sel_one(x: int) -> int:
        return 1

    sels = {
        "FAST-one": _sel_one,
        "FAST-log": _sel_log,
        "FAST-sqrt": _sel_sqrt,
        "FAST-all": _sel_all,
    }
    if algo not in sels:
        raise ValueError(f"Unsupported FAST algo: {algo}")

    _, _, order = _fast.fast_(
        str(input_file),
        selsize=sels[algo],
        r=r,
        b=b,
        bbox=True,
        k=k,
        memory=True,
        B=budget,
    )
    return order


def pytest_addoption(parser: pytest.Parser) -> None:
    group = parser.getgroup("fast-tcp")
    group.addoption(
        "--fast-tcp", action="store_true", help="Enable FAST TCP ordering of tests"
    )
    group.addoption(
        "--fast-tcp-algo",
        action="store",
        default="FAST-pw",
        help="FAST variant: FAST-pw | FAST-one | FAST-log | FAST-sqrt | FAST-all",
    )
    group.addoption(
        "--fast-tcp-r", action="store", type=int, default=1, help="FAST rows (r)"
    )
    group.addoption(
        "--fast-tcp-b", action="store", type=int, default=10, help="FAST bands (b)"
    )
    group.addoption(
        "--fast-tcp-k",
        action="store",
        type=int,
        default=5,
        help="k-shingle size for black-box mode",
    )
    group.addoption(
        "--fast-tcp-budget",
        action="store",
        type=int,
        default=0,
        help="Budget B (0=all tests)",
    )


def pytest_collection_modifyitems(
    session: pytest.Session, config: pytest.Config, items: List[pytest.Item]
) -> None:
    if not items or not config.getoption("--fast-tcp"):
        return

    algo = config.getoption("--fast-tcp-algo")
    r = int(config.getoption("--fast-tcp-r"))
    b = int(config.getoption("--fast-tcp-b"))
    k = int(config.getoption("--fast-tcp-k"))
    budget = int(config.getoption("--fast-tcp-budget"))

    nodeids = [it.nodeid for it in items]
    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        input_file = _write_blackbox_input(nodeids, tmp_dir)
        order = _run_fast_blackbox(
            input_file=input_file, algo=algo, r=r, b=b, k=k, budget=budget
        )

    # Map back to pytest items (FAST uses 1-based indices in input order)
    index_to_item = {i + 1: it for i, it in enumerate(items)}
    new_items = [index_to_item[i] for i in order if i in index_to_item]

    # If a budget excluded some, append the rest in original order
    if len(new_items) < len(items):
        seen = set(new_items)
        for it in items:
            if it not in seen:
                new_items.append(it)

    items[:] = new_items
