import tempfile
from pathlib import Path
from typing import List

import pytest

from fast_tcp import fast


def _write_blackbox_input(nodeids: List[str], tmp_dir: Path, k: int) -> Path:
    # Each line represents a test case. We'll use the nodeid string as the sequence to shingle.
    input_file = tmp_dir / "suite-bbox.txt"
    with open(input_file, "w") as f:
        for nodeid in nodeids:
            f.write(nodeid + "\n")
    return input_file


def _run_fast_blackbox(
    input_file: Path, algo: str, r: int, b: int, k: int, budget: int
) -> List[int]:
    # Map algorithm variant to callable
    if algo == "FAST-pw":
        _, _, order = fast.fast_pw(
            str(input_file), r=r, b=b, bbox=True, k=k, memory=True, B=budget
        )
        return order
    elif algo in {"FAST-one", "FAST-log", "FAST-sqrt", "FAST-all"}:

        def sel_all(x: int) -> int:
            return x

        def sel_sqrt(x: int) -> int:
            import math

            return int(math.sqrt(x)) + 1

        def sel_log(x: int) -> int:
            import math

            return int(math.log(x, 2)) + 1

        def sel_one(x: int) -> int:
            return 1

        sels = {
            "FAST-one": sel_one,
            "FAST-log": sel_log,
            "FAST-sqrt": sel_sqrt,
            "FAST-all": sel_all,
        }
        _, _, order = fast.fast_(
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
    else:
        raise ValueError(f"Unsupported FAST algo: {algo}")


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
        help="k-shingle size for BB",
    )
    group.addoption(
        "--fast-tcp-budget",
        action="store",
        type=int,
        default=0,
        help="Budget B (0=all)",
    )


def pytest_collection_modifyitems(
    session: pytest.Session, config: pytest.Config, items: List[pytest.Item]
) -> None:
    if not config.getoption("--fast-tcp") or not items:
        return

    algo = config.getoption("--fast-tcp-algo")
    r = int(config.getoption("--fast-tcp-r"))
    b = int(config.getoption("--fast-tcp-b"))
    k = int(config.getoption("--fast-tcp-k"))
    budget = int(config.getoption("--fast-tcp-budget"))

    # Build the input file with nodeids for black-box prioritization
    nodeids = [it.nodeid for it in items]
    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        input_file = _write_blackbox_input(nodeids, tmp_dir, k)

        # Compute prioritized order (list of 1-based indices)
        order = _run_fast_blackbox(input_file, algo, r, b, k, budget)

    # Map 1-based indices back to the pytest item list order
    # FAST returns positions by input order, so we map them: 1 -> nodeids[0], etc.
    index_to_item = {i + 1: it for i, it in enumerate(items)}

    # Keep only those returned by FAST; if budget < len(items), some will be omitted
    new_items = [index_to_item[i] for i in order if i in index_to_item]

    # If budget omitted some, append the rest in original order to still run all tests
    if len(new_items) < len(items):
        seen = set(new_items)
        for it in items:
            if it not in seen:
                new_items.append(it)

    items[:] = new_items
