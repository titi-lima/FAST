#!/usr/bin/env python3
import argparse
import glob
import os
import pickle
from pathlib import Path


def find_latest_pickle(out_dir: Path) -> Path:
    # Match both with or without nested "prioritized" depending on producer
    patterns = [
        str(out_dir / "*" / "FAST-*-bbox-*.pickle"),
        str(out_dir / "*" / "prioritized" / "FAST-*-bbox-*.pickle"),
    ]
    candidates = []
    for pat in patterns:
        candidates.extend(glob.glob(pat))
    if not candidates:
        raise FileNotFoundError(f"No prioritized pickle found under {out_dir}")
    candidates.sort(key=os.path.getmtime, reverse=True)
    return Path(candidates[0])


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Map prioritized IDs to JUnit method selectors"
    )
    ap.add_argument(
        "--out-dir",
        required=True,
        help="Directory where CLI wrote results (e.g., <example>/out)",
    )
    ap.add_argument(
        "--selectors",
        required=True,
        help="Path to selectors.txt produced by tests-map.py",
    )
    ap.add_argument(
        "--out", required=True, help="Output path for prioritized-selectors.txt"
    )
    args = ap.parse_args()

    out_dir = Path(args.out_dir).resolve()
    selectors_path = Path(args.selectors).resolve()
    out_path = Path(args.out).resolve()

    pickle_path = find_latest_pickle(out_dir)
    order = pickle.load(open(pickle_path, "rb"))
    selectors = selectors_path.read_text(encoding="utf-8").splitlines()

    # IDs are 1-based
    prioritized = [selectors[i - 1] for i in order if 1 <= i <= len(selectors)]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(prioritized) + "\n", encoding="utf-8")
    print(f"Mapped {len(prioritized)} selectors from {pickle_path} -> {out_path}")


if __name__ == "__main__":
    main()
