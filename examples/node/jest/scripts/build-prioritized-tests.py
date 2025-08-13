#!/usr/bin/env python3
import argparse
import glob
import pickle
from pathlib import Path


def find_latest_pickle(out_dir: Path) -> Path:
    # Match both with or without nested "prioritized" depending on producer
    patterns = [
        str(out_dir / "*" / "FAST-*-bbox-*.pickle"),
        str(out_dir / "*" / "prioritized" / "FAST-*-bbox-*.pickle"),
        str(out_dir / "FAST-*-bbox-*.pickle"),
    ]
    candidates: list[str] = []
    for pat in patterns:
        candidates.extend(glob.glob(pat))
    if not candidates:
        raise SystemExit(f"No prioritized pickle found under {out_dir}")
    candidates.sort(key=lambda p: Path(p).stat().st_mtime, reverse=True)
    return Path(candidates[0])


def main() -> None:
    ap = argparse.ArgumentParser(description="Map prioritized IDs to Jest (file\tname)")
    ap.add_argument(
        "--out-dir", required=True, help="Directory where CLI wrote results (.fast/out)"
    )
    ap.add_argument(
        "--tests", required=True, help="Path to test-names.tsv produced by tests-map.py"
    )
    ap.add_argument(
        "--out", required=True, help="Output path for prioritized-tests.tsv"
    )
    args = ap.parse_args()

    out_dir = Path(args.out_dir).resolve()
    names_path = Path(args.tests).resolve()
    out_path = Path(args.out).resolve()

    pickle_path = find_latest_pickle(out_dir)
    order = pickle.load(open(pickle_path, "rb"))
    pairs = [
        tuple(line.split("\t", 1))
        for line in names_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]

    # IDs are 1-based
    prioritized_pairs = [pairs[i - 1] for i in order if 1 <= i <= len(pairs)]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        "\n".join(f"{f}\t{n}" for f, n in prioritized_pairs) + "\n", encoding="utf-8"
    )
    print(f"Mapped {len(prioritized_pairs)} tests from {pickle_path} -> {out_path}")


if __name__ == "__main__":
    main()
