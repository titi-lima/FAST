#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


TEST_NAME_RE = re.compile(
    r"(?:it|test)\s*(?:\.(?:only|skip))?\s*\(\s*([\'\"`])(?P<name>.*?)\1\s*,",
    re.S,
)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def discover_test_files(root: Path) -> list[Path]:
    files: list[Path] = []
    for pattern in (
        "**/*.test.js",
        "**/*.test.ts",
        "**/*.spec.js",
        "**/*.spec.ts",
        "**/*.test.mjs",
        "**/*.spec.mjs",
        "**/*.test.cjs",
        "**/*.spec.cjs",
        "**/*.test.jsx",
        "**/*.spec.jsx",
        "**/*.test.tsx",
        "**/*.spec.tsx",
    ):
        files.extend(root.glob(pattern))
    filtered = []
    for p in files:
        rp = p.resolve()
        if any(part == "node_modules" for part in rp.parts):
            continue
        if any(part in {"dist", "build", ".fast"} for part in rp.parts):
            continue
        filtered.append(rp)
    return sorted(set(filtered))


def extract_test_names(file_path: Path) -> list[str]:
    src = read_text(file_path)
    return [m.group("name") for m in TEST_NAME_RE.finditer(src)]


def build_test_map(tests_dir: Path) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for test_file in discover_test_files(tests_dir):
        names = extract_test_names(test_file)
        for name in names:
            pairs.append((str(test_file.resolve()), name))
    return pairs


def main() -> None:
    ap = argparse.ArgumentParser(description="Emit Vitest test map (file\tname)")
    ap.add_argument(
        "--tests-dir", required=True, help="Directory to scan (example root)"
    )
    ap.add_argument(
        "--out", required=True, help="Output TSV with file and full test name"
    )
    args = ap.parse_args()

    tests_dir = Path(args.tests_dir).resolve()
    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    pairs = build_test_map(tests_dir)
    lines = [f"{file}\t{name}" for file, name in pairs]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {len(pairs)} test entries to {out_path}")


if __name__ == "__main__":
    main()

