#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def strip_comments(go_src: str) -> str:
    no_block = re.sub(r"/\*.*?\*/", " ", go_src, flags=re.S)
    no_line = re.sub(r"//.*", " ", no_block)
    return no_line


def extract_test_names(go_src: str) -> list[str]:
    text = strip_comments(go_src)
    names: list[str] = []
    for m in re.finditer(
        r"\bfunc\s+(Test[A-Za-z0-9_]+)\s*\(\s*t\s*\*testing\.T\s*\)", text
    ):
        names.append(m.group(1))
    return names


def build_test_list(tests_dir: Path) -> list[str]:
    test_names: list[str] = []
    for go_file in sorted(tests_dir.rglob("*_test.go")):
        src = read(go_file)
        for name in extract_test_names(src):
            test_names.append(name)
    return test_names


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Emit Go test function names in discovery order"
    )
    ap.add_argument(
        "--tests-dir", required=True, help="Path to module dir to scan for *_test.go"
    )
    ap.add_argument("--out", required=True, help="Output path for test-names.txt")
    args = ap.parse_args()

    tests_dir = Path(args.tests_dir)
    names = build_test_list(tests_dir)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(names) + "\n", encoding="utf-8")
    print(f"Wrote {len(names)} test names to {out_path}")


if __name__ == "__main__":
    main()

