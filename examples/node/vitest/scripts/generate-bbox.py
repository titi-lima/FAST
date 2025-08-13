#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


TEST_DEF_RE = re.compile(
    r"(?P<kind>it|test)\s*(?:\.(?:only|skip))?\s*\(\s*(?P<q>['\"`])(?P<name>.*?)(?P=q)\s*,",
    re.S,
)


def find_test_body(file_src: str, test_name: str) -> str | None:
    for m in TEST_DEF_RE.finditer(file_src):
        if m.group("name") != test_name:
            continue
        idx = m.end()
        func_start = file_src.find("{", idx)
        if func_start == -1:
            continue
        depth = 0
        body_chars: list[str] = []
        for ch in file_src[func_start:]:
            body_chars.append(ch)
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth <= 0:
                    break
        body = "".join(body_chars)
        if body:
            return body
    return None


def tokenize(text: str) -> str:
    tokens = re.findall(r"[A-Za-z0-9_]+", text)
    return " ".join(tokens)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Generate bbox tokens per Vitest test using a tests-map TSV (file\tname)"
        )
    )
    ap.add_argument(
        "--tests-dir", required=True, help="Directory to scan (example root)"
    )
    ap.add_argument(
        "--tests-map", required=True, help="Path to test-names.tsv from tests-map.py"
    )
    ap.add_argument(
        "--out", required=True, help="Output bbox path, e.g., .fast/in/vitest-bbox.txt"
    )
    args = ap.parse_args()

    tests_dir = Path(args.tests_dir).resolve()
    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    lines = Path(args.tests_map).read_text(encoding="utf-8").splitlines()

    cases: list[str] = []
    for line in lines:
        if not line.strip():
            continue
        try:
            file_str, test_name = line.split("\t", 1)
        except ValueError:
            file_str, test_name = "", line.strip()
        file_path = (
            (tests_dir / Path(file_str)).resolve()
            if not Path(file_str).is_absolute()
            else Path(file_str)
        )
        src = read_text(file_path)
        body = find_test_body(src, test_name)
        if body is None:
            body = src
        cases.append(tokenize(body))

    if not cases:
        raise SystemExit("No tests discovered to generate bbox tokens.")

    out_path.write_text("\n".join(cases) + "\n", encoding="utf-8")
    print(f"Wrote {len(cases)} test cases to {out_path}")


if __name__ == "__main__":
    main()

