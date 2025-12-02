#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


def read_text(path: Path) -> str:
    with open(path, "r", encoding="utf-8") as f:
        return f.read()


def strip_comments(go_src: str) -> str:
    no_block = re.sub(r"/\*.*?\*/", " ", go_src, flags=re.S)
    no_line = re.sub(r"//.*", " ", no_block)
    return no_line


def extract_test_functions(go_src: str) -> list[str]:
    text = strip_comments(go_src)
    tests: list[str] = []
    for m in re.finditer(r"\bfunc\s+(Test[A-Za-z0-9_]+)\s*\(.*?\)\s*\{", text):
        start = m.end() - 1
        depth = 0
        body_chars: list[str] = []
        for ch in text[start:]:
            body_chars.append(ch)
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth <= 0:
                    break
        body = "".join(body_chars)
        if body:
            tests.append(body)
    return tests


def tokenize(text: str) -> str:
    tokens = re.findall(r"[A-Za-z0-9_]+", text)
    return " ".join(tokens)


def collect_tests(tests_dir: Path) -> list[str]:
    cases: list[str] = []
    for path in tests_dir.rglob("*_test.go"):
        src = read_text(path)
        bodies = extract_test_functions(src)
        if bodies:
            for body in bodies:
                cases.append(tokenize(body))
        else:
            cases.append(tokenize(src))
    return cases


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate bbox input from Go tests")
    ap.add_argument("--tests-dir", required=True, help="Path with *_test.go files")
    ap.add_argument(
        "--out", required=True, help="Output .txt path (e.g., .fast/in/gotest-bbox.txt)"
    )
    args = ap.parse_args()

    tests_dir = Path(args.tests_dir).resolve()
    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cases = collect_tests(tests_dir)
    if not cases:
        raise SystemExit(f"No Go test cases found in {tests_dir}")

    with open(out_path, "w", encoding="utf-8") as f:
        for line in cases:
            f.write(line + "\n")

    print(f"Wrote {len(cases)} test cases to {out_path}")


if __name__ == "__main__":
    main()



