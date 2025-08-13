#!/usr/bin/env python3
import argparse
import os
import re
from pathlib import Path


def read_text(path: Path) -> str:
    with open(path, "r", encoding="utf-8") as f:
        return f.read()


def strip_comments(java_src: str) -> str:
    # Remove block comments and line comments
    no_block = re.sub(r"/\*.*?\*/", " ", java_src, flags=re.S)
    no_line = re.sub(r"//.*", " ", no_block)
    return no_line


def extract_test_methods(java_src: str) -> list[str]:
    text = strip_comments(java_src)
    tests: list[str] = []

    # Split by occurrences of @Test (JUnit 5). Keep text after each annotation.
    parts = re.split(r"@Test\b", text)
    if len(parts) <= 1:
        return tests

    # For every part after the first, try to capture the first method body
    for part in parts[1:]:
        # Find the method signature and body starting brace
        m = re.search(r"\)\s*\{", part)
        if not m:
            # Fallback: attempt to find first opening brace
            brace_idx = part.find("{")
            if brace_idx == -1:
                continue
            start = brace_idx
        else:
            start = m.end() - 1  # position at '{'

        # Extract balanced braces to get the method body
        depth = 0
        body_chars: list[str] = []
        for ch in part[start:]:
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
    # Basic tokenization: words and numbers
    tokens = re.findall(r"[A-Za-z0-9_]+", text)
    return " ".join(tokens)


def collect_tests(tests_dir: Path) -> list[str]:
    cases: list[str] = []
    for path in tests_dir.rglob("*.java"):
        src = read_text(path)
        methods = extract_test_methods(src)
        if methods:
            for method_src in methods:
                cases.append(tokenize(method_src))
        else:
            # Fallback: tokenize entire file as a single test case representation
            cases.append(tokenize(src))
    return cases


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate black-box input for FAST TCP from JUnit tests"
    )
    parser.add_argument(
        "--tests-dir", required=True, help="Path to src/test/java directory"
    )
    parser.add_argument(
        "--out", required=True, help="Output file path, e.g., data/ant-sample-bbox.txt"
    )
    args = parser.parse_args()

    tests_dir = Path(args.tests_dir).resolve()
    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cases = collect_tests(tests_dir)
    if not cases:
        raise SystemExit(f"No test cases found in {tests_dir}")

    with open(out_path, "w", encoding="utf-8") as f:
        for line in cases:
            f.write(line + "\n")

    print(f"Wrote {len(cases)} test cases to {out_path}")


if __name__ == "__main__":
    main()
