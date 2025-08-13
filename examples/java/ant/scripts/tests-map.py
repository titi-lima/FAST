#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def extract_package(java_src: str) -> str:
    m = re.search(r"\bpackage\s+([\w\.]+)\s*;", java_src)
    return m.group(1) if m else ""


def extract_class(java_src: str) -> str:
    # take the first public class name as the test container
    m = re.search(r"\bclass\s+([A-Za-z0-9_]+)\b", java_src)
    return m.group(1) if m else "UnknownTestClass"


def strip_comments(java_src: str) -> str:
    no_block = re.sub(r"/\*.*?\*/", " ", java_src, flags=re.S)
    no_line = re.sub(r"//.*", " ", no_block)
    return no_line


def extract_test_methods(java_src: str) -> list[str]:
    text = strip_comments(java_src)
    methods: list[str] = []
    parts = re.split(r"@Test\b", text)
    for chunk in parts[1:]:
        # find the next method name after @Test
        m = re.search(r"\b([A-Za-z0-9_]+)\s*\(.*?\)\s*\{", chunk)
        if m:
            methods.append(m.group(1))
    return methods


def build_selectors(tests_dir: Path) -> list[str]:
    selectors: list[str] = []
    for java_file in sorted(tests_dir.rglob("*.java")):
        src = read(java_file)
        pkg = extract_package(src)
        cls = extract_class(src)
        fqcn = f"{pkg}.{cls}" if pkg else cls
        for method in extract_test_methods(src):
            selectors.append(f"{fqcn}#{method}")
    return selectors


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Emit JUnit5 method selectors in discovery order"
    )
    ap.add_argument("--tests-dir", required=True, help="Path to src/test/java")
    ap.add_argument(
        "--out", required=True, help="Output file to write selectors (one per line)"
    )
    args = ap.parse_args()

    tests_dir = Path(args.tests_dir)
    selectors = build_selectors(tests_dir)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(selectors) + "\n", encoding="utf-8")
    print(f"Wrote {len(selectors)} selectors to {out_path}")


if __name__ == "__main__":
    main()
