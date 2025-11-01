#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def _run(cmd: list[str]) -> int:
    return subprocess.call(cmd)


def ensure_junit_console(example_dir: Path) -> None:
    lib_dir = example_dir / "lib"
    lib_dir.mkdir(parents=True, exist_ok=True)
    jar = lib_dir / "junit-platform-console-standalone.jar"
    if jar.exists():
        return
    script = example_dir / "scripts" / "get-junit.sh"
    if script.exists():
        subprocess.check_call(["bash", str(script)])


def ant_build(example_dir: Path) -> None:
    build_xml = example_dir / "build.xml"
    subprocess.check_call(["ant", "-f", str(build_xml), "clean", "compile"])


def generate_bbox(tests_dir: Path, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    _gen_bbox_impl.generate_bbox_from_junit_tests(
        tests_dir=tests_dir, out_path=out_path
    )


def prioritize_bbox(
    data_dir: Path, out_dir: Path, algo: str, repetitions: int, debug: bool = False
) -> None:
    cmd = [
        sys.executable,
        "-m",
        "fast_tcp.cli",
        "--test-dir",
        str(data_dir),
        "--algo",
        algo,
        "--entity",
        "bbox",
        "--repetitions",
        str(repetitions),
        "--file-naming",
        "entity-suffix",
        "--output-dir",
        str(out_dir),
    ]
    if debug:
        cmd.append("--debug")
    subprocess.check_call(cmd)


def build_selectors(tests_dir: Path, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    selectors = _tests_map_impl.build_selectors(tests_dir)
    out_path.write_text("\n".join(selectors) + "\n", encoding="utf-8")


def map_prioritized_to_selectors(
    out_dir: Path, selectors_path: Path, out_path: Path
) -> None:
    import glob
    import pickle

    patterns = [
        str(out_dir / "*" / "FAST-*-bbox-*.pickle"),
        str(out_dir / "*" / "prioritized" / "FAST-*-bbox-*.pickle"),
    ]
    candidates: list[str] = []
    for pat in patterns:
        candidates.extend(glob.glob(pat))
    if not candidates:
        raise FileNotFoundError(f"No prioritized pickle found under {out_dir}")
    candidates.sort(key=os.path.getmtime, reverse=True)
    pickle_path = Path(candidates[0])

    order = pickle.load(open(pickle_path, "rb"))
    selectors = selectors_path.read_text(encoding="utf-8").splitlines()
    prioritized = [selectors[i - 1] for i in order if 1 <= i <= len(selectors)]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(prioritized) + "\n", encoding="utf-8")


def run_junit_in_order(example_dir: Path, prioritized_selectors: Path) -> int:
    lib_jar = example_dir / "lib" / "junit-platform-console-standalone.jar"
    main_cp = example_dir / "build" / "classes" / "main"
    test_cp = example_dir / "build" / "classes" / "test"
    classpath = f"{main_cp}:{test_cp}:{lib_jar}"
    args: list[str] = []
    for line in prioritized_selectors.read_text(encoding="utf-8").splitlines():
        if line.strip():
            args.extend(["--select-method", line.strip()])
    cmd = [
        "java",
        "-cp",
        classpath,
        "org.junit.platform.console.ConsoleLauncher",
        "--fail-if-no-tests",
        "--class-path",
        f"{main_cp}:{test_cp}",
        *args,
    ]
    return _run(cmd)


def main() -> int:
    parser = argparse.ArgumentParser(description="FAST TCP Ant integration helper")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_macro = sub.add_parser("macro", help="Run full pipeline for an Ant project")
    p_macro.add_argument("--project-dir", required=True)
    p_macro.add_argument("--algo", default="FAST-pw")
    p_macro.add_argument("--repetitions", type=int, default=3)
    p_macro.add_argument(
        "--debug",
        action="store_true",
        help="Print debug information including timing for preparation and prioritization",
    )

    args = parser.parse_args()

    if args.cmd == "macro":
        project_dir = Path(args.project_dir).resolve()
        ensure_junit_console(project_dir)
        ant_build(project_dir)

        tests_dir = project_dir / "src" / "test" / "java"
        data_dir = project_dir / ".fast" / "in"
        out_dir = project_dir / ".fast" / "out"
        data_dir.mkdir(parents=True, exist_ok=True)
        out_dir.mkdir(parents=True, exist_ok=True)

        bbox_path = data_dir / f"{project_dir.name}-bbox.txt"
        selectors_path = data_dir / "selectors.txt"
        prioritized_selectors_path = data_dir / "prioritized-selectors.txt"

        build_selectors(tests_dir, selectors_path)
        generate_bbox(tests_dir, bbox_path)
        prioritize_bbox(
            data_dir, out_dir, args.algo, args.repetitions, debug=args.debug
        )
        map_prioritized_to_selectors(
            out_dir, selectors_path, prioritized_selectors_path
        )
        return run_junit_in_order(project_dir, prioritized_selectors_path)

    return 0


# Lightweight reuse of example scripts without importing their CLIs directly
class _gen_bbox_impl:
    @staticmethod
    def generate_bbox_from_junit_tests(*, tests_dir: Path, out_path: Path) -> None:
        import re

        def read_text(path: Path) -> str:
            return path.read_text(encoding="utf-8")

        def strip_comments(java_src: str) -> str:
            no_block = re.sub(r"/\*.*?\*/", " ", java_src, flags=re.S)
            no_line = re.sub(r"//.*", " ", no_block)
            return no_line

        def extract_test_methods(java_src: str) -> list[str]:
            text = strip_comments(java_src)
            tests: list[str] = []
            parts = re.split(r"@Test\b", text)
            if len(parts) <= 1:
                return tests
            for part in parts[1:]:
                import re as _re

                m = _re.search(r"\)\s*\{", part)
                if not m:
                    brace_idx = part.find("{")
                    if brace_idx == -1:
                        continue
                    start = brace_idx
                else:
                    start = m.end() - 1
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
            import re as _re

            tokens = _re.findall(r"[A-Za-z0-9_]+", text)
            return " ".join(tokens)

        cases: list[str] = []
        for path in tests_dir.rglob("*.java"):
            src = read_text(path)
            methods = extract_test_methods(src)
            if methods:
                for method_src in methods:
                    cases.append(tokenize(method_src))
            else:
                cases.append(tokenize(src))

        out_path.write_text("\n".join(cases) + "\n", encoding="utf-8")


class _tests_map_impl:
    @staticmethod
    def build_selectors(tests_dir: Path) -> list[str]:
        import re

        def read(path: Path) -> str:
            return path.read_text(encoding="utf-8")

        def extract_package(java_src: str) -> str:
            m = re.search(r"\bpackage\s+([\w\.]+)\s*;", java_src)
            return m.group(1) if m else ""

        def extract_class(java_src: str) -> str:
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
                m = re.search(r"\b([A-Za-z0-9_]+)\s*\(.*?\)\s*\{", chunk)
                if m:
                    methods.append(m.group(1))
            return methods

        selectors: list[str] = []
        for java_file in sorted(tests_dir.rglob("*.java")):
            src = read(java_file)
            pkg = extract_package(src)
            cls = extract_class(src)
            fqcn = f"{pkg}.{cls}" if pkg else cls
            for method in extract_test_methods(src):
                selectors.append(f"{fqcn}#{method}")
        return selectors


if __name__ == "__main__":
    raise SystemExit(main())
