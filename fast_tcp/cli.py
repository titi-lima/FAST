"""
Command-line interface for FAST: Scalable Similarity-based Test Case Prioritization
"""

import argparse
import sys
import os
import glob
import shutil
import tempfile
import re
from pathlib import Path
from typing import List

from . import prioritize


def create_parser():
    """Create the main argument parser."""
    parser = argparse.ArgumentParser(
        description="Test Case Prioritization Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  tcp-prioritize --test-dir /path/to/tests --algo FAST-pw --entity bbox
  tcp-prioritize --test-dir ./my_tests --algo FAST-pw --entity function --repetitions 10
  tcp-prioritize --test-dir /tests --algo STR --entity bbox --output-dir ./results
        """,
    )

    # Main prioritization arguments (no subcommands needed)
    parser.add_argument(
        "--test-dir",
        required=True,
        help="Directory containing test files to prioritize",
    )
    parser.add_argument(
        "--algo",
        required=True,
        choices=[
            "FAST-pw",
            "FAST-one",
            "FAST-log",
            "FAST-sqrt",
            "FAST-all",
        ],
        help="Algorithm used for prioritization",
    )
    parser.add_argument(
        "--entity",
        required=True,
        choices=["bbox"],
        help="BB prioritization",
    )
    parser.add_argument(
        "--repetitions",
        type=int,
        default=1,
        help="Number of prioritizations to compute (default: 1)",
    )
    parser.add_argument(
        "--output-dir",
        help="Output directory for results (default: current directory + 'output')",
    )
    parser.add_argument(
        "--pattern",
        default="*.txt",
        help="File pattern to match test files (default: *.txt)",
    )
    parser.add_argument(
        "--file-naming",
        choices=["auto", "entity-suffix", "custom"],
        default="auto",
        help="How to identify entity type from filenames (default: auto)",
    )

    return parser


def create_init_parser():
    parser = argparse.ArgumentParser(
        prog="fast-tcp init",
        description="Initialize FAST TCP integration",
    )
    sub = parser.add_subparsers(dest="tool", required=True)

    p_ant = sub.add_parser("ant", help="Initialize Ant project integration")
    p_ant.add_argument("--project-dir", default=".", help="Path to Ant project root")
    p_ant.add_argument("--algo", default="FAST-pw")
    p_ant.add_argument("--repetitions", type=int, default=3)

    p_pytest = sub.add_parser("pytest", help="Initialize Pytest project integration")
    p_pytest.add_argument(
        "--project-dir", default=".", help="Path to Pytest project root"
    )
    p_pytest.add_argument(
        "--enable",
        dest="enable",
        action="store_true",
        help="Enable FAST TCP by default via addopts",
    )
    p_pytest.add_argument(
        "--no-enable", dest="enable", action="store_false", help="Do not modify addopts"
    )
    p_pytest.set_defaults(enable=True)
    p_pytest.add_argument(
        "--algo", default="FAST-pw", help="FAST variant for default config"
    )
    p_pytest.add_argument(
        "--r", type=int, default=1, help="FAST rows (r) for default config"
    )
    p_pytest.add_argument(
        "--b", type=int, default=10, help="FAST bands (b) for default config"
    )
    p_pytest.add_argument(
        "--k", type=int, default=5, help="k-shingle size for default config"
    )
    p_pytest.add_argument(
        "--budget", type=int, default=0, help="Budget B for default config (0=all)"
    )

    p_vitest = sub.add_parser("vitest", help="Initialize Vitest project integration")
    p_vitest.add_argument(
        "--project-dir", default=".", help="Path to Vitest project root"
    )
    p_vitest.add_argument("--algo", default="FAST-pw", help="FAST variant to use")
    p_vitest.add_argument("--repetitions", type=int, default=3)
    p_vitest.add_argument(
        "--no-scripts",
        action="store_true",
        help="Do not modify package.json scripts; only copy .fast helpers",
    )

    p_gotest = sub.add_parser(
        "gotest", help="Initialize Go (testing) project integration"
    )
    p_gotest.add_argument(
        "--project-dir", default=".", help="Path to Go module/project root"
    )
    p_gotest.add_argument("--algo", default="FAST-pw", help="FAST variant to use")
    p_gotest.add_argument("--repetitions", type=int, default=3)

    return parser


def _safe_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _init_ant(project_dir: Path, algo: str, repetitions: int) -> int:
    build_xml = project_dir / "build.xml"
    if not build_xml.exists():
        print(f"No build.xml found under {project_dir}")
        return 2

    original = build_xml.read_text(encoding="utf-8")
    updated = original

    def replace_or_insert_target(text: str, target_xml: str) -> str:
        # Match any number of our comment lines preceding the fast-tcp target, plus the target itself
        pattern = r"(?:\s*<!--\s*FAST TCP initialization target\s*-->\s*)*<target\s+name=\"fast-tcp\"[\s\S]*?</target>"
        if re.search(pattern, text):
            return re.sub(pattern, target_xml.strip(), text, count=1)
        return text.replace("\n</project>", f"\n{target_xml}\n</project>")

    # Copy packaged macros to .fast/fast-tcp.xml
    macros_src = Path(__file__).parent / "integrations" / "ant" / "fast-tcp.xml"
    macros_dst = project_dir / ".fast" / "fast-tcp.xml"
    macros_dst.parent.mkdir(parents=True, exist_ok=True)
    macros_dst.write_text(macros_src.read_text(encoding="utf-8"), encoding="utf-8")

    # Ensure import present near top
    if "fast-tcp.xml" not in updated:
        lines = updated.splitlines()
        for i, line in enumerate(lines):
            if line.strip().startswith("<project"):
                insert_idx = i + 1
                break
        else:
            insert_idx = 1
        lines.insert(insert_idx, '    <import file=".fast/fast-tcp.xml"/>')
        updated = "\n".join(lines) + ("\n" if updated.endswith("\n") else "")

    # Ensure/replace fast-tcp target
    new_target = (
        "    <!-- FAST TCP initialization target -->\n"
        '    <target name="fast-tcp">\n'
        f'        <fast-tcp-prioritize-and-run projectDir="${{basedir}}" algo="{algo}" repetitions="{repetitions}"/>\n'
        "    </target>"
    )
    updated = replace_or_insert_target(updated, new_target)

    if updated != original:
        backup = project_dir / ".fast" / "build.xml.bak"
        backup.parent.mkdir(parents=True, exist_ok=True)
        backup.write_text(original, encoding="utf-8")
        build_xml.write_text(updated, encoding="utf-8")
        print(f"Updated {build_xml} (backup at {backup})")
    else:
        print("No changes needed; FAST TCP already initialized.")

    print("Done. Run: ant fast-tcp")
    return 0


def _pytest_ini_update_addopts(existing: str, new_tokens: List[str]) -> str:
    import re as _re

    # Normalize accidental one-liner: "[pytest] addopts = ..." -> two lines
    existing = _re.sub(
        r"^\[pytest\][ \t]*addopts\s*=\s*",
        "[pytest]\naddopts = ",
        existing,
        flags=_re.M,
    )

    lines = existing.splitlines(keepends=True)
    if not lines:
        lines = ["[pytest]\n"]

    # Find the [pytest] header
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "[pytest]":
            header_idx = i
            break

    addopts_line = "addopts = " + " ".join(new_tokens) + "\n"

    if header_idx is None:
        # Append a new section at the end
        if lines and not lines[-1].endswith("\n"):
            lines[-1] = lines[-1] + "\n"
        lines.extend(["[pytest]\n", addopts_line])
        result = "".join(lines)
        if not result.endswith("\n"):
            result += "\n"
        return result

    # Ensure header ends with a newline
    if not lines[header_idx].endswith("\n"):
        lines[header_idx] = lines[header_idx].rstrip("\n") + "\n"

    # Determine the end of the [pytest] section
    j = header_idx + 1
    while j < len(lines) and not lines[j].lstrip().startswith("["):
        j += 1

    # Locate addopts within the section
    addopts_idx = None
    for k in range(header_idx + 1, j):
        if _re.match(r"^[ \t]*addopts\s*=", lines[k]):  # do not consume newlines
            addopts_idx = k
            break

    if addopts_idx is None:
        lines.insert(header_idx + 1, addopts_line)
    else:
        current = _re.sub(r"^\s*addopts\s*=\s*", "", lines[addopts_idx]).strip()
        tokens = current.split()
        for t in new_tokens:
            if t not in tokens:
                tokens.append(t)
        lines[addopts_idx] = "addopts = " + " ".join(tokens) + "\n"

    result = "".join(lines)
    if not result.endswith("\n"):
        result += "\n"
    return result


def _init_pytest(
    project_dir: Path, *, enable: bool, algo: str, r: int, b: int, k: int, budget: int
) -> int:
    cfg_path = project_dir / "pytest.ini"
    if not cfg_path.exists():
        cfg_path.parent.mkdir(parents=True, exist_ok=True)
        cfg_text = "[pytest]\n"
        if enable:
            addopts = [
                "--fast-tcp",
                "--fast-tcp-algo",
                algo,
                "--fast-tcp-r",
                str(r),
                "--fast-tcp-b",
                str(b),
                "--fast-tcp-k",
                str(k),
                "--fast-tcp-budget",
                str(budget),
            ]
            cfg_text += "addopts = " + " ".join(addopts) + "\n"
        _safe_write(cfg_path, cfg_text)
        print(f"Created {cfg_path}")
    else:
        text = cfg_path.read_text(encoding="utf-8")
        if enable:
            tokens = [
                "--fast-tcp",
                "--fast-tcp-algo",
                algo,
                "--fast-tcp-r",
                str(r),
                "--fast-tcp-b",
                str(b),
                "--fast-tcp-k",
                str(k),
                "--fast-tcp-budget",
                str(budget),
            ]
            updated = _pytest_ini_update_addopts(text, tokens)
            if updated != text:
                cfg_path.write_text(updated, encoding="utf-8")
                print(f"Updated {cfg_path}")
            else:
                print("No changes needed; Pytest already configured for FAST TCP.")
        else:
            print(
                "Pytest plugin is auto-discovered; run pytest with --fast-tcp to enable."
            )

    print("Done. Run: pytest")
    return 0


def _init_vitest(
    project_dir: Path, *, algo: str, repetitions: int, no_scripts: bool
) -> int:
    fast_dir = project_dir / ".fast"
    in_dir = fast_dir / "in"
    out_dir = fast_dir / "out"
    scripts_dir = fast_dir / "tools" / "vitest"
    scripts_dir.mkdir(parents=True, exist_ok=True)
    in_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Materialize helper scripts from packaged templates
    templates_dir = (
        Path(__file__).resolve().parent / "integrations" / "vitest" / "templates"
    )
    mapping = {
        templates_dir / "tests-map.cjs": scripts_dir / "tests-map.cjs",
        templates_dir / "generate-bbox.cjs": scripts_dir / "generate-bbox.cjs",
        templates_dir / "run-fast.cjs": scripts_dir / "run-fast.cjs",
    }
    for src, dst in mapping.items():
        if not src.exists():
            print(f"Warning: missing template {src}")
            continue
        dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")

    # Copy USAGE.md template
    usage_src = templates_dir / "USAGE.md"
    usage_dst = fast_dir / "USAGE.md"
    if usage_src.exists():
        usage_dst.write_text(usage_src.read_text(encoding="utf-8"), encoding="utf-8")

    # Optionally update package.json scripts
    if not no_scripts:
        pkg = project_dir / "package.json"
        if pkg.exists():
            import json

            try:
                data = json.loads(pkg.read_text(encoding="utf-8"))
            except Exception:
                data = {}
            scripts = data.get("scripts", {})
            scripts.setdefault("test", "vitest")
            scripts["test:fast"] = "node .fast/tools/vitest/run-fast.cjs"
            data["scripts"] = scripts
            pkg.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
            print("Updated package.json scripts: added 'test:fast'")
        else:
            print("No package.json found; skipped scripts update.")

    print("Done. Run: npm run test:fast")
    return 0


def _init_gotest(project_dir: Path, *, algo: str, repetitions: int) -> int:
    fast_dir = project_dir / ".fast"
    in_dir = fast_dir / "in"
    out_dir = fast_dir / "out"
    scripts_dir = fast_dir / "tools" / "gotest"
    scripts_dir.mkdir(parents=True, exist_ok=True)
    in_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    templates_dir = (
        Path(__file__).resolve().parent / "integrations" / "gotest" / "templates"
    )

    mapping = {
        templates_dir / "tests-map.py": scripts_dir / "tests-map.py",
        templates_dir / "generate-bbox.py": scripts_dir / "generate-bbox.py",
        templates_dir
        / "build-prioritized-tests.py": scripts_dir
        / "build-prioritized-tests.py",
        templates_dir / "run-fast.sh": scripts_dir / "run-fast.sh",
    }
    for src, dst in mapping.items():
        if not src.exists():
            print(f"Warning: missing template {src}")
            continue
        text = src.read_text(encoding="utf-8")
        # Inject defaults for algo/repetitions into runner script
        if dst.name == "run-fast.sh":
            text = re.sub(r"^ALG=\".*?\"$", f'ALG="{algo}"', text, flags=re.M)
            text = re.sub(r"^REPS=\"?\d+\"?$", f"REPS={repetitions}", text, flags=re.M)
        dst.write_text(text, encoding="utf-8")
        if dst.suffix == ".sh":
            try:
                os.chmod(dst, 0o755)
            except Exception:
                pass

    usage_src = templates_dir / "USAGE.md"
    usage_dst = fast_dir / "USAGE.md"
    if usage_src.exists():
        usage_dst.write_text(usage_src.read_text(encoding="utf-8"), encoding="utf-8")

    print("Done. Run: bash .fast/tools/gotest/run-fast.sh")
    return 0


def run_init(argv: List[str]) -> int:
    parser = create_init_parser()
    args = parser.parse_args(argv)
    tool = args.tool
    if tool == "ant":
        project_dir = Path(args.project_dir).resolve()
        return _init_ant(project_dir, algo=args.algo, repetitions=args.repetitions)
    if tool == "pytest":
        project_dir = Path(args.project_dir).resolve()
        return _init_pytest(
            project_dir,
            enable=args.enable,
            algo=args.algo,
            r=args.r,
            b=args.b,
            k=args.k,
            budget=args.budget,
        )
    if tool == "vitest":
        project_dir = Path(args.project_dir).resolve()
        return _init_vitest(
            project_dir,
            algo=args.algo,
            repetitions=args.repetitions,
            no_scripts=args.no_scripts,
        )
    if tool == "gotest":
        project_dir = Path(args.project_dir).resolve()
        return _init_gotest(project_dir, algo=args.algo, repetitions=args.repetitions)
    print(f"Unknown init tool: {tool}")
    return 2


def discover_test_files(test_dir, pattern="*.txt", file_naming="auto", entity=None):
    """
    Discover test files in a directory based on patterns and naming conventions.

    Args:
        test_dir (str): Directory to search for test files
        pattern (str): File pattern to match (default: *.txt)
        file_naming (str): How to identify entity type from filenames
        entity (str): Required entity type for filtering

    Returns:
        dict: Mapping of entity types to file paths
    """
    test_dir = Path(test_dir)
    if not test_dir.exists():
        raise FileNotFoundError(f"Test directory not found: {test_dir}")

    # Find all files matching the pattern
    pattern_path = test_dir / pattern
    all_files = glob.glob(str(pattern_path))

    if not all_files:
        raise ValueError(f"No files found matching pattern {pattern} in {test_dir}")

    discovered_files = {}

    if file_naming == "auto":
        # Try to automatically detect entity type from filename
        for file_path in all_files:
            filename = Path(file_path).stem.lower()

            # Check for entity keywords in filename
            if any(keyword in filename for keyword in ["bbox", "black", "bb"]):
                discovered_files["bbox"] = file_path
            elif any(keyword in filename for keyword in ["function", "func", "fn"]):
                discovered_files["function"] = file_path
            elif any(keyword in filename for keyword in ["branch", "br"]):
                discovered_files["branch"] = file_path
            elif any(keyword in filename for keyword in ["line", "ln"]):
                discovered_files["line"] = file_path
            else:
                # If no entity type detected, try to infer from content
                entity_type = infer_entity_from_content(file_path)
                if entity_type:
                    discovered_files[entity_type] = file_path

    elif file_naming == "entity-suffix":
        # Expect files named like: dataset-entity.txt
        for file_path in all_files:
            filename = Path(file_path).stem
            if "-" in filename:
                parts = filename.split("-")
                potential_entity = parts[-1].lower()
                if potential_entity in ["bbox", "function", "branch", "line"]:
                    discovered_files[potential_entity] = file_path

    elif file_naming == "custom":
        # For custom naming, use the first file found for the specified entity
        if entity and all_files:
            discovered_files[entity] = all_files[0]

    # Filter by requested entity if specified
    if entity and entity in discovered_files:
        return {entity: discovered_files[entity]}
    elif entity and entity not in discovered_files:
        raise ValueError(f"No file found for entity type '{entity}' in {test_dir}")

    return discovered_files


def infer_entity_from_content(file_path):
    """
    Infer entity type from file content patterns.

    Args:
        file_path (str): Path to the test file

    Returns:
        str: Inferred entity type or None
    """
    try:
        with open(file_path, "r") as f:
            first_line = f.readline().strip()

        # Black-box files typically contain command-line like strings
        if (
            first_line.startswith("-P[")
            or "error" in first_line
            or any(c in first_line for c in ["-F[", "-C["])
        ):
            return "bbox"

        # White-box files contain space-separated numbers
        if first_line and all(part.isdigit() for part in first_line.split()):
            # For white-box, we'll default to 'function' if we can't be more specific
            return "function"

    except Exception:
        pass

    return None


def run_prioritization(args):
    """Run test case prioritization."""
    try:
        # Discover test files in the specified directory
        print(f"Discovering test files in: {args.test_dir}")
        discovered_files = discover_test_files(
            args.test_dir,
            pattern=args.pattern,
            file_naming=args.file_naming,
            entity=args.entity,
        )

        if not discovered_files:
            print(f"No test files found for entity '{args.entity}' in {args.test_dir}")
            return

        print(f"Found files: {discovered_files}")

        # Set up output directory
        if args.output_dir:
            output_dir = Path(args.output_dir)
        else:
            output_dir = Path.cwd() / "output"

        output_dir.mkdir(parents=True, exist_ok=True)

        # Create a temporary directory structure that mimics the expected input format
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create input directory structure
            temp_input_dir = temp_path / "input"
            temp_output_dir = temp_path / "output"
            temp_input_dir.mkdir(parents=True, exist_ok=True)
            temp_output_dir.mkdir(parents=True, exist_ok=True)

            # Copy discovered file to expected location
            entity_file = discovered_files[args.entity]
            # Create a custom dataset name based on the directory
            dir_name = Path(args.test_dir).name
            synthetic_dataset = f"{dir_name}_v1"
            dataset_dir = temp_input_dir / synthetic_dataset
            dataset_dir.mkdir(parents=True, exist_ok=True)

            # Copy the file with the expected naming convention
            target_file = dataset_dir / f"{dir_name}-{args.entity}.txt"
            shutil.copy2(entity_file, target_file)

            # Change to temporary directory to run prioritization
            original_cwd = os.getcwd()
            os.chdir(temp_path)

            try:
                # Mock sys.argv for the prioritize module
                sys.argv = [
                    "prioritize.py",
                    synthetic_dataset,
                    args.entity,
                    args.algo,
                    str(args.repetitions),
                ]

                print(
                    f"Running {args.algo} prioritization on {args.entity} with {args.repetitions} repetitions..."
                )

                # Import and run the prioritize functionality
                # Call the main function from prioritize module
                prioritize.main()

                # Copy results back to the desired output directory
                prioritized_dir = temp_output_dir / synthetic_dataset / "prioritized"
                if prioritized_dir.exists():
                    final_output_dir = output_dir / synthetic_dataset
                    if final_output_dir.exists():
                        shutil.rmtree(final_output_dir)
                    shutil.copytree(prioritized_dir, final_output_dir)
                    print(f"Results saved to: {final_output_dir}")
                else:
                    print("Warning: No prioritized results found")

            except SystemExit:
                # Handle sys.exit() calls in the original module
                pass
            finally:
                os.chdir(original_cwd)

    except Exception as e:
        print(f"Error running prioritization: {e}")
        sys.exit(1)


def main():
    """Main CLI entry point."""
    # Support a top-level 'init' subcommand without breaking existing flags-only usage
    if len(sys.argv) > 1 and sys.argv[1] == "init":
        return run_init(sys.argv[2:])

    parser = create_parser()
    args = parser.parse_args()
    run_prioritization(args)
    return 0


if __name__ == "__main__":
    main()
