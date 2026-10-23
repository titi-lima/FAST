"""Create disposable, revision-checked project checkouts and frozen environments."""
import argparse
import shutil
import subprocess
from pathlib import Path
from common import DEFAULT_OUTPUT, HERE, PROJECTS, ROOT, environment_python


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--projects", nargs="+", choices=PROJECTS, default=list(PROJECTS))
    args = parser.parse_args()
    output = args.output_dir.resolve()
    uv = shutil.which("uv")
    if not uv:
        parser.error("Install uv and put it on PATH before preparing the experiments")
    for name in args.projects:
        spec = PROJECTS[name]
        checkout = output / "checkouts" / name
        if not checkout.exists():
            checkout.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(["git", "clone", "--depth", "1", "--branch", spec["ref"], spec["url"], str(checkout)], check=True)
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=checkout, text=True).strip()
        dirty = subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"], cwd=checkout, text=True)
        if commit != spec["commit"] or dirty:
            raise SystemExit(f"{name}: checkout differs from the pinned revision; use a new --output-dir")
        python = environment_python(output, name)
        if not python.exists():
            subprocess.run([uv, "venv", str(python.parent.parent), "--python", "3.11.13"], check=True)
        subprocess.run([uv, "pip", "install", "--python", str(python),
                        "-r", str(HERE / "requirements" / f"{name}.txt"),
                        "--no-deps", "-e", str(ROOT), "-e", str(checkout)], check=True)
        subprocess.run([uv, "pip", "check", "--python", str(python)], check=True)
        print(f"Prepared {name} at {commit}", flush=True)


if __name__ == "__main__":
    main()
