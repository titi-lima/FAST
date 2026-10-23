"""Run pinned OSS suites; persist test identities and all outcomes before timing claims."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import time

from common import DEFAULT_OUTPUT, HERE, PROJECTS, ROOT, environment_python, validate_execution

OUT = DEFAULT_OUTPUT
PLUGINS = {name: spec["plugins"] for name, spec in PROJECTS.items()}


def execute(name, condition, label, collect=False):
    project = OUT / "checkouts" / name
    py = environment_python(OUT, name)
    dest = OUT / "runs" / name
    dest.mkdir(parents=True, exist_ok=True)
    audit = dest / f"{label}.json"
    env = os.environ.copy()
    env.update(PYTHONPATH=str(Path(__file__).parent), PYTEST_DISABLE_PLUGIN_AUTOLOAD="1",
               FAST_AUDIT_FILE=str(audit), PYTHONHASHSEED="0", TZ="UTC",
               HYPOTHESIS_STORAGE_DIRECTORY=str(OUT / "hypothesis" / name))
    env["PATH"] = str(py.parent) + os.pathsep + env.get("PATH", "")
    cmd = [str(py), "-m", "pytest", "-p", "fast_tcp.integrations.pytest.plugin",
           "-p", "fast_verification_audit"]
    for plugin in PLUGINS[name]:
        cmd += ["-p", plugin]
    config = project / ("pyproject.toml" if name in ("rich", "attrs") else "setup.cfg")
    cmd += ["-c", str(config), "--rootdir", str(project), "-o", "addopts=", "tests", "-q", "--tb=short"]
    if condition == "fast":
        cmd += ["--fast-tcp"]
    if collect:
        cmd += ["--collect-only"]
    start = time.perf_counter()
    with (dest / f"{label}.stdout.txt").open("w") as stdout, (dest / f"{label}.stderr.txt").open("w") as stderr:
        try:
            proc = subprocess.run(cmd, cwd=project, env=env, stdout=stdout, stderr=stderr, timeout=600)
            code = proc.returncode
        except subprocess.TimeoutExpired:
            code = 124
    elapsed = time.perf_counter() - start
    data = json.loads(audit.read_text()) if audit.exists() else {}
    from common import outcomes as classify_outcomes
    outcomes = classify_outcomes(data)
    collected = data.get("collected", [])
    complete = set(collected) == set(outcomes) and len(collected) == len(outcomes)
    try:
        validate_execution(data, enabled=condition == "fast")
        valid = code == 0
    except ValueError:
        valid = False
    if collect:
        valid = code == 0 and bool(collected) and not data.get("collection_errors")
    meta = {"project": name, "condition": condition, "label": label, "command": cmd,
            "wall_s": elapsed, "returncode": code, "valid": valid,
            "complete": complete, "collected": len(collected), "executed": len(outcomes),
            "outcomes": dict(Counter(outcomes.values())), "outcome_by_id": outcomes,
            "fast": data.get("fast")}
    (dest / f"{label}.meta.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps({k: v for k, v in meta.items() if k not in ("command", "outcome_by_id")}), flush=True)
    return meta


def benchmark(name, repetitions):
    project = OUT / "checkouts" / name
    if (OUT / "runs" / name).exists():
        raise SystemExit(f"Existing results for {name}; use a new --output-dir to preserve them")
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=project, text=True).strip()
    dirty = subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"], cwd=project, text=True)
    if commit != PROJECTS[name]["commit"] or dirty:
        raise SystemExit(f"{name}: wrong project revision; run prepare_projects.py first")
    # Only this harness's disposable checkout cache is reset. No source changes.
    fast_dir = project / ".fast"
    if fast_dir.exists():
        shutil.rmtree(fast_dir)
    # Initialize snapshot storage without changing the project's pytest config.
    py = environment_python(OUT, name)
    subprocess.run([str(py), "-c", "from fast_tcp.snapshot_cache import initialize_snapshot_cache; import sys; initialize_snapshot_cache(sys.argv[1])", str(project)], check=True)
    collected = execute(name, "baseline", "collection", collect=True)
    if not collected["valid"]:
        return {"project": name, "valid": False, "reason": "collection failed"}
    baseline = execute(name, "baseline", "screen-baseline")
    cold = execute(name, "fast", "cold")
    if not baseline["valid"] or not cold["valid"] or baseline["outcome_by_id"] != cold["outcome_by_id"]:
        return {"project": name, "valid": False, "reason": "screening or outcome equivalence failed"}
    records = []
    for rep in range(1, repetitions + 1):
        conditions = ["baseline", "fast"] if rep % 2 else ["fast", "baseline"]
        for condition in conditions:
            record = execute(name, condition, f"{condition}-{rep}")
            if not record["valid"] or record["outcome_by_id"] != baseline["outcome_by_id"]:
                return {"project": name, "valid": False, "reason": "incomplete run or outcome mismatch"}
            records.append(record)
    base = [r["wall_s"] for r in records if r["condition"] == "baseline"]
    fast = [r["wall_s"] for r in records if r["condition"] == "fast"]
    b, f = statistics.mean(base), statistics.mean(fast)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=project, text=True).strip()
    freeze = subprocess.check_output([str(py), "-c",
        "import importlib.metadata as m; print('\\n'.join(sorted(d.metadata['Name'] + '==' + d.version for d in m.distributions())))"], text=True)
    (OUT / "runs" / name / "environment.txt").write_text(freeze)
    return {"project": name, "valid": True, "commit": commit, "collected": baseline["collected"],
            "outcomes": baseline["outcomes"], "repetitions": repetitions,
            "verified_executions": 2 + 2 * repetitions,
            "baseline_mean_s": b, "baseline_std_s": statistics.stdev(base) if len(base)>1 else 0,
            "fast_mean_s": f, "fast_std_s": statistics.stdev(fast) if len(fast)>1 else 0,
            "wall_difference_s": f-b, "wall_difference_percent": 100*(f-b)/b,
            "cold_wall_s": cold["wall_s"], "warm_fast_driver_mean_s": statistics.mean(r["fast"]["driver_s"] for r in records if r["condition"]=="fast"),
            "warm_prioritization_mean_s": statistics.mean(r["fast"]["prioritization_s"] for r in records if r["condition"]=="fast")}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--projects", nargs="+", choices=PROJECTS, default=list(PROJECTS))
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--repetitions", type=int, default=5)
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error("--repetitions must be positive")
    OUT = args.output_dir.resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "verified-summary.json"
    result = json.loads(path.read_text()) if path.exists() else {"platform": platform.platform(), "projects": []}
    result["tool_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/"fast_tcp").rglob("*.py"))}
    result["tool_commit"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    result["tool_dirty"] = bool(subprocess.check_output(["git", "status", "--porcelain", "--", "fast_tcp"], cwd=ROOT, text=True).strip())
    result["hardware"] = {"platform": platform.platform(), "machine": platform.machine(), "processor": platform.processor(), "logical_cpus": os.cpu_count()}
    result["harness_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob("*.py"))}
    for name in args.projects:
        print("VERIFY_PROJECT " + name, flush=True)
        summary = benchmark(name, args.repetitions)
        result["projects"] = [r for r in result["projects"] if r["project"] != name] + [summary]
        path.write_text(json.dumps(result, indent=2))
        print("PROJECT_RESULT " + json.dumps(summary), flush=True)

    if not all(project["valid"] for project in result["projects"]):
        raise SystemExit(1)
