"""Recompute published statistics and verify complete execution from raw reports."""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import statistics
import zipfile
from common import HERE, PROJECTS, validate_execution


def require(condition, message):
    if not condition:
        raise ValueError(message)


def same(actual, expected, label):
    require(math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-10),
            f"{label}: {actual} != {expected}")


def audit(read, names, summary, synthetic=None):
    total = 0
    require(bool(summary["projects"]), "No project results")
    for project in summary["projects"]:
        name = project["project"]
        require(project["valid"] and project["commit"] == PROJECTS[name]["commit"], f"{name}: invalid project/revision")
        repetitions = project["repetitions"]
        require(repetitions >= 1, f"{name}: no repetitions")
        labels = ["screen-baseline", "cold"] + [f"{c}-{i}" for i in range(1, repetitions + 1) for c in ("baseline", "fast")]
        found = {p for p in names if p.startswith(f"runs/{name}/") and p.endswith(".meta.json") and not p.endswith("/collection.meta.json")}
        require(found == {f"runs/{name}/{label}.meta.json" for label in labels}, f"{name}: unexpected or missing runs")
        collection = read(f"runs/{name}/collection.json")
        collection_meta = read(f"runs/{name}/collection.meta.json")
        require(collection_meta["returncode"] == 0 and collection["exitstatus"] == 0 and not collection["collection_errors"], f"{name}: collection failed")
        expected = None
        records = {}
        for label in labels:
            raw = read(f"runs/{name}/{label}.json")
            meta = read(f"runs/{name}/{label}.meta.json")
            enabled = label == "cold" or label.startswith("fast-")
            actual = validate_execution(raw, enabled=enabled)
            require(meta["returncode"] == 0 and meta["valid"] and meta["wall_s"] > 0, f"{name}/{label}: failed process")
            require(set(raw["collected"]) == set(collection["collected"]), f"{name}/{label}: collection differs")
            require(actual == meta["outcome_by_id"] and dict(Counter(actual.values())) == meta["outcomes"], f"{name}/{label}: metadata differs from reports")
            require(len(actual) == meta["executed"] == meta["collected"], f"{name}/{label}: wrong counts")
            if expected is None:
                expected = actual
            require(actual == expected, f"{name}/{label}: test outcomes differ across conditions")
            if enabled:
                require(raw["fast"] == meta["fast"], f"{name}/{label}: FAST timing differs")
                require(all(math.isfinite(v) and v >= 0 for v in raw["fast"].values()), f"{name}/{label}: invalid timings")
            else:
                require(raw["fast"] is None, f"{name}/{label}: baseline unexpectedly used FAST")
            records[label] = meta
        require(len(expected) == project["collected"] and dict(Counter(expected.values())) == project["outcomes"], f"{name}: summary counts differ")
        require(project["verified_executions"] == len(labels), f"{name}: summary run count differs")
        for condition in ("baseline", "fast"):
            values = [records[f"{condition}-{i}"]["wall_s"] for i in range(1, repetitions + 1)]
            same(statistics.mean(values), project[f"{condition}_mean_s"], f"{name}/{condition}: mean")
            same(statistics.stdev(values) if len(values) > 1 else 0, project[f"{condition}_std_s"], f"{name}/{condition}: standard deviation")
        for component, field in (("driver_s", "warm_fast_driver_mean_s"), ("prioritization_s", "warm_prioritization_mean_s")):
            same(statistics.mean(records[f"fast-{i}"]["fast"][component] for i in range(1, repetitions + 1)), project[field], f"{name}/{field}")
        same(records["cold"]["wall_s"], project["cold_wall_s"], f"{name}: cold run")
        difference = project["fast_mean_s"] - project["baseline_mean_s"]
        same(difference, project["wall_difference_s"], f"{name}: wall difference")
        same(100 * difference / project["baseline_mean_s"], project["wall_difference_percent"], f"{name}: wall percentage")
        total += len(labels)
        print(f"{name}: {len(labels)} complete matching executions, {len(expected)} tests; statistics verified")
    if synthetic is not None:
        audit_synthetic(read, synthetic)
    print(f"Verified {total} complete real-project executions")


def audit_synthetic(read, synthetic, require_snapshots=False):
    require(len(synthetic) == 20, "Synthetic experiment must have 20 steps")
    warmup = read("synthetic/warmup.audit.json")
    validate_execution(warmup, enabled=True)
    require(len(warmup["collected"]) == 50, "Wrong synthetic warmup size")
    for step, row in enumerate(synthetic, 1):
        baseline = read(f"synthetic/{step:02d}-baseline.audit.json")
        fast = read(f"synthetic/{step:02d}-fast.audit.json")
        b, f = validate_execution(baseline), validate_execution(fast, enabled=True)
        require(b == f and set(b.values()) == {"passed"}, f"Synthetic step {step}: outcomes differ")
        require(row["iteration"] == step and row["tests"] == len(f) == 50 + 5 * step - 5 * (step // 10), f"Synthetic step {step}: wrong suite size")
        if require_snapshots:
            expected = {node.split("::", 1)[0] for node in f}
            observed = {path for path in fast.get("snapshot_files", []) if path.startswith("tests/") and path.endswith(".py")}
            require(expected == observed, f"Synthetic step {step}: snapshot files differ from tests")
        for key, value in fast["fast"].items():
            same(value, row[key], f"Synthetic step {step}/{key}")
    print("Synthetic: 20 baseline/FAST pairs plus warmup; all tests passed, recorded component timings verified")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, help="Audit a new run instead of the bundled paper results")
    args = parser.parse_args()
    if args.results_dir:
        root = args.results_dir.resolve()
        read = lambda name: json.loads((root / name).read_text())
        synthetic = read("synthetic/measurements.json") if (root / "synthetic/measurements.json").exists() else None
        if (root / "verified-summary.json").exists():
            audit(read, {p.relative_to(root).as_posix() for p in root.rglob("*.json")}, read("verified-summary.json"), synthetic)
        else:
            require(synthetic is not None, "No experiment results found")
        if synthetic is not None:
            audit_synthetic(read, synthetic, require_snapshots=True)
    else:
        data = HERE / "data"
        for name, digest in json.loads((data / "SHA256SUMS.json").read_text()).items():
            require(hashlib.sha256((data / name).read_bytes()).hexdigest() == digest, f"Checksum mismatch: {name}")
        with zipfile.ZipFile(data / "full-run-reports.zip") as archive:
            read = lambda name: json.loads(archive.read(name))
            summary = json.loads((data / "verified-summary.json").read_text())
            synthetic = json.loads((data / "synthetic-measurements.json").read_text())
            require(summary == read("verified-summary.json"), "Bundled summaries differ")
            require(synthetic == read("synthetic/measurements.json"), "Bundled synthetic measurements differ")
            audit(read, set(archive.namelist()), summary, synthetic)
        corrected = data / "corrected-synthetic"
        with zipfile.ZipFile(corrected / "reports.zip") as archive:
            print("Corrected synthetic run (isolated snapshot, retained separately from historical evidence):")
            audit_synthetic(lambda name: json.loads(archive.read(name)), json.loads((corrected / "measurements.json").read_text()), require_snapshots=True)


if __name__ == "__main__":
    main()
