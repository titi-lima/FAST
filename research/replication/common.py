"""Paths, project revisions, and outcome checks shared by the experiment CLIs."""
from pathlib import Path
import json
import os
import math
import hashlib
import platform
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_OUTPUT = ROOT / "output/research/replication"
PROJECTS = json.loads((HERE / "projects.json").read_text())


def environment_python(output, name):
    return output / "envs" / name / ("Scripts/python.exe" if os.name == "nt" else "bin/python")


def provenance():
    return {
        "platform": platform.platform(), "logical_cpus": os.cpu_count(),
        "tool_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "tool_sha256": {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in sorted((ROOT / "fast_tcp").rglob("*.py"))},
        "harness_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob("*.py"))},
    }


def outcomes(data):
    result = {}
    for node, reports in data.get("reports", {}).items():
        setup, call, teardown = (reports.get(phase, {}) for phase in ("setup", "call", "teardown"))
        if setup.get("outcome") == "failed" or teardown.get("outcome") == "failed":
            value = "error"
        elif setup.get("xfail") is not None:
            value = "xfailed"
        elif setup.get("outcome") == "skipped":
            value = "skipped"
        elif call.get("xfail") is not None:
            value = "xfailed" if call.get("outcome") == "skipped" else "xpassed"
        else:
            value = call.get("outcome", "missing")
        result[node] = value
    return result


def validate_execution(data, enabled=False):
    actual = outcomes(data)
    collected = data.get("collected", [])
    if not collected or len(set(collected)) != len(collected) or set(collected) != set(actual) or data.get("testscollected") != len(collected):
        raise ValueError("Collected tests and executed test identities do not match")
    if data.get("exitstatus") != 0 or data.get("collection_errors") or data.get("testsfailed"):
        raise ValueError("Test execution or collection failed")
    for node, value in actual.items():
        phases = data["reports"][node]
        if value in ("error", "failed", "missing") or not {"setup", "teardown"} <= phases.keys():
            raise ValueError(f"Incomplete or unsuccessful execution: {node}: {value}")
    if enabled and not isinstance(data.get("fast"), dict):
        raise ValueError("FAST did not run")
    if enabled:
        for key in ("partition_s", "preparation_s", "prioritization_s", "driver_s"):
            value = data["fast"].get(key)
            if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ValueError(f"Missing or invalid FAST timing: {key}")
    return actual
