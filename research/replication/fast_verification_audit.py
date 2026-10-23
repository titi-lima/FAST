"""Read-only pytest observer used identically in both benchmark conditions."""
import json
import os
from pathlib import Path
import sys
import time
import subprocess

state = {"collected": [], "reports": {}, "collection_errors": [], "fast": None}


def pytest_configure(config):
    from fast_tcp.integrations.pytest import plugin
    original = plugin._run_fast_blackbox
    def measured(**kwargs):
        start = time.perf_counter()
        result = original(**kwargs)
        state["fast"] = {
            "partition_s": result[0], "preparation_s": result[1],
            "prioritization_s": result[2], "driver_s": time.perf_counter()-start,
        }
        return result
    plugin._run_fast_blackbox = measured


def pytest_collection_finish(session):
    state["collected"] = [item.nodeid for item in session.items]


def pytest_collectreport(report):
    if report.failed:
        state["collection_errors"].append(report.nodeid)


def pytest_runtest_logreport(report):
    state["reports"].setdefault(report.nodeid, {})[report.when] = {
        "outcome": report.outcome, "duration": report.duration,
        "xfail": getattr(report, "wasxfail", None),
    }


def pytest_sessionfinish(session, exitstatus):
    if state["fast"] is not None:
        repo = Path(session.config.rootpath) / ".fast" / "snapshot.git"
        state["snapshot_files"] = subprocess.check_output(
            ["git", "--git-dir", str(repo), "ls-tree", "-r", "--name-only", "refs/snapshots/latest"],
            text=True).splitlines()
    state.update(exitstatus=int(exitstatus), python=sys.version,
                 testsfailed=session.testsfailed, testscollected=session.testscollected)
    Path(os.environ["FAST_AUDIT_FILE"]).write_text(json.dumps(state, indent=2))
