"""Read-only execution evidence for the terminal recording (no ordering hooks)."""
import hashlib
import json
from pathlib import Path

state = {"collected": [], "reports": {}}


def pytest_collection_finish(session):
    state["collected"] = [item.nodeid for item in session.items]


def pytest_runtest_logreport(report):
    state["reports"].setdefault(report.nodeid, {})[report.when] = report.outcome


def pytest_sessionfinish(session, exitstatus):
    root = Path(session.config.rootpath)
    output = root / ".fast" / "last-demo-run.json"
    output.parent.mkdir(exist_ok=True)
    state.update(exitstatus=int(exitstatus),
                 source_sha256={str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted(root.rglob("*.py")) if ".fast" not in p.parts})
    output.write_text(json.dumps(state, indent=2))
