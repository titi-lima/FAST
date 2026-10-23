"""Measure a real 20-step mutating pytest workload, replacing the simulated trace."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from common import DEFAULT_OUTPUT, environment_python, validate_execution, provenance
OUT = DEFAULT_OUTPUT / "synthetic"
PYTHON = environment_python(DEFAULT_OUTPUT, "rich")


def run(project, label, enabled):
    audit = OUT / f"{label}.audit.json"
    env = os.environ.copy()
    env.update(PYTHONPATH=str(Path(__file__).parent), PYTEST_DISABLE_PLUGIN_AUTOLOAD="1",
               FAST_AUDIT_FILE=str(audit), PYTHONHASHSEED="0")
    cmd = [str(PYTHON),"-m","pytest","-p","fast_tcp.integrations.pytest.plugin",
           "-p","fast_verification_audit","-o","addopts=","tests","-q"]
    if enabled: cmd += ["--fast-tcp"]
    start=time.perf_counter()
    result=subprocess.run(cmd,cwd=project,env=env,text=True,capture_output=True,timeout=30)
    elapsed=time.perf_counter()-start
    (OUT/f"{label}.stdout.txt").write_text(result.stdout)
    (OUT/f"{label}.stderr.txt").write_text(result.stderr)
    data=json.loads(audit.read_text())
    if result.returncode:
        raise RuntimeError(f"{label} failed; see {OUT}")
    validate_execution(data, enabled=enabled)
    if enabled:
        expected = {p.relative_to(project).as_posix() for p in (project / "tests").glob("*.py")}
        observed = {p for p in data.get("snapshot_files", []) if p.startswith("tests/") and p.endswith(".py")}
        if expected != observed:
            raise RuntimeError(f"{label}: snapshot does not contain the generated test files")
    if any(r.get("call", {}).get("outcome") != "passed" for r in data["reports"].values()):
        raise RuntimeError(f"{label}: synthetic tests did not all pass")
    return data,elapsed


def main():
    global OUT, PYTHON
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--python", type=Path, help="Prepared Python with fast-tcp and pytest (defaults to rich environment)")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    OUT = output / "synthetic"
    PYTHON = args.python.absolute() if args.python else environment_python(output, "rich")
    OUT.mkdir(parents=True,exist_ok=True)
    project=OUT/"project"
    if project.exists():raise SystemExit("Preserve or move an earlier synthetic run before rerunning.")
    (OUT / "metadata.json").write_text(json.dumps(provenance(), indent=2))
    freeze = subprocess.check_output([str(PYTHON), "-c", "import importlib.metadata as m; print('\\n'.join(sorted(d.metadata['Name'] + '==' + d.version for d in m.distributions())))"], text=True)
    (OUT / "environment.txt").write_text(freeze)
    tests=project/"tests";tests.mkdir(parents=True)
    # Isolate the workload from the parent repository's ignored output/ path.
    subprocess.run(["git", "init", "--quiet", str(project)], check=True)
    def write_file(path, index):
        path.write_text("\n".join(f"def test_case_{i}():\n    assert {index}+{i} >= 0\n" for i in range(5)))
    for i in range(10):write_file(tests/f"test_group_{i:02d}.py",i)
    (project/"pytest.ini").write_text("[pytest]\n")
    run(project,"warmup",True)
    rows=[]
    for i in range(1,21):
        write_file(tests/f"test_zz_extra_{i:02d}.py",100+i)
        modified=tests/f"test_group_{i%10:02d}.py"
        with modified.open("a") as f:f.write(f"\n# mutation {i}\n")
        if i in (10,20):(tests/f"test_zz_extra_{i-9:02d}.py").unlink()
        baseline,bt=run(project,f"{i:02d}-baseline",False)
        fast,ft=run(project,f"{i:02d}-fast",True)
        assert set(baseline["collected"])==set(fast["collected"])
        row={"iteration":i,"tests":len(fast["collected"]),"baseline_wall_s":bt,
             "fast_wall_s":ft,**fast["fast"]}
        rows.append(row)
        (OUT/"measurements.json").write_text(json.dumps(rows,indent=2))
        print(json.dumps(row),flush=True)
    assert rows[-1]["tests"]==140
    print("SYNTHETIC_VERIFIED 20 mutation steps; every execution passed; final suite 140 tests",flush=True)


if __name__=="__main__":main()
