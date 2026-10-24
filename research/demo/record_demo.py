"""Capture genuine shell/Vim actions for a three-minute pytest demonstration."""
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

import pexpect

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = Path(os.environ.get("FAST_DEMO_OUTPUT", ROOT / "output/video/three-minute"))
PROJECT = Path(os.environ.get("FAST_DEMO_PROJECT", "/tmp/fast-three-minute-demo/project"))
VENV = Path(os.environ.get("FAST_DEMO_VENV", sys.prefix))
PACE = float(os.environ.get("FAST_DEMO_PACE", "1.18"))
ANSI = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")
READY = "\x1b]777;ready\x07"


def seed_project():
    PROJECT.mkdir(parents=True, exist_ok=True)
    if any(PROJECT.iterdir()):
        raise SystemExit("Use a fresh FAST_DEMO_PROJECT directory")
    subprocess.run(["git", "init", "--quiet", str(PROJECT)], check=True)
    tests = PROJECT / "tests"
    tests.mkdir()
    (PROJECT / "pricing.py").write_text(
        "def subtotal(prices):\n    return round(sum(prices), 2)\n\n\n"
        "def shipping_cost(amount, express=False):\n"
        "    if express:\n        return 12\n"
        "    return 0 if amount >= 50 else 5\n\n\n"
        "def tax(amount, rate):\n    return round(amount * rate / 100, 2)\n\n\n"
        "def discounted_price(price, percent):\n    return price * (1 - percent / 10)\n\n\n"
        "def total_price(amount, rate=0, express=False):\n"
        "    return amount + tax(amount, rate) + shipping_cost(amount, express)\n")
    (PROJECT / "conftest.py").write_text(
        "import time\nimport pytest\n\n\n"
        "@pytest.fixture(autouse=True)\n"
        "def visible_progress():\n"
        "    # A short demo delay makes test execution visible.\n"
        "    time.sleep(0.03)\n")
    cases = {
        "subtotal": ("prices, expected", [([], 0), ([1], 1), ([1, 2], 3), ([2.5, 3.5], 6),
            ([10, 20], 30), ([0, 5], 5), ([1, 1, 1], 3), ([50, 50], 100), ([.1, .2], .3), ([100], 100)],
            "subtotal", "subtotal(prices)"),
        "shipping": ("amount, express, expected", [(0, False, 5), (10, False, 5), (49, False, 5),
            (50, False, 0), (100, False, 0), (0, True, 12), (10, True, 12), (49, True, 12),
            (50, True, 12), (100, True, 12)], "shipping_cost", "shipping_cost(amount, express)"),
        "tax": ("amount, rate, expected", [(0, 20, 0), (100, 0, 0), (100, 5, 5), (100, 10, 10),
            (100, 20, 20), (50, 10, 5), (50, 20, 10), (200, 5, 10), (25, 20, 5), (10, 10, 1)],
            "tax", "tax(amount, rate)"),
        "totals": ("amount, rate, express, expected", [(0, 0, False, 5), (10, 0, False, 15),
            (50, 0, False, 50), (100, 0, False, 100), (100, 20, False, 120),
            (50, 10, False, 55), (100, 0, True, 112), (100, 20, True, 132),
            (10, 10, True, 23), (0, 0, True, 12)], "total_price", "total_price(amount, rate, express)"),
        "zero_discount": ("price, expected", [(p, p) for p in (0, 5, 10, 25, 50, 75, 100, 250, 500, 1000)],
            "discounted_price", "discounted_price(price, 0)"),
    }
    for name, (params, rows, function, call) in cases.items():
        text = f"import pytest\nfrom pricing import {function}\n\n\n@pytest.mark.parametrize(\"{params}\", [\n"
        text += "".join(f"    {row!r},\n" for row in rows)
        text += f"])\ndef test_{name}({params}):\n    assert {call} == expected\n"
        (tests / f"test_{name}.py").write_text(text)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    seed_project()
    env = os.environ.copy()
    env.update(PATH=str(VENV / "bin") + os.pathsep + env.get("PATH", ""),
               TERM="xterm-256color", PS1="demo$ ", PS2="> ", COLUMNS="110", LINES="32",
               PROMPT_COMMAND="printf '\\033]777;ready\\007'", BASH_SILENCE_DEPRECATION_WARNING="1",
               PYTHONHASHSEED="0", PYTEST_ADDOPTS="", PY_COLORS="1",
               PYTHONPATH=str(HERE) + os.pathsep + str(PROJECT), PYTEST_PLUGINS="demo_observer",
               VIMINIT="set nocompatible number noswapfile nobackup nowritebackup viminfo= laststatus=0 | syntax on")
    env.pop("PYTEST_DISABLE_PLUGIN_AUTOLOAD", None)
    started = time.monotonic()
    cast = (OUT / "fast-tcp-demo.cast").open("w")
    cast.write(json.dumps({"version": 2, "width": 110, "height": 32,
                           "title": "fast-tcp pytest demonstration", "timestamp": int(time.time())}) + "\n")

    class Capture:
        def write(self, data):
            cast.write(json.dumps([round(time.monotonic() - started, 4), "o", data]) + "\n")
            cast.flush()
        def flush(self):
            cast.flush()

    terminal = pexpect.spawn("/bin/bash", ["--noprofile", "--norc", "-o", "pipefail", "-i"],
                            cwd=str(PROJECT), env=env, encoding="utf-8", dimensions=(32, 110), timeout=60)
    terminal.logfile_read = Capture()
    terminal.delaybeforesend = 0
    terminal.expect_exact(READY)
    terminal.expect_exact("demo$ ")
    evidence = {"runs": {}, "scenes": [], "fixture_sleep_per_test_s": .03, "initial_tests": 50}

    def drain():
        try:
            while True:
                terminal.read_nonblocking(16384, timeout=0)
        except pexpect.TIMEOUT:
            pass

    def pause(seconds):
        end = time.monotonic() + seconds * PACE
        while time.monotonic() < end:
            time.sleep(min(.035, max(0, end - time.monotonic())))
            drain()
        drain()

    def type_text(text, delay=.065):
        for character in text:
            terminal.send(character)
            pause(delay)

    def scene(name):
        evidence["scenes"].append({"scene": name, "start": round(time.monotonic() - started, 3)})
        print("SCENE " + name, flush=True)
        terminal.sendcontrol("l")
        terminal.expect_exact("demo$ ")
        pause(.15)

    def command(label, text, hold=1.6, total=None, passed=None, failed=0, first_prefix=None):
        type_text(text)
        terminal.sendline("")
        terminal.expect_exact(READY)
        captured = terminal.before
        terminal.expect_exact("demo$ ")
        clean = ANSI.sub("", captured)
        record = {"command": text, "output": captured}
        (OUT / f"{label}.txt").write_text(captured.replace("\r\n", "\n").replace("\r", "\n"))
        if total is not None:
            audit = json.loads((PROJECT / ".fast/last-demo-run.json").read_text())
            assert len(audit["collected"]) == total, (label, audit)
            outcomes = [report.get("call") for report in audit["reports"].values()]
            assert outcomes.count("passed") == passed and outcomes.count("failed") == failed, (label, outcomes)
            assert audit["exitstatus"] == (1 if failed else 0), label
            if not failed:
                assert len(outcomes) == total and set(outcomes) == {"passed"}, label
            if first_prefix:
                assert all(node.startswith(first_prefix) for node in audit["collected"][:10]), label
            record["audit"] = audit
            (OUT / f"{label}.audit.json").write_text(json.dumps(audit, indent=2))
        evidence["runs"][label] = record
        pause(hold)
        return clean

    def vim_open(path):
        type_text("vim " + path)
        terminal.sendline("")
        pause(.9)

    def vim_save():
        terminal.send("\x1b")
        pause(.35)
        terminal.send(":wq\r")
        terminal.expect_exact(READY)
        terminal.expect_exact("demo$ ")

    def trace(label, total, passed, partition=None):
        clean = command(label, "pytest -q --fast-tcp-debug",
                        total=total, passed=passed)
        if partition:
            assert partition in clean, (label, clean)
        return clean

    scene("pricing project")
    command("files", "ls tests")
    command("pricing", "cat pricing.py", hold=2)
    command("test_example", "cat tests/test_shipping.py", hold=2)
    command("demo_fixture", "cat conftest.py")

    scene("initialize pytest")
    command("initialize", "fast-tcp init pytest")
    command("config", "cat pytest.ini", hold=2)

    scene("cold cache and warm cache")
    trace("cold", 50, 50, "50 new / 0 unchanged / 0 deleted")
    trace("warm", 50, 50, "0 new / 50 unchanged / 0 deleted")
    command("cache", "ls .fast")
    command("snapshot", "git --git-dir=.fast/snapshot.git ls-tree -r --name-only refs/snapshots/latest")

    scene("add regression tests")
    vim_open("tests/test_zz_discount.py")
    terminal.send("i")
    code = ("from pricing import discounted_price\n\n"
            "def test_twenty_percent_discount():\n    assert discounted_price(100, 20) == 80\n\n"
            "def test_ten_percent_discount():\n    assert discounted_price(100, 10) == 90\n")
    type_text(code)
    pause(1.7)
    vim_save()
    assert (PROJECT / "tests/test_zz_discount.py").read_text().rstrip() == code.rstrip()

    scene("normal pytest order")
    normal = command("normal", "pytest -o addopts= -x -q", hold=2, total=52, passed=50, failed=1)
    assert "1 failed, 50 passed" in normal
    scene("new tests first")
    command("fast", "pytest -x -v", hold=2, total=52, passed=0, failed=1)
    assert evidence["runs"]["fast"]["audit"]["collected"][0].startswith("tests/test_zz_discount.py")

    scene("fix the discount calculation")
    vim_open("pricing.py")
    terminal.send(r"/percent \/ 10" + "\r")
    pause(.6)
    terminal.send("$hh")
    terminal.send("cw")
    type_text("100", .18)
    terminal.send("\x1b")
    pause(1.5)
    vim_save()
    assert "percent / 100)" in (PROJECT / "pricing.py").read_text()
    command("fixed", "pytest -q", total=52, passed=52)

    scene("change an existing tax test")
    vim_open("tests/test_tax.py")
    terminal.send("/25, 20\r")
    pause(.6)
    terminal.send("cw")
    type_text("24.99", .18)
    terminal.send("\x1b")
    pause(1.5)
    vim_save()
    assert "(24.99, 20, 5)" in (PROJECT / "tests/test_tax.py").read_text()
    command("modified_first", "pytest -v", total=52, passed=52,
            first_prefix="tests/test_tax.py")

    scene("consolidate the discount tests")
    vim_open("tests/test_zz_discounts.py")
    terminal.send("i")
    merged = ("import pytest\nfrom pricing import discounted_price\n\n"
              "@pytest.mark.parametrize(\"price, percent, expected\", [\n"
              + "".join(f"    ({price}, 0, {price}),\n" for price in (0, 5, 10, 25, 50, 75, 100, 250, 500, 1000))
              + "    (100, 20, 80),\n    (100, 10, 90),\n])\n"
              "def test_discount(price, percent, expected):\n"
              "    assert discounted_price(price, percent) == expected\n")
    type_text(merged)
    pause(1.8)
    vim_save()
    assert (PROJECT / "tests/test_zz_discounts.py").read_text().rstrip() == merged.rstrip()
    command("replacement_tests", "pytest -q", total=64, passed=64)
    command("remove_old_tests", "rm tests/test_zero_discount.py tests/test_zz_discount.py")
    command("refactored_files", "ls tests")
    trace("refactored", 52, 52, "12 new / 40 unchanged / 12 deleted")

    scene("select FAST-one")
    command("one", "pytest --fast-tcp-algo FAST-one -q", total=52, passed=52)
    scene("unchanged suite reuses the cache")
    trace("final_warm", 52, 52, "0 new / 52 unchanged / 0 deleted")

    duration = round(time.monotonic() - started, 3)
    terminal.logfile_read = None
    cast.write(json.dumps([duration, "o", ""]) + "\n")
    cast.close()
    terminal.close(force=True)
    evidence.update(verified=True, duration_s=duration, final_tests=52,
                    normal_failure_position=51, fast_failure_position=1,
                    note="Genuine shell/Vim capture. Constructed pricing example; not performance evidence.")
    (OUT / "evidence.json").write_text(json.dumps(evidence, indent=2))
    print(json.dumps({key: value for key, value in evidence.items() if key != "runs"}), flush=True)


if __name__ == "__main__":
    main()
