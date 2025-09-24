#!/usr/bin/env node
const { spawnSync } = require("child_process");
const fs = require("fs");
const path = require("path");

const ROOT = process.cwd();
const FAST_DIR = path.join(ROOT, ".fast");
const TOOLS = path.join(FAST_DIR, "tools", "vitest");
const IN_DIR = path.join(FAST_DIR, "in");
const OUT_DIR = path.join(FAST_DIR, "out");

function run(cmd, args, opts = {}) {
  const res = spawnSync(cmd, args, { stdio: "inherit", ...opts });
  if (res.status !== 0) process.exit(res.status || 1);
}

function ensureDirs() {
  fs.mkdirSync(IN_DIR, { recursive: true });
  fs.mkdirSync(OUT_DIR, { recursive: true });
}

function whichPython() {
  const candidates = [
    path.join(ROOT, "virtual", "bin", "python"),
    "python3",
    "python",
  ];
  for (const p of candidates) {
    const res = spawnSync(p, ["--version"]);
    if (res.status === 0) return p;
  }
  console.error("Python not found. Please install Python 3.");
  process.exit(1);
}

function escapeRegex(str) {
  return str.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}

function main() {
  ensureDirs();
  const PY = whichPython();

  // 1) Discover tests
  run("node", [
    path.join(TOOLS, "tests-map.cjs"),
    "--tests-dir",
    ROOT,
    "--out",
    path.join(IN_DIR, "test-names.tsv"),
  ]);

  // 2) Generate bbox
  run("node", [
    path.join(TOOLS, "generate-bbox.cjs"),
    "--tests-dir",
    ROOT,
    "--tests-map",
    path.join(IN_DIR, "test-names.tsv"),
    "--out",
    path.join(IN_DIR, "vitest-bbox.txt"),
  ]);

  // 3) Prioritize via Python CLI
  run(PY, [
    "-m",
    "fast_tcp.cli",
    "--test-dir",
    IN_DIR,
    "--algo",
    process.env.FAST_TCP_ALGO || "FAST-pw",
    "--entity",
    "bbox",
    "--repetitions",
    process.env.FAST_TCP_REPS || "3",
    "--file-naming",
    "entity-suffix",
    "--output-dir",
    OUT_DIR,
  ]);

  // 4) Map prioritized IDs to (file, name) and run Vitest
  const names = fs
    .readFileSync(path.join(IN_DIR, "test-names.tsv"), "utf8")
    .split(/\r?\n/)
    .filter(Boolean)
    .map((line) => {
      const [file, name] = line.split("\t");
      return { file, name };
    });

  function findLatestPickleDir(rootDir) {
    const datasets = fs
      .readdirSync(rootDir)
      .map((d) => path.join(rootDir, d))
      .filter((p) => fs.statSync(p).isDirectory());
    let latest = null;
    for (const ds of datasets) {
      const pr = path.join(ds, "prioritized");
      if (fs.existsSync(pr) && fs.statSync(pr).isDirectory()) {
        if (!latest || fs.statSync(pr).mtimeMs > fs.statSync(latest).mtimeMs)
          latest = pr;
      }
    }
    return latest;
  }

  const prDir = findLatestPickleDir(OUT_DIR) || OUT_DIR;
  const pkls = fs
    .readdirSync(prDir)
    .filter(
      (f) => f.startsWith("FAST-") && f.includes("-bbox-") && f.endsWith(".tsv")
    );
  let orderIdx = [];
  if (pkls.length) {
    const tsvPath = path.join(
      prDir,
      pkls.sort(
        (a, b) =>
          fs.statSync(path.join(prDir, b)).mtimeMs -
          fs.statSync(path.join(prDir, a)).mtimeMs
      )[0]
    );
    const lines = fs
      .readFileSync(tsvPath, "utf8")
      .split(/\r?\n/)
      .filter(Boolean);
    if (lines[0].includes("\t")) {
      orderIdx = names.map((_, i) => i + 1);
    } else {
      orderIdx = lines
        .map((x) => parseInt(x, 10))
        .filter((n) => Number.isFinite(n));
    }
  } else {
    orderIdx = names.map((_, i) => i + 1);
  }

  const prioritized = orderIdx.map((i) => names[i - 1]).filter(Boolean);

  for (const { file, name } of prioritized) {
    const pat = escapeRegex(name);
    console.log(`==> ${name} (${file})`);
    run("npx", ["--yes", "vitest", "run", file, "-t", pat, "--silent"]);
  }
}

if (require.main === module) main();
