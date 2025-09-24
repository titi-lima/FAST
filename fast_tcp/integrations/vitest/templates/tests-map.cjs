#!/usr/bin/env node
const fs = require("fs");
const path = require("path");

const TEST_NAME_RE =
  /(?:it|test)\s*(?:\.(?:only|skip))?\s*\(\s*([\'\"`])(?<name>.*?)\1\s*,/gs;

function readText(p) {
  return fs.readFileSync(p, "utf8");
}

function parseArgs(argv) {
  const out = {};
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    if (a.startsWith("--")) {
      const key = a.slice(2);
      const val =
        argv[i + 1] && !argv[i + 1].startsWith("--") ? argv[++i] : true;
      out[key] = val;
    }
  }
  return out;
}

const TEST_EXTS = new Set([".js", ".ts", ".mjs", ".cjs", ".jsx", ".tsx"]);

function isTestFile(filePath) {
  const base = path.basename(filePath);
  const ext = path.extname(base);
  if (!TEST_EXTS.has(ext)) return false;
  return /\.(test|spec)\.[^.]+$/.test(base);
}

function discoverTestFiles(root) {
  const out = [];
  function walk(dir) {
    const entries = fs.readdirSync(dir, { withFileTypes: true });
    for (const e of entries) {
      const p = path.join(dir, e.name);
      if (e.isDirectory()) {
        const name = e.name;
        if (
          name === "node_modules" ||
          name === "dist" ||
          name === "build" ||
          name === ".fast"
        )
          continue;
        walk(p);
      } else if (e.isFile()) {
        if (isTestFile(p)) out.push(path.resolve(p));
      }
    }
  }
  walk(root);
  return Array.from(new Set(out)).sort();
}

function extractTestNames(filePath) {
  const src = readText(filePath);
  const names = [];
  for (const m of src.matchAll(TEST_NAME_RE)) {
    names.push(m.groups.name);
  }
  return names;
}

function main() {
  const args = parseArgs(process.argv.slice(2));
  const testsDir = path.resolve(args["tests-dir"] || ".");
  const outPath = path.resolve(args["out"]);
  if (!outPath) {
    console.error("--out is required");
    process.exit(2);
  }
  fs.mkdirSync(path.dirname(outPath), { recursive: true });
  const pairs = [];
  for (const file of discoverTestFiles(testsDir)) {
    for (const name of extractTestNames(file)) {
      pairs.push(`${file}\t${name}`);
    }
  }
  fs.writeFileSync(outPath, pairs.join("\n") + "\n", "utf8");
  console.log(`Wrote ${pairs.length} test entries to ${outPath}`);
}

if (require.main === module) main();
