#!/usr/bin/env node
const fs = require("fs");
const path = require("path");

const TEST_DEF_RE =
  /(?<kind>it|test)\s*(?:\.(?:only|skip))?\s*\(\s*(?<q>['\"`])(?<name>.*?)(?<q2>['\"`])\s*,/gs;

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

function findTestBody(fileSrc, testName) {
  for (const m of fileSrc.matchAll(TEST_DEF_RE)) {
    if (m.groups.name !== testName) continue;
    const idx = m.index + m[0].length;
    const funcStart = fileSrc.indexOf("{", idx);
    if (funcStart === -1) continue;
    let depth = 0;
    const chars = [];
    for (let i = funcStart; i < fileSrc.length; i++) {
      const ch = fileSrc[i];
      chars.push(ch);
      if (ch === "{") depth++;
      else if (ch === "}") {
        depth--;
        if (depth <= 0) break;
      }
    }
    const body = chars.join("");
    if (body) return body;
  }
  return null;
}

function tokenize(text) {
  const tokens = text.match(/[A-Za-z0-9_]+/g) || [];
  return tokens.join(" ");
}

function main() {
  const args = parseArgs(process.argv.slice(2));
  const testsDir = path.resolve(args["tests-dir"] || ".");
  const testsMap = path.resolve(args["tests-map"]);
  const outPath = path.resolve(args["out"]);
  if (!testsMap || !outPath) {
    console.error("--tests-map and --out are required");
    process.exit(2);
  }
  fs.mkdirSync(path.dirname(outPath), { recursive: true });
  const lines = fs
    .readFileSync(testsMap, "utf8")
    .split(/\r?\n/)
    .filter(Boolean);
  const cases = [];
  for (const line of lines) {
    const [file, name] = line.split("\t");
    const absFile = path.isAbsolute(file) ? file : path.resolve(testsDir, file);
    const src = readText(absFile);
    const body = findTestBody(src, name) || src;
    cases.push(tokenize(body));
  }
  if (cases.length === 0) {
    console.error("No tests discovered to generate bbox tokens.");
    process.exit(1);
  }
  fs.writeFileSync(outPath, cases.join("\n") + "\n", "utf8");
  console.log(`Wrote ${cases.length} test cases to ${outPath}`);
}

if (require.main === module) main();
