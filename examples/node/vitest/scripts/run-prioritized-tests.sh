#!/usr/bin/env bash
set -euo pipefail

# End-to-end: install deps, generate bbox, prioritize, map, then run Vitest in order

ROOT_DIR="$(cd "$(dirname "$0")"/../.. && pwd)"   # examples/node
EX_DIR="$ROOT_DIR/vitest"                               # this example root
REPO_DIR="$(cd "$ROOT_DIR"/../.. && pwd)"          # repo root
DATA_DIR="$EX_DIR/.fast/in"
OUT_DIR="$EX_DIR/.fast/out"

ALG="FAST-pw"
REPS=3

mkdir -p "$DATA_DIR" "$OUT_DIR"

PY="$REPO_DIR/virtual/bin/python"
if [ ! -x "$PY" ]; then
  PY="python3"
fi

echo "[1/6] Install Node dev deps"
(cd "$EX_DIR" && npm install --silent) | cat

echo "[2/6] Discover tests (file + fullName)"
"$PY" "$EX_DIR/scripts/tests-map.py" \
  --tests-dir "$EX_DIR" \
  --out "$DATA_DIR/test-names.tsv"

echo "[3/6] Generate bbox tokens from tests"
"$PY" "$EX_DIR/scripts/generate-bbox.py" \
  --tests-dir "$EX_DIR" \
  --tests-map "$DATA_DIR/test-names.tsv" \
  --out "$DATA_DIR/vitest-bbox.txt"

echo "[4/6] Prioritize with FAST TCP ($ALG, reps=$REPS)"
"$PY" -m fast_tcp.cli \
  --test-dir "$DATA_DIR" \
  --algo "$ALG" \
  --entity bbox \
  --repetitions "$REPS" \
  --file-naming entity-suffix \
  --output-dir "$OUT_DIR" | cat

echo "[5/6] Map prioritized IDs to Vitest (file + name)"
"$PY" "$EX_DIR/scripts/build-prioritized-tests.py" \
  --out-dir "$OUT_DIR" \
  --tests "$DATA_DIR/test-names.tsv" \
  --out "$DATA_DIR/prioritized-tests.tsv"

echo "[6/6] Run Vitest tests in prioritized order"
while IFS=$'\t' read -r file name; do
  [ -z "$file" ] && continue
  PATTERN=$("$PY" -c 'import re,sys; print(re.escape(sys.argv[1]))' "$name")
  echo "==> $name ($file)"
  (cd "$EX_DIR" && npx --yes vitest run "$file" -t "$PATTERN" --silent) | cat || true
done < "$DATA_DIR/prioritized-tests.tsv"

echo "Done."


